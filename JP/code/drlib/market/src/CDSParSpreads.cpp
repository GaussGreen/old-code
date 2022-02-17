//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CDSParSpreads.cpp
//
//   Description : Holds the current par spreads for CDSs
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 30, 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CDSPARSPREADS_CPP
#include "edginc/CDSParSpreads.hpp"
#include "edginc/DefaultRates.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/PrepayCurve.hpp"
#include "edginc/ZeroSysEqn.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/Addin.hpp"
#include "edginc/NonPricingModel.hpp"



DRLIB_BEGIN_NAMESPACE

CDSParSpreads::~CDSParSpreads(){}

CDSParSpreads::CDSParSpreads(): 
    CDSParSpreadsBase(TYPE), hasDefaulted(false)
{}

/** A constructor for convenience to be called by other QLib classs */
CDSParSpreads::CDSParSpreads(string                iname,
                             DateTime              ivalueDate,
                             int                   ispotOffset,
                             double                iparRecovery,
                             int                   iparSwapFreq,
                             ExpiryArraySP         iexpiries,
                             CDoubleArraySP        ispreads,
                             CDoubleArraySP        iupfronts,
                             bool                  iparAccrueFee,
                             string                iparDCC,       // day count convention for par CDS
                             string                iparBDC,       // bad day convention for par CDS
                             HolidayConstSP        iparHols,      // holidays for par CDS
                             YieldCurveConstSP     idiscount,     // corresponding discount curve
                             DecretionCurveConstSP iprepay        // prepay curve
                             )
    :
    hasDefaulted(false), CDSParSpreadsBase(iname,
                                           ivalueDate,
                                           ispotOffset,
                                           iparRecovery,
                                           iparSwapFreq,
                                           iexpiries,
                                           ispreads,
                                           iupfronts,
                                           iparAccrueFee,
                                           iparDCC,
                                           iparBDC,  
                                           iparHols, 
                                           idiscount,
                                           iprepay)
{}

/** adjust dates */
class AdjustDatesFunc : public ZeroSysEqnFunc
{
    const CDSParSpreads& parSpreadCurve;
    DefaultRatesSP       defaultRates;
    DateTimeArray        endDates;
    DateTimeArray&       adjEndDates;
    int                  count;
    string               ivl;

public:
    AdjustDatesFunc(
        const CDSParSpreads& pParSpreadCurve, 
        DefaultRatesSP       pDefaultRates,
        DateTimeArray&       pAdjEndDates)
        : parSpreadCurve(pParSpreadCurve),
          defaultRates(pDefaultRates),
          endDates(pParSpreadCurve.getExpiries()->size()),
          adjEndDates(pAdjEndDates)
    {
        MaturityPeriod period(parSpreadCurve.getSwapFrequency());
        period.decompose(count, ivl);

        // pre-calculate unadjusted end dates
        for (int i = 0 ; i < endDates.size() ; i++)
        {
            endDates[i] = (*pParSpreadCurve.getExpiries())[i]->toDate(pParSpreadCurve.getEffDate());
        }
    }

    ~AdjustDatesFunc()
    {
    }

    void operator()(int n, double* x, double* errors)
    {
        static const string method = "AdjustDatesFunc::operator()";

        // validate current guess
        for (int i = 0 ; i < n ; i++)
        {
            if (x[i] < 0.0)
            {
                string msg = Format::toString(
                    "Root solver implies negative spread (%.2f%%) for date %s",
                    1e2 * x[i], adjEndDates[i].toString().c_str());
                throw ModelException(method, msg);
            }
        }

        // recreate clean spread curve
        DoubleArray adjSpreads(n);
        for (int k = 0 ; k < adjSpreads.size() ; k++)
            adjSpreads[k] = x[k];

        defaultRates = DefaultRatesSP(
            new CDSHelper::CParSpreadDefaultRates(defaultRates->getValueDate(), adjEndDates, adjSpreads));

        // recalculate errors
        for (int j = 0 ; j < n ; j++)
        {
            DoubleArrayConstSP prices = parSpreadCurve.getUpfrontFees();
            double price = prices.get() ? (*prices)[j] : 0.0;
            double recalculated = calcParSpread(price, endDates[j]);
            errors[j] = recalculated - (*parSpreadCurve.getParSpreads())[j];
        }
    }

    // calculate par spread from current clean spread curve
    double calcParSpread(double price, const DateTime& maturity)
    {
        // clean spread => par spread
        DateTime baseDate = parSpreadCurve.getEffDate();
        DateTimeArraySP paymentDates(SwapTool::paymentDates(baseDate, maturity, count, ivl, true));

        CDSPricer pricer(
            1.0,                            // fee, 1.0 so fee value = annuity
            *paymentDates,
            defaultRates,
            1.0,                            // notional
            parSpreadCurve.getRecovery(),
            parSpreadCurve.getAccrueFee(),
            DayCountConventionSP(parSpreadCurve.dayCountConv()),
            parSpreadCurve.getValueDate(),
            parSpreadCurve.getEffDate(),    // protection start date
            maturity,                       // protection end date
            parSpreadCurve.getEffDate(),    // accrued start date
            parSpreadCurve.getYieldCurve(),
            parSpreadCurve.getPrepayCurve()
            );

        return pricer.impliedParSpread(price);
    }
};

CDSParSpreadsSP CDSParSpreads::adjustDates(string name, ExpiryArraySP expiries) const
{
    static const string method = "CDSParSpreads::adjustDates";

    if (expiries->size() != getExpiries()->size())
    {
        string msg = Format::toString(
            "number of adjusted expiries (%d) must match number of expiries (%d)",
            expiries->size(), getExpiries()->size());
        throw ModelException(method, msg);
    }

    /*
     * 1. Bootstrap the Cds curve using the market rates.
     * 2. Imply the adjSpreads from this curve and use this as the initial
     *    guess.
     * 3. Solve for adjSpreads such that the curve bootstrapped from this
     *    will return the market rates.
     */
    DateTime baseDate = getEffDate();

    /*
     * Instead of going through the par spreads we will use the clean
     * spreads internally. This should enable us to skip a whole bunch
     * of bootstrap routines. The guess spreads will be expressed as
     * continuously compounded rates.
     *
     * This relies on knowing that the dates for the clean spread curve
     * are the same dates as the adjusted spread curve.
     */
    int size = expiries->size();
    DoubleArray guess_spreads(size);
    DateTimeArray adjEndDates(size);
    CleanSpreadCurveSP marketCreditCurve = CDSHelper::getCleanSpreadCurve(*this, getValueDate(), true);
    for (int i = 0 ; i < size ; i++)
    {
        adjEndDates[i] = (*expiries)[i]->toDate(baseDate);
        guess_spreads[i] = marketCreditCurve->getCleanSpread(adjEndDates[i], CompoundBasis::CONTINUOUS);
    }

    // perform numerical estimation
    AdjustDatesFunc func(*this, defaultRates(), adjEndDates);
    ZeroSysEqnSolver solver(1.0e-10, 100);
    solver.solve(func, guess_spreads);

    // convert back from clean spreads to par spreads
    DoubleArraySP parSpreads(new DoubleArray(adjEndDates.size()));
    for (int j = 0 ; j < size ; j++)
    {
        (*parSpreads)[j] = func.calcParSpread(0.0, adjEndDates[j]);

        // do sanity check on resultant par spread
        if ((*parSpreads)[j] < 0.0)
        {
            string msg = Format::toString("Negative spread %f%% found for date %s",
                (*parSpreads)[j] * 100.0,
                adjEndDates[j].toString().c_str());
            throw ModelException(method, msg);
        }
    }

    // create result par spread curve
    DoubleArrayConstSP prices = getUpfrontFees();
    return CDSParSpreadsSP(new CDSParSpreads(name,
                      DateTime(getValueDate()),
                      getSpotOffset(),
                      getRecovery(),
                      getSwapFrequency(),
                      expiries,
                      parSpreads,
                      DoubleArraySP(prices.get() ? new DoubleArray(*prices) : NULL),
                      isFeeAccrued(),
                      dayCountConv()->toString(),
                      getBadDayConvention()->toString(),
                      getHolidays(),
                      getYieldCurve(),
                      getPrepayCurveObj()));
}


/** Get today and the par spreads curve from the market data cache */
void CDSParSpreads::getMarket(const IModel* model, const MarketData *market){
    static const string method = "CDSParSpreads::getMarket";
        
    // Call getMarket in the parent class
    CDSParSpreadsBase::getMarket(model, market);
    
    // Specific checks for this class
    if (hasDefaulted) {
        if (defaultDate.empty()) {
            // if the name has defaulted, defaultDate is mandatory
            throw ModelException(method,
                                 "No default date specified for defaulted name " +
                                 name);
        } else {
            // check if default date is in the past
            if (defaultDate > valueDate) {
                throw ModelException(method,
                                     "Default date in the future"
                                     " (Default date is " +
                                     defaultDate.toString() +
                                     ", Value date is " +
                                     valueDate.toString() + ")" );
            }    
        }
    }
}


/** Returns name identifying CDS Par Curve for CREDIT_DEFAULT_SENS */
string CDSParSpreads::sensName(CreditDefaultSensBase* sens) const {
    return name;
}

/** Shifts the object using given shift */
bool CDSParSpreads::sensShift(CreditDefaultSensBase* sens) {
    static const string method = "CDSParSpreads::sensShift(CreditDefaultSensBase* sens)";
    if (!defaulted()) {
        // We override the default date and recovery only if the name has not already defaulted
        // otherwise, we keep all the historical values unchanged
        defaultDate = sens->getCreditEventDate(valueDate);
        hasDefaulted = true;
        parRecovery = sens->getOverriddenRecovery(parRecovery);
        // Reset the "local" cleanSpreads cache
        cleanSpreads.reset();
    }
    return false;
}

/** Returns true if the name has defaulted */
bool CDSParSpreads::defaulted() const{
    return hasDefaulted;
}
    
/** Returns the date that this name defaulted. If a default has not
    happened an 'empty' DateTime is returned */
const DateTime& CDSParSpreads::getDefaultDate() const{
    return defaultDate;
}
                                      

void CDSParSpreads::acceptWrapperNameCollector(const CDSParSpreads* parSpreads, 
                                               WrapperNameCollector* collector)
{
    collector->addName(parSpreads->getName());
}


/** Invoked when this class is 'loaded' */
class CDSParSpreadsHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CDSParSpreads, clazz);
        SUPERCLASS(CDSParSpreadsBase);
        IMPLEMENTS(CreditDefaultSensBase::IShift);
        IMPLEMENTS(IDiscountCurveRisky);
        EMPTY_SHELL_METHOD(defaultCDSParSpreads);
        FIELD(hasDefaulted, "has a default happened");
        FIELD_MAKE_OPTIONAL(hasDefaulted); // for backwards compatibility
        FIELD(defaultDate, "If a default has happened, when");
        FIELD_MAKE_OPTIONAL(defaultDate);
 
        ClassSetAcceptMethod(CDSParSpreads::acceptWrapperNameCollector);
    }

    static IObject* defaultCDSParSpreads(){
        return new CDSParSpreads();
    }
};


CClassConstSP const CDSParSpreads::TYPE = CClass::registerClassLoadMethod(
    "CDSParSpreads", typeid(CDSParSpreads), CDSParSpreadsHelper::load);

DEFINE_TEMPLATE_TYPE(CDSParSpreadsWrapper);


class AdjustCdsRatesAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    CDSParSpreadsSP parSpreads;
    ExpiryArraySP   expiries;
    MarketDataSP    market;

    static IObject* defaultAdjustCdsRatesAddin()
    {
        return new AdjustCdsRatesAddin();
    }


    static IObjectSP adjustDates(AdjustCdsRatesAddin* params)
    {
        return params->adjust();
    }

    IObjectSP adjust()
    {
        static const string method = "AdjustCdsRatesddin::adjust";

        try
        {
            // create a non pricing model in order to get the market data
            NonPricingModel dummyModel;
            // must call getMarket() to ensure a prepay curve exists
            parSpreads->getMarket(&dummyModel, market.get());
            return IObjectSP(parSpreads->adjustDates("adjusted", expiries)->getParSpreads().clone());
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }


    /** for reflection */
    AdjustCdsRatesAddin():  CObject(TYPE), market(new MarketData()) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        REGISTER(AdjustCdsRatesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAdjustCdsRatesAddin);
        FIELD       (parSpreads, "par spread curve");
        FIELD       (expiries,   "adjusted maturity dates");
        FIELD       (market,     "Market data cache");
        FIELD_MAKE_OPTIONAL(market);

        Addin::registerClassObjectMethod(
            "ADJUST_PAR_SPREADS",
            Addin::MARKET,
            "adjust market CDS rates for an alternative set of dates",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*) adjustDates);
    }
};

CClassConstSP const AdjustCdsRatesAddin::TYPE = CClass::registerClassLoadMethod(
    "AdjustCdsRatesAddin", typeid(AdjustCdsRatesAddin), load);
   

DRLIB_END_NAMESPACE
