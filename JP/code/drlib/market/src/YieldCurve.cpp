//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YieldCurve.cpp
//
//   Description : yield curve interface
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_YIELDCURVE_CPP
#include "edginc/DeterministicYieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/RegressionTestMode.hpp"
#include "edginc/FourPlusIZeroCurve.hpp"
#include "edginc/UntweakableYC.hpp"
#include "edginc/MDFUtil.hpp"
#include <set>


DRLIB_BEGIN_NAMESPACE
////// IYieldCurve /////////

double IYieldCurve::parSwapRate(
                   const DateTime&           startDate, 
                   const DateTime&           endDate,
                   const MaturityPeriod&     period, 
                   const DayCountConvention& dcc,
                   const Stub&               stubType,
                   bool                      stubAtEnd,
                   const BadDayConvention&   accBadDayConv,
                   const BadDayConvention&   payBadDayConv,
                   const Holiday&            holidays) const 
{
    return SwapTool::parSwapRate(*this, startDate, endDate, period, dcc, stubType, stubAtEnd, accBadDayConv, payBadDayConv, holidays);
}


double IYieldCurve::parSwapRate(const DateTime& maturity) const
{
    throw ModelException("IYieldCurve::parSwapRate", "parSwapRate() not implemented for this curve implementation type");
}

// ALIB: tcurve.c#674 (GtoCurveDatesAndRates)
CashFlowArraySP IYieldCurve::getRatesAndDates(
    const ExpiryArray&         startPeriods, 
    const MaturityPeriodArray& addIntervals, 
    const DateTimeArray*       requiredDates) const
{
    static const string method = "IYieldCurve::getRatesAndDates";

    try
    {
        /*
         * Step one - get the required date list
         */
        set<DateTime> reqdDates;
        if (requiredDates != NULL)
        {
            reqdDates.insert(requiredDates->begin(), requiredDates->end());
        }
        else
        {
            DateTimeArray tmpDates = zeroDates();
            reqdDates.insert(tmpDates.begin(), tmpDates.end());
        }

        if (reqdDates.size() == 0)
            throw ModelException(method, "No dates requested");

        /*
         * Step two - get the curve date list by also adding the additional points.
         * 
         * We need to get the additional date interval for a given date, and then 
         * add this interval to the date and see if it falls before the next date 
         * in the date list.  If it does, then we must add a new date.
         */
        if (startPeriods.size() != addIntervals.size())
        {
            string msg = Format::toString(
                "start periods size (%d) and add intervals size (%d) must be the same",
                startPeriods.size(),
                addIntervals.size());
            throw ModelException(method, msg);
        }

        ExpiryArray startIntervals(startPeriods);
        MaturityPeriodArray gapIntervals(addIntervals);
        int numIntervals = startPeriods.size();

        set<DateTime> dates;
        DateTime baseDate = getSpotDate();

        DateTime lastDate = baseDate;
        if (numIntervals > 0 && lastDate < gapIntervals[0]->toDate(baseDate))
        {
            numIntervals++;
            startIntervals.insert(startIntervals.begin(), ExpirySP(new MaturityPeriod("0D")));
            gapIntervals.insert(gapIntervals.begin(), MaturityPeriodSP(new MaturityPeriod("0D")));
        }

        int intervalIdx = 0;
        for (set<DateTime>::iterator i = reqdDates.begin() ; i != reqdDates.end() ; i++)
        {
            // always add required date
            dates.insert(*i);

            if (numIntervals > intervalIdx)
            {
                // 0D can be provided as an addInterval which will not add any points
                if (!gapIntervals[intervalIdx]->isZeroLength())
                {
                    // add points in interval [last date, this date] determined by
                    // intervals array.
                    DateTime addDate = gapIntervals[intervalIdx]->toDate(lastDate);

                    while (addDate < *i)
                    {
                        // add addDate and increment it
                        dates.insert(addDate);
                        
                        // see whether to increment intervalIdx
                        while (intervalIdx < numIntervals - 1
                            && addDate >= startIntervals[intervalIdx+1]->toDate(baseDate))
                            intervalIdx++;

                        addDate = gapIntervals[intervalIdx]->toDate(addDate);
                    }
                }

                // see whether to increment intervalIdx
                while (intervalIdx < numIntervals - 1
                    && *i >= startIntervals[intervalIdx+1]->toDate(baseDate))
                    intervalIdx++;
            }

            lastDate = *i;
        }

        /*
         * Step three - now we have all the dates we need to return.
         */
        CashFlowArraySP datesAndRates(new CashFlowArray(dates.size()));
        int j = 0;
        for (set<DateTime>::iterator k = dates.begin() ; k != dates.end() ; k++)
        {
            DateTime date = *k;

            if (date == baseDate)
            {
                // not as arbitrary as seems!
                date = date.rollDate(1);
            }

            (*datesAndRates)[j].date = *k;
            (*datesAndRates)[j].amount = zero(date);
            j++;
        }

        return datesAndRates;
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}


IYieldCurve::~IYieldCurve(){}
void IYieldCurve::load(CClassSP& clazz){
    REGISTER_INTERFACE(IYieldCurve, clazz);
    EXTENDS(IMarketFactor);
    EXTENDS(IDiscountCurve);
}

CClassConstSP const IYieldCurve::TYPE = CClass::registerInterfaceLoadMethod(
    "IYieldCurve", typeid(IYieldCurve), load);

////// IDeterministicYieldCurve /////////
IDeterministicYieldCurve::~IDeterministicYieldCurve(){}
void IDeterministicYieldCurve::load(CClassSP& clazz){
    REGISTER_INTERFACE(IDeterministicYieldCurve, clazz);
    EXTENDS(IYieldCurve);
}

CClassConstSP const IDeterministicYieldCurve::TYPE =
CClass::registerInterfaceLoadMethod(
    "IDeterministicYieldCurve", typeid(IDeterministicYieldCurve), load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(IDeterministicYieldCurveWrapper);

const string YieldCurve::MMRT_RATE   = "M";
const string YieldCurve::MMRT_RATE2  = "MM";
const string YieldCurve::SWAP_RATE   = "S";
const string YieldCurve::SWAP_RATE2  = "SW";
const string YieldCurve::FUTURE_RATE = "F";
const string YieldCurve::TURN_RATE   = "A";
const string YieldCurve::BRAZIL_FUTURE_PRICE = "BF";

YieldCurve::~YieldCurve() {
    // empty
}


/** Calculates a par swap rate for the swap starting at startDate,
 * maturing at maturityDate, with fixed day count convention 
 * fixedDayCountConv and fixed payments occuring at time 
 * intervals defined by interval. In other words, the routine calculates 
 * the fixed rate such that the present value of the fixed and 
 * floating sides are equal. The floating side is assumed to be at par.
 */
// DEFAULT IMPLEMENTATION
double YieldCurve::couponRate(
    const DateTime&           startDate,       // (I) Date instrument begins at
    const DateTime&           maturityDate,    // (I) Date instrument matures at
    const MaturityPeriod&     interval,        // (I) Time between payments 
    bool                      stubAtEnd,
    const DayCountConvention* dcc) const {
    
    static const string method = "YieldCurve::couponRate";
    try {
        int    count;
        string period;

        interval.decompose(count, period);

        CashFlowArray cfl = SwapTool::cashflows(startDate,
                                                maturityDate,
                                                stubAtEnd,
                                                1.0,
                                                count,
                                                period,
                                                dcc);
        
        if (cfl.empty()) {
            throw ModelException(method, "no cashflows between " + 
                                 startDate.toString() + " and " +
                                 maturityDate.toString());
        }

        // don't want final principal repayment
        cfl.back().amount -= 1.0;

        const DateTime& spotDate = getSpotDate();
        auto_ptr<IYieldCurve::IKey> key(logOfDiscFactorKey());

        // Get present value of 1 at startDate
        double startDatePV = exp(key->calc(spotDate, startDate));

        // Compute coupon for the given zero-coupon rates
        double couponsPV = 0.0;
        for (int i = 0; i < cfl.size()-1; i++) {
            couponsPV += cfl[i].amount * exp(key->calc(spotDate, cfl[i].date));
        }
        double lastPV = exp(key->calc(spotDate, cfl.back().date));
        couponsPV += cfl.back().amount * lastPV;

        if (!Maths::isPositive(couponsPV)) {
            throw ModelException(method, "coupons with rate=1.0 value to <= 0.0");
        }

        return (startDatePV - lastPV)/couponsPV;  // This is the coupon rate
    }
    catch (exception&e ) {
        throw ModelException(e, method);
    }
}


/** Optimized equalTo for performance, used for caching only */
// DEFAULT IMPLEMENTATION - call Object's equalTo
bool YieldCurve::zeroCurveEquals(const IYieldCurve* yc2) const {
    return equalTo(yc2);
}


/** Optimized hashCode for performance, used for caching only */
// DEFAULT IMPLEMENTATION - call Object's hashCode
int YieldCurve::zeroCurveHash() const{
    return hashCode();
}


/** Interpolate a futures rate between two dates 
 * see CompoundBasis for basis values
 */
// DEFAULT IMPLEMENTATION
double YieldCurve::future(const DateTime&           lodate, 
                          const DateTime&           hidate,
                          const DayCountConvention* dcc,
                          int                       basis,
                          double                    irVol) const {
    
    static const string method = "YieldCurve::future";
    try {
        double   fwdRate = fwd(lodate, hidate, dcc, basis);
        DateTime valueDate = getSpotDate();

        double futRate = RateConversion::forwardToFuture(valueDate,
                                                         fwdRate,
                                                         irVol,
                                                         lodate,
                                                         hidate,
                                                         dcc);       
        return futRate;
    }
    catch (exception&e ) {
        throw ModelException(e, method);
    }
}


/** Calculates present value to baseDate of supplied cash flows. Cash
    flows on or before baseDate (eg value date) are ignored. No
    ordering of the cashflows is assumed */
double YieldCurve::pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const
{
    return pv(cashFlows, baseDate, true);
}

/** Calculates present value to baseDate of supplied cash flows. Cash
    flows on or before baseDate (eg value date) are ignored depending
    on the optional parameter ignoreCashFlowsOnBaseDate. No
    ordering of the cashflows is assumed */
double YieldCurve::pv(const CashFlowArray& cashFlows,
                      const DateTime& baseDate,
                      const bool ignoreCashFlowsOnBaseDate) const
{
    double sum = 0.0;
    for (int i = 0; i < cashFlows.size(); i++){
        if (cashFlows[i].date.isGreater(baseDate) ||
            (!ignoreCashFlowsOnBaseDate && (cashFlows[i].date == baseDate)))
        {
            sum += pv(baseDate, cashFlows[i].date) * cashFlows[i].amount;
        }
    }
    return sum;
}


/** Just does fwd(settles(start), rateMat(settles(start), end, bdc),
    dcc, basis) */
double YieldCurve::fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       int                       basis) const{
    DateTime start(settles(refixDate));
    DateTime end(rateMaturity(start, rateMat, bdc));
    return fwd(start, end, dcc, basis);
}

// apparently this is a sensible default behaviour regardless of the value
// of useAssetRecovery
double YieldCurve::riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional,
                           bool            useAssetRecovery,
                           double          assetRecovery) const {
    return riskyPV(lodate, hidate, cashFlow, recoveryNotional);
}


class DefaultLogOfDiscFactorKey: public YieldCurve::IKey{
public:
    DefaultLogOfDiscFactorKey(const YieldCurve* yc):yc(yc){}
    /** Returns the log of the discount factor between the two dates */
    virtual double calc(const DateTime&  loDate,
                        const DateTime&  hiDate){
        return (yc->fwd(loDate, hiDate, &dcc, CompoundBasis::CONTINUOUS) *
                -Actual365F::yearFraction(loDate, hiDate));
    }
private:
    const YieldCurve* yc;
    Actual365F        dcc;
};

/** Returns a key used to optimise repeated calculations of
    discount factors - see above. The default implementation has no
    performance improvements. */
YieldCurve::IKey* YieldCurve::logOfDiscFactorKey() const{
    return new DefaultLogOfDiscFactorKey(this);
}
   
/** Invoked when this class is 'loaded' */
void YieldCurve::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(YieldCurve, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IYieldCurve);
}

YieldCurve::YieldCurve(CClassConstSP clazz): MarketObject(clazz){}

CClassConstSP const YieldCurve::TYPE = CClass::registerClassLoadMethod(
    "YieldCurve", typeid(YieldCurve), load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(YieldCurveWrapper);

// initialise type for array of YieldCurve
DEFINE_TEMPLATE_TYPE(YieldCurveArray);

DEFINE_TEMPLATE_TYPE(YieldCurveWrapperArray);

class PvAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the yield curve and dates to 
        get pv factor between
    */
    YieldCurveSP    yc;
    DateTimeArraySP dates;     
    bool            fromFirstDate;
	bool            useEstimatingCurve;

    /** the 'addin function' - builds array of correct type */
    static double pv(PvAddin* params){
        static const string routine = "PvAddin::pv";
        try{
			params->yc->setProjectionCurve(params->useEstimatingCurve);

            if (params->dates->size() == 0){
                throw ModelException(routine, "no dates passed in");
            }
            if (params->dates->size() == 1){
                return params->yc->pv((*params->dates)[0]);
            } 
            else if (params->dates->size() == 2){
                return params->yc->pv((*params->dates)[0],
                                      (*params->dates)[1]);
            } 
            else {
                throw ModelException(routine, "Too many dates supplied");
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    static IObjectSP pvCurve(PvAddin* params){
        static const string routine = "PvAddin::pvCurve";
        try {
            DateTime today;
            int i;

            params->yc->setProjectionCurve(params->useEstimatingCurve);

            int outSize = params->dates->size();

            if (params->fromFirstDate) {
                outSize--;
            }
            else {
                today = params->yc->getToday();
            }            

            DoubleArraySP output(new DoubleArray(outSize));

            // Calculate the output using the pv function
            for (i = 0; i < outSize; i++) {
                if (params->fromFirstDate) {
                    (*output)[i] = params->yc->pv((*params->dates)[0],(*params->dates)[i+1]);
                }
                else {
                    (*output)[i] = params->yc->pv((*params->dates)[i]);
                }
            }

            if (RegressionTestMode::isOn())
            {
                // Re-calculate the output using the logOfDiscFactorKey function.
                // This is, effectively, performing the same calculation as above
                // in a different way, and comparing the results. If the performance
                //  of this add-in is an issue, this should be reworked.
                double output2;
                auto_ptr<IYieldCurve::IKey> key(params->yc->logOfDiscFactorKey());
                for (i = 0; i < outSize; i++) {
                    if (params->fromFirstDate) {
                        output2 = exp(key->calc((*params->dates)[0], (*params->dates)[i+1]));
                    }
                    else {
                        output2 = exp(key->calc(today, (*params->dates)[i]));
                    }
                    
                    #define ABS_TOL    1e-14
                    #define PRECISION  "%.15f"
                    if ((output2 - (*output)[i] > ABS_TOL) || 
                        (output2 - (*output)[i] < -ABS_TOL)) 
                    {
                        throw ModelException(routine,
                                            "Calculating PV for date '" +
                                            (*params->dates)[i+(params->fromFirstDate ? 1 : 0)].toString() +
                                            "' using 'pv' and 'logOfDiscFactorKey' differ " +
                                            "by more than the maximum tolerance (" + 
                                            Format::toString(PRECISION, ABS_TOL) + "): " +
                                            Format::toString(PRECISION, (*output)[i]) + " vs " + 
                                            Format::toString(PRECISION, output2));
                    }
                }
            }

            return IObjectSP(output);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    PvAddin():  CObject(TYPE), fromFirstDate(false), useEstimatingCurve(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PvAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPvAddin);
        FIELD(yc, "Yield curve");
        FIELD(dates, "Date(s) for pv factor");
        FIELD(fromFirstDate, "used in PV_FACTORS only");
		FIELD(useEstimatingCurve, "use estimating curve");
        FIELD_MAKE_OPTIONAL(fromFirstDate);
		FIELD_MAKE_OPTIONAL(useEstimatingCurve);

        Addin::registerInstanceDoubleMethod("PV_FACTOR",
                                            Addin::RISK,
                                            "Returns a discount factor",
                                            TYPE,
                                            (Addin::DoubleMethod*)pv);

        Addin::registerClassObjectMethod("PV_FACTORS",
                                         Addin::RISK,
                                         "Returns a discount curve",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)pvCurve);
    }

    static IObject* defaultPvAddin(){
        return new PvAddin();
    }
    
};

CClassConstSP const PvAddin::TYPE = CClass::registerClassLoadMethod(
    "PvAddin", typeid(PvAddin), load);



class CouponRateAddin: public CObject{
    static CClassConstSP const TYPE;

    YieldCurveSP         yc;
    DateTime             startDate;     
    DateTime             maturityDate;   
    string               interval;
    bool                 stubAtEnd;
    DayCountConventionSP dcc;
    MarketDataConstSP    marketData;

    /** the 'addin function' */
    static double coupon(CouponRateAddin* params){
        static const string routine = "CouponRateAddin::coupon";
        try {
            MaturityPeriod       interval(params->interval);
            
	if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->dcc, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            }     
            return params->yc->couponRate(params->startDate,
                                          params->maturityDate,
                                          interval,
                                          params->stubAtEnd,
                                          params->dcc.get());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    CouponRateAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CouponRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCouponRateAddin);
        FIELD(yc, "Yield curve");
        FIELD(startDate, "startDate");
        FIELD(maturityDate, "maturityDate");
        FIELD(interval, "interval");
        FIELD(stubAtEnd, "stubAtEnd");
        FIELD(dcc, "dcc");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);
        Addin::registerInstanceDoubleMethod("COUPON_RATE",
                                            Addin::RISK,
                                            "Returns a coupon rate",
                                            TYPE,
                                            (Addin::DoubleMethod*)coupon);
    }

    static IObject* defaultCouponRateAddin(){
        return new CouponRateAddin();
    }
    
};

CClassConstSP const CouponRateAddin::TYPE = CClass::registerClassLoadMethod(
    "CouponRateAddin", typeid(CouponRateAddin), load);

class fwdRateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    YieldCurveSP         yc;
    DateTime             lodate;            
    DateTime             hidate;
    DayCountConventionSP dcc;
    int                  compoundBasis;
    bool                 useEstimatingCurve;
    MarketDataConstSP    marketData;

    // Return the rate
    static double forward(fwdRateAddin* params){
        static const string routine = "fwdRateAddin::forward";
        try{
            if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->dcc, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            }                        

	    params->yc->setProjectionCurve(params->useEstimatingCurve);
            return params->yc->fwd(params->lodate,
                                   params->hidate,
                                   params->dcc.get(),
                                   params->compoundBasis);
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }
    

    fwdRateAddin(): CObject(TYPE), useEstimatingCurve(false) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(fwdRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultfwdRateAddin);
        FIELD(yc, "yield curve");
        FIELD(lodate, "low date");
        FIELD(hidate, "high date");
        FIELD(dcc, "dcc");
        FIELD(compoundBasis, "compounding basis");
        FIELD(useEstimatingCurve, "use estimating curve");
        FIELD_MAKE_OPTIONAL(useEstimatingCurve);
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);
        Addin::registerInstanceDoubleMethod("FWD_RATE",
                                            Addin::RISK,
                                            "Returns a fwd rate",
                                            TYPE,
                                            (Addin::DoubleMethod*)forward);
    }
    
    static IObject* defaultfwdRateAddin(){
        return new fwdRateAddin();
    }
 
};

CClassConstSP const fwdRateAddin::TYPE = CClass::registerClassLoadMethod("fwdRateAddin", 
                                                                         typeid(fwdRateAddin), load);




class FutureRateAddin: public CObject{
    static CClassConstSP const TYPE;

    YieldCurveSP         yc;
    DateTime             lodate;     
    DateTime             hidate;   
    DayCountConventionSP dcc;
    int                  basis;
    double               irvol;
    MarketDataConstSP    marketData;

    /** the 'addin function' */
    static double future(FutureRateAddin* params){
        static const string routine = "FutureRateAddin::future";
        try {
	    if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->dcc, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            }  

            return params->yc->future(params->lodate,
                                      params->hidate,
				      params->dcc.get(),
                                      params->basis,
                                      params->irvol);

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    FutureRateAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FutureRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFutureRateAddin);
        FIELD(yc, "Yield curve");
        FIELD(lodate, "lodate");
        FIELD(hidate, "hidate");
        FIELD(dcc, "dcc");
        FIELD(basis, "basis");
        FIELD(irvol, "irvol");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);
        Addin::registerInstanceDoubleMethod("FUTURE_RATE",
                                            Addin::RISK,
                                            "Returns a future rate",
                                            TYPE,
                                            (Addin::DoubleMethod*)future);
    }

    static IObject* defaultFutureRateAddin(){
        return new FutureRateAddin();
    }
    
};

CClassConstSP const FutureRateAddin::TYPE = CClass::registerClassLoadMethod(
    "FutureRateAddin", typeid(FutureRateAddin), load);

// addin for PV factors with and without ccy basis
class PvFactor: public CObject{
    static CClassConstSP const TYPE;

    string            name;
    DateTimeArray     dates;     
    MarketDataSP      market;
    IModelSP          model;
    bool              useProjectionCurve;
    bool              fromFirstDate;

    /** the 'addin function' - builds array of correct type */
    static IObjectSP pv(PvFactor* params){
        static const string routine = "PvFactor::pv";
        try{
            if (params->dates.empty()){
                throw ModelException(routine, "no dates passed in");
            }
            
            YieldCurveWrapper yc(params->name);

             // clone the model
            IModelSP mdl(params->model.get());
            // and get the market data fetcher
            MarketDataFetcherSP mdf = mdl->getMDF();
            // and allow the use of currency basis
            MDFUtil::setUseCurrencyBasis(*mdf, true);

            yc.getData(mdl.get(), params->market.get());
            
            yc->setProjectionCurve(params->useProjectionCurve);

            int numDates = params->dates.size();
            if (params->fromFirstDate)
            {
                if (numDates < 2)
                {
                    throw ModelException(routine, "must supply at least 2 dates with fromFirstDate");
                }
                else
                {
                    numDates--;
                }
            }

            DoubleArraySP output(new DoubleArray(numDates));

            for (int i = 0; i < numDates; i++) {
                if (params->fromFirstDate)
                {
                    (*output)[i] = yc->pv(params->dates[0], params->dates[i+1]);
                }
                else
                {
                    (*output)[i] = yc->pv(params->dates[i]);
                }
            }

            return IObjectSP(output); 
        } catch (exception& e){
            throw ModelException(e, routine);
        }
   }

    /** for reflection */
    PvFactor():  CObject(TYPE), useProjectionCurve(false), fromFirstDate(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PvFactor, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPvFactor);
        FIELD(name, "Yield curve name");
        FIELD(dates, "Date(s) for pv factor");
        FIELD(market, "market");
        FIELD(model, "model");
        FIELD(useProjectionCurve, "use projection curve?");
        FIELD(fromFirstDate, "if >1 date, pv to first date, not value date?");
        FIELD_MAKE_OPTIONAL(fromFirstDate);
        Addin::registerClassObjectMethod("PV_FACTOR_CCY_BASIS",
                                         Addin::RISK,
                                         "Returns a discount factor",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)pv);

    }

    static IObject* defaultPvFactor(){
        return new PvFactor();
    }
    
};

CClassConstSP const PvFactor::TYPE = CClass::registerClassLoadMethod(
    "PvFactor", typeid(PvFactor), load);


// addin for discount factor from spot date to today
class SpotDatePV: public CObject, virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    string            name;
    MarketDataSP      market;

    // EdrAction version of addin
    IObjectSP run() {
        return CDoubleSP(CDouble::create(discFactor(this)));
    }

    /** the 'addin function' */
    static double discFactor(SpotDatePV* params){
        static const string method = "SpotDatePV::discFactor";
        try{

            YieldCurveWrapper yc(params->name);
            NonPricingModel dummyModel;
      
            // get the default market data fetcher
            MarketDataFetcherSP mdf = dummyModel.getMDF();
            // and allow the use of currency basis
            MDFUtil::setUseCurrencyBasis(*mdf, true);

            yc.getData(&dummyModel, params->market.get());
            
            double pvFactor =
                yc->pv(yc->getSpotDate());

            return pvFactor;

        } catch (exception& e){
            throw ModelException(e, method);
        }
   }

    /** for reflection */
    SpotDatePV():  CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotDatePV, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultSpotDatePV);
        FIELD(name, "Yield curve name");
        FIELD(market, "market");
        Addin::registerClassDoubleMethod("SPOTDATE_PV",
                                         Addin::RISK,
                                         "Returns a discount factor from the spot date to today",
                                         TYPE,
                                         (Addin::DoubleMethod*)discFactor);

    }

    static IObject* defaultSpotDatePV(){
        return new SpotDatePV();
    }
};

CClassConstSP const SpotDatePV::TYPE = CClass::registerClassLoadMethod(
    "SpotDatePV", typeid(SpotDatePV), load);


// construct yield curve from date & rates

class YieldCurveAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

private:
    DateTime             today;
    DateTime             baseDate;
    DateTimeArray        dates;
    DoubleArray          rates;
    DayCountConventionSP dcc;
    int                  basis;
    MarketDataConstSP    marketData;


    static IObject* defaultYieldCurveAddin()
    {
        return new YieldCurveAddin();
    }


    YieldCurveAddin() : CObject(TYPE)
    {
    }


    IObjectSP createYieldCurve()
    {
        static const string method = "YieldCurveAddin::createYieldCurve";

        try 
        {
            if (dates.size() != rates.size())
            {
                string msg = Format::toString(
                    "number of dates (%d) fiffer from number of rates (%d)", 
                    dates.size(), rates.size());
                throw ModelException(method, msg);
            }
	    
	    
            //if (dcc.empty())
	    if (!dcc.get())
            {
                throw ModelException(method, "no day count convention specified");
            }
	    
	    if (marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                marketData->fetchForNonMarketObjects(dcc, 
						     IModelConstSP(new NonPricingModel()),
						     "Business252");
            } 
            FourPlusIZeroCurveSP zc(new FourPlusIZeroCurve(baseDate, dcc, basis, dates, rates));
            IObject* object = new UntweakableYC(today, baseDate, *zc, zc.get());
            return IObjectSP(object);
        }
        catch (exception &e) 
        {
            throw ModelException(e, method);
        }
    }


    static void load(CClassSP& clazz)
    {
        REGISTER(YieldCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultYieldCurveAddin);

        FIELD(today,    "today");
        FIELD(baseDate, "date all rates start at");
        FIELD(dates,    "dates");
        FIELD(rates,    "rates");
        FIELD(dcc,      "day count convention, usually Act/365F");
        FIELD(basis,    "rate type");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerObjectMethod(
            "ZERO_CURVE_MAKE",
            Addin::MARKET,
            "Create yield curve from dates and rates",
            true,
            Addin::returnHandle,
            &YieldCurveAddin::createYieldCurve);
    }
};


CClassConstSP const YieldCurveAddin::TYPE = 
CClass::registerClassLoadMethod("YieldCurveAddin", typeid(YieldCurveAddin), load);


// Add-in for zero dates
class ZeroDatesAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

private:
    YieldCurveConstSP yc;

    ZeroDatesAddin(): CObject(TYPE) {}

    static IObject* defaultZeroDatesAddin()
    {
        return new ZeroDatesAddin();
    }

    IObjectSP zeroDates()
    {
        return IObjectSP(yc->zeroDates().clone());
    }

    static void load(CClassSP& clazz)
    {
        REGISTER(ZeroDatesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZeroDatesAddin);

        FIELD(yc, "Yield curve");

        Addin::registerObjectMethod("ZERO_DATES",
                                    Addin::RISK,
                                    "Returns the zero dates",
                                    false,
                                    Addin::expandSimple,
                                    &ZeroDatesAddin::zeroDates);
    }
};


CClassConstSP const ZeroDatesAddin::TYPE = CClass::registerClassLoadMethod(
    "ZeroDatesAddin", typeid(ZeroDatesAddin), load);


// Add-in for zero dates and rates
class ZeroDatesAndRatesAddin : public CObject
{
public:
    static CClassConstSP const TYPE;

private:
    YieldCurveConstSP yieldCurve;
    DateTimeArraySP   requiredDates;
    ExpiryArraySP     startPeriods;
    ExpiryArraySP     addIntervals; // must be maturity periods
    bool              useEstimatingCurve;

    ZeroDatesAndRatesAddin(): CObject(TYPE), useEstimatingCurve(false) {}

    static IObject* defaultZeroDatesAndRatesAddin()
    {
        return new ZeroDatesAndRatesAddin();
    }

    ObjectArraySP zeroDatesAndRates()
    {
        CashFlowArraySP cfl;
        yieldCurve->setProjectionCurve(useEstimatingCurve);

        if (requiredDates.get() == NULL && startPeriods.get() == NULL && addIntervals.get() == NULL)
        {
            cfl = yieldCurve->getRatesAndDates();
        }
        else
        {
            if (startPeriods.get() == NULL || addIntervals.get() == NULL)
            {
                throw ModelException("ZeroDatesAndRatesAddin::zeroDatesAndRates",
                    "startPeriods and addIntervals must be provided");
            }

            // QLib add-in takes ExpiryArray, but really it should be an array of tenors
            MaturityPeriodArray intervals;
            for (int i = 0 ; i < addIntervals->size() ; i++)
                intervals.push_back(MaturityPeriodSP(dynamic_cast<MaturityPeriod*>((*addIntervals)[i].get())));

            cfl = yieldCurve->getRatesAndDates(*startPeriods, intervals, requiredDates.get());
        }

        // bundle results in a form that expands to 2 Excel columns
        ObjectArraySP results(new ObjectArray(0));
        results->push_back(IObjectSP(CashFlow::dates(*cfl).clone()));
        results->push_back(CashFlow::amounts(*cfl));
        return results;
    }

    static void load(CClassSP& clazz)
    {
        REGISTER(ZeroDatesAndRatesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZeroDatesAndRatesAddin);

        FIELD       (yieldCurve,         "Yield curve");
        FIELD       (startPeriods,       "Start periods");
        FIELD       (addIntervals,       "Add intervals");
        FIELD       (requiredDates,      "Required dates");
		FIELD(useEstimatingCurve, "use estimating curve");

        FIELD_MAKE_OPTIONAL(requiredDates);
        FIELD_MAKE_OPTIONAL(startPeriods);
        FIELD_MAKE_OPTIONAL(addIntervals);
        FIELD_MAKE_OPTIONAL(useEstimatingCurve);

        Addin::registerObjectMethod("ZERO_CURVE_RATES",
                                    Addin::RISK,
                                    "Returns the zero dates",
                                    false,
                                    Addin::expandMulti,
                                    &ZeroDatesAndRatesAddin::zeroDatesAndRates);
    }
};


CClassConstSP const ZeroDatesAndRatesAddin::TYPE = CClass::registerClassLoadMethod(
    "ZeroDatesAndRatesAddin", typeid(ZeroDatesAndRatesAddin), load);


// Util: curve discount factor
class YCDiscountFactorAddin: public CObject
{
public:
    static CClassConstSP const TYPE;
    YieldCurveSP        yieldCrv;
    DateTime            date;
    YCDiscountFactorAddin(): CObject(TYPE) {}
    static IObject* defaultYCDiscountFactorAddin()
    {
        return new YCDiscountFactorAddin();
    }
    double discountFactor()
    {
        return yieldCrv->pv(date);
    }
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(YCDiscountFactorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultYCDiscountFactorAddin);
        FIELD(yieldCrv, "yield curve");
        FIELD(date, "date");
        Addin::registerDoubleMethod("YC_DISCOUNT_FACTOR",
            Addin::UTILITIES,
            "Returns discount factor from zero curve",
            &YCDiscountFactorAddin::discountFactor);                                    
    }
};
CClassConstSP const YCDiscountFactorAddin::TYPE=CClass::registerClassLoadMethod(
    "YCDiscountFactorAddin", typeid(YCDiscountFactorAddin), load);


// Util: curve forward rate
class YCForwardRateAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    YieldCurveSP         yieldCrv;
    DateTime             dateStart;
    DateTime             dateEnd;
    DayCountConventionSP dcc;
    string               basis;
    MarketDataConstSP    marketData;


    YCForwardRateAddin(): CObject(TYPE) {}

    static IObject* defaultYCForwardRateAddin()
    {
        return new YCForwardRateAddin();
    }

    double forwardRate()
    {
        if (marketData.get())
        {
                // May need to fetch holidays for Business/252 DCC
                marketData->fetchForNonMarketObjects(dcc, 
						     IModelConstSP(new NonPricingModel()),
						     "Business252");
        }    
        return yieldCrv->fwd(dateStart, dateEnd, dcc.get(), CompoundBasis::toInt(basis));
    }

    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(YCForwardRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultYCForwardRateAddin);
        FIELD(yieldCrv, "yield curve");
        FIELD(dateStart, "start date");
        FIELD(dateEnd, "end date");
        FIELD(dcc, "day count convention");
        FIELD(basis, "compound basis");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerDoubleMethod("YC_FORWARD_RATE",
            Addin::UTILITIES,
            "Returns the forward rate from curve",
            &YCForwardRateAddin::forwardRate);                                    
    }

};

CClassConstSP const YCForwardRateAddin::TYPE=CClass::registerClassLoadMethod(
    "YCForwardRateAddin", typeid(YCForwardRateAddin), load);


// Util: curve par swap rate
class YCParSwapRateAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    YieldCurveSP         yieldCrv;
    DateTime             startDate;
    DateTime             endDate;
    string               period;
    DayCountConventionSP dcc;
    string               stubRule;
    bool                 stubAtEnd;
    string               accBadDayConv;
    string               payBadDayConv;
    HolidaySP            holidays;
    MarketDataConstSP    marketData;

    YCParSwapRateAddin(): CObject(TYPE) {}

    static IObject* defaultYCParSwapRateAddin()
    {
        return new YCParSwapRateAddin();
    }

    double parSwapRate()
    {
        if (marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                marketData->fetchForNonMarketObjects(dcc, 
						     IModelConstSP(new NonPricingModel()),
						     "Business252");
            }  
        MaturityPeriod _period(period);
        StubSP _stub(StubFactory::make(stubRule));
        BadDayConventionSP _accBadDayConv_ptr (BadDayConventionFactory::make(accBadDayConv));
        BadDayConventionSP _payBadDayConv_ptr (BadDayConventionFactory::make(payBadDayConv));

	double rtn = yieldCrv->parSwapRate(startDate, endDate, _period, *dcc.get(), *_stub.get(), 
                                stubAtEnd, *_accBadDayConv_ptr.get(), *_payBadDayConv_ptr.get(), *holidays.get());
        return rtn;
    }

    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(YCParSwapRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultYCParSwapRateAddin);
        FIELD(yieldCrv, "yield curve");
        FIELD(startDate, "start date");
        FIELD(endDate,   "end date");
        FIELD(period, "payment period");
        FIELD(dcc, "day count convention");
        FIELD(stubRule, "stub rule");
        FIELD(stubAtEnd, "stub at end");
        FIELD(accBadDayConv, "accrual bad day convention");
        FIELD(payBadDayConv, "pay bad day convention");
        FIELD(holidays, "holidays");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerDoubleMethod("YC_PAR_SWAP_RATE",
            Addin::UTILITIES,
            "Returns the forward rate from curve",
            &YCParSwapRateAddin::parSwapRate);                                    
    }

};

CClassConstSP const YCParSwapRateAddin::TYPE=CClass::registerClassLoadMethod(
    "YCParSwapRateAddin", typeid(YCParSwapRateAddin), load);


DRLIB_END_NAMESPACE
