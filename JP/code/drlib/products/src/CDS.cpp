//----------------------------------------------------------------------------
//
//   Filename    : CDS.cpp
//
//   Description : General Credit Default Swap
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDS.hpp"
#include "edginc/ICreditEventOverride.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/CreditIndex.hpp"
#include "edginc/Settlement.hpp"

DRLIB_BEGIN_NAMESPACE

class CDS::CDSOutputs {
public:
    CDSOutputs() : ctgPv(0.0), feePv(0.0), accInt(0.0)
    {}

    // Returns the price as ctgPv - feePv -  accInt
    // Convention is long protection, matching that of CDO
    const double price() const {
        return ctgPv - feePv - accInt;
    }

    // fields
    double ctgPv;
    double feePv;
    double accInt;
};


/******************************************************************************
 *                           CDS::PriceCache
 ******************************************************************************/

/** PriceCache is used to price and cache the price of the underlying 
 * CDSs in a CIS. Note this class is tightly linked to the model since
 * it controls how the pricing gets computed. Actually, it should be 
 * moved into the model...*/
class CDS::PriceCache : public CObject {
public:
    static CClassConstSP const TYPE;

    /** Constructor - Note that the AdjustedCDSParSpreadsSP will be
     * deep-copied, ie, the underlying curve to adjust and the actual
     * (index basis) adjustment will be copied here. */ 
    PriceCache(AdjustedCDSParSpreadsSP     adjCDS,
               double                      namesWeight,
               ICreditEventOverrideNameSP  nameOverride);

    /** Notification that some internal fields have changed */
    void fieldsUpdated(const CFieldArray& fields);

    /** Obtains the price of the curve, computing if required */
    double getPrice(const CDS* const cis,
                    IForwardRatePricerSP model);  

    /** Return the known cash flows for this name's CDS, taking potential
     * defaults into account */
    CashFlowArraySP knownCashFlows(const CDS* const cis,
                                   IForwardRatePricerSP model) const;

private:
    /** Prices a defaulted curve */
    double priceDefaultedCurve(const CDS* const cis,
                               IForwardRatePricerSP model);

    /** For reflection - this class needs to be registered to received
     * notifications via fieldsUpdated */
    PriceCache();
    static void load(CClassSP& clazz);
    static IObject* defaultPriceCache();

    // Fields
    bool   isPriceValid; // Whether the price cached is up-to-date or not
    double price;        // Cached price of the adjCDS curve 
    AdjustedCDSParSpreadsSP     adjCDS;       // The curve whose price we cache
    double                      adjCDSWeight; // Weight of the curve in the index
    ICreditEventOverrideNameSP  nameOverride; // Default parameters override
};


/* Constructor - Note that the AdjustedCDSParSpreads will be
 * deep-copied, ie, the underlying curve to adjust and the actual
 * (index basis) adjustment will be copied here. */
CDS::PriceCache::PriceCache(
        AdjustedCDSParSpreadsSP     copyCDS,
        double                      namesWeight,
        ICreditEventOverrideNameSP  nameOverride) :
    CObject(TYPE),
    isPriceValid(false),
    adjCDS(copyCDS.clone()), // Clone the SP
    adjCDSWeight(namesWeight),
    nameOverride(nameOverride)
{}


CDS::PriceCache::PriceCache() : 
    CObject(TYPE), isPriceValid(false), adjCDSWeight(0.0)
{}


IObject* CDS::PriceCache::defaultPriceCache() {
    return new PriceCache();
}


/* Notification that some internal fields have changed -> Reset the 
 * isPriceValid flag since the price is no longer valid */
void CDS::PriceCache::fieldsUpdated(const CFieldArray& fields) {
    isPriceValid = false;
}


/* Obtains the price of the underlying curve. Computes the price if 
 * required (using instrument data) and caches the result */
double CDS::PriceCache::getPrice(const CDS* const cis,
                                 IForwardRatePricerSP model) {
    static const string method = "CDS::PriceCache::getPrice";

    if (!isPriceValid) {
        if (adjCDSWeight == 0.0) {
            price = 0.0;  // Zero weight implies zero price
        }
        else {
            // Price assuming that all the CIS's notional is used to price this 
            // name's CDS. Will be scaled at the end before assigning it to "price"
            double totalPrice = 0.0;
        
            if (adjCDS->defaulted()) {
                totalPrice = priceDefaultedCurve(cis, model);
            }
            else { // The CDS has not defaulted
                DateTime valueDate = cis->getValueDate();
                IDecretionCurveConstSP prepay = adjCDS->getPrepayCurve();
                IBadDayAdjusterSP idxAdjuster = (cis->getParSpreadsWrapper()).getSP();
                CDS::CDSOutputs dummyCDSOutputs;
                totalPrice = cis->getPV(valueDate, valueDate,
                                        *(adjCDS.get()), prepay, model, 
                                        idxAdjuster, dummyCDSOutputs);
            }

            // Up to here, we have assumed the whole CIS's notional was for this
            // name's CDS. Scale the result using this name's weight.
            price = totalPrice * adjCDSWeight;
        }

        // We have just computed the price so it is now valid
        isPriceValid = true;
    }
    return price;
}


/** Price a defaulted curve */
double CDS::PriceCache::priceDefaultedCurve(const CDS* const cis,
                                            IForwardRatePricerSP model)
{
//     static const string method("CDS::PriceCache::priceDefaultedCurve");

    const DateTime& curveDefaultDate = adjCDS->getDefaultDate();

    ICreditFeeLegSP fLeg            = cis->getFeeLeg();
    ICreditContingentLegSP cLeg     = cis->getContingentLeg();
    const DateTime& valueDate       = cis->getValueDate();
    IBadDayAdjusterSP idxAdjuster   = cis->getParSpreadsWrapper().getSP();
    IDiscountCurveConstSP  discount = cis->getYieldCurveWrapper().getSP();
    IDecretionCurveConstSP prepay   = adjCDS->getPrepayCurve();

    double ctgLegPv = 0.0;
    if (!!cLeg) {
        ctgLegPv = cLeg->getContingentLegDefaultedPV(
            valueDate, 
            *(adjCDS.get()),
            discount,
            prepay,
            curveDefaultDate,
            idxAdjuster,
            false, // This is CIS so do not allow including payments on valueDate
            nameOverride,
            cis->triggerDelay,
            cis->defaultToSettlementDelay,
            cis->lastTriggerDate);
    }

    double feeLegPv = 0.0;
    if (!!fLeg) {
        feeLegPv = fLeg->getFeeLegDefaultedPV(
            valueDate, 
            curveDefaultDate,
            !cLeg ? DateTime() : cLeg->firstObservationStartDate(),
            !cLeg ? DateTime() : cLeg->lastObservationEndDate(),
            false, // This is CIS so do not allow including payments on valueDate
            discount,
            prepay,
            model,
            idxAdjuster,
            nameOverride,
            cis->triggerDelay,
            cis->defaultToSettlementDelay,
            cis->lastTriggerDate);
    }

    // Fee Leg has the right sign here (negative if paying fees) so
    // just add the legs' pvs
    const double value = ctgLegPv + feeLegPv;
    return value;
}

/** Returns all known cash flows */
CashFlowArraySP CDS::PriceCache::knownCashFlows(
    const CDS* const cis,
    IForwardRatePricerSP model) const 
{
    CashFlowArraySP knownCashFlows;

    // Risky payments depend on whether a default happened or not
    if (adjCDS->defaulted()) {
        IBadDayAdjusterSP idxAdjuster = (cis->getParSpreadsWrapper()).getSP();
        // Get the fee leg cashflows - pass an emtpy DateTime to indicate
        // we are interested in all cashflows, ie, no cashflows should be
        // excluded
        knownCashFlows = cis->feeLeg->generateKnownCashFlowsGivenDefault(
            cis->getValueDate(),
            adjCDS->getDefaultDate(),
            DateTime(), // exclude no cashflows
            !cis->contingentLeg ? DateTime() : 
                                  cis->contingentLeg->firstObservationStartDate(),
            !cis->contingentLeg ? DateTime() : 
                                  cis->contingentLeg->lastObservationEndDate(),
            false, // This is CIS so do not allow including payments on valueDate
            model,
            nameOverride,
            idxAdjuster,
            cis->triggerDelay,
            cis->defaultToSettlementDelay,
            cis->lastTriggerDate);

        //getFeeLegCashFlowsGivenDefault(cis, DateTime(), model);
    }
    else {
        // Get the fee leg to produce the cash flows for this non-defaulted name
        knownCashFlows = cis->feeLeg->generateKnownCashFlows(model);
    }

    for (int i=0; i < knownCashFlows->size(); ++i) {
        (*knownCashFlows)[i].amount *= adjCDSWeight;
    }

    // Aggregate payments happening on the same date
    CashFlow::aggregate(*knownCashFlows);

    return knownCashFlows;
}


/** Invoked when Class is 'loaded' */
void CDS::PriceCache::load(CClassSP& clazz){
    // clazz->setPublic(); - Do NOT make visible to EAS/spreadsheet
    REGISTER(CDS::PriceCache, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultPriceCache);

    FIELD(isPriceValid, "Whether the cached price is or not valid");
    FIELD(price,        "The cached price");
    FIELD(adjCDSWeight, "The weight of this curve in the index");
    FIELD(adjCDS,       "The curve whose price is cached here");
    FIELD(nameOverride, "Override for this name's default parameters");
}

CClassConstSP const CDS::PriceCache::TYPE = 
    CClass::registerClassLoadMethod("CDS::PriceCache", 
                                    typeid(CDS::PriceCache), 
                                    load);

// Array has to have its own type
// Work around bug in msvc7
typedef CDS::PriceCacheArray CDSPriceCacheArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("CDS::PriceCacheArray", CDSPriceCacheArray);

//--------------------------------------------------------------------------------

//private class
class CDSClosedFormCDSPS: public ClosedFormCDSPS::IProduct
{
public:
    CDSClosedFormCDSPS(const CDS* imnt) : imnt(imnt)
    {}

    void price(ClosedFormCDSPS* model,
               Control*         control, 
               CResults*        results) const
    {
        //get the forward rate pricer from the model
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp) {
            throw ModelException("CredDefSwapClosedFormCDSPS::price",
                                 "Model must implement IHasForwardRatePricer");
        }
        IForwardRatePricerSP ifrp = ihfrp->getForwardRatePricer();

        //delegate back to the instrument
        imnt->priceSingleName(results, control, ifrp);
    }

private:
    const CDS* imnt; // a reference
};

//--------------------------------------------------------------------------------

/** Private class - Handles the interaction between model and instrument */
class CDSClosedFormCDSBasket: public ClosedFormCDSBasket::IProduct
{
public:
    CDSClosedFormCDSBasket(const CDS* imnt): imnt(imnt)
    {}

    void price(ClosedFormCDSBasket* model,
               Control*         control, 
               CResults*        results) const
    {
        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp) {
            throw ModelException("CreditIndexSwapClosedFormCDSBasket::price",
                                 "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

        imnt->priceMultiName(results, control, frModel);
    }

private:
    const CDS* imnt;  // a reference
};

//--------------------------------------------------------------------------------

//------------
// CDS methods
//------------

CDS::~CDS(){}
// Store the name of the adjCurvesCache field
const string CDS::adjCurvesCacheFieldName = "adjCurvesCache";

/** Main pricing method */
void CDS::priceSingleName(CResults* results, 
                          Control* control, 
                          IForwardRatePricerSP model) const
{
    /* see if there is a default */
    if(underlying->defaulted()) {
        const DateTime& defaultDate = underlying->getDefaultDate();
        double defaultedPrice = getPVGivenDefault(valueDate, 
                                                  defaultDate, 
                                                  *(underlying.getSP()), 
                                                  model);
        results->storePrice(defaultedPrice, discount->getCcy());        
    }
    else {
        /* Adjust the settlement date for holidays */
        DateTime adjSettlementDate = valueDate; //badDayAdjust(valueDate);

        /* Calculate the pv as of valueDate conditional on on default to adjSettlementDate */
        CDSOutputs myOutput;

        IDecretionCurveConstSP prepay = underlying->getPrepayCurve();
        IBadDayAdjusterConstSP bda(IBadDayAdjusterConstSP::attachToRef(this));
        double price = getPV(valueDate, adjSettlementDate,
                             *(underlying.getSP()), prepay,
                             model, bda, myOutput);

        results->storePrice(price, discount->getCcy());

        //--------------------------
        // Now handle OutputRequests
        //--------------------------

        if (control && control->isPricing()) {
            addOutputRequests(control, results, model, false, myOutput);
        }
    }
}

void CDS::priceMultiName(CResults* results, 
                         Control* control, 
                         IForwardRatePricerSP model) const
{
    static const string method = "CDS::priceMultiName";

    //ensure that the underlying is a multi-name form
    //...currently this means just a CreditIndex
    CreditIndexConstSP creditIndex;
    try {
        creditIndex = CreditIndexConstSP::dynamicCast(underlying.getSP());
    }
    catch (exception&) {
        throw ModelException(method,
            "Pricing via ClosedFormCDSBasket requires that the underlying "
            "is a CreditIndex.");
    }

    //----------------------------------
    // Firsly set up the adjusted names
    // and the associated price caches
    // - note that we need to do this only once
    // - and therefore we just test for >0 elements in the
    // - already present in the cache
    //----------------------------------
    if (adjCurvesCache->size() == 0)
    {
        // Need to populate adjCurvesCache with the portfolio of 
        // names in the index, basis-adjusted.
        AdjustedCDSParSpreadsArrayConstSP adjCDSArray = 
            creditIndex->getAllBasisAdjustedCDSParSpreads();

        // Allocate the adjCurvesCache array, and store each adjusted curve
        // from adjCDSArray in it
        int numNamesInPortfolio = adjCDSArray->size();
        adjCurvesCache = PriceCacheArraySP(new PriceCacheArray(numNamesInPortfolio));

        // If there is no override in the instrument, try getting the  override 
        // from the creditIndex object (which may also be "NULL")
        if (!creditEventOverride) {
            creditEventOverride = creditIndex->getCreditEventOverride();
        }

        // Note the adjCurvesCache class is tightly linked to the model, since
        // it controls how the pricing gets computed. Actually, it should be 
        // moved into the model...
        ICreditEventOverrideNameSP defOverride;

        // Get the names and weights of all names in the portfolio
        CStringArray adjCDSName(numNamesInPortfolio);
        CDoubleArray adjCDSWeight(numNamesInPortfolio);
        double cumulativeWeight = 0.0;
        for (int i=0; i < numNamesInPortfolio; ++i) {
            adjCDSName[i] = (*adjCDSArray)[i]->getName();
            adjCDSWeight[i] = creditIndex->getNamesWeight(adjCDSName[i]);
            cumulativeWeight += adjCDSWeight[i];
        }

        // For improved performance and clarity, duplicate the loop in the if and 
        // else, depending on whether we need to worry about credit event override
        // or not.
        // Note that each name's weight is scaled up (or down) so that the sum of all 
        // names covered for protection adds up to 1.0. This is in order to allow
        // names to be present in the credit index but not covered for protection
        if (!creditEventOverride) {
            // No override -> initialise the adjCDSArray with NULL override
            for (int i=0; i < numNamesInPortfolio; ++i) {
                (*adjCurvesCache)[i] = PriceCacheSP(
                    new PriceCache((*adjCDSArray)[i], 
                                adjCDSWeight[i]/cumulativeWeight, 
                                defOverride));
            }
        }
        else {
            for (int i=0; i < numNamesInPortfolio; ++i) {
                // Get the override for each name
                defOverride = creditEventOverride->getOverrideForName(adjCDSName[i]);

                (*adjCurvesCache)[i] = PriceCacheSP(
                    new PriceCache((*adjCDSArray)[i], 
                                    adjCDSWeight[i]/cumulativeWeight, 
                                    defOverride));
            }
        }
    }

    //-------------------
    // Now do the pricing
    //-------------------
    try {
        double price = 0.0;
        
        int numNamesInPortfolio = adjCurvesCache->size();
        // Price each CDS in adjCurvesCache and add their prices up
        for (int i=0; i < numNamesInPortfolio; ++i) {
            double thisSwapPv = (*adjCurvesCache)[i]->getPrice(this, model);
            price += thisSwapPv;
        }
        
        results->storePrice(price, discount->getCcy());

        // If pricing, produce the output requests too (ie, not if we
        // are tweaking)
        if (control && control->isPricing()) {
            addOutputRequests(control, results, model, true, CDSOutputs()); // CDSOuputs not used here
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CDS::addOutputRequests(Control*               control,
                            CResults*              results,
                            IForwardRatePricerSP   model,
                            bool                   isMultiName,
                            const CDSOutputs&      myOutput) const
{
    static const string method = "CDS::addOutputRequests";

    // supporting information
    DateTimeArraySP              feeDates;
    DateTimeArraySP              riskyFeeDates;
    AbstractCashFlowArrayConstSP fees;
    if (!!feeLeg) {
        feeDates = feeLeg->getCashFlowDates();
        riskyFeeDates = feeLeg->getRiskyCashFlowDates();
        fees = feeLeg->getCashFlows(model);
    }
    else {
        feeDates.reset(new DateTimeArray(0));
        riskyFeeDates.reset(new DateTimeArray(0));
        fees.reset(new AbstractCashFlowArray(0));
    }

    DateTime accruedStart;
    DateTime protStart;
    DateTime protEnd;
    double   recoveryRate;
    if (!!contingentLeg) {
        accruedStart = contingentLeg->firstObservationStartDate();
        protStart = valueDate.max(contingentLeg->firstObservationStartDate());
        protEnd   = contingentLeg->lastObservationEndDate();
        recoveryRate = contingentLeg->getRecovery(*(underlying.getSP()));
    }
    else {
        //zero length protection period
        accruedStart = valueDate;
        protStart = valueDate;
        protEnd   = valueDate;
        recoveryRate = 1.0; //100%
    }

    IDecretionCurveConstSP prepay = underlying->getPrepayCurve();

    DefaultRatesSP  psDefRates       = underlying->defaultRates();
    CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();
    IObjectSP       currentSpread    = underlying->getCurrentSpreadOrUntweakable(
                                                                    underlying->spotDate(valueDate),
                                                                    getMaturity());
    string          ccyName          = discount->getName();
    CashFlowArraySP feePayments      = AbstractCashFlow::asCashFlows(fees, model);
    DateTime        effDate          = underlying->spotDate(valueDate);

    ICreditEventOverrideNameSP nameOverride;
    if (!!creditEventOverride) {
        nameOverride = creditEventOverride->getOverrideForName(underlying->getName());
    }

    // IMPLIED_CDS_SPREAD and CDS_RISKY_DURATION
    OutputRequest* requestImpliedSpread = 
        control->requestsOutput(OutputRequest::IMPLIED_CDS_SPREAD);

    OutputRequest* requestDuration = 
        control->requestsOutput(OutputRequest::CDS_RISKY_DURATION);

    if (requestImpliedSpread || requestDuration)
    {
        try
        {
            //for now, just use risky fee dates
            //longer term, we should probably take conventions off the
            //curve to build a par instrument.
            CDSPricer impliedSpreadPricer = CDSPricer(
                1.0, // coupon rate
                *(riskyFeeDates.get()),
                psDefRates,
                -1.0, // unit notional
                recoveryRate,
                payAccruedFee,
                accrualDCC,
                valueDate,
                protStart,
                protEnd,
                accruedStart,
                discount.getSP(),
                prepay);

            double parSpread = 0.0;
            double riskyDuration = 0.0;
                
            // do only one call to impliedParSpreadAndDuration(...)
            impliedSpreadPricer.impliedParSpreadAndDuration(parSpread, 
                                                            riskyDuration);

            if (requestImpliedSpread) {
                results->storeRequestResult(requestImpliedSpread, 
                                            parSpread);
            }
            if (requestDuration) {
                results->storeRequestResult(requestDuration, 
                                            riskyDuration);
            }
        }
        catch (exception)
        {
            if (requestImpliedSpread)
            {
                results->storeNotApplicable(requestImpliedSpread);
            }
            if (requestDuration)
            {
                results->storeNotApplicable(requestDuration);
            }
        }
    }

    // CLEAN_DEFAULT_SPREAD_CURVE
    OutputRequest* request = 
        control->requestsOutput(OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE);

    if (request && cleanSpreadCurve.get()) {
        string name = underlying->getName();
        OutputNameConstSP defaultOutput(new OutputName(name));
        IObjectSP    cflows(new CashFlowList(cleanSpreadCurve.get()));
        
        results->storeRequestResult(request, cflows, defaultOutput);
    }

    // PAYMENT_DATES
    request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
    if (request) {
        OutputRequestUtil::recordPaymentDates(control,results,feeDates.get()); 
    }

    // KNOWN_CASHFLOWS, actually input cashflows
    request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
    if (request) {
        if (isMultiName)
        {
            CashFlowArraySP cfls = knownCashFlows(model);
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfls.get()); 
        }
        else
        {
            CashFlowArraySP cfls = knownCashFlowsWithPrepayment(model);
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfls.get()); 
        }
    }

    // IMPLIED_DEFAULT_PROBABILITY at maturity
    request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
    if (request)
    {
        try
        {
            // set the rolling effective date - for bootstrapping
            double defaultProb = 1.0 - psDefRates->calcDefaultPV(effDate, 
                                                                protEnd);
            results->storeRequestResult(request, defaultProb);
        }
        catch (exception&)
        {
            results->storeNotApplicable(request);
        }
    }

    // THETA_ACCRUED_INTEREST
    request = control->requestsOutput(OutputRequest::THETA_ACCRUED_INTEREST);
    if (request) {
        double accIntToday = getAccruedInterest(valueDate, model);
        accIntToday *= -1; // Traditionally interest is reported as > 0 if paid (!)

        DateTime aiAccrueStartDate;
        DateTime aiAccrueEndDate;
        DateTime aiPaymentDate;
        double   aiAmount;

        //get information about the fee in scope
        feeLeg->getActiveFee(valueDate, //withRespectTo
                             accruedStart, //earliestAccrualDate
                             model,
                             aiAccrueStartDate,
                             aiAccrueEndDate,
                             aiPaymentDate,
                             aiAmount);

        HolidaySP hols(Holiday::weekendsOnly());
        BadDayConventionSP badDay(BadDayConventionFactory::make("Following"));
        DateTime tomorrow = valueDate.rollDate(1);
        tomorrow = badDay->adjust(tomorrow, hols.get());

        double accIntTomorrow = getAccruedInterest(tomorrow, model);
        accIntTomorrow *= -1; // Traditionally interest is reported as > 0 if paid (!)

        // Account for the case where today and tomorrow are in different fee periods
        double thetaAI = 0.0;
        if ((valueDate < aiPaymentDate && aiPaymentDate <= tomorrow) &&
            !(accruedStart >= aiPaymentDate) && 
            !(valueDate <= accruedStart && accruedStart <= aiPaymentDate))
        {
            thetaAI = (accIntTomorrow + aiAmount) - accIntToday;
        }
        else {
            thetaAI = accIntTomorrow - accIntToday;
        }

        results->storeRequestResult(request, thetaAI);
    } 

    // ACCRUED_INTEREST
    request = control->requestsOutput(OutputRequest::ACCRUED_INTEREST);
    if (request) {
        if (!feeLeg) {
            results->storeNotApplicable(request);
        }
        else {
            // do not include payments on valueDate if multiName
            double accInt = getAccruedInterest(valueDate, model);
            accInt *= -1; // Traditionally interest is reported as > 0 if paid (!)

            results->storeRequestResult(request, accInt);
        }
    } 

    // RECOVERY_VALUE
    request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
    if (request) {
        if (!contingentLeg) {
            results->storeNotApplicable(request);
        }
        else {
            double defaultValue = contingentLeg->recoveredValue(
                valueDate,*(underlying.getSP()), prepay,
                IDiscountCurveRisky::RECOVER_R);
            results->storeRequestResult(request, defaultValue);
        }
    } 

    // CURRENT_SPREAD
    request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
    if (request) {
        results->storeRequestResult(request, currentSpread);
    }

    // IND_CDS_PAR_SPREAD
    // This one is the same as CURRENT_SPREAD
    // but goes in its own packet and is qualified by curve name
    request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
    if (request) {
        results->storeRequestResult(
            request, 
            currentSpread,
            OutputNameConstSP(new OutputName(underlying->getName())));
    }

    // PAR_SPREAD_CURVE
    request = control->requestsOutput(OutputRequest::PAR_SPREAD_CURVE);
    if (request)
    {
        results->storeRequestResult(
            request,
            ExpiryResult::createExpiryResultArray(
                *underlying->getParSpreadsExpiries().get(),
                *underlying->getParSpreads().get()),
            OutputNameSP(new OutputName(underlying->getName())));
    }

    // FEE_LEG_RISKLESS_FV - fee leg fees discounted against the riskless curve
    // as opposed to the fairvalue of riskless fees in the fee leg.
    request = control->requestsOutput(OutputRequest::FEE_LEG_RISKLESS_FV);
    if (request) {
        double factor = prepay->getFactor(valueDate);
        double balance = prepay->pv(valueDate);  
        double risklessPV = discount->pv(*feePayments, valueDate) * factor * balance;
        results->storeRequestResult(request, risklessPV);
    }

    // INDEX_BASIS
    request = control->requestsOutput(OutputRequest::INDEX_BASIS);
    if (request)
    {
        if (isMultiName)
        {
            try
            {
                CreditIndexConstSP creditIndex = CreditIndexConstSP::dynamicCast(underlying.getSP());
                CreditIndexBasisConstSP indexBasis = creditIndex->getIndexBasis();
                const DoubleArrayConstSP basis     = indexBasis->getBasis();
                const ExpiryArrayConstSP expiries  = indexBasis->getBasisExpiries();

                ExpiryResultArraySP outputBasis    = 
                    ExpiryResult::createExpiryResultArray(*expiries, *basis);

                OutputRequestUtil::recordExpiryResultArray(
                    control, 
                    results, 
                    outputBasis.get(),
                    OutputRequest::INDEX_BASIS);
            }
            catch (exception&)
            {
                results->storeNotApplicable(request);
            }
        }
        else
        {
            results->storeNotApplicable(request);
        }
    }

    request = control->requestsOutput(OutputRequest::FEE_LEG_FV);
    if (request) {
        if (isMultiName) {
            results->storeNotApplicable(request);
        }
        else {
            results->storeRequestResult(request, myOutput.feePv);
        }
    }

    request = control->requestsOutput(OutputRequest::CONTINGENT_LEG_FV);
    if (request) {
        if (isMultiName) {
            results->storeNotApplicable(request);
        }
        else {
            results->storeRequestResult(request, myOutput.ctgPv);
        }
    }

    request = control->requestsOutput(OutputRequest::ACCRUED_INT_FV);
    if (request) {
        if (isMultiName) {
            results->storeNotApplicable(request);
        }
        else {
            results->storeRequestResult(request, myOutput.accInt);
        }
    }
}

bool CDS::paysAccrued() const
{
    return payAccruedFee;
}

//--------------------
// CInstrument methods
//--------------------

/** Called immediately after object constructed */
void CDS::validatePop2Object()
{
    static const string method = "CDS::validatePop2Object";

    try
    {
        if (payAccruedFee)
        {
            if (accrualDayCountConvention != "")
            {
                //build from the specified string
                //- if not specified, getMarket will
                //- create from the underlying
                accrualDCC.reset(
                    DayCountConventionFactory::make(
                        accrualDayCountConvention));
            }
        }

        if (settlementBadDayConvention.empty())
        {
            //build a "None" flavour of bad day adjuster
            settlementBDC = BadDayConventionSP(
                BadDayConventionFactory::make("None"));
        }
        else
        {
            // build from the supplied string
            settlementBDC = BadDayConventionSP(
                BadDayConventionFactory::make(settlementBadDayConvention));
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** Allow the instrument to retrieve its market data */
void CDS::GetMarket(const IModel*            model,
	                const CMarketDataSP      market)
{
    const static string method = "CDS::getMarket";
    try {
        market->GetReferenceDate(valueDate);

        // Get the interest rate discount curve
        discount.getData(model, market);

        // Get the credit underlying
        //underlying.getData(model, market);
        ICDSParSpreads::getMarketData(model, 
                                      market.get(),
                                      discount.getName(),
                                      underlying);

        if (feeLeg.get()) {
            feeLeg->getMarket(model, market.get());
        }

        if (!accrualDCC) {
            // default to par curve DCC
            accrualDCC = DayCountConventionSP(underlying->dayCountConv());
        }

        // Get settlement holidays
        if (settlementHols.isEmpty()) {
            // default to no holidays (so business days = calendar days)
            settlementHols.setObject(MarketObjectSP(Holiday::noHolidays()));
        } 
        else {  
            settlementHols.getData(model, market);
        }

        // If there is a credit event override, let it get its market data
        if (!!creditEventOverride) {
            creditEventOverride->getMarket(model, market.get());
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** called after market data has been retrieved */
void CDS::Validate() {
    const static string method = "CDS::Validate";
    try {
        if (!!contingentLeg) {
            DateTime protEndDate = contingentLeg->lastObservationEndDate();

            if (lastTriggerDate.empty()) {
                // lastTriggerDate is optional and may not be present. If so,
                // set it to the maturity date
                lastTriggerDate = protEndDate;
            }
            else if (lastTriggerDate < protEndDate) {
                throw ModelException(method,
                                        "lastTriggerDate (" + 
                                        lastTriggerDate.toString() + 
                                        ") cannot be before protection end date (" +
                                        protEndDate.toString() + ").");
            }
        }

        //we lose this validation now that we are mixing single & multi
        //name capability...if the wrong name is used, we'll just get
        //an empty sp back for the per-name params, rather than failing
        //I've left this in so that it doesnt get "forgotten" since
        //we may need to do something tighter at a later date
        //if (!!creditEventOverride &&
        //    (creditEventOverride->getName() != underlying.getName()))
        //{
        //    throw ModelException(method,
        //                            "A creditEventOverride is present, but its "
        //                            "name (" + creditEventOverride->getName() +
        //                            ") is not the same as the underlying curve's "
        //                            "name (" + underlying.getName() + ").");
        //}

        if ((!!triggerDelay && !defaultToSettlementDelay) ||
            (!triggerDelay && !!defaultToSettlementDelay))
        {
            throw ModelException(method,
                                    "triggerDelay and defaultToSettlementDelay must "
                                    "both be provided or both be removed.");
        }

        if (!!triggerDelay) { // both delays are provided
            if (defaultToSettlementDelay->intValue() < triggerDelay->intValue()) {
                throw ModelException(method,
                                        "The defaultToSettlementDelay (" +
                                        Format::toString(defaultToSettlementDelay->intValue()) +
                                        ") cannot be smaller than the triggerDelay(" +
                                        Format::toString(triggerDelay->intValue()) + 
                                        ").");
            }

            if (triggerDelay->intValue() < 0) {
                throw ModelException(method,
                                        "The triggerDelay is negative (" +
                                        Format::toString(triggerDelay->intValue()) +
                                        "), and negative delays are not accepted.");
            }
            if(defaultToSettlementDelay->intValue() < 0) {
                throw ModelException(method,
                                        "The defaultToSettlementDelay is negative (" +
                                        Format::toString(defaultToSettlementDelay->intValue()) +
                                        "), and negative delays are not accepted.");
            }
        }

        //if the underlyer is a credit index
        //then we restrict it to be of the same currency
        //as the discount curve
        //...nb if we are pricing with ClosedFormCDSPS then
        //...this is not a requirement, but we dont have the
        //...model here
        //the cast throws an exception if it fails
        bool doTest = true;
        CreditIndexSP ci;
        try
        {
            ci = CreditIndexSP::dynamicCast(underlying.getSP());
        }
        catch (exception&)
        {
            //swallow the exception
            doTest = false;
        }

        if (doTest && (discount->getCcy() != ci->getCcy()))
        {
            throw ModelException(method,
                                    "In the current implementation the " 
                                    "currency of the CIS instrument (" + 
                                    discount->getCcy() + 
                                    ") must be the same as the currency of " 
                                    "the index being swapped (" +
                                    ci->getCcy() + ").");
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/* Notification that (some) underlying fields have changed */
void CDS::fieldsUpdated(const CFieldArray& fields)
{
    // Since fields in the instrument are used to price each of the
    // elements in the adjCurvesCache array, notify them that 
    // something has changed. 
    // Of course, unless the field changing is the actual adjCurvesCache,
    // in which case the appropriate element(s) of the array have been
    // notified already
    unsigned int i;
    for (i=0; i<fields.size(); ++i) {
        if (fields[i]->getName() != adjCurvesCacheFieldName) {
            // Some field(s) other that the caches has changed, so need to 
            // update all elements of the adjCurvesCache array
            int numNamesInPortfolio = adjCurvesCache->size();
            for (int j=0; j < numNamesInPortfolio; ++j) {
                (*adjCurvesCache)[j]->fieldsUpdated(fields);
            }
            break; // We are done
        }
    }
}

/** Returns the value date (aka today) the instrument is currently pricing for */
DateTime CDS::getValueDate() const
{
    return valueDate;
}

/** Returns the name of the instrument's discount currency */
string CDS::discountYieldCurveName() const
{
    return discount.getName();
}

//--------------
// ICDS methods
//--------------

/** Return the fee leg */
ICreditFeeLegSP CDS::getFeeLeg() const
{
    return feeLeg;
}

/** Return the contingent leg */
ICreditContingentLegSP CDS::getContingentLeg() const
{
    return contingentLeg;
}

/** Return Accrued Interest (cash amount, not percentage) for settlement on 
    settlementDate */
double CDS::getAccruedInterest(const DateTime&      valuationDate,
                               IForwardRatePricerSP model) const
{
    double ai = !feeLeg ? 0.0 : feeLeg->getAccruedInterest(valuationDate, model);
    return ai;
}


// A solaris opt bug requires adding some funky stuff for it to stop 
// crashing in CDS::getPVGivenDefault in the presence of exceptions (test 
// cdsnewinp/cds_override_wrongDefaultDate_error.xml was crashing)
#ifdef sun
static void dummy(int a, ...){}
#else
static void dummy(int a){}
#endif

/**
  * Returns the present value of the instrument at defaultDate,
  * conditional on default at defaultDate.
  */
double CDS::getPVGivenDefault(const DateTime&            valuationDate,
                              const DateTime&            defaultDate,
				              const IDiscountCurveRisky& crv,
                              IForwardRatePricerSP       model) const
{
    static const string method = "CDS::getPVGivenDefault";

    IDecretionCurveConstSP prepay = underlying->getPrepayCurve();
    IBadDayAdjusterConstSP bdAdjuster(IBadDayAdjusterConstSP::attachToRef(this));
    YieldCurveConstSP discountCurve = discount.getSP();
    ICreditEventOverrideNameSP nameOverride;
    if (!!creditEventOverride) {
        nameOverride = creditEventOverride->getOverrideForName(underlying->getName());
    }

    double feeLegValue = 0.0;
    try {
        // problems with solaris.opt - put this statement in a separate
        // try catch block and add a call to the static dummy method.
        if (!!feeLeg) {
            try {
                // problems with solaris.opt - put this statement in a separate
                // try catch block and add a call to the static dummy method.
                dummy(0);
                feeLegValue = feeLeg->getFeeLegDefaultedPV(
                    valuationDate, defaultDate, 
                    !contingentLeg ? DateTime() : 
                                     contingentLeg->firstObservationStartDate(),
                    !contingentLeg ? DateTime() : 
                                     contingentLeg->lastObservationEndDate(),
                    true, // Pricing a single name, so allow including payments on valueDate
                    discountCurve, prepay, model,
                    bdAdjuster,
                    nameOverride, triggerDelay,
                    defaultToSettlementDelay, lastTriggerDate);
            }
            catch (exception& e) {
                // problems with solaris.opt - add a call to the static dummy method
                dummy(0);
                throw ModelException(e, method, "Failed to price the fee leg");
            }
        }
    }
    catch (exception&) {
        // problems with solaris.opt - add a call to the static dummy method.
        dummy(0);
        throw ;
    }

    double contingentLegValue = 0.0;
    if (!!contingentLeg)
    {
        try {
            // problems with solaris.opt - put this statement in a separate
            // try catch block and add a call to the static dummy method.
            dummy(0);
            contingentLegValue = contingentLeg->getContingentLegDefaultedPV(
                valuationDate,
                crv, discountCurve, prepay,
                defaultDate, bdAdjuster,
                true, // Pricing a single name, so allow including payments on valueDate
                nameOverride, triggerDelay,
                defaultToSettlementDelay, lastTriggerDate);
        }
        catch (exception& e) {
            // problems with solaris.opt - add a call to the static dummy method
            dummy(0);
            throw ModelException(e, method, "Failed to price the ctg leg");
        }
    }

    // Fee Leg has the right sign here (negative if paying fees) so
    // just add the legs' pvs
    double price = contingentLegValue + feeLegValue;
   
    return price;
}

/** Determines if the instrument is perpetual. */
bool CDS::hasFiniteMaturity() const
{
    //no concept of perpetuals
    return true;
}

/** Returns the maturity of the instrument */
DateTime CDS::getMaturity() const {
    // max of { last fee leg pay date , last protection end date }

    IBadDayAdjusterConstSP bda(IBadDayAdjusterConstSP::attachToRef(this));

    DateTime lastFeeDate = (!feeLeg) ? valueDate : feeLeg->getLastPayDate();
    DateTime lastCtgDate = (!contingentLeg) ? valueDate : 
                                              contingentLeg->lastPayDate(bda);

    return lastFeeDate.max(lastCtgDate);
}

/**
 * Returns the earliest of the first accrual period, start of the contingent leg,
 * the first cash flow, etc.
 */
DateTime CDS::getStartDate() const
{
    static const string method = "CDS::getStartDate";

    DateTime firstFeeDate = valueDate;
    if (!!feeLeg)
    {
        //get first accrual start date
        AccrualPeriodArrayConstSP accPrds = feeLeg->getAccrualPeriods();
        AccrualPeriodSP accPrd = AccrualPeriodSP::constCast(accPrds->front());
        firstFeeDate = accPrd->startDate();
    }

    DateTime firstCtgDate = (!contingentLeg) ? valueDate : contingentLeg->firstObservationStartDate();

    return firstFeeDate.min(firstCtgDate);
}

/**
 * Get the base (initial) notional of the instrument.
 * Positive notional = long risk/short protection.
 */
double CDS::getNotional() const
{
    static const string method = "CDS::getNotional";
    
    //we can return a single value here if and only if the legs have a
    //single value, and they are the same

    double notional;

    if ((!feeLeg) && (!contingentLeg))
    {
        //no legs, so return a notional of 0.0
        notional = 0.0;
    }
    else if ((!feeLeg) && (!!contingentLeg))
    {
        //contingent leg only
        notional = contingentLeg->getContingentLegNotional();
    }
    else if ((!!feeLeg) && (!contingentLeg))
    {
        //fee leg only
        notional = feeLeg->getFeeLegNotional();
    }
    else
    {
        double feeLegNotional = (!feeLeg) ? 0.0 : feeLeg->getFeeLegNotional();
        double ctgLegNotional = (!contingentLeg) ? 0.0 : contingentLeg->getContingentLegNotional();

        if (feeLegNotional != ctgLegNotional)
        {
            throw ModelException(method,
                "Legs have different notionals");
        }
        else
        {
            //doesnt matter which
            notional = feeLegNotional;
        }
    }

    return notional;
}

/**
 * Set the base (initial) notional of the instrument.
 * Positive notional = short risk/long protection.
 */
void CDS::setNotional(double newNotional)
{
    static const string method = "CDS::setNotional";

    if (!!feeLeg)
    {
        feeLeg->setFeeLegNotional(newNotional);
    }

    if (!!contingentLeg)
    {
        contingentLeg->setContingentLegNotional(newNotional);
    }
}

/** Returns the day count convention used for accruals */
DayCountConventionSP CDS::getAccrualDcc() const
{
    return accrualDCC;
}

CashFlowArraySP CDS::getInstrumentCashFlows(IForwardRatePricerSP model) const
{
    //get flows from the fee leg
    AbstractCashFlowArrayConstSP feeLegCfls = feeLeg->getCashFlows(model);

    //return the conversion from abstract to cashflow
    return AbstractCashFlow::asCashFlows(feeLegCfls, model);
}

/** Returns the wrapper to the yield curve */
YieldCurveWrapper CDS::getYieldCurveWrapper() const
{
    return discount;
}

/** Returns the wrapper describing the credit curve */
ICDSParSpreadsWrapper CDS::getParSpreadsWrapper() const
{
    return underlying;
}

/**This has been overridden to set getPV = getFeeLegPV + getContingentLegPV.
   If you override this, you should still ensure that this is the case!*/
double CDS::getPV(const DateTime&              valuationDate,
                  const DateTime&              settlementDate,
                  const IDiscountCurveRisky&   crv,
                  const IDecretionCurveConstSP prepay,
                  IForwardRatePricerSP         model,
                  IBadDayAdjusterConstSP       bda,
                  CDSOutputs&                  myOutput) const
{
    //pricing this way (via individual legs) gives rise to differences versus CredDefSwap
    // 1) Contingent leg : the integration timeline misses out the fee payment dates
    //    which has a small effect when the risk free curve is not built with flat forwards
    // 2) Fee leg : CDSPricer negated riskless cashflows. In this implementation
    //    we leave them as specified in the inputs (migrated test cases had their riskless cashflows
    //    manipulated)

    //legs may be optional in the parent structure
    if (!!feeLeg)
    {
        DateTime earliestRiskyDate = (!contingentLeg) ? valuationDate
                                                      : contingentLeg->firstObservationStartDate();
        DateTime latestRiskyDate = (!contingentLeg) ? feeLeg->getLastPayDate()
                                                    : contingentLeg->lastObservationEndDate();

        myOutput.feePv = feeLeg->getFeeLegPV(valuationDate, valuationDate,
                                             earliestRiskyDate, latestRiskyDate,
                                             *(discount.getSP()), crv, prepay,
                                             false, accrualDCC, false, model); //do not include accrued
        if (payAccruedFee) {
            myOutput.accInt = feeLeg->getFeeLegAI(valuationDate, settlementDate,
                                                  earliestRiskyDate, latestRiskyDate,
                                                  accrualDCC, crv,
                                                  discount.getSP(), prepay, model);
        }
    }

    if (!!contingentLeg)
    {
        myOutput.ctgPv = contingentLeg->getContingentLegPV(valuationDate, settlementDate, crv, bda);
    }

    return myOutput.price();
}

double CDS::getPV(const DateTime&              valuationDate,
                  const IDiscountCurveRisky&   crv,
                  const IDecretionCurveConstSP prepay,
                  IForwardRatePricerSP         model,
                  IBadDayAdjusterConstSP       bda,
                  CDSOutputs&                  myOutput) const
{
    return getPV(valuationDate, valuationDate, crv, prepay, model, bda, myOutput);
}

//------------------------
// IBadDayAdjuster methods
//------------------------

/** Returns "date" bad day adjusted using the bad day convention
    * and holidays in this object */
DateTime CDS::badDayAdjust(const DateTime& date) const
{
    return settlementBDC->adjust(date, settlementHols.get());
}

/** Add a number of business days to a date */
DateTime CDS::addBusinessDays(const DateTime& from, int busDays) const
{
    return settlementHols->addBusinessDays(from, busDays);
}

//--------------------
//Theta::Shift methods
//--------------------

bool CDS::sensShift(Theta* shift)
{
    valueDate = shift->rollDate(valueDate);
    return true;
}

//--------------------------------------
// ClosedFormCDSPS::IIntoProduct methods
//--------------------------------------

/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPS::IProduct* CDS::createProduct(ClosedFormCDSPS* model) const
{
    return new CDSClosedFormCDSPS(this);
}

//------------------------------------------
// ClosedFormCDSBasket::IIntoProduct methods
//------------------------------------------

/** Implementation of ClosedFormCDSBasket::IntoProduct interface */
ClosedFormCDSBasket::IProduct* CDS::createProduct(ClosedFormCDSBasket* model) const
{
    return new CDSClosedFormCDSBasket(this);
}

//----------------------
// CDS methods (private)
//----------------------
CDS::CDS() : CInstrument(TYPE)
{
    //initialise price cache
    adjCurvesCache = PriceCacheArraySP(new PriceCacheArray(0));
}

void CDS::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CDS, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(ICDS);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
    IMPLEMENTS(ClosedFormCDSBasket::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(valueDate,                  "valuation date");
    FIELD(payAccruedFee,              "is accrued interest paid upon default");
    FIELD(accrualDayCountConvention,  "required when payAccruedFee = True");
    FIELD(feeLeg,                     "the annuity leg");
    FIELD(contingentLeg,              "the protection leg");
    FIELD(underlying,                 "risky asset");
    FIELD(discount,                   "risk-free discount curve");

    FIELD(creditEventOverride,        "Override for credit event parameters");
    FIELD(triggerDelay,               "Delay, in days, between default and "
                                      "eventDeterminationDate - used to "
                                      "estimate eventDeterminationDate "
                                      "before it is known.");
    FIELD(defaultToSettlementDelay,   "Delay between credit event and settlement.");
    FIELD(lastTriggerDate,            "Last date when a default occurred during "
                                      "the protection period can be triggered.");
    FIELD(settlementHols,             "Holidays for the tranche. Used to "
                                      "adjust the delays regarding defaults' "
                                      "settlements. Default: weekends only");
    FIELD(settlementBadDayConvention, "used in conjunction with settlement holidays");

    FIELD(accrualDCC,                 "transient form of accrualDayCountConvention");
    FIELD(settlementBDC,              "transient form of settlementBadDayConvention");

    FIELD(adjCurvesCache,             "transient field with the cache for the "
                                      "curves' prices.");

    FIELD_MAKE_OPTIONAL(valueDate);     //will be returned from market cache if necessary
    FIELD_MAKE_OPTIONAL(feeLeg);        //allows structure without annuity leg
    FIELD_MAKE_OPTIONAL(contingentLeg); //allows structure without protection leg
    FIELD_MAKE_OPTIONAL(creditEventOverride);
    FIELD_MAKE_OPTIONAL(triggerDelay);
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelay);
    FIELD_MAKE_OPTIONAL(lastTriggerDate);
    FIELD_MAKE_OPTIONAL(settlementHols);
    FIELD_MAKE_OPTIONAL(settlementBadDayConvention);

    FIELD_MAKE_TRANSIENT(accrualDCC);
    FIELD_MAKE_TRANSIENT(settlementBDC);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(adjCurvesCache);
}

IObject* CDS::defaultConstructor()
{
    return new CDS();
}

//Multi name known cashflows
CashFlowArraySP CDS::knownCashFlows(IForwardRatePricerSP model) const
{
    CashFlowArraySP cashFlow;
    int numNamesInPortfolio = adjCurvesCache->size();

    if (numNamesInPortfolio) { // ie, there is at least 1 element
        // Initialise cashFlow with the cash flows of the 1st name
        cashFlow = (*adjCurvesCache)[0]->knownCashFlows(this, model);
    }

    // Obtain the known cashflows of the remaining names, and merge them
    CashFlowArraySP nextNamesCashFlows;
    for (int i=1; i < numNamesInPortfolio; ++i) {
        nextNamesCashFlows = (*adjCurvesCache)[i]->knownCashFlows(this, model);
        cashFlow = CashFlow::merge(cashFlow, nextNamesCashFlows);
    }
    return cashFlow;
}

//Fee leg scaled by pre-payment
CashFlowArraySP CDS::knownCashFlowsWithPrepayment(IForwardRatePricerSP model) const
{
    static const string method = "CDS::knownCashFlowsWithPrepayment";

    try
    {
        CashFlowArraySP cashFlow = CashFlowArraySP(new CashFlowArray(0));

        if (!!feeLeg)
        {
            IDecretionCurveConstSP psPrepay = underlying->getPrepayCurve();
            double factor = psPrepay->getFactor(valueDate);

            AbstractCashFlowArrayConstSP allcfls = feeLeg->getCashFlows(model);
            for (int i=0; i<allcfls->size(); i++)
            {
                DateTime date = (*allcfls)[i]->getPayDate();
                double amount = (*allcfls)[i]->getAmount(model);

                //only risky cashflows get scaled by the pre-payment curve
                //which is consistent with the pricing
                double balance = 1.0;
                double useFactor = 1.0;
                if (!((*allcfls)[i]->isRiskFree()))
                {
                    balance = psPrepay->pv(date.rollDate(-1)); 
                    useFactor = factor;
                }
                CashFlow cf(psPrepay->getSettlement()->settles(date),
                            amount * useFactor * balance);

                cashFlow->push_back(cf);
            }
        }

        return cashFlow;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const CDS::TYPE = 
CClass::registerClassLoadMethod("CDS", typeid(CDS), load);

//for linking in
bool CDSLoad()
{
    return (CDS::TYPE != 0);
}

DRLIB_END_NAMESPACE
