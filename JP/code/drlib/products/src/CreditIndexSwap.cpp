//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditIndexSwap.cpp
//
//   Description : Credit default Index Swap: Instrument equivalent to
//                 a CredDefSwap, but where the underlier is an index 
//                 (CDX or iTraxx style), rather than a single name.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditIndexSwap.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/AdjustedCDSParSpreads.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/CreditIndex.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/ICreditEventOverride.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


/******************************************************************************
 *                           CreditIndexSwap::PriceCache
 ******************************************************************************/

/** PriceCache is used to price and cache the price of the underlying 
 * CDSs in a CIS. Note this class is tightly linked to the model since
 * it controls how the pricing gets computed. Actually, it should be 
 * moved into the model...*/
class CreditIndexSwap::PriceCache : public CObject {
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
    double getPrice(const CreditIndexSwap* const cis,
                    IForwardRatePricerSP model);  

    /** Return the known cash flows for this name's CDS, taking potential
     * defaults into account */
    CashFlowArraySP knownCashFlows(const CreditIndexSwap* const cis,
                                   IForwardRatePricerSP model) const;

private:
    /** Prices a defaulted curve */
    double priceDefaultedCurve(const CreditIndexSwap* const cis,
                               IForwardRatePricerSP model);

    /** Prices the fee leg of a defaulted name using the credit event override 
        if available */
    double priceFeeLegGivenDefault(const CreditIndexSwap* const cis,
                                   IForwardRatePricerSP model) const;

    /** Prices the contingent leg of a defaulted CDS, using the credit event
     * override if available */
    double priceContingentLegGivenDefault(const CreditIndexSwap* const cis) const;

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
CreditIndexSwap::PriceCache::PriceCache(
        AdjustedCDSParSpreadsSP     copyCDS,
        double                      namesWeight,
        ICreditEventOverrideNameSP  nameOverride) :
    CObject(TYPE),
    isPriceValid(false),
    adjCDS(copyCDS.clone()), // Clone the SP
    adjCDSWeight(namesWeight),
    nameOverride(nameOverride)
{}


CreditIndexSwap::PriceCache::PriceCache() : 
    CObject(TYPE), isPriceValid(false), adjCDSWeight(0.0)
{}


IObject* CreditIndexSwap::PriceCache::defaultPriceCache() {
    return new PriceCache();
}


/* Notification that some internal fields have changed -> Reset the 
 * isPriceValid flag since the price is no longer valid */
void CreditIndexSwap::PriceCache::fieldsUpdated(const CFieldArray& fields) {
    isPriceValid = false;
}


/* Obtains the price of the underlying curve. Computes the price if 
 * required (using instrument data) and caches the result */
double CreditIndexSwap::PriceCache::getPrice(const CreditIndexSwap* const cis,
                                             IForwardRatePricerSP model) {
    static const string method = "CreditIndexSwap::PriceCache::getPrice";

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
                CashFlowArraySP riskyFeePayments    = cis->feeLeg->getRiskyCashFlows(model);
                CashFlowArraySP risklessFeePayments = cis->feeLeg->getRisklessCashFlows(model);
                DefaultRatesSP defaultRates         = adjCDS->defaultRates();
                double recoveryRate                 = adjCDS->getRecovery();

                // JLH - We know this is WRONG: accrual start and end dates are NOT
                // necessarily protectionStartDate and protectionEndDate. Even worse 
                // is the fact that CDSPricer does not take a proper feeLeg but a 
                // CashFlowArray... Do not worry about this yet.
                CDSPricer cdsPricer(riskyFeePayments,
                                    risklessFeePayments,
                                    defaultRates,
                                    cis->notional,
                                    recoveryRate,
                                    cis->payAccruedFee,
                                    cis->swapDcc,
                                    cis->valueDate,
                                    cis->protectionStartDate,
                                    cis->protectionEndDate,   // Protection and accrued end
                                    cis->protectionStartDate, // Accrual start date
                                    cis->discount.getSP());
            
                // CAUTION: The risky fees are obtained by calling 
                // feeLeg->getxxxxCashFlows, which return feePayments with the same
                // sign as the notional, which is wrong. Work around this by changing
                // the sign of the riskless payments here.
                // Note: the same has been done in CDSPricer.cpp
                totalPrice = -1.0 * cdsPricer.price();
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


/** Prices the fee leg of a defaulted name using the credit event override info 
    if available */
double CreditIndexSwap::PriceCache::priceFeeLegGivenDefault(
    const CreditIndexSwap* const cis,
    IForwardRatePricerSP model) const
{
    const DateTime& defaultDate   = adjCDS->getDefaultDate();
    IDecretionCurveConstSP prepay = adjCDS->getPrepayCurve();

    double feeLegPv = 0.0;
    if (!!cis->feeLeg) {
        feeLegPv = cis->feeLeg->getFeeLegDefaultedPV(
            cis->valueDate,
            defaultDate,
            cis->protectionStartDate,
            cis->protectionEndDate,
            false, // do not include payments on valueDate
            cis->discount.getSP(),
            prepay,
            model,
            cis->creditIndex.getSP(),  // IBadDayAdjusterSP
            nameOverride,
            cis->triggerDelay,
            cis->defaultToSettlementDelay,
            cis->lastTriggerDate);
    }

    return feeLegPv;
}


/** Prices the contingent leg of a defaulted CDS, using the credit event
 * override if available */
double CreditIndexSwap::PriceCache::priceContingentLegGivenDefault(
    const CreditIndexSwap* const cis) const 
{
    const DateTime& curveDefaultDate = adjCDS->getDefaultDate();

    CashFlowArraySP contingentCashFlows(new CashFlowArray(0));
    if ((curveDefaultDate < cis->protectionStartDate) ||
        (curveDefaultDate > cis->protectionEndDate)) // Both are included
    {
        // Default happened outside the protection period, so
        // contingentCashFlows is empty
    }
    else {
        const double recoveryRate = adjCDS->getRecovery();
        if (!nameOverride) {
            // There is no override, so generate the one cashflow "manually".
            DateTime defaultPaymentDate;
            if (!cis->defaultToSettlementDelay) {
                // No adjustment is done here.
                defaultPaymentDate = curveDefaultDate;
            }
            else {
                // Note this means both triggerDelay and defaultToSettlementDelay
                // are not NULL.
                defaultPaymentDate = 
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        curveDefaultDate,
                        cis->defaultToSettlementDelay,
                        cis->valueDate,
                        cis->valueDate,
                        cis->lastTriggerDate,
                        cis->creditIndex.getSP());  // IBadDayAdjusterSP
            }

            // Generate the one cashflow "manually"
            if ((defaultPaymentDate > cis->valueDate) &&
                (curveDefaultDate <= cis->lastTriggerDate))
            {
                // Notional < 0 means short protection
                // Notional > 0 means long protection
                contingentCashFlows->push_back(
                    CashFlow(defaultPaymentDate,
                             cis->notional * (1.0 - recoveryRate)));
            }
            else {
                // Since this contingent payment would have been paid before
                // excludeCashFlowsBeforeDate, or the default was triggered 
                // too late, return an empty cash flow array
            }
        }
        else {
            // Use the override
            CtgLegLossPerDefaultArraySP contingentReductions =
                nameOverride->historicContingentLegLosses(
                    cis->notional,
                    recoveryRate,
                    cis->lastTriggerDate,
                    curveDefaultDate,
                    cis->creditIndex.getSP());  // IBadDayAdjusterSP

            CashFlowArraySP allCashFlows = 
                CtgLegLossPerDefault::getReductions(contingentReductions);

            // Insert all future cashflows
            int numCashFlows = allCashFlows->size();
            contingentCashFlows->reserve(numCashFlows);
            for (int i=0; i < numCashFlows; ++i) {
                CashFlow cf((*allCashFlows)[i]); // for ease
                if (cf.date > cis->valueDate) {
                    contingentCashFlows->push_back(cf);
                }
            }
        }
    }

    double contingentLegValue = 0.0;
    auto_ptr<IYieldCurve::IKey> key(cis->discount->logOfDiscFactorKey());
    // pv feeCashFlows
    for (int i=0; i < contingentCashFlows->size(); ++i) {
        CashFlow cf = (*contingentCashFlows)[i]; // For ease
        const double pv = exp(key->calc(cis->valueDate, cf.date));
        contingentLegValue += pv * cf.amount;
    }

    IDecretionCurveConstSP psPrepay = adjCDS->getPrepayCurve();
    const double factor  = psPrepay->getFactor(cis->valueDate);
    const double balance = psPrepay->pv(curveDefaultDate);

    contingentLegValue *= factor * balance;

    return contingentLegValue;
}



/** Price a defaulted curve */
double CreditIndexSwap::PriceCache::priceDefaultedCurve(
    const CreditIndexSwap* const cis,
    IForwardRatePricerSP model)
{
    static const string method(
        "CreditIndexSwap::PriceCache::priceDefaultedCurve");

    double feeLegValue = priceFeeLegGivenDefault(cis, model);
    double contingentLegValue = priceContingentLegGivenDefault(cis);

    double price = contingentLegValue + feeLegValue;
    return price;
}


//------------------------------------------------------------------------------
//                    Methods related to output requests
//------------------------------------------------------------------------------

/** Returns all known cash flows */
CashFlowArraySP CreditIndexSwap::PriceCache::knownCashFlows(
    const CreditIndexSwap* const cis,
    IForwardRatePricerSP model) const 
{
    CashFlowArraySP knownCashFlows;

    // Risky payments depend on whether a default happened or not
    if (adjCDS->defaulted()) {
        // Get the fee leg cashflows - pass an emtpy DateTime to indicate
        // we are interested in all cashflows, ie, no cashflows should be
        // excluded
        knownCashFlows = cis->feeLeg->generateKnownCashFlowsGivenDefault(
            cis->getValueDate(),
            adjCDS->getDefaultDate(),
            DateTime(), // exclude no cashflows
            cis->protectionStartDate,
            cis->protectionEndDate,
            false, // This is CIS so do not allow including payments on valueDate
            model,
            nameOverride,
            cis->creditIndex.getSP(), // IBadDayAdjuster
            cis->triggerDelay,
            cis->defaultToSettlementDelay,
            cis->lastTriggerDate);
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
void CreditIndexSwap::PriceCache::load(CClassSP& clazz){
    // clazz->setPublic(); - Do NOT make visible to EAS/spreadsheet
    REGISTER(CreditIndexSwap::PriceCache, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultPriceCache);

    FIELD(isPriceValid, "Whether the cached price is or not valid");
    FIELD(price,        "The cached price");
    FIELD(adjCDSWeight, "The weight of this curve in the index");
    FIELD       (adjCDS,       "The curve whose price is cached here");
    FIELD       (nameOverride, "Override for this name's default parameters");
}


CClassConstSP const CreditIndexSwap::PriceCache::TYPE = 
    CClass::registerClassLoadMethod("CreditIndexSwap::PriceCache", 
                                    typeid(CreditIndexSwap::PriceCache), 
                                    load);

// Array has to have its own type
// Work around bug in msvc7
typedef CreditIndexSwap::PriceCacheArray CreditIndexSwapPriceCacheArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("CreditIndexSwap::PriceCacheArray", CreditIndexSwapPriceCacheArray);


/******************************************************************************
 *                              CreditIndexSwap definitions
 ******************************************************************************/
/** Credit Index Swap (CIS) is an instrument used to buy/sell protection
 * on an index (ie, in all the component names of the index).
 * The swap is priced as a basket of CDS, where the underlying names in 
 * the CDSs are the (index-basis adjusted) underlying names in the index.
 * The index basis computation and adjustment is performed within the 
 * creditIndex field of the CIS.
 * The adjCurvesCache is used to cache the price of each (index-basis 
 * adjusted) CDS */
CreditIndexSwap::CreditIndexSwap() : CInstrument(TYPE)
{}

CreditIndexSwap::~CreditIndexSwap()
{}


/** Allow the instrument to retrieve its market data */
void CreditIndexSwap::GetMarket(const IModel* model, 
                                const CMarketDataSP market) 
{
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
    creditIndex.getData(model, market);

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
            if(!!defOverride) {
                defOverride->getMarket(model, market.get());
            }
            (*adjCurvesCache)[i] = PriceCacheSP(
                 new PriceCache((*adjCDSArray)[i], 
                                adjCDSWeight[i]/cumulativeWeight, 
                                defOverride));
        }
    }
}


/** Called immediately after object constructed */
void CreditIndexSwap::validatePop2Object() {
    static const string method = "CreditIndexSwap::validatePop2Object";
    if (protectionStartDate >= protectionEndDate) { 
        // JLH - should be > rather than >= but we currenty do not count
        // days "properly", and dcc->year(date, date) returns 0.
        throw ModelException(method,
                             "protectionStartDate (" + 
                             protectionStartDate.toString() + 
                             ") should be before protectionEndDate (" +
                             protectionEndDate.toString() + ").");
    }
    swapDcc.reset(DayCountConventionFactory::make(dcc));

    if (lastTriggerDate.empty()) {
        // lastTriggerDate is optional and may not be present. If so,
        // set it to the protectionEndDate
        lastTriggerDate = protectionEndDate;
    }
    else if (lastTriggerDate < protectionEndDate) {
        throw ModelException(method,
                             "lastTriggerDate (" + 
                             lastTriggerDate.toString() + 
                             ") cannot be before protectionEndDate (" +
                             protectionEndDate.toString() + ").");
    }

    if ((!!triggerDelay && !defaultToSettlementDelay) ||
        (!triggerDelay && !!defaultToSettlementDelay))
    {
        throw ModelException(method,
                             "triggerDelay and defaultToSettlementDelay are "
                             "optional but both of them have to be provided "
                             "if one is provided - alternativeley, remove "
                             "both of them.");
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

}


/** Called after market data has been retrieved */
void CreditIndexSwap::Validate() {
    static const string method = "CreditIndexSwap::validate";
    if (discount->getCcy() != creditIndex->getCcy()) {
        throw ModelException(method,
                             "In the current implementation the " 
                             "currency of the CIS instrument (" + 
                             discount->getCcy() + 
                             ") must be the same as the currency of " 
                             "the index being swapped (" +
                             creditIndex->getCcy() + ").");
    }
}


/* Notification that (some) underlying fields have changed */
void CreditIndexSwap::fieldsUpdated(const CFieldArray& fields) {
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


/** Returns the value date (aka today) the instrument is currently
    pricing for */
DateTime CreditIndexSwap::getValueDate() const {
    return valueDate;
}


/** Returns the name of the instrument's discount currency. */
string CreditIndexSwap::discountYieldCurveName() const {
    return discount.getName();
}


/** When to stop tweaking */
DateTime CreditIndexSwap::endDate(const Sensitivity* sensControl) const {
    // Get last fee payment date
    const DateTime& feeLastDate  = feeLeg->getLastPayDate();
    const DateTime& sensLastDate = creditIndex->stopTweaking(feeLastDate);

    if (!sensControl) {
        return sensLastDate;
    }
    // JLH - Just like in CredDefSwaps, this is a mess:
    // Need to route endDate through model not through instrument
    IModel* model = sensControl->getModel();
    ClosedFormCDSBasket* cdsModel = dynamic_cast<ClosedFormCDSBasket*>(model);
    return cdsModel ? cdsModel->endDate(sensControl, this, sensLastDate): 
                      sensLastDate;
}


void CreditIndexSwap::price(CResults* results, Control* control, IForwardRatePricerSP model) const {
    static const string method = "CreditIndexSwap::price";

    try {
        double price = 0.0;
        
        int numNamesInPortfolio = adjCurvesCache->size();
        // Price each CDS in adjCurvesCache and add their prices up
        for (int i=0; i < numNamesInPortfolio; ++i) {
            double thisSwapPrice = (*adjCurvesCache)[i]->getPrice(this, model);
            price += thisSwapPrice;
        }
        
        results->storePrice(price, discount->getCcy());

        // If pricing, produce the output requests too (ie, not if we
        // are tweaking)
        if (control && control->isPricing()) {
            addOutputRequests(control, results, model);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


//------------------------------------------------------------------------------
//                    Methods related to output requests
//------------------------------------------------------------------------------

/** Deals with the CIS output requests */
void CreditIndexSwap::addOutputRequests(Control* control, 
                                        Results* results,
                                        IForwardRatePricerSP model) const 
{
    static const string method = "CreditIndexSwap::addOutputRequests";

    try {
        OutputRequest* request;

        // INDEX_BASIS
        request = control->requestsOutput(OutputRequest::INDEX_BASIS);
        if (request) {
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

        // PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(paymentDates(model));
            OutputRequestUtil::recordPaymentDates(control, results, dates.get());
        }

        // KNOWN_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArraySP cfl(knownCashFlows(model));
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get());
        }
     }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}


/** When do payments occur */
DateTimeArraySP CreditIndexSwap::paymentDates(IForwardRatePricerSP model) const {
    AbstractCashFlowArrayConstSP cashFlows = feeLeg->getCashFlows(model);
    int numberOfCashFlows = cashFlows->size();
    DateTimeArraySP paydates(new DateTimeArray(numberOfCashFlows));

    for (int i=0; i < numberOfCashFlows; ++i) {
        (*paydates)[i] = (*cashFlows)[i]->getPayDate();
    }
    return paydates;
}


/** Returns all known cash flows.
 * Caution: Due to the use of CashFlow::merge, this may be slow.
 * Since it is only called once per pricing (not for greeks), this has 
 * not been optimized - if needed, it could follow a similar logic to
 * DateTime::merge, which takes a vector of arrays and merges them
 * in one call */
CashFlowArraySP CreditIndexSwap::knownCashFlows(IForwardRatePricerSP model) const {
    CashFlowArraySP cashFlow;
    int numNamesInPortfolio = adjCurvesCache->size();

    // Obtain the known cashflows of the remaining names, and merge them
    CashFlowArraySP nextNamesCashFlows;
    for (int i=0; i < numNamesInPortfolio; ++i) {
        nextNamesCashFlows = (*adjCurvesCache)[i]->knownCashFlows(this, model);
        cashFlow = CashFlow::merge(cashFlow, nextNamesCashFlows);
    }
    return cashFlow;
}


/** Implementation of Theta::Shift */
bool CreditIndexSwap::sensShift(Theta* shift) {
    valueDate = shift->rollDate(valueDate);
    return true; // Shift components also
}


void CreditIndexSwap::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditIndexSwap, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ClosedFormCDSBasket::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultCreditIndexSwap);

    FIELD(protectionStartDate, "Protection and accrual start date");
    FIELD(protectionEndDate,   "Protection end date");
    FIELD(payAccruedFee,       "Make accrued payment on default?");
    FIELD(dcc,                 "Day count conv for accrual periods");
    FIELD(notional,            "Original CIS notional before any defaults. "
                               "Notional > 0 means long protection.");
    FIELD(discount,            "Discount curve");
    FIELD(creditIndex,         "Index to swap");
    FIELD(feeLeg,              "Fee leg of the swap");
    FIELD(adjCurvesCache,      "Transient field with the cache for the "
                               "curves' prices.");
    FIELD(swapDcc,             "Day Count Convention");
    FIELD(creditEventOverride, "Override for credit event parameters");
    FIELD(valueDate,           "Valuation Date");
    FIELD(lastTriggerDate,     "Last date when a default occurred during "
                               "the protection period can be triggered.");
    FIELD(triggerDelay,        "Delay, in days, between default and "
                               "eventDeterminationDate - used to "
                               "estimate eventDeterminationDate "
                               "before it is known.");
    FIELD(defaultToSettlementDelay,"Delay between credit event and "
                                   "settlement.");

    FIELD_MAKE_OPTIONAL(creditEventOverride);
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD_MAKE_OPTIONAL(lastTriggerDate);
    FIELD_MAKE_OPTIONAL(triggerDelay);
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelay);
    
    FIELD_MAKE_TRANSIENT(swapDcc);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(adjCurvesCache);
    //FIELD_MAKE_NONTWEAKABLE(creditIndex); // We are pricing as a basket of CDS
                                          // so do not tweak the creditIndex. 
                                          // The names' CDSs are refered to from 
                                          // adjCurvesCache and tweaked there

    // Verify the name adjCurvesCache field - putting this check in this 
    // static method will ensure we catch issues as soon as the app is loaded.
    if (clazz->hasDeclaredField(adjCurvesCacheFieldName) == CFieldConstSP()) {
        throw ModelException("CreditIndexSwap::load",
                             "Internal error: field " + 
                             adjCurvesCacheFieldName +
                             "has changed name. Need to update "
                             "adjCurvesCacheFieldName.");
    }
}


IObject* CreditIndexSwap::defaultCreditIndexSwap() {
    return new CreditIndexSwap();
}


/*****************************************************************************
 *                   Model - instrument interaction
 *****************************************************************************/

/** Private class - Handles the interaction between model and instrument */
class CreditIndexSwapClosedFormCDSBasket: public ClosedFormCDSBasket::IProduct {
private:
    const CreditIndexSwap* cf;  // A reference to the actual CIS instrument

public:
    CreditIndexSwapClosedFormCDSBasket(const CreditIndexSwap* cf): cf(cf)
    {}

    void price(ClosedFormCDSBasket* model,
               Control*         control, 
               CResults*        results) const
    {
        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp)
        {
            throw ModelException("CreditIndexSwapClosedFormCDSBasket::price",
                "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

        cf->price(results, control, frModel);
    }
};


/** Implementation of ClosedFormCDSBasket::IntoProduct interface */
ClosedFormCDSBasket::IProduct* CreditIndexSwap::createProduct(
    ClosedFormCDSBasket* model) const 
{
    return new CreditIndexSwapClosedFormCDSBasket(this);
}


CClassConstSP const CreditIndexSwap::TYPE = CClass::registerClassLoadMethod(
    "CreditIndexSwap", typeid(CreditIndexSwap), load);

// Store the name of the adjCurvesCache field
const string CreditIndexSwap::adjCurvesCacheFieldName = "adjCurvesCache";


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool CreditIndexSwapLoad() {
    return (CreditIndexSwap::TYPE != 0);
}


DRLIB_END_NAMESPACE
