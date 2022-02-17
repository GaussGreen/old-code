#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/MCProductClient.hpp"
#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

MCProductClient::KnownCashFlowsHelper::KnownCashFlowsHelper(): 
    knownFlows(new CashFlowArray(0)){};
 
void MCProductClient::KnownCashFlowsHelper::addFlow(const DateTime& when,
                                                    double          howMuch) {
    // take care to maintain date order
    CashFlow cf(when, howMuch);
    CashFlowArraySP cfa(new CashFlowArray(0));
    cfa->push_back(cf);
    // merge needs 2 existing CFAs - so start carefully
    if (knownFlows.get()) {
        knownFlows = CashFlow::merge(cfa, knownFlows);
    } else {
        knownFlows = cfa;
    }
};
    
const CashFlowArray* MCProductClient::KnownCashFlowsHelper::getFlows() {
    return knownFlows.get();
};

/** 'Full' constructor for single factor payoffs */
MCProductClient::MCProductClient(
    const CAsset*               asset,
    const DateTime&             today,
    const YieldCurve*           discount,
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcPastValues, 
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
IMCProduct(asset, today, discount, refLevel, simSeries, mcPastValues, instSettle, matDate),
knownCashFlows(new KnownCashFlowsHelper()){}

/** 'Full' constructor for single factor payoffs, including
    MCSqrtAnnualQuadVar-type past values */
MCProductClient::MCProductClient(
    const CAsset*               asset,        // single factor
    const DateTime&             today,        // value date
    const YieldCurve*           discount,     // for pv'ing
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcPastValues,                  // SVGenSpot-type historic values
    const IPastValuesConstSP&   mcSqrtAnnualQuadVarPastValues, // MCSqrtAnnualQuadVar-type historic values
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
IMCProduct(asset, today, discount, refLevel, simSeries, mcPastValues, instSettle, matDate),
knownCashFlows(new KnownCashFlowsHelper()),
mcSqrtAnnualQuadVarPastValues(mcSqrtAnnualQuadVarPastValues){}    

/** 'Full' constructor for multi-factor payoffs */
MCProductClient::MCProductClient(
    const IMultiMarketFactors*  mFactors,
    const DateTime&             today,
    const YieldCurve*           discount,
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcPastValues, 
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
IMCProduct(mFactors, today, discount, refLevel, simSeries, mcPastValues, instSettle, matDate),
knownCashFlows(new KnownCashFlowsHelper()){}    

/** 'Full' constructor for multi-factor payoffs, including
    MCSqrtAnnualQuadVar-type past values*/
MCProductClient::MCProductClient(
    const IMultiMarketFactors*  mFactors, // typically a MultiAsset
    const DateTime&             today,    // value date
    const YieldCurve*           discount, // for pv'ing
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcSpotPastValues,              // SVGenSpot-type historic values
    const IPastValuesConstSP&   mcSqrtAnnualQuadVarPastValues, // MCSqrtAnnualQuadVar-type historic values
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
IMCProduct(mFactors, today, discount, refLevel, simSeries, mcPastValues, instSettle, matDate),
knownCashFlows(new KnownCashFlowsHelper()),
mcSqrtAnnualQuadVarPastValues(mcSqrtAnnualQuadVarPastValues){}    

MCProductClient::MCProductClient(
    const IMultiMarketFactors*  mFactors, // typically a MultiAsset
    const DateTime&             today,    // value date
    const YieldCurve*           discount): // for pv'ing
IMCProduct(mFactors, today, discount), 
knownCashFlows(new KnownCashFlowsHelper()) {};

MCProductClient::MCProductClient(
    const DateTime&             today,    // value date
    const YieldCurve*           discount, // for pv'ing
    const SimSeriesConstSP&     simSeries):
IMCProduct(today, discount, simSeries), 
knownCashFlows(new KnownCashFlowsHelper()) {};

MCProductClient::~MCProductClient() {
    delete knownCashFlows;
}

void MCProductClient::validate() {
    const static string routine("MCProductClient::validate");
    // some idiot checking
    if (!simSeries){
        throw ModelException(routine, "NULL sim series");
    }
    /* validate numAssets is consistent */
    if (NbAssets != simSeries->getNumAssets()){
        throw ModelException(routine, "MultiAsset has "+
                             Format::toString(NbAssets)+
                             " but Sim Series has "+
                             Format::toString(simSeries->getNumAssets()));
    }
    if (mcPastValues.get() && NbAssets != mcPastValues->getNumAssets()){
        throw ModelException(routine, "MultiAsset has "+
                             Format::toString(NbAssets)+" assets"
                             " but Past Values has data for"+
                             Format::toString(mcPastValues->getNumAssets())+
                             " assets");
    }
}        

/** Returns the SVGenSpot-type PastValues (ie historic samples) */
const IPastValues* MCProductClient::getSVGenSpotPastValues() const{
    return IMCProduct::getMCPastValues();
}

/** Returns the MCSqrtAnnualQuadVar-type PastValues (ie historic samples) */
const IPastValues* MCProductClient::getMCSqrtAnnualQuadVarPastValues() const{
    if (!mcSqrtAnnualQuadVarPastValues){
        throw ModelException("MCProductClient::getMCSqrtAnnualQuadVarPastValues",
                             "internal error: mcSqrtAnnualQuadVarPastValues is null");
    }
    return mcSqrtAnnualQuadVarPastValues.get();
}

double MCProductClient::pvFromPaymentDate() const {
    if (paymentDate.empty()) {
        throw ModelException("MCProductClient::pvFromPaymentDate",
                             "payment date has not been set");
    }
    if (paymentDate.isGreater(Today)){
        return 1.0;
    }
    return 0.0; // cashflow is in the past
}

// This overrides the default provided by IMCProduct
void MCProductClient::recordEvents(Control* control,
                                   IMCPrices&  prices,
                                   Results* results) {
    static const string method("MCProductClient::recordEvents");
    try {
        // see if product handles this first
        IHandlePaymentEvents* handler =
            dynamic_cast<IHandlePaymentEvents*>(this);
  
        if (handler) {
            // this method/interface may need work
            handler->recordEvents(control, results);
        }
        else {
            // built-in implementation for single cashflow payoffs
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                if (paymentDate.empty()) {
                    throw ModelException(method, 
                                         "payment date has not been set");
                }

                DateTimeArray date(1, paymentDate);
                OutputRequestUtil::recordPaymentDates(control,results,&date);
            }

            // This needs to be specialised for SV-compliant products
            // since the prices class contains already PVd quantities.
            // As far as I can see, the payoff needs to do some work
            // so I've provided a helper class. This allows a more general
            // solution since the product may place multiple flows in
            // this helper class, whereas the non-SV default implementation
            // assumed a single flow at maturity.
            request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished() &&
                !simSeries->getLastDate().isGreater(Today) &&
                !settlement->isPhysical()) {

                const CashFlowArray* kcf = knownCashFlows->getFlows();
                if (kcf) {
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            kcf); 
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
