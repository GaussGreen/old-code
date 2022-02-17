#ifndef EDR_MCPRODUCTCLIENT_HPP
#define EDR_MCPRODUCTCLIENT_HPP

#include "edginc/config.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/MCProduct.hpp"

DRLIB_BEGIN_NAMESPACE


/** Products that are users of state variables derive from this. */ 
class MCARLO_DLL MCProductClient: public IMCProduct,
                       virtual public IStateVariableClient {
public:
    /** Returns the SVGenSpot-type PastValues (ie historic samples) */
    const IPastValues* getSVGenSpotPastValues() const;

    /** Returns the MCSqrtAnnualQuadVar-type PastValues (ie historic samples) */
    const IPastValues* getMCSqrtAnnualQuadVarPastValues() const;

    /** Overrides method assuming that payoffs will be doing
        the discounting on their own */
    virtual double pvFromPaymentDate() const;

    virtual void recordEvents(Control* control,
                              IMCPrices&  prices,
                              Results* results);

    // Helper class for the output request KNOWN_CASHFLOWS
    // Could have done this via an interface, but we need to support KCFs
    // and there is a well-defined minimum of functionality.
    // It's not clear to me if or how derived classes may wish to extend this,
    // so for now I leave a naive implementation. The SPI has a more general
    // form but that will be forced to change in an SV-world (since pv-ing will
    // be done inside the simulation).
    class MCARLO_DLL KnownCashFlowsHelper {
    public:
        KnownCashFlowsHelper();
        
        void addFlow(const DateTime& when,
                     double          howMuch);
        
        const CashFlowArray* getFlows();
        
    private:
        CashFlowArraySP   knownFlows;
    };

    /** Returns true if there are future simulation dates. 
        Overridden default for SV since then the IR dates 
        are relevant for the decision. Should possibly examine
        all the SVs but for now this should suffice. */
    virtual bool hasFuture() const {
        return getSimSeries()->getLastDate().isGreater(getToday());
// Something like this needed, but not clear quite how to make it work.
//         if (paymentDate.empty()) {
//             throw ModelException("MCProductClient::hasFuture",
//                                  "paymentDate not yet set!");
//         }
//         return paymentDate.isGreater(getToday());
    };

protected:
    /** 'Full' constructor for single factor payoffs */
    MCProductClient(const CAsset*               asset,        // single factor
                    const DateTime&             today,        // value date
                    const YieldCurve*           discount,     // for pv'ing
                    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
                    const SimSeriesConstSP&     simSeries,    // sim dates
                    const IPastValuesConstSP&   mcPastValues, // historic values
                    const InstrumentSettlement* instSettle,   // used for pay date
                    const DateTime&             matDate);     // used for pay date

    /** 'Full' constructor for single factor payoffs, including
        MCSqrtAnnualQuadVar-type past values */
    MCProductClient(const CAsset*               asset,        // single factor
                    const DateTime&             today,        // value date
                    const YieldCurve*           discount,     // for pv'ing
                    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
                    const SimSeriesConstSP&     simSeries,    // sim dates
                    const IPastValuesConstSP&   mcPastValues,                  // SVGenSpot-type historic values
                    const IPastValuesConstSP&   mcSqrtAnnualQuadVarPastValues, // MCSqrtAnnualQuadVar-type historic values
                    const InstrumentSettlement* instSettle,   // used for pay date
                    const DateTime&             matDate);     // used for pay date

    /** 'Full' constructor for multi-factor payoffs */
    MCProductClient(const IMultiMarketFactors*  mFactors, // typically a MultiAsset
                    const DateTime&             today,    // value date
                    const YieldCurve*           discount, // for pv'ing
                    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
                    const SimSeriesConstSP&     simSeries,    // sim dates
                    const IPastValuesConstSP&   mcPastValues, // historic values
                    const InstrumentSettlement* instSettle,   // used for pay date
                    const DateTime&             matDate);     // used for pay date

    /** 'Full' constructor for multi-factor payoffs, including
        MCSqrtAnnualQuadVar-type past values */
    MCProductClient(const IMultiMarketFactors*  mFactors, // typically a MultiAsset
                    const DateTime&             today,    // value date
                    const YieldCurve*           discount, // for pv'ing
                    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
                    const SimSeriesConstSP&     simSeries,    // sim dates
                    const IPastValuesConstSP&   mcSpotPastValues,              // SVGenSpot-type historic values
                    const IPastValuesConstSP&   mcSqrtAnnualQuadVarPastValues, // MCSqrtAnnualQuadVar-type historic values
                    const InstrumentSettlement* instSettle,   // used for pay date
                    const DateTime&             matDate);     // used for pay date

    /** Simpathl constructor */
    MCProductClient(const IMultiMarketFactors* mFactors,
                    const DateTime& today,
                    const YieldCurve* discount);

    /** Lean ctor */
    MCProductClient(const DateTime& today,
                    const YieldCurve* discount,
                    const SimSeriesConstSP&     simSeries);

    virtual ~MCProductClient();

    virtual void validate();

    // this is how the SV-compliant products support KCF
    // Protected scope so derived classes can use it freely.
    KnownCashFlowsHelper* knownCashFlows; 

private:
    IPastValuesConstSP mcSqrtAnnualQuadVarPastValues; // holds historic samples
};


/** Products that work with both path and slice based Pricing  
    must implement the IMCStatelessProductClient */
class MCARLO_DLL IMCStatelessProductClient : public virtual IHistoricalContextGenerator {
public:
    /** Product should return the list of indices in the timeline it is 
    interested in.  The statelessPayoff function will only be called
    on the marked dates. */
    virtual vector<int> finalize( 
        const DateTimeArray& timeline ) = 0;

    /** The payoff function should use the current date index to determine 
        what action it should take.  Typically, it will update the 
        historical context on each day until expiry at which time
        it will insert the PV of the product into the prices structure.
        The currentDateIdx is with respect to the instrument's perspective
        (as specified by the finalize method) rather than the overall 
        timeline's perspective. */
    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices ) = 0;

    /* return array of past dates, used for stateless past Pricing */
    virtual DateTimeArray getPastDates() = 0; // xxx should be made pure virtual

};

DRLIB_END_NAMESPACE

#endif // EDR_MCPRODUCTCLIENT_HPP
