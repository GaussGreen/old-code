#ifndef EDR_MCPRODUCT_HPP
#define EDR_MCPRODUCT_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/MCPathGenerator.hpp"
#include "edginc/PastValues.hpp"
#include "edginc/MCPrices.hpp" // Why is MCPricesSP being passed to handlePast?
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

//class IMCPrices;
//class IMCPricesSP;
class IMultiFactors;
class IMultiMarketFactors;
class SensControl;
class Control;
class EventResults;
class MonteCarlo;
class IMCQuickGreeks;
class IMCQuickXGamma;

class MCARLO_DLL IMCProduct {
public:
    /** These ints are used in a bitwise manner in
        createOrigPrices to indicate what quick greek
        calculations are needed */
    static const int QUICK_GREEKS_BIT;
    static const int QUICK_X_GAMMA_BIT;
    static const int CACHE_PRODUCT_BIT;
    // friend class MonteCarlo;

    /** Invoked once the 'past' has been completed (invoked once only
        per pricing). Default implementation does nothing */
    virtual void       finaliseSimulationDates();

    /** Create IMCPrices object for first pricing call. The mode
        parameter indicates which 'quick greek' calculations are
        required ie it indicates whether the additional information
        typically needed by quick greeks should be stored or not. The
        parameter is a bit wise integer using.  QUICK_X_GAMMA_BIT and
        QUICK_GREEKS_BIT */
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode);

    /** What derived classes must implement. Will be called possibly once for
        the past (if there is any past dates) and then again for each
        simulation. Note must support all dates being historic */
    virtual void       payoff(const MCPathGenerator* pathGen,
                              IMCPrices&               prices) = 0;

    /* Default implementation of this just does
       "payment in future"? discount->pv(Today, paymentDate): 0.0.
       Implementations that override this default must take into account
       whether or not the cashflow is historic. They should also override
       endDate() */
    virtual double     pvFromPaymentDate() const;

    // methods to allow easy interaction with path generators

    //// Returns the multiFactors - if market data is consistent with
    //// a IMultiFactors
    const IMultiFactors* getMultiFactors() const;

    //// Returns the MultiMarketFactors object
    const IMultiMarketFactors* getMultiMarketFactors() const;

    //// Returns the num of assets (in the payoff view of the world)
    int getNumAssets() const;

    //// Returns the value date
    const DateTime& getToday() const;

    //// Returns the SimSeries
    const SimSeries* getSimSeries() const;

    /** Returns true if there are future simulation dates.
        Useful in context of the payoff call for any past dates */
    virtual bool hasFuture() const {
        return getSimSeries()->getLastDate().isGreater(getToday());
    };

    //// Returns the RefLevel
    const IRefLevel* getRefLevel() const;

    //// Returns the discount YieldCurve
    const YieldCurve* getDiscount() const;

    //// Returns the MC PastValues (ie historic samples)
    const IPastValues* getMCPastValues() const;

    //// Returns MAX(simulation start date, today)
    const DateTime& getEffectiveSimStartDate() const;

    /** the date to stop tweaking for the supplied tweak. The
        default implementation uses the last simulation date and the
        payment date. Derived classes overriding pvFromPaymentDate() need to
        override this method too */
    virtual DateTime endDate(const Sensitivity*  sensControl) const;

    /** invoked after final simulated path is run. Default does nothing.
        Allows derived classes to store debug information for example */
    virtual void recordExtraOutput(Control*      control,
                                   Results*      results,
                                   const IMCPrices& prices) const;

    // collects events which should be sitting there
    // called after past is run
    virtual void retrieveEvents(EventResults* events) const;

    // go get the events based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen,
                           EventResults* events,
                           const DateTime& eventDate);

     /** Returns an array of 'forward start dates' (ie cliquet start dates)
        which should be used for volatilty interpolation. The default method
        here uses the IMCProductLN interface if the product implements it
        otherwise it fails (so this method would need to be overridden) */
    virtual DateTimeArray getVolStartDates(
        const MCPathGenerator* pathGenerator) const;

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    void setSimStartDateToFirstRefDate();

    virtual ~IMCProduct();

    /** This method is called every time the path generator is changed
        (which is, at the moment, when the past path generator is
        created, and then when the future path generator is
        created. Default does nothing  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    /** return all state variable interfaces */
    virtual void getSvSet(vector<IStateVariableSP>&) {};

    // runs the past without storing results at all - just wraps a call to handlePast()
    // may need to keep hold of path generator until finished with to avoid reading freed memeory
    MCPathGeneratorSP runPast(MonteCarlo* mcarlo);

    // runs the MC simulation
    virtual MCPathGeneratorSP price(MonteCarlo*       mcarlo,
				                    Control*          control,
                                    Results*          results);

    virtual void validate();
    /** Returns null if not supported. Implemented via dynamic cast. Allows
        per instrument choice */
    virtual IMCQuickGreeks* quickGreeksSupported();
    /** Returns null if not supported. Implemented via dynamic cast. Allows
        per instrument choice */
    virtual IMCQuickXGamma* quickXGammaSupported();

    virtual void finalize() {}

protected:

    /** 'Full' constructor for single factor payoffs */
    IMCProduct(const CAsset*               asset,        // single factor
              const DateTime&             today,        // value date
              const YieldCurve*           discount,     // for pv'ing
              const IRefLevelConstSP&     refLevel,     // how to 'avg in'
              const SimSeriesConstSP&     simSeries,    // sim dates
              const IPastValuesConstSP&   mcPastValues, // historic values
              const InstrumentSettlement* instSettle,   // used for pay date
              const DateTime&             matDate);     // used for pay date

    /** 'Full' constructor for multi-factor payoffs */
    IMCProduct(const IMultiMarketFactors*  mFactors, // typically a MultiAsset
              const DateTime&             today,    // value date
              const YieldCurve*           discount, // for pv'ing
              const IRefLevelConstSP&     refLevel,     // how to 'avg in'
              const SimSeriesConstSP&     simSeries,    // sim dates
              const IPastValuesConstSP&   mcPastValues, // historic values
              const InstrumentSettlement* instSettle,   // used for pay date
              const DateTime&             matDate);     // used for pay date

    /** Simpathl constructor */
    IMCProduct(const IMultiMarketFactors*  mFactors, // typically a MultiAsset
              const DateTime&             today,    // value date
              const YieldCurve*           discount);

    /** lean ctor */
    IMCProduct(const DateTime&             today,    // value date
               const YieldCurve*           discount,
               const SimSeriesConstSP&     simSeries);

    //// sets the payment date to given value
    void setPaymentDate(const DateTime& payDate);

    /** sets the payment date calculated from supplied InstrumentSettlementSP
        and maturity date */
    void setPaymentDate(const InstrumentSettlement* instSettle,
                        const DateTime&             matDate);

    /** sets the RefLevel obejct (does not take copy) */
    void setRefLevel(const IRefLevelConstSP& refLevel);

    /** set the SimSeries object (does not take copy) */
    void setSimSeries(const SimSeriesConstSP& simSeries);

    /** set the PastValues object (does not take copy) */
    void setMCPastValues(const IPastValuesConstSP& mcPastValues);

    /** Returns true upon construction of IMCProduct and is switched to false
        once the future path generator has been created. */
    bool doingPast() const{ return doingThePast;} // in-line for perf

    // default implementation to record events (pay dates, cashflows etc)
    // will try and use any product's implementation first
    // Protected & virtual to allow MCProductClient to provide alternative
    // default implementation
    virtual void recordEvents(Control* control,
                              IMCPrices&  prices,
                              Results* results);

    // handles the past
    MCPathGeneratorSP handlePast(MonteCarlo*       mcarlo,
                                   Control*        control,
                                   IMCPricesSP&       prices,
                                   bool            needPrices);

    void postResults(CControl*     control,
                     IMCPrices&       prices,
                     CResults*     results);

    //// turns Asset into IMultiFactors
    static IMultiFactorsConstSP convertToMultiFactors(const Asset* asset);

    ///  Protected  Fields //////////////
	bool                    doingThePast;
    DateTime                Today;      // value date
    int                     NbAssets;   // number of assets (not factors)
    const YieldCurve*       discount;   // discount curve

    InstrumentSettlementConstSP settlement;     // instruments settlement

    /// fields - accessed via methods for flexibility
    mutable IMultiFactorsConstSP mAsset;  /* This is the equity view of assets
                                             - may be null */
    IMultiMarketFactorsConstSP mfAsset;// This is the IMCProduct's view of assets
    IRefLevelConstSP     refLevel;    // 'averaging' in
    SimSeriesConstSP     simSeries;   // dates to simulate on
    IPastValuesConstSP   mcPastValues; // holds historic samples
    DateTime             paymentDate;  // Used for PV
};
DECLARE_REF_COUNT(IMCProduct);
DRLIB_END_NAMESPACE

#endif // EDR_MCPRODUCT_HPP
