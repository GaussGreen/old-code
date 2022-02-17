//----------------------------------------------------------------------------
//
//   Filename    : CDS.hpp
//
//   Description : General Credit Default Swap
//
//----------------------------------------------------------------------------

#ifndef QLIB_CDS_HPP
#define QLIB_CDS_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormCDSBasket.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICreditEventOverride)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)

/** General Credit Default Swap instrument */

class PRODUCTS_DLL CDS: public CInstrument,
                        virtual public ICDS,
                        virtual public IBadDayAdjuster,
                        virtual public Theta::Shift, 
                        virtual public ClosedFormCDSPS::IIntoProduct,    //single name swap
                        virtual public ClosedFormCDSBasket::IIntoProduct //index swap
{
public:
    class CDSOutputs;

    //------------
    // CDS methods
    //------------
    static CClassConstSP const TYPE;

    virtual ~CDS();

    /** Main pricing methods */
    void priceSingleName(CResults* results, Control* control, IForwardRatePricerSP model) const;
    void priceMultiName(CResults* results, Control* control, IForwardRatePricerSP model) const;

    //Does the instrument pay accrued on default
    bool paysAccrued() const;

    //--------------------
    // CInstrument methods
    //--------------------

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Allow the instrument to retrieve its market data */
    virtual void GetMarket(const IModel*            model,
	                       const CMarketDataSP      market);

    /** called after market data has been retrieved */
    virtual void Validate();

    /** Notification that (some) underlying fields have changed */
    void fieldsUpdated(const CFieldArray& fields);

    /** Returns the value date (aka today) the instrument is currently pricing for */
    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency */
    virtual string discountYieldCurveName() const;

    //--------------
    // ICDS methods
    //--------------

    /** Return the fee leg */
    virtual ICreditFeeLegSP getFeeLeg() const;

    /** Return the contingent leg */
    virtual ICreditContingentLegSP getContingentLeg() const;

    /** Return Accrued Interest (cash amount, not percentage) for settlement 
        on settlementDate */
    virtual double getAccruedInterest(const DateTime&      settlementDate,
                                      IForwardRatePricerSP model) const;

    /**
     * Returns the present value of the instrument at valuationDate,
     * conditional on default at defaultDate.
     */
    virtual double getPVGivenDefault(const DateTime&            valuationDate, 
                                     const DateTime&            defaultDate,
				                     const IDiscountCurveRisky& crv,
                                     IForwardRatePricerSP       model) const;

    /** Determines if the instrument is perpetual. */
    virtual bool hasFiniteMaturity() const;

    /** Returns the maturity of the instrument */
    virtual DateTime getMaturity() const;

    /**
     * Returns the earliest of the first accrual period, start of the contingent leg,
     * the first cash flow, etc.
     */
    virtual DateTime getStartDate() const;

    /**
     * Get the base (initial) notional of the instrument.
     * Positive notional = long risk/short protection.
     * If the notional is varying, this will return the base notional.
     */
    virtual double getNotional() const;

    /**
     * Set the base (initial) notional of the instrument.
     * Positive notional = short risk/long protection.
     */
    virtual void setNotional(double newNotional);

    /** Returns the day count convention used for accruals */
    virtual DayCountConventionSP getAccrualDcc() const;

    virtual CashFlowArraySP getInstrumentCashFlows(IForwardRatePricerSP model) const;

    /** Returns the wrapper to the yield curve */
    virtual YieldCurveWrapper getYieldCurveWrapper() const;

    /** Returns the wrapper describing the credit curve */
    virtual ICDSParSpreadsWrapper getParSpreadsWrapper() const;


    /**
     *  Compute the PV of a vanilla cds with the given legs given no default before
     *  valuation date.
     *  Overrides ICDS implementation
     */
    virtual double getPV(const DateTime&              valuationDate,
                         const IDiscountCurveRisky&   crv,
                         const IDecretionCurveConstSP prepay,
                         IForwardRatePricerSP         model,
                         IBadDayAdjusterConstSP       bda,
                         CDSOutputs&                  myOutput) const;

    /**
     * Compute the PV of a vanilla cds at valuationDate for unconditional settlement and
     * payment on paymentDate.  Conditional on no new defaults before valuationDate if
     * valuationDate is in the future.
     * Overrides ICDS implementation
     */
    virtual double getPV(const DateTime&              valuationDate,
                         const DateTime&              settlementDate,
                         const IDiscountCurveRisky&   crv,
                         const IDecretionCurveConstSP prepay,
                         IForwardRatePricerSP         model,
                         IBadDayAdjusterConstSP       bda,
                         CDSOutputs&                  myOutput) const;

    void addOutputRequests(Control*               control,
                           CResults*              results,
                           IForwardRatePricerSP   model,
                           bool                   isMultiName,
                           const CDSOutputs&      myOutput) const;

    //------------------------
    // IBadDayAdjuster methods
    //------------------------

    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;

    //--------------------
    //Theta::Shift methods
    //--------------------

    virtual bool sensShift(Theta* shift);

    //--------------------------------------
    // ClosedFormCDSPS::IIntoProduct methods
    //--------------------------------------

    /** Implementation of ClosedFormCDSPS::IntoProduct interface */
    virtual ClosedFormCDSPS::IProduct* createProduct(ClosedFormCDSPS* model) const;

    //------------------------------------------
    // ClosedFormCDSBasket::IIntoProduct methods
    //------------------------------------------

    /** Implementation of ClosedFormCDSBasket::IntoProduct interface */
    virtual ClosedFormCDSBasket::IProduct* createProduct(ClosedFormCDSBasket* model) const;

    // PriceCache should be private (to avoid other classes using it) but
    // PriceCacheArray needs to be public to be able to register the type.
    class PriceCache;
    DECLARE(PriceCache);

private:
    CDS();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //Multi name known cashflows
    CashFlowArraySP knownCashFlows(IForwardRatePricerSP model) const;
    //Single name known cashflows scaled by prepayment
    CashFlowArraySP knownCashFlowsWithPrepayment(IForwardRatePricerSP model) const;

    // Used to store the name of the adjCurvesCache field (in case it changes)
    static const string adjCurvesCacheFieldName;

    //-------
    // Fields
    //-------
    DateTime                       valueDate;
    bool                           payAccruedFee;
    string                         accrualDayCountConvention;

    ICreditFeeLegSP                feeLeg;
    ICreditContingentLegSP         contingentLeg;

    ICDSParSpreadsWrapper          underlying;
    YieldCurveWrapper              discount;

    mutable ICreditEventOverrideSP creditEventOverride;        // Override for credit event parameters
    CIntSP                         triggerDelay;               // Delay between default and eventDeterminationDate, in days
    CIntSP                         defaultToSettlementDelay;   // Delay between credit event and settlement
    DateTime                       lastTriggerDate;            // Last date when a default occurred during 
                                                           // the protection period can be triggered
    HolidayWrapper                 settlementHols;             // Holidays, used to settle payments
    string                         settlementBadDayConvention; // used in conjunction with holidays 

    DayCountConventionSP           accrualDCC;                 // used if payAccruedFee = true 
    BadDayConventionSP             settlementBDC;              // transient; created from settlementBadDayConvention

    mutable PriceCacheArraySP      adjCurvesCache;             // cache for curves prices used when pricing as a basket
};

DECLARE(CDS)

DRLIB_END_NAMESPACE
#endif

