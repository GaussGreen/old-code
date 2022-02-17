//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CredDefSwap.hpp
//
//   Description : Credit default swap
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 8, 2001
//
//
//----------------------------------------------------------------------------

#ifndef CREDDEFSWAP_HPP
#define CREDDEFSWAP_HPP
#include "edginc/Instrument.hpp"
#include "edginc/Generic1FactorCredit.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/ClosedFormCDSPSandFA.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CreditSupport.hpp"
#include "edginc/ICDS.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(IDiscountCurveRisky)
FORWARD_DECLARE_WRAPPER(ICreditEventOverrideName);

// For defining the gradient of par spreads with respect to clean spreads
typedef struct _TEqSpreadDataResult
{
    CDoubleMatrixSP  dSi_dlj;
    CDoubleArraySP   rRate;      // interest rate ( r )
    CDoubleArraySP   lRate;      // default rate ( lambda )
    CDoubleArraySP   Si;         // par spreads
    CDoubleArraySP   Smax;       // max par spreads
}TEqSpreadDataResult;

/** Credit default swap instrument - fixed amounts at known future dates */

class PRODUCTS_DLL CredDefSwap:
    public Generic1FactorCredit,
    virtual public ClosedFormCDSPS::IIntoProduct,
    virtual public ClosedFormFA::IIntoProduct,
    virtual public ClosedFormCDSPSandFA::IIntoProduct,
    virtual public LastSensDate,
    virtual public DeltaToCredit::IShift,
    virtual public ObjectIteration::IOverride,
    virtual public IGetMarket,
    virtual public CreditSupport::Interface,
    virtual public ICDS,
    virtual public ICDSConvention,
    virtual public ICreditFeeLeg,          // we allow the CDS to behave as a fee due to its structure
    virtual public IFixedRateCreditFeeLeg, // and we make it a fixed rate variety
    virtual public ICreditContingentLeg,   // similarly allow ctg leg behaviour
    virtual public IBadDayAdjuster
{
public:
    static CClassConstSP const TYPE;
    friend class CredDefSwapHelper;
    friend class CredDefSwapClosedFormCDSPS;
    friend class CredDefSwapClosedFormFA;
    friend class CredDefSwapClosedFormCDSPSandFA;
    friend class QPCDSParSpreads;
    friend class QPCDSFirmAsset;
    friend class QPCDSPSandFA;
    friend class CDSCreditSupport;

    static const double CDS_ADJUSTMENT;

    //--------------------
    // CredDefSwap methods
    //--------------------
    virtual ~CredDefSwap();

    virtual void validatePop2Object();

    virtual void Validate();

    //--------------
    // ICDS methods
    //--------------

    /** Return the fee leg */
    virtual ICreditFeeLegSP getFeeLeg() const;

    /** Return the contingent leg */
    virtual ICreditContingentLegSP getContingentLeg() const;

    /*=========================================================================
     * Methods required by ICDSWithConvention
     * These are dreadfully hacked, and should only be used in money-at-risk
     * production settings with extreme caution, after having read the code
     * carefully and satisfied yourself that you can live with the consequences
     * ABSOLUTELY NOT FOR CASUAL CONSUMPTION - YOU HAVE BEEN WARNED!
     * Charles Morcom 18 January 2006
     *=======================================================================*/
    virtual YieldCurveWrapper getYieldCurveWrapper() const;
    virtual ICDSParSpreadsWrapper getParSpreadsWrapper() const;

    /**Return the accrued interest due given settlement at a specified date
     * This is signed, and scaled by the notional and the fee (i.e. actual value).*/
    virtual double getAccruedInterest(const DateTime&      settlementDate,
                                      IForwardRatePricerSP model) const;

    /**Returns the present value of the CDS at valuationDate, conditional on default at defaultDate.
     * For a vanilla CDS, this would be (1-R)+accrued interest discounted risk-free
     * for any applicable recovery payment delay.
     * How should this be defined if the credit is known to have defaulted already?
     * Or if the instrument itself has been overridden as defaulted? Not sure about this
     * yet. This should probably give value for _total_ default, so may need more
     * methods to handle multiply defaultable instruments like index CDS.*/
    virtual double getPVGivenDefault(
        const DateTime& valuationDate, 
        const DateTime& defaultDate,
        const IDiscountCurveRisky& crv,
        IForwardRatePricerSP model) const;

    /**Returns true, since CredDefSwap is never perpetual*/
    virtual bool hasFiniteMaturity() const;

    /**The end date of the CDS - either the last fee payment date, or the end of
       the protection period, whichever is later */
    virtual DateTime getMaturity() const;

    /**Returns earlier of first accrual date, or protection start date*/
    virtual DateTime getStartDate() const;

    virtual double getNotional() const;

    virtual void setNotional(double newNotional);

    virtual DayCountConventionSP getAccrualDcc() const;

    virtual CashFlowArraySP getInstrumentCashFlows(
        IForwardRatePricerSP            model) const;

    /*=========================================================================
     * END ICDSWithConvention methods
     *=======================================================================*/

    //-----------------------
    // ICDSConvention methods
    //-----------------------

    /**Make a new CDS with the same defining characteristics as this one, but with
     * different start and end dates. Assumes protection and accruals both start and
     * end on the same date. */
    virtual ICDSSP generateCDS(const DateTime& startDate,
                               const DateTime& endDate,
                               double feeRate) const;

    //----------------------
    // ICreditFeeLeg methods
    //----------------------

    /** Return the value of this leg */
    virtual double price(
        const DateTime&             today,
        const DateTime&             valDateCF,
        const IDiscountCurveRiskySP effectiveCurve,
        const YieldCurveWrapper&    discount,
        const double&               lowStrike,
        const double&               highStrike,
        const double&               outstandingNotional,
        const CashFlowArray&        pastTrancheLosses,
        double&                     riskyDurationTotal,  // (O) notional weighted risky duration
        double&                     riskyNotionalsMean,  // (O) mean value of notional per period
        double&                     risklessCFPV,        // (O) fair value of riskless payments
        bool                        computeExtra,        // true: populate arrays
        DoubleArray&                debugUnitPrice,      // price for each leg unit
        DoubleArray&                debugUnitHistPrice,  // price for each leg unit due to historical default
        CashFlowArraySP             rebatePayments,
        BoolArrayConstSP            payAsYouGoArray,
        IntArrayConstSP             numDelayDaysArray,
        DateTimeArrayConstSP        startDates,
        DateTimeArrayConstSP        endDates,
        DateTimeArrayConstSP        paymentDates,
        IBadDayAdjusterConstSP      bda,
        IForwardRatePricerSP        model) const;

    /** Return all cash flow dates */
    virtual DateTimeArraySP getCashFlowDates() const;

    /** Return risky cash flow dates */
    virtual DateTimeArraySP getRiskyCashFlowDates() const;

    /** Return riskfree cash flow dates */
    virtual DateTimeArraySP getRisklessCashFlowDates() const;

    /** Returns the accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getAccrualPeriods() const;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRiskyAccrualPeriods() const;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRisklessAccrualPeriods() const;

    /** Returns all cash flows */
    virtual AbstractCashFlowArrayConstSP getCashFlows(IForwardRatePricerSP model) const;

    /** Returns risky cash flows only */
    virtual CashFlowArraySP getRiskyCashFlows(IForwardRatePricerSP model) const;

    /** Returns risk free cash flows only */
    virtual CashFlowArraySP getRisklessCashFlows(IForwardRatePricerSP model) const;

    /** Returns risky notional dates */
    virtual DateTimeArraySP getRiskyNotionalDates(IForwardRatePricerSP model) const;

	/** Returns risky coupon notional types */
	virtual CouponNotionalTypesArraySP getRiskyCouponNotionalTypes() const;

	/** Returns risky observation dates */
	virtual DateTimeArraySP getRiskyObservationDates() const;
	
	/** Return known cash flows corresponding to a CDO tranche */
    virtual CashFlowArraySP generateKnownCashFlows(
         const DateTime       today,
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         const double         pastTrancheLoss,
         IForwardRatePricerSP model);

    /** Return known cash flows corresponding to a non-defaulted CDS */
    virtual CashFlowArraySP generateKnownCashFlows(IForwardRatePricerSP model) const;

    /** Estimates the known cash flows (so they are not really "known")
     * corresponding to a CDO tranche - takes into account estimated losses
     * in the future */
    virtual CashFlowArraySP estimateKnownCashFlows(
         const double        initialTrancheSize,
         const DateTimeArray pastTrancheLossDates,
         const DoubleArray   pastTrancheLosses,
         IForwardRatePricerSP model);

    /** Returns the known cash flows corresponding to a defaulted CDS
     * (taking accrued payments into consideration), that happen after a specific
     * date (typically used to exclude cashflows in the past, already paid).
     * If excludePaymentsBeforeDate is empty, all cashflows will be returned.
     *
     * CAUTION:
     * It does not aggregate cashflows happening on the same dates into one.
     * If this is required (e.g., to output these cashflows for MiddleOffice use)
     * you should call CashFlow::agregate(...) on the resulting CashFlowArray -
     * It is not done here for performance reasons */
    virtual CashFlowArraySP generateKnownCashFlowsGivenDefault(
        const DateTime&            valueDate,
        const DateTime&            defaultDate,
        const DateTime&            excludePaymentsBeforeDate,
        const DateTime&            protectionStartDate,
        const DateTime&            protectionEndDate,
        const bool                 allowIncludingTodaysPayments,
        IForwardRatePricerSP       model,
        ICreditEventOverrideNameSP creditEventOverride,
        IBadDayAdjusterConstSP     badDayAdjuster,
        CIntSP                     triggerDelay,
        CIntSP                     defaultToSettlementDelay,
        const DateTime&            lastTriggerDate) const;

    /** Returns the pay date which is last in terms of time */
    virtual DateTime getLastPayDate() const;

    /** Returns the observation date which is last in terms of time */
    virtual DateTime getLastObservationDate() const;

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate) const;

    /** Feedback method for getting information about a fee cashflow */
    virtual void getActiveFee(
        const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
        const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start 
                                                  //     dates, and we are interested in the first fee
        IForwardRatePricerSP model,               // (I) for calculating the amount
        DateTime&            accrueStartDate,     // (O) when the fee starts accruing
        DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
        DateTime&            paymentDate,         // (O) when the fee is paid
        double&              amount) const;       // (O) the cashflow amount

    /** Price this fee leg assuming it corresponds to a defaulted CDS */
    virtual double pv(const DateTime&         valueDate,
                      const DateTime&         defaultDeterminationDate,
                      const DateTime&         accrualPaymentDate,
                      const YieldCurveConstSP discount,
                      const bool              computeAccrual,
                      IForwardRatePricerSP    model) const;

    /** Compute PV of fee leg at valuation date. */
    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               IForwardRatePricerSP         model) const;

    /**
    * Compute PV of fee leg at valuationDate, but for unconditional payment at
    * paymentDate, conditional on no new defaults before valuationDate.
    */
    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               IForwardRatePricerSP         model) const;

    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               bool                         defaultValueOnly,
                               IForwardRatePricerSP         model) const;

    /** Returns the accrued interest */
    virtual double getFeeLegAI(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestAccrualStart,
                               const DateTime&              latestAccrualEnd,
                               const DayCountConventionSP   dcc, //allows an override to be specified
                               const IDiscountCurveRisky&   crv,
                               const IDiscountCurveConstSP  discount,
                               const IDecretionCurveConstSP prepay,
                               IForwardRatePricerSP         model) const;

    /** Prices the leg sufferring a default, under the (optional) credit event */
    virtual double getFeeLegDefaultedPV(
        const DateTime&              valuationDate,
        const DateTime&              defaultDate,
        const DateTime&              protectionStartDate,
        const DateTime&              protectionEndDate,
        const bool                   allowIncludingTodaysPayments,
        const IDiscountCurveConstSP  discount,
        const IDecretionCurveConstSP prepay,
        IForwardRatePricerSP         model,
        IBadDayAdjusterConstSP       badDayAdjuster,
        ICreditEventOverrideNameSP   creditEventOverride,
        CIntSP                       triggerDelay,
        CIntSP                       defaultToSettlementDelay,
        DateTime                     lastTriggerDate) const;

    /** Returns the leg notional */
    virtual double getFeeLegNotional() const;

    /** Returns the leg notional */
    virtual void setFeeLegNotional(double newNotional);

	/** creates the SV Generator
		Note this method throws exception - ASarma
	*/
	virtual ICreditFeeLegSVGenSP createSVGen(
		const DateTime& valueDate,
		const DateTimeArray& modifiedFeeLegObservationDates,
		const DateTimeArray& productTimeline, //product timeline
		const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
		double lossConfigNotional
		) const;

    //-------------------------------
    // IFixedRateCreditFeeLeg methods
    //-------------------------------

    /**Get the fee rate of the fee leg.*/
    virtual double getRate() const;

    /**Change the fee rate of the fee leg. Note that, if the rate is currently zero,
       this will throw an exception.*/
    virtual void setRate(double newRate);

    //-----------------------------
    // ICreditContingentLeg methods
    //-----------------------------

    /** Price this contingent leg. */
    virtual double price(
        double                       initialNotional,     // highStrike - lowStrike
        double                       outstandingNotional, // initialTrancheSize - pastTrancheLoss
        const DateTime&              today,
        const DateTime&              valDateCF, // to be scrapped
        const IDiscountCurveRiskySP  effectiveCurve,
        const CashFlowArray&         pastTrancheLosses,
        const BoolArray&             payPastTrancheLosses,
        bool                         computeDebugPrices, // true: populate arrays below
        DoubleArray&                 debugUnitPrice, // price for each leg unit
        DoubleArray&                 debugUnitHistPrice, // price for each leg unit due to historical default
        IBadDayAdjusterConstSP       bda) const;

    /** Returns the earliest observation start date */
    virtual DateTime firstObservationStartDate() const;

    /** Returns the last pay date */
    virtual DateTime lastPayDate(IBadDayAdjusterConstSP bda) const;

    /** Returns the last observation date */
    virtual DateTime lastObservationEndDate() const;

    /** When to stop tweaking */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate,
                                    IBadDayAdjusterConstSP bda) const;

    virtual CashFlowArraySP generateKnownCashFlows(
        const DateTime&        today,
        double                 initialTrancheSize,
        const CashFlowArray&   pastTrancheLosses,
        const BoolArray&       payPastTrancheLosses,
        IBadDayAdjusterConstSP bda) const;

    /** CAUTION: HACK!
     * This method is for the benefit of the fee leg only: it returns
     * payment details required in the fee leg (because they should
     * have been added to the CDO instrument in the first place, but
     * now it is too late to change them). Although it would be possible
     * to add these fields to the fee leg also, it has been estimated
     * that this would cause confussion among the library users
     * (specially the start/end arrays) and instead they will be
     * obtained from the contingent leg and passed to the fee leg as
     * required.
     * All parameters are really outputs of the method - the contents of
     * the smart pointers passed in will be discarded */
    virtual void getPaymentInformation(BoolArraySP&     payAsYouGoArray,
                                       IntArraySP&      numDelayDaysArray,
                                       DateTimeArraySP& startDate,
                                       DateTimeArraySP& endDate,
                                       DateTimeArraySP& paymentDate) const;

    /* Checks if the input date is covered for protection, i.e., falls in
     * one of the observation periods of this leg */
    virtual bool isDateCoveredForProtection(const DateTime& date) const;

    /**Compute PV of contingent leg at valuation date, with instrument default settlement and
       payment date behaviour.*/
    virtual double getContingentLegPV(const DateTime&            valuationDate,
                                      const IDiscountCurveRisky& crv,
                                      IBadDayAdjusterConstSP     bda) const;

    /**Compute PV of contingent leg at valuationDate, but for unconditional payment at
       paymentDate, conditional on no new defaults before valuationDate.*/
    virtual double getContingentLegPV(const DateTime&            valuationDate,
                                      const DateTime&            paymentDate,
                                      const IDiscountCurveRisky& crv,
                                      IBadDayAdjusterConstSP     bda) const;

    /** Return the recovery rate used
        crv should be the underlying if the leg does not support an override */
    virtual double getRecovery(const IDiscountCurveRisky& crv) const;

    /** Returns the amount that would be recovered upon default */
    virtual double recoveredValue(const DateTime& valueDate,
                                  const IDiscountCurveRisky& crv,
                                  const IDecretionCurveConstSP prepay,
                                  const IDiscountCurveRisky::RecoveryType recType) const;

    /** Returns the amount that would be recovered upon default,
        using the supplied recovery rate */
    virtual double recoveredValue(const DateTime& valueDate,
                                  const IDiscountCurveRisky& crv,
                                  const IDecretionCurveConstSP prepay,
                                  const double recoveryToUse,
                                  const IDiscountCurveRisky::RecoveryType recType) const;

    /** Prices the leg sufferring a default, under the (optional) credit event */
    virtual double getContingentLegDefaultedPV(
        const DateTime&              valuationDate,
        const IDiscountCurveRisky&   crv,
        const IDiscountCurveConstSP  discount,
        const IDecretionCurveConstSP prepay,
        const DateTime&              defaultDate,
        IBadDayAdjusterConstSP       badDayAdjuster,
        const bool                   allowIncludingTodaysPayments,
        ICreditEventOverrideNameSP   creditEventOverride,
        CIntSP                       triggerDelay,
        CIntSP                       defaultToSettlementDelay,
        DateTime                     lastTriggerDate) const;

    /** Returns the leg notional */
    virtual double getContingentLegNotional() const;

    /** Sets the leg notional */
    virtual void setContingentLegNotional(double newNotional);

    //-----------------------------
    // existing CredDefSwap methods
    //-----------------------------

    /** Returns all known cash flows */
    CashFlowArraySP knownCashFlows()const;

    /** when do payments occur ? */
    DateTimeArraySP paymentDates() const;

    /** used to override CREDIT_SPREAD_RHO on cdsParSpreads
        ObjectIteration::IOverride interface */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** instrument specific market data handling */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Implementation of ClosedFormCDSPS::IntoProduct interface */
    virtual ClosedFormCDSPS::IProduct* createProduct(ClosedFormCDSPS* model) const;

    /** Implementation of ClosedFormFA::IntoProduct interface */
    virtual ClosedFormFA::IProduct* createProduct(ClosedFormFA* model) const;

    /** Implementation of ClosedFormCDSPSandFA::IntoProduct interface */
    virtual ClosedFormCDSPSandFA::IProduct* createProduct(ClosedFormCDSPSandFA* model) const;

    /** Implementation of CreditSupport::Interface interface */
    virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(DeltaToCredit* shift) const;

    /** Shifts the object using given shift */
    virtual bool sensShift(DeltaToCredit* shift);

    /** convert liquidity spread rho into par spread equivalent */
    virtual double calcParSpreadEquivalent(const double& lsSpreadRho) const;

    //++++++++++++++++++++++++++++++++++++++++
    //  IBadDayAdjuster methods
    //
    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const;
    //
    //  IBadDayAdjuster methods
    //------------------------------------------


    /** converts a CDS par spread curve into a clean spread curve. This method lives here, because it uses the CDS pricing
        heavily. The CDS code implements its own default rate classes, which I chose not to move into the market directory,
        as the CDS pricing code and the default rate classes interfere heavily with each other - instead, a new class
        CleanSpreadCurve is born (ASe) */
    static CleanSpreadCurveSP getCleanSpreadCurve(CDSParSpreadsSP     parSpreadCurve,
                                                  YieldCurveSP        discountCurve,
                                                  const DateTime&     valueDate,
                                                  const bool          returnFwdRates);

    CredDefSwap(const DateTime&        valueDate,
                const bool&            oneContract,
                const double&          notional,
                const string&          ccyTreatment,
                InstrumentSettlementSP instSettle,
                InstrumentSettlementSP premiumSettle,
                CAssetWrapper          asset,
                ICDSParSpreadsWrapper  cdsParSpreads,
                YieldCurveWrapper      discount,
                const DateTime&        swapEffectiveDate,
                const DateTime&        protectionEndDate,
                CashFlowArraySP        feePayments,
                const bool&            useSwapRecovery,
                const double&          swapRecovery,
                const bool&            payAccruedFee,
                const string&          dcc,
                const string&          bdc,
                HolidayWrapper         settlementHols);

    // REDUCED INTERFACE CLASS
    class PRODUCTS_DLL QuickPricer: public Generic1FactorCredit,
                       virtual public ClosedFormCDSPS::IIntoProduct,
                       virtual public ClosedFormFA::IIntoProduct,
                       virtual public ClosedFormCDSPSandFA::IIntoProduct,
                       virtual public LastSensDate,
                       virtual public IGetMarket,
                       public CreditSupport::Interface
    {
    public:
        static CClassConstSP const TYPE;
        friend class CDSQuickPricerHelper;

        void validatePop2Object();

        // Generic1FactorCredit interface
        virtual void Validate();
        virtual void GetMarket(const IModel* model, const CMarketDataSP market);
        virtual DateTime getValueDate() const;

        // IIntoProduct interfaces
        virtual ClosedFormCDSPS::IProduct* createProduct(ClosedFormCDSPS* model) const;
        virtual ClosedFormFA::IProduct* createProduct(ClosedFormFA* model) const;
        virtual ClosedFormCDSPSandFA::IProduct* createProduct(ClosedFormCDSPSandFA* model) const;

        // CreditSupport interface
        virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

        // LastSensDate interface
        virtual DateTime endDate(const Sensitivity* sensControl) const;

        // IGetMarket interface
        virtual void getMarket(const IModel* model, const MarketData* market);

        // These methods are currently only used internally but are public in CredDefSwap
        // so included to maintain consistency
        CashFlowArraySP knownCashFlows() const;
        DateTimeArraySP paymentDates() const;

//         QuickPricer(const DateTime&         valueDate,
//                     const bool&             oneContract,
//                     const double&           notional,
//                     const string&           ccyTreatment,
//                     InstrumentSettlementSP  instSettle,
//                     InstrumentSettlementSP  premiumSettle,
//                     CAssetWrapper           asset,
//                     ICDSParSpreadsWrapper   cdsParSpreads,
//                     YieldCurveWrapper       discount,
//                     const DateTime&         swapEffectiveDate,
//                     const bool&             useSwapRecovery,
//                     const double&           swapRecovery,
//                     const bool&             payAccruedFee,
//                     const string&           dcc,
//                     const string&           bdc,
//                     const double&           dealSpread,
//                     const int&              frequency,
//                     const DateTime&         maturityDate);

    private:
        QuickPricer();
        QuickPricer(const QuickPricer& rhs);            // not implemented
        QuickPricer& operator=(const QuickPricer& rhs); // not implemented

        // Inputs
        //double                  notional;
        DateTime                swapEffectiveDate;
        bool                    useSwapRecovery;
        double                  swapRecovery;       // overidden by par recovery rate if useSwapRecovery=FALSE
        bool                    payAccruedFee;      // whether to make accrued payment on default
        string                  dcc;                // DCC for accrual periods
        string                  bdc;
        double                  dealSpread;         // CDS spread
        int                     frequency;          // annual frequency of CDS payments
        DateTime                maturityDate;
        smartPtr<CredDefSwap>   cds;
        HolidayWrapper          settlementHols;     // Holidays, used to settle payments
    };

protected:
    // for inheritance - only done for 'flow' version of CredDefSwap which is
    // called CreditDefaultSwap
    CredDefSwap(CClassConstSP clazz);

    /** Validates that the dcc field is not empty (it is flagged as
        optional) and that the bdc is not supplied (it is to be
        removed) */
    void validateDCCSuppliedAndNoBDC() const;

private:
    CredDefSwap();
    CredDefSwap(const CredDefSwap& rhs);
    CredDefSwap& operator=(const CredDefSwap& rhs);
    void priceParSpreads(CResults* results, Control* control, IForwardRatePricerSP model) const;
    void priceFirmAsset(CResults* results, Control* control, IForwardRatePricerSP model) const;

    /** Returns the CDS price when the underlying name has defaulted */
    double getPriceIfDefaulted(Control* control,
        IForwardRatePricerSP model) const;

    /**Return the accrued interest due given settlement at a specified date
     * This is signed, and scaled by the notional and the fee (i.e. actual value).*/
    double getAccruedInterest(const DateTime& settlementDate,
                              IForwardRatePricerSP model,
                              int& feeIdx) const;

    void addRequests(Control*             control,
                     CResults*            results,
                     CashFlowArraySP      cleanSpreadCurve,
                     IObjectSP            currentSpread,
                     IForwardRatePricerSP model) const;

    /** Obtains the fee leg cashflows using the credit event override if
        available.
        Depending on the value of the pricingAccrue flag, two things can be
        configured:
        - Riskless fees are only included if pricingAccrue is false (ie,
          when pricing a defaulted fee leg).
        - Cashflows on "pretendedValueDate" may get included in the valuation
          for backwards compatibility under certain conditions, but only if
          pricingAccrue is false.
        The "pretendedDefaultDate" and "pretendedValueDate" are
        required to support the (THETA_)ACCRUED_INTEREST output requests */
    double priceFeeLegGivenDefault(const DateTime& pretendedDefaultDate,
                                   const DateTime& pretendedValueDate,
                                   int& feeIdx,
                                   IForwardRatePricerSP model,
                                   const bool pricingAccrue) const;

    /** Prices the contingent leg of a defaulted CDS, using the credit event
     * override if available */
    double priceContingentLegGivenDefault() const;

    // Calculate the theta accrued interest - diff between accrued interest tommorrow and today
    double calcThetaAccruedInterest(IForwardRatePricerSP model)const;

    // Converts basis of default rate from continuous to annualised
    CashFlowArraySP annualiseDefaultRates(const CashFlowArray& defaultRates)const;

    double getRecovery() const;
    DateTime getProtectionEndDate() const;

    // Inputs
    DateTime                swapEffectiveDate; // aka protection start date
    DateTime                protectionEndDate; /* optional, use last fee date
                                                  if not specified */

	double                  feeRate; /* optional for backwards compatibility- set to 1.0 if not used */
    mutable CashFlowArraySP feePayments;
    bool                    useSwapRecovery;
    double                  swapRecovery;
    bool                    payAccruedFee;
    string                  dcc;
    string                  bdc;
    mutable CashFlowArraySP adjFeePayments;
    CashFlowArraySP         payAsYouGo;      //for ABCDS, optional, riskless discounting

    ICreditEventOverrideNameSP creditEventOverride; // Override for credit event parameters
    CIntSP triggerDelay; // Delay between default and eventDeterminationDate, in days
    CIntSP defaultToSettlementDelay; // Delay between credit event and settlement

    DateTime lastTriggerDate; // Last date when a default occurred during
                              // the protection period can be triggered

    HolidayWrapper settlementHols; // Holidays, used to settle payments

    // Transient
    mutable DayCountConventionSP swpAccrualDCC;          // used if payAccruedFee = true
    mutable BadDayConventionSP   BDC;                    // bad day conv for BM adjustment
    bool                         createE2Csensitivities; // whether to create debt/equity sensitivities
    mutable bool                 e2cBasePriceCalculated;
    mutable double               e2cBasePrice;
};


typedef smartConstPtr<CredDefSwap>  CredDefSwapConstSP;
typedef smartPtr<CredDefSwap>       CredDefSwapSP;

DRLIB_END_NAMESPACE
#endif




