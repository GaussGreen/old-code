//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : VanillaCreditContingentLeg.hpp
//
//   Description : Contingent Leg for a Vanilla Credit Instrument
//
//   Author      : Doris Morris
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef VANILLA_CREDIT_CONTINGENT_LEG_HPP
#define VANILLA_CREDIT_CONTINGENT_LEG_HPP


#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/Theta.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Settlement);
FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(DayCountConvention);
FORWARD_DECLARE_WRAPPER(Holiday);

/** Base Contingent Leg class for a Vanilla Credit Instrument*/
class PRODUCTS_DLL VanillaCreditContingentLeg : public CObject,
                                                public virtual ICreditContingentLeg,
                                                public virtual ICreditContingentLegGenerator,
                                                public virtual Theta::Shift, 
                                                public virtual IGetMarket
{
public:

    //-----------------------------------
    // VanillaCreditContingentLeg methods
    //-----------------------------------

    static CClassConstSP const TYPE;

    VanillaCreditContingentLeg(const  DateTime& valueDate,
                               const  DateTime& startDate,
                               const  DateTime& endDate,
                               double baseNotional,
                               double recoveryRate,
                               double delay,
                               const  YieldCurveWrapper discountCrv,
                               const  BadDayConventionSP pmtBdc);

    virtual ~VanillaCreditContingentLeg();

    /**Set the base notional of an instrument. This is here because some derivatives
       need to be able to control the notional of their underlyings to guarantee
       that the valuations make sense. In general the instrument will be loaded with
       the correct notional, so this method will not be commonly used.*/
    void setNotional(double newNotional);

    //----------------
    // CObject methods
    //----------------

    virtual void validatePop2Object();

    //-----------------------------
    // ICreditContingentLeg methods
    //-----------------------------

    /** Price this contingent leg. */
    virtual double price(
        double                       initialNotional,     // highStrike - lowStrike
        double                       outstandingNotional, /* initialTrancheSize -
                                                             pastTrancheLoss */
        const DateTime&              today,
        const DateTime&              valDateCF, // to be scrapped
        const IDiscountCurveRiskySP  effectiveCurve,
        const CashFlowArray&         pastTrancheLosses,
        const BoolArray&             payPastTrancheLosses,
        bool                         computeDebugPrices, /* true: populate arrays 
                                                            below */
        DoubleArray&                 debugUnitPrice, // price for each leg unit
        DoubleArray&                 debugUnitHistPrice, /* price for each leg unit
                                                          * due to historical default */
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

    //--------------------------------------
    // ICreditContingentLegGenerator methods
    //--------------------------------------

    virtual ICreditContingentLegSP generateCreditContingentLeg(const DateTime& startDate,
                                                               const DateTime& endDate) const;

    //---------------------
    // Theta::Shift methods
    //---------------------

    virtual bool sensShift(Theta* shift);

    //-------------------
    // IGetMarket methods
    //-------------------

    virtual void getMarket(const IModel* model, const MarketData* market);



    ///****************** ICreditVanillaInstrument ******************/
    ///**Accrued interest (cash amount, not percentage) for settlement on settlementDate.*/
    //virtual double getAccruedInterest(const DateTime& settlementDate) const;

    ///**Calculates PV of the instrument, given a valuation date. 
    //   Settlement is whatever the instrument determines it to be (e.g. T+1 CDS; T+3 bond, etc.)
    //   relative to the valuation date. Value should be conditional on no new defaults 
    //   before valuationDate.*/
    //virtual double getPV(const DateTime&                 valuationDate, 
    //                     const IDiscountCurveRisky&      crv) const;

    ///**Calculates the PV of the instrument on a given valuationDate for unconditional settlement
    //   and payment on paymentDate. Used to calculate forward PVs and differentiate between
    //   conditional and unconditional settlement. Value should be conditional on no new defaults
    //   before valuationDate if valuationDate is in the future.*/
    //virtual double getPV(const DateTime&                 valuationDate, 
    //                     const DateTime&                 paymentDate, 
    //                     const IDiscountCurveRisky&      crv) const;

    ///**Calculates the PV of an instrument on defaultDate, given that default has just occured.
    //   Note that the curve is needed in case there is a delay between default and recovery
    //   payment. This method assumes a complete default, so should be interpreted with care
    //   when the instrument/curve has multiple names.*/
    //virtual double getPVGivenDefault(const DateTime&                 defaultDate,
    //                                 const IDiscountCurveRisky&      crv) const;

    ///** Is perReturns true if the instrument has a finite maturity. This will be false if
    //    it is perpetual.*/
    //virtual bool hasFiniteMaturity() const;

    ///**Returns the maturity of the leg. */
    //virtual DateTime getMaturity() const;

    ///**Returns the earliest possible "interesting" date: this will be the earliest of the
    //   start of the first accrual period, the start of the contingent leg, the first
    //   cash-flow, etc.*/
    //virtual DateTime getStartDate() const;

    ///** Return the base notional. */
    //virtual double getNotional() const;


    ///* Get cashflows */
    //virtual CashFlowArraySP getInstrumentCashFlows() const;

    ///*************************************************************/

    ///***************** IInstrument ***********************************/
    //virtual void GetMarket(const IModel*, const CMarketDataSP);

    ///** Called once before the initial pricing */
    //virtual void Validate();

    ///** override a control shift (eg for delta on trees) - may return
    //    null to use original. */
    //virtual CSensControl* AlterControl( const IModel*          modelParams,
    //                                    const CSensControl*    sensControl) const;

    ///** Returns the value date (aka today) the instrument is currently
    //    pricing for */
    //virtual DateTime getValueDate() const;

    ///** price a dead instrument until settlement - exercised, expired,
    //    knocked out etc.  returns true if it is dead (and priced), false
    //    if it is not dead */
    //virtual bool priceDeadInstrument(CControl* control,
    //                                 CResults* results) const;

    ///** Returns the name of the instrument's discount currency. */
    //virtual string discountYieldCurveName() const;

    ///*************************************************************/

private:

    // For reflection
    static void load (CClassSP& clazz);
    static IObject* defaultConstructor();
    VanillaCreditContingentLeg();


    // To create a VanillaCreditFeeLeg from a small set of input data
    VanillaCreditContingentLeg(const DateTime& startDate,
                               const DateTime& endDate);

    //-------
    // Fields
    //-------

    mutable double       baseNotional;          // initial notional
    DateTime             valueDate;             // Valuation date
    DateTime             effectiveDate;         // aka protection start date
    DateTime             maturityDate;          // Maturity date aka protection end date
    string               recoveryTypeString;    // Recovery type 1-R
    mutable double       recoveryRate;          // recovery rate

    double               delay;                 // aka payment offset for pyramid

    BadDayConventionSP   pmtBdc;                // Payment Bad Day Convention
    HolidayWrapper       pmtHol;                // Payment holiday

    YieldCurveWrapper    discount;              // discount curve

    /* converted data fields from input */
    IDiscountCurveRisky::RecoveryType recovery; // recoveryType -> recovery
};


typedef smartPtr<VanillaCreditContingentLeg> VanillaCreditContingentLegSP;
typedef smartConstPtr<VanillaCreditContingentLeg> VanillaCreditContingentLegConstSP;


DRLIB_END_NAMESPACE

#endif // VANILLA_CREDIT_CONTINGENT_LEG_HPP
