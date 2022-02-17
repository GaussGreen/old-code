//----------------------------------------------------------------------------
//
//   Filename    : BasicCreditContingentLeg.hpp
//
//   Description : Single period contingent leg
//
//----------------------------------------------------------------------------

#ifndef QLIB_BASICCREDITCONTINGENTLEG_HPP
#define QLIB_BASICCREDITCONTINGENTLEG_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

/** Single period contingent leg */
class PRODUCTS_DLL BasicCreditContingentLeg : public CObject,
                                              public virtual ICreditContingentLeg,
                                              public virtual ICreditContingentLegGenerator
{
public:

    //---------------------------------
    // BasicCreditContingentLeg methods
    //---------------------------------

    static CClassConstSP const TYPE;

    BasicCreditContingentLeg(const  DateTime& protectionStartDate,
                             const  DateTime& protectionEndDate);

    virtual ~BasicCreditContingentLeg();

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
    virtual double getContingentLegPV(const DateTime&             valuationDate, 
                                      const IDiscountCurveRisky&  crv,
                                      IBadDayAdjusterConstSP      bda) const;
    
    /**Compute PV of contingent leg at valuationDate, but for unconditional payment at 
       paymentDate, conditional on no new defaults before valuationDate.*/
    virtual double getContingentLegPV(const DateTime&            valuationDate, 
                                      const DateTime&            paymentDate, 
                                      const IDiscountCurveRisky& crv,
                                      IBadDayAdjusterConstSP     bda) const;

    /** return the recovery rate to be used */
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

private:

    // For reflection
    static void load (CClassSP& clazz);
    static IObject* defaultConstructor();
    BasicCreditContingentLeg();

    //-------
    // Fields
    //-------

    double    notional;        //positive => long protection
    DateTime  protectionStartDate;
    DateTime  protectionEndDate;
    bool      useSwapRecovery; // allow overriding of the underlying recovery rate
    CDoubleSP swapRecovery;    // the recovery rate to use if useSwapRecovery is False
};

DECLARE(BasicCreditContingentLeg)

DRLIB_END_NAMESPACE

#endif

