//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITCONTINGENTLEGBASE_HPP
#define QR_CREDITCONTINGENTLEGBASE_HPP

#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/IProtectionProvider.hpp"

DRLIB_BEGIN_NAMESPACE

/** Common base class for sharing implementation across 
    credit contingent legs */
class MARKET_DLL CreditContingentLegBase: public CObject,
                                          public virtual ICreditContingentLeg,
                                          public virtual IProtectionProvider
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    virtual ~CreditContingentLegBase();

    /** Called immediately after object constructed. This method uses
        the protected methods payingAsYouGo and numDelayDays - derived
        classes must do any necessary validation before this method is
        called */
    virtual void validatePop2Object();

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
        const DateTime&      today,
        double               initialTrancheSize,
        const CashFlowArray& pastTrancheLosses,
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
    void getPaymentInformation(BoolArraySP&     payAsYouGoArray,
                               IntArraySP&      numDelayDaysArray,
                               DateTimeArraySP& startDate,
                               DateTimeArraySP& endDate,
                               DateTimeArraySP& paymentDate) const;

    //returns a reference to the notionals of different periods
    virtual const DoubleArray&  notionals()  const;

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
    //  IProtectionProvider methods
    //-----------------------------

    /* Checks if the input date is covered for protection, i.e., falls in
     * one of the observation periods of this leg */
    virtual bool isDateCoveredForProtection(const DateTime& date) const;

protected:
    //// how many periods make up this contingent leg
    int numPeriods() const;

    /** Returns the number of delay days for payment for the specified
        period given by the supplied integer when 'paying as you
        go' (otherwise not defined). The integer must be within range */
    virtual int numDelayDays(int periodIdx) const = 0;

    /** Returns whether for the specified period we are 'paying as you go' */
    virtual bool payingAsYouGo(int periodIdx) const = 0;

    //// constructor for derived classes
    CreditContingentLegBase(CClassConstSP clazz);

    /** Explicit constructor for derived classes */
    CreditContingentLegBase(
        CClassConstSP clazz,
        const DateTimeArray& obsStartDate,
        const DateTimeArray& obsEndDate,
        const DateTimeArray& payDate,
        const DoubleArray& notional);


private:
    void calculatePayAsYouGo(
        int                         period,
        double                      initialNotional,     // highStrike - lowStrike
        double                      outstandingNotional, /* initialTrancheSize -
                                                            pastTrancheLoss */
        const DateTime&             today,
        const DateTime&             valDateCF, // to be scrapped
        const IDiscountCurveRiskySP effectiveCurve,
        const CashFlowArray&        pastTrancheLosses,
        const BoolArray&            payPastTrancheLosses,
        double&                     unitPrice, // price for each leg unit
        double&                     unitHistPrice, /* price for each leg unit
                                                    * due to historical default*/
        IBadDayAdjusterConstSP      bda) const;


    void calculatePayAtTheEnd(
        int                         period,
        double                      initialNotional,     // highStrike - lowStrike
        double                      outstandingNotional, /* initialTrancheSize -
                                                            pastTrancheLoss */
        const DateTime&             today,
        const DateTime&             valDateCF, // to be scrapped
        const IDiscountCurveRiskySP effectiveCurve,
        const CashFlowArray&        pastTrancheLosses,
        const BoolArray&            payPastTrancheLosses,
        double&                     unitPrice, // price for each leg unit
        double&                     unitHistPrice) const; /* price for each leg
                                                             unit due to
                                                             historical default*/

    static void load(CClassSP& clazz);

    //-------
    // FIELDS
    //-------
    DateTimeArray obsStartDate;/* observation start date */
    DateTimeArray obsEndDate;  /* observation end date */
    DateTimeArray payDate;     /* paydate */
    DoubleArray   notional;    /* notional for each unit */
};

DRLIB_END_NAMESPACE
#endif
