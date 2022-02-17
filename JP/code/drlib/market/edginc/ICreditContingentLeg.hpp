//----------------------------------------------------------------------------
//
//   File        : ICreditContingentLeg.hpp
//
//   Description : A general interface describing the functionality
//                 of a contingent (protection) leg.
//                 Usable by all forms of credit instrument.
//
//----------------------------------------------------------------------------

#ifndef QR_ICREDITCONTINGENTLEG_HPP
#define QR_ICREDITCONTINGENTLEG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE(FlatFwdZeroCurve);
FORWARD_DECLARE(IDiscountCurveRisky);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IDecretionCurve);
FORWARD_DECLARE(ICreditEventOverrideName);

/** Interface describing a contingent leg. Basically this pays
    out on default(s) */
class MARKET_DLL ICreditContingentLeg: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    ICreditContingentLeg();
    virtual ~ICreditContingentLeg();

    /** Price this contingent leg. */
    virtual double price(
        double                      initialNotional,     // highStrike - lowStrike
        double                      outstandingNotional, // initialTrancheSize - pastTrancheLoss
        const DateTime&             today,
        const DateTime&             valDateCF, // to be scrapped
        const IDiscountCurveRiskySP effectiveCurve,
        const CashFlowArray&        pastTrancheLosses,
        const BoolArray&            payPastTrancheLosses,
        bool                        computeDebugPrices, // true: populate arrays below
        DoubleArray&                debugUnitPrice, // price for each leg unit
        DoubleArray&                debugUnitHistPrice, // price for each leg unit due to historical default
        IBadDayAdjusterConstSP      bda) const = 0;

    /** Returns the earliest observation start date */
    virtual DateTime firstObservationStartDate() const = 0;

    /** Returns the last pay date */
    virtual DateTime lastPayDate(IBadDayAdjusterConstSP bda) const = 0;
    
    /** Returns the last observation date */
    virtual DateTime lastObservationEndDate() const = 0;
    
    /** When to stop tweaking */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate,
                                    IBadDayAdjusterConstSP bda) const = 0;
    
    virtual CashFlowArraySP generateKnownCashFlows(
        const DateTime&        today,
        double                 initialTrancheSize,
        const CashFlowArray&   pastTrancheLosses,
        const BoolArray&       payPastTrancheLosses,
        IBadDayAdjusterConstSP bda) const = 0;

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
                                       DateTimeArraySP& paymentDate) const = 0;

    /* Checks if the input date is covered for protection, i.e., falls in
     * one of the observation periods of this leg */
    virtual bool isDateCoveredForProtection(const DateTime& date) const = 0;

    /**Compute PV of contingent leg at valuation date, with instrument default settlement and
       payment date behaviour.*/
    virtual double getContingentLegPV(const DateTime&            valuationDate, 
                                      const IDiscountCurveRisky& crv,
                                      IBadDayAdjusterConstSP     bda) const = 0;
    
    /**Compute PV of contingent leg at valuationDate, but for unconditional payment at 
       paymentDate, conditional on no new defaults before valuationDate.*/
    virtual double getContingentLegPV(const DateTime&            valuationDate, 
                                      const DateTime&            paymentDate, 
                                      const IDiscountCurveRisky& crv,
                                      IBadDayAdjusterConstSP     bda) const = 0;

    /** Return the recovery rate used
        crv should be the underlying if the leg does not support an override */
    virtual double getRecovery(const IDiscountCurveRisky& crv) const = 0;

    /** Returns the amount that would be recovered upon default
        recType allows for various recovery strategies */
    virtual double recoveredValue(const DateTime& valueDate,
                                  const IDiscountCurveRisky& crv,
                                  const IDecretionCurveConstSP prepay,
                                  const IDiscountCurveRisky::RecoveryType recType) const = 0;

    /** Returns the amount that would be recovered upon default,
        using the supplied recovery rate */
    virtual double recoveredValue(const DateTime& valueDate,
                                  const IDiscountCurveRisky& crv,
                                  const IDecretionCurveConstSP prepay,
                                  const double recoveryToUse,
                                  const IDiscountCurveRisky::RecoveryType recType) const = 0;

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
        DateTime                     lastTriggerDate) const = 0;
    
    /** Returns the leg notional */
    virtual double getContingentLegNotional() const = 0;

    /** Sets the leg notional */
    virtual void setContingentLegNotional(double newNotional) = 0;

private:
    static void load(CClassSP& clazz);
};

// Support for smart pointers
typedef smartPtr<ICreditContingentLeg> ICreditContingentLegSP;
typedef smartConstPtr<ICreditContingentLeg> ICreditContingentLegConstSP;

/** Describes something (like a CDS or bond) which has a contingent leg which is
  * publicly accessible.
  */
class MARKET_DLL IHasCreditContingentLeg : public virtual IObject {
public:

    static CClassConstSP const TYPE;

    IHasCreditContingentLeg();
    virtual ~IHasCreditContingentLeg();

    /** Return the contingent leg */
    virtual ICreditContingentLegSP getContingentLeg() const = 0;

private:
    static void load(CClassSP& clazz);
};

/** A credit contingent leg generator is an object which is capable of generating a contingent
  * leg, if you give it new start and end dates. */
class MARKET_DLL ICreditContingentLegGenerator : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    ICreditContingentLegGenerator();
    virtual ~ICreditContingentLegGenerator();

    /**Generate a new contingent leg with specified start and end dates whose other
       characteristics are defined by this. "startDate" should be the start of protection;
       "endDate" should be the end of protection.*/
    virtual ICreditContingentLegSP generateCreditContingentLeg(
        const DateTime& startDate,
        const DateTime& endDate
        ) const = 0;

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
