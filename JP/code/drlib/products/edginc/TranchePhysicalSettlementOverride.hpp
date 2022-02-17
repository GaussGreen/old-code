//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TranchePhysicalSettlementOverride.hpp
//
//   Description : Class used by CDO in order to override a name's default-
//                 related parameters, when the default is settled physically.
//
//   Author      : Jose Hilera
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHEPHYSICALSETTLEMENTOVERRIDE_HPP
#define QLIB_TRANCHEPHYSICALSETTLEMENTOVERRIDE_HPP

#include "edginc/TrancheCreditEventOverride.hpp"
#include "edginc/NoticeOfPhysicalSettlement.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(DeliveryDetails);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);

/** Class used by CDO in order to override a name's default-related parameters,
 * when the default is settled physically. */
class TranchePhysicalSettlementOverride: public TrancheCreditEventOverride { 
public:
    static CClassConstSP const TYPE;
    virtual ~TranchePhysicalSettlementOverride();

    /** Public constructor, with all fields passed in */
    TranchePhysicalSettlementOverride(
        CIntSP triggerDelay,
        DateTime eventDeterminationDate,
        CDoubleSP recovery,
        DateTime valueDate,
        NoticeOfPhysicalSettlementSP NoPS,
        DeliveryDetailsArraySP deliveries,
        CBoolSP settleDeliveriesSeparately,
        CIntSP defaultToCalculationDelay,
        const DateTime& calculationDate,
        const DateTime& cutoffDate,
        CBoolSP moreDeliveriesPending,
        double notionalFractionSettlingPhysically);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the contingent leg loss cashflows corresponding to a defaulted 
     * name. This object contains the credit event parameters */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

    /** Returns the fee leg loss cashflows corresponding to a defaulted name,
     * assuming that the default settles in one date.
     * This object contains the credit event parameters. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

    /* Returns the overall recovery rate accross all settlements, such that
     * total loss = overall defaulted notional * (1-overall recovery rate) */
    virtual double getOverallRecoveryRate(double recoveryRate) const;

    /** Returns the percentage of the notional that has been triggered 
     * for protection */
    virtual double getOverallDefaultedNotionalFraction() const;

private:
    // For reflection
    static void load (CClassSP& clazz);
    TranchePhysicalSettlementOverride();
    static IObject* defaultTranchePhysicalSettlementOverride();

    /** Returns the overall loss resulting from the settlement of this name */
    double getOverallLoss(const double notional,
                          const double recoveryRate) const;

    /** Returns the notional specified in the NoPS, or 
     * notionalFractionSettlingPhysically if no NoPS have been provided */
    double getTotalExpectedNotional() const;

    /** Returns the overall notional fraction scheduled for delivery */
    double getNotionalScheduledForDelivery() const;

    double getUnscheduledNotionalFraction() const;
    
    /** Estimates the calculation date for the outstanding notional with no
     * delivery information */
    DateTime commonCalculationDate(const DateTime& creditEventDate,
                                   const DateTime& lastTriggerDate,
                                   IBadDayAdjusterConstSP bda) const;

    // Fields                      ***WARNING***
    //
    // Any fields changed/added/removed here should also be updated
    // in TranchePartCashSettlementOverride, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly
    /** Specification of the amounts being settled */
    NoticeOfPhysicalSettlementSP NoPS; 

    /** Array with the details of all deliveries */
    DeliveryDetailsArraySP deliveries;

    /** If true, each delivery will have its own calculation and settlement 
     * date; otherwise all deliveries will share a common calculation and 
     * settlement date. */
    bool settleDeliveriesSeparately;

    /** An estimate for the delay between a credit event and the corresponding 
     * calculation date. Used to estimate the calculation date either for all 
     * obligations (if field settleDeliveriesSeparately is false) or for the 
     * remaining expected obligations with no delivery information yet 
     * available (if "Settle deliveries separately" is true). */
    CIntSP defaultToCalculationDelay;

    /** Date when the recovery rate is determined if NOT settling deliveries
     * separately*/
    DateTime calculationDate;

    /** Indicates the last date when deliveries are accepted. If the "Recover
     * notional" flag is set, any outstanding notional not delivered on or 
     * before this date will be deemed recovered.
     * If this date is empty and "more deliveries pending" is true (see below), 
     * outstanding deliveries will be expected indefinitely. */
    DateTime cutoffDate;
    
    /** Indicates whether more deliveries are to be expected before the 
     * cut-off date. If field "Settle deliveries separately" is false no more
     * deliveries will be expected after the calculation date, regardless of 
     * the value of this field. If "Settle deliveries separately" is true, 
     * this flag determines whether the settlement of this name's default will 
     * be complete when the last scheduled delivery settles, or more unscheduled 
     * deliveries are to be expected. 
     * No deliveries will be deemed pending after cut-off date regardless of the 
     * value of this flag. */
    bool moreDeliveriesPending;

    // Transient field
    /** Fraction of the notional expected to settle physically */
    double notionalFractionSettlingPhysically;
};

DECLARE (TranchePhysicalSettlementOverride);

DRLIB_END_NAMESPACE

#endif
