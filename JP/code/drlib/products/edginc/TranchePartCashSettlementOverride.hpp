//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TranchePartCashSettlementOverride.hpp
//
//   Description : Class used by CDO in order to override a name's default-
//                 related parameters, when the default is settled physically
//                 with a part cash settlement
//
//   Author      : Jose Hilera
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHEPARTCASHSETTLEMENTOVERRIDE_HPP
#define QLIB_TRANCHEPARTCASHSETTLEMENTOVERRIDE_HPP

#include "edginc/TrancheCreditEventOverride.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(DeliveryDetails);
FORWARD_DECLARE(TranchePhysicalSettlementOverride);
FORWARD_DECLARE(TrancheCashSettlementOverride);
FORWARD_DECLARE(NoticeOfPhysicalSettlement);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);

/** Class used by CDO in order to override a name's default-related parameters,
 * when the default is settled physically with a part cash settlement. */
class TranchePartCashSettlementOverride: public TrancheCreditEventOverride { 
public:
    static CClassConstSP const TYPE;
    virtual ~TranchePartCashSettlementOverride();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the contingent leg loss cashflows corresponding to a defaulted 
     * name. */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

    /** Returns the fee leg loss cashflows corresponding to a defaulted name,
     * assuming that the default settles in one date. */
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
    TranchePartCashSettlementOverride();
    static IObject* defaultTranchePartCashSettlementOverride();

    // Fields

    // ... for physical settlement
    /** Specification of the amounts being settled */
    NoticeOfPhysicalSettlementSP NoPS; 
    /** Array with the details of all deliveries */
    DeliveryDetailsArraySP deliveries;
    /** If true, each delivery will have its own calculation and settlement 
     * date; otherwise all deliveries will share a common calculation and 
     * settlement date. */
    CBoolSP settleDeliveriesSeparately;
    /** An estimate for the delay between a credit event and the corresponding 
     * calculation date. Used to estimate the calculation date either for all 
     * obligations (if field settleDeliveriesSeparately is false) or for the 
     * remaining expected obligations with no delivery information yet 
     * available (if "Settle deliveries separately" is true). */
    CIntSP physicalDefaultToCalculationDelay;
    /** Date when the recovery rate is determined if NOT settling deliveries
     * separately*/
    DateTime physicalCalculationDate;
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
    CBoolSP moreDeliveriesPending;

    // ... for cash settlement
    /** % of the notional being triggered for (cash-settled) protection, 
     * with 100% = (<contingent notional in the observation period when 
     * name defaulted> / < tranche width in percentage>) * 
     * (<name's notional> / <portfolio notional>) */
    CDoubleSP cashDefaultedNotionalFraction;
    /** Delay between credit event date and calculation date */
    CIntSP cashDefaultToCalculationDelay;
    /** Date when the recovery rate is determined */
    DateTime cashCalculationDate;
    /** % of the notional being considered for cash-settlement.
     * The fraction of the notional between cashNotionalFraction and 
     * defaultedNotionalFraction (ie, the difference between what was 
     * expected to settle in cash and what has actually been triggered for 
     * cash-protection) will be deemed recovered. */
    double cashNotionalFraction;
    /** Recovery rate override. Optional, hence CDoubleSP. If not
     * present will use the "recovery" field (in TrancheCreditEventOverride)
     * for both the Physical and Cash settlements */
    CDoubleSP cashRecovery;


    // Transient fields
    // Why produce these two overrides internally with the input data, rather 
    // than having them as an input? because both overrides derive from 
    // TrancheCreditEventOverride and therefore require the fields 
    // declared there - forcing the user to introduce this data twice (or three
    // times, since it is also present in this class) is anoying and 
    // error-prone. The drawback is that we need to update the fields in this
    // class whenever the fields in TranchePhysicalSettlementOverride or 
    // TrancheCashSettlementOverride change
    TranchePhysicalSettlementOverrideSP physicalSettlement;
    TrancheCashSettlementOverrideSP     cashSettlement;
};

DECLARE (TranchePartCashSettlementOverride);

DRLIB_END_NAMESPACE

#endif
