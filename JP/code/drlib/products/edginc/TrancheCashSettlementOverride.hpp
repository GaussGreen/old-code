//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TrancheCashSettlementOverride.hpp
//
//   Description : Class used by CDO in order 
//                 to override a name's default-related parameters, when
//                 the default is settled in cash.
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHECASHSETTLEMENTOVERRIDE_HPP
#define QLIB_TRANCHECASHSETTLEMENTOVERRIDE_HPP

#include "edginc/TrancheCreditEventOverride.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);

/** Class used by credit instruments (e.g., CIS) in order to override
 * a name's default-related parameters, when the default is settled
 * in cash. */
class TrancheCashSettlementOverride: public TrancheCreditEventOverride {
public:
    static CClassConstSP const TYPE;
    virtual ~TrancheCashSettlementOverride();


    /** Public constructor, with all fields passed in */
    TrancheCashSettlementOverride(CIntSP triggerDelay,
                                  DateTime eventDeterminationDate,
                                  CDoubleSP recovery,
                                  DateTime valueDate,
                                  CDoubleSP defaultedNotionalFractionSP,
                                  CIntSP defaultToCalculationDelay,
                                  DateTime calculationDate,
                                  double notionalFraction);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();
    
    /** Returns the contingent leg loss cashflows corresponding to a defaulted 
     * name, assuming that the default settles in one date. */
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
    /** Returns the date of the calculation date based on the information in 
     * this override - it can be an estimate or the actual date if known */
    DateTime getCalculationDate(const DateTime& creditEventDate,
                                const DateTime& lastTriggerDate,
                                IBadDayAdjusterConstSP bda) const;

    /* Returns the loss amount associated to this name's default */
    double getLossAmount(const double nameNotional,
                         const double recoveryRate) const;

    // For reflection
    static void load (CClassSP& clazz);
    TrancheCashSettlementOverride();
    static IObject* defaultTrancheCashSettlementOverride();

    // Fields
    // WARNING: any fields changed/added/removed here should also be updated
    // in TranchePartCashSettlementOverride, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly

    /** % of the notional being cash-settled, with 
     * 100% = (<contingent notional in the observation period when name 
     * defaulted> / < tranche width in percentage>) * 
     * (<name's notional> / <portfolio notional>) */
    double defaultedNotionalFraction;

    /** Delay between credit event date and calculation date */
    CIntSP defaultToCalculationDelay;

    /** Date when the recovery rate is determined */
    DateTime calculationDate;

    // Transient
    /** % of the notional being considered for cash-settlement. It will typically
     * be 100% unless this override is part of a physical settlement with part
     * cash.
     * The fraction of the notional between notionalFraction and 
     * defaultedNotionalFraction (ie, the difference between what has 
     * been triggered for protection and what could have been, ie, the fraction 
     * that has NOT been triggered for protection) will be deemed recovered. */
    double notionalFraction;
};

DECLARE (TrancheCashSettlementOverride);

DRLIB_END_NAMESPACE

#endif
