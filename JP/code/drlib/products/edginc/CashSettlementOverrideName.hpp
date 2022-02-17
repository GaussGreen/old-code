//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CashSettlementOverrideName.hpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled in cash.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CASHSETTLEMENTOVERRIDENAME_HPP
#define QLIB_CASHSETTLEMENTOVERRIDENAME_HPP

#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);

/** Class used by credit instruments (e.g., CIS) in order to override
 * a name's default-related parameters, when the default is settled
 * in cash. */
class PRODUCTS_DLL CashSettlementOverrideName: public DetailedCreditEventOverrideName
{
public:
    static CClassConstSP const TYPE;
    virtual ~CashSettlementOverrideName();

    /** Public constructor, with all fields passed in */
    CashSettlementOverrideName(ICDSParSpreadsWrapper name,
                               DateTime eventDeterminationDate,
                               CIntSP triggerDelay,
                               CDoubleSP recovery,
                               double notionalFraction,
                               CIntSP defaultToSettlementDelay,
                               DateTime calculationDate,
                               int calculationToSettlementDelay);

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

    /* Returns the overall recovery rate accross all settlements, such that
     * total loss = overall defaulted notional * (1-overall recovery rate) */
    virtual double getOverallRecoveryRate(double recoveryRate) const;

    /** Returns the percentage of the notional that has been triggered 
     * for protection */
    virtual double getOverallDefaultedNotionalFraction() const;

   /** Returns the accrual pay date. For cash settlements this is the same
    * as the settlement date. */
   virtual const DateTime getAccrualPayDate(const DateTime& creditEventDate,
                                            const DateTime& lastTriggerDate,
                                            const DateTime& valueDate,
                                            IBadDayAdjusterConstSP bda) const;

private:
    /** Returns the date of the cash settlement based on the information in 
     * this override - it can be an estimate or the actual date if known */
    const DateTime getSettlementDate(
        const DateTime& creditEventDate,
        const DateTime& lastTriggerDate,
        const DateTime& valueDate,
        IBadDayAdjusterConstSP bda) const;

    // For reflection
    static void load (CClassSP& clazz);
    CashSettlementOverrideName();
    static IObject* defaultCashSettlementOverrideName();

    // Fields
    // WARNING: any fields changed/added/removed here should also be updated
    // in PartCashSettlementOverrideName, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly
    double notionalFraction;          // % of the notional being cash-settled
    CIntSP defaultToSettlementDelay;  // Delay between credit event and settlement
    DateTime calculationDate;         // Date when the recovery rate is determined
    int calculationToSettlementDelay; // Delay between calculation and settlement
};

DECLARE (CashSettlementOverrideName);

DRLIB_END_NAMESPACE

#endif
