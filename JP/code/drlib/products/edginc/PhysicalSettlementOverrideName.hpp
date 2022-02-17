//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Filename    : PhysicalSettlementOverrideName.hpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled physically.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_PHYSICALSETTLEMENTOVERRIDENAME_HPP
#define QLIB_PHYSICALSETTLEMENTOVERRIDENAME_HPP

#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/NoticeOfPhysicalSettlement.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(DeliveryDetails);
FORWARD_DECLARE(IBadDayAdjuster);

/** Class used by credit instruments (e.g., CIS) in order to override
 * a name's default-related parameters, when the default is settled
 * physically. */
class PhysicalSettlementOverrideName: public DetailedCreditEventOverrideName {
public:
    static CClassConstSP const TYPE;
    virtual ~PhysicalSettlementOverrideName();

    /** Public constructor, with all fields passed in */
    PhysicalSettlementOverrideName(ICDSParSpreadsWrapper name,
                                   DateTime eventDeterminationDate,
                                   CIntSP triggerDelay,
                                   CDoubleSP recovery,
                                   NoticeOfPhysicalSettlementSP NoPS,
                                   CIntSP defaultToSettlementDelay,
                                   DeliveryDetailsArraySP deliveries,
                                   double totalExpectedNotional);

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

    /** Returns the accrual pay date. For physical settlements this is the date
        of the first delivery */
    virtual const DateTime getAccrualPayDate(const DateTime& creditEventDate,
                                             const DateTime& lastTriggerDate,
                                             const DateTime& valueDate,
                                             IBadDayAdjusterConstSP bda) const;

private:
    /** Returns the overall notional fraction we have deliveries scheduled for */
    double getNotionalScheduledForDelivery() const;

    /* Returns the fraction of the notional which is expected to be settled
     * but still has no settlement (delivery) date */
    double getUnscheduledNotional() const;

    /** Returns the estimated settlement date based on the 
     * defaultToSettlementDelay */
    const DateTime estimateSettlementDate(const DateTime& creditEventDate,
                                          const DateTime& lastTriggerDate,
                                          const DateTime& valueDate,
                                          IBadDayAdjusterConstSP bda) const;

    /** Returns the overall loss resulting from the settlement of this name */
    double getOverallLoss(const double notional,
                          const double recoveryRate) const;

    // For reflection
    static void load (CClassSP& clazz);
    PhysicalSettlementOverrideName();
    static IObject* defaultPhysicalSettlementOverrideName();

    // Fields
    // WARNING: any fields changed/added/removed here should also be updated
    // in PartCashSettlementOverrideName, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly
    NoticeOfPhysicalSettlementSP NoPS; // Specification of the amounts being settled
    CIntSP defaultToSettlementDelay;   // Delay between credit event and settlement
    DeliveryDetailsArraySP deliveries; // Array with the details of all deliveries

    // Transient field
    double totalExpectedNotional; // Fraction of the notional expected to settle physically
};

DECLARE(PhysicalSettlementOverrideName);


DRLIB_END_NAMESPACE

#endif
