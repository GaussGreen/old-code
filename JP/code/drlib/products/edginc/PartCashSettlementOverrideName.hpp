//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : PartCashSettlementOverrideName.hpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled physically but there is a partial
//                 cash settlement.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_PARTCASHSETTLEMENTOVERRIDENAME_HPP
#define QLIB_PARTCASHSETTLEMENTOVERRIDENAME_HPP

#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(NoticeOfPhysicalSettlement);
FORWARD_DECLARE(DeliveryDetails);
FORWARD_DECLARE(PhysicalSettlementOverrideName);
FORWARD_DECLARE(CashSettlementOverrideName);


/** Class used by credit instruments (e.g., CIS) in order to override
 * a name's default-related parameters, when the default is settled
 * physically but there is a partial cash settlement. */
class PRODUCTS_DLL PartCashSettlementOverrideName : public DetailedCreditEventOverrideName 
{
public:
    static CClassConstSP const TYPE;
    virtual ~PartCashSettlementOverrideName();

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
    // For reflection
    static void load (CClassSP& clazz);
    PartCashSettlementOverrideName();
    static IObject* defaultPartCashSettlementOverrideName();

    // Fields
    // ... for physical settlement
    NoticeOfPhysicalSettlementSP NoPS;  // Specification of the amounts being settled
    CIntSP defaultToSettlementDelayPhys;// Delay between credit event and settlement
    DeliveryDetailsArraySP deliveries;  // Array with the details of all deliveries

    // ... for the cash part
    double notionalFractionCash;        // % of the notional being cash-settled
    CIntSP defaultToSettlementDelayCash;// Delay between credit event and settlement
    DateTime calculationDate;           // Date when the recovery rate is determined
    int calculationToSettlementDelay;   // Delay between calculation and settlement


    // Transient fields
    // Why produce these two overrides internally with the input data, rather 
    // than having them as an input? because both overrides derive from 
    // DetailedCreditEventOverrideName and therefore require the fields 
    // declared there - forcing the user to introduce this data twice (or three
    // times, since it is also present in this class) is anoying and 
    // error-prone. The drawback is that we need to update the fields in this
    // class whenever the fields in PhysicalSettlementOverrideName or 
    // CashSettlementOverrideName change
    PhysicalSettlementOverrideNameSP physicalSettlement;
    CashSettlementOverrideNameSP     cashSettlement;
};

DECLARE (PartCashSettlementOverrideName);

DRLIB_END_NAMESPACE

#endif
