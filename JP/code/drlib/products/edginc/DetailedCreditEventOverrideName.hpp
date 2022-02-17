//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : DetailedCreditEventOverrideName.hpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters. Used
//                 as a base class for other overrides.
//                 CAUTION: This class contains overrides for some of the credit 
//                 event parameters of a single name - However, it is incomplete 
//                 in the sense that the parameters being overriden are not 
//                 present here and need to be passed in as arguments whenever 
//                 required. These parameters (e.g., creditEventDate and 
//                 lastTriggerDate in methods getFeeLegCashFlows and 
//                 getContingentLegCashFlows) should therefore be consistent 
//                 with the override represented by instances of this class.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_DETAILEDCREDITEVENTOVERRIDENAME_HPP
#define QLIB_DETAILEDCREDITEVENTOVERRIDENAME_HPP

#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IForwardRatePricer);

/** Class used by credit instruments (e.g., CIS) in order to override
 * a name's default-related parameters. Used as a base class for other
 * overrides.
 * CAUTION: This class contains overrides for some of the credit event
 * parameters of a single name - However, it is incomplete in the sense
 * that the parameters being overriden are not present here and need to be
 * passed in as arguments whenever required. These parameters (e.g., 
 * creditEventDate and lastTriggerDate in methods getFeeLegCashFlows and 
 * getContingentLegCashFlows) should therefore be consistent with the override
 * represented by instances of this class */
class PRODUCTS_DLL DetailedCreditEventOverrideName :
    public CObject,
    public virtual ICreditEventOverrideName
{

public:
    static CClassConstSP const TYPE;
    virtual ~DetailedCreditEventOverrideName();

    /** Name of the curve associated to this override */
    string getName() const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

     /** GetMarket implementation*/
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the fee leg loss cashflows corresponding to a defaulted name,
     * assuming that the default settles in one date. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const;

   /** Returns the date when the accrual and (1-R) are paid */
   virtual const DateTime getAccrualPayDate(const DateTime& creditEventDate,
                                            const DateTime& lastTriggerDate,
                                            const DateTime& valueDate,
                                            IBadDayAdjusterConstSP bda) const = 0;

protected:
    DetailedCreditEventOverrideName(CClassConstSP clazz);

    /** Constructor with all fields passed in */
    DetailedCreditEventOverrideName(CClassConstSP clazz,
                                    ICDSParSpreadsWrapper name,
                                    DateTime eventDeterminationDate,
                                    CIntSP triggerDelay,
                                    CDoubleSP recovery);

    /** Returns the event determination date with the information in this 
     * override - it can be an estimate or the actual eventDeterminationDate,
     * if set. If there is no valid eventDeterminationDate (ie, the default
     * cannot be triggered) an empty DateTime is returned */
    const DateTime getEventDeterminationDate(
        const DateTime& creditEventDate,
        const DateTime& lastTriggerDate,
        IBadDayAdjusterConstSP bda) const;

    /* Returns the overriden recovery rate */
    double getOverridenRecoveryRate(double recoveryRate) const;

    // Fields
    // WARNING: any fields changed/added/removed here should also be updated
    // in PartCashSettlementOverrideName, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly
    ICDSParSpreadsWrapper name;       // Name of the curve the override applies to
    DateTime eventDeterminationDate;  // When fees and accrual stop
    CIntSP triggerDelay; // Delay between default and eventDeterminationDate, in days
    CDoubleSP recovery;               // Optional, hence CDoubleSP
    DateTime valueDate;

private:
    // For reflection
    static void load (CClassSP& clazz);
    DetailedCreditEventOverrideName();
};

DECLARE (DetailedCreditEventOverrideName);

DRLIB_END_NAMESPACE

#endif
