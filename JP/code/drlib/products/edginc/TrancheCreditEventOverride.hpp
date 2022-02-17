//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : TrancheCreditEventOverride.hpp
//
//   Description : Class used by CDO in order 
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
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHECREDITEVENTOVERRIDE_HPP
#define QLIB_TRANCHECREDITEVENTOVERRIDE_HPP


#include "edginc/Atomic.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/ITrancheCreditEventOverride.hpp" 
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);

/** Class used by credit instruments (e.g., CDO) in order to override
 * a name's default-related parameters. Used as a base class for other
 * overrides.
 * CAUTION: This class contains overrides for some of the credit event
 * parameters of a single name - However, it is incomplete in the sense
 * that the parameters being overriden are not present here and need to be
 * passed in as arguments whenever required. These parameters (e.g., 
 * creditEventDate and lastTriggerDate in methods getFeeLegCashFlows and 
 * getContingentLegCashFlows) should therefore be consistent with the override
 * represented by instances of this class */
class TrancheCreditEventOverride: public CObject,
                                  public virtual ITrancheCreditEventOverride
{
public:
    static CClassConstSP const TYPE;
    virtual ~TrancheCreditEventOverride();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** GetMarket implementation*/
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the event determination date with the information in this 
     * override - it can be an estimate or the actual eventDeterminationDate,
     * if set. If there is no valid eventDeterminationDate (ie, the default
     * cannot be triggered) an empty DateTime is returned */
    virtual const DateTime getEventDeterminationDate(
        const DateTime& creditEventDate,
        const DateTime& lastTriggerDate,
        IBadDayAdjusterConstSP bda) const;

protected:
    TrancheCreditEventOverride(CClassConstSP clazz);

    /** Constructor with all fields passed in */
    TrancheCreditEventOverride(CClassConstSP clazz,
                               CIntSP triggerDelay,
                               const DateTime& eventDeterminationDate,
                               CDoubleSP recovery,
                               const DateTime& valueDate);

    // Fields
    // WARNING: any fields changed/added/removed here should also be updated
    // in TranchePartCashSettlementOverride, and the constructor in this class 
    // that takes all fields passed in should also be updated accordingly
    /** Delay between default and eventDeterminationDate, in days */
    CIntSP triggerDelay;

    /** Date when fees and accrual stop */
    DateTime eventDeterminationDate;  

    /** Recovery rate override. Optional, hence CDoubleSP */
    CDoubleSP recovery;

    /** Value date */
    DateTime valueDate;

private:
    // For reflection
    static void load (CClassSP& clazz);
    TrancheCreditEventOverride();
};

DECLARE (TrancheCreditEventOverride);

DRLIB_END_NAMESPACE

#endif
