//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : ITrancheCreditEventOverride.hpp
//
//   Description : Class defining the interface used by credit instruments 
//                 (e.g., CDO) in order to override at instrument level a 
//                 name's default-related parameters.
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ITRANCHECREDITEVENTOVERRIDE_HPP
#define QLIB_ITRANCHECREDITEVENTOVERRIDE_HPP

#include "edginc/Atomic.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/AccrualPeriod.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

/* 
   Pyramid is currently not using the event handling functionality as it was
   designed to work. Whenever there is a default they keep rolling the market
   credit event date until the calculation date is reached, and at that time
   they mark the trade as triggered and set the calculation date. 
   However, if any trades on that name still have not been triggered, they 
   keep rolling the market credit event date. This results in defaults being 
   triggered before the actual market default date and needs to be accomodated
   in QLib - QLIB_REMOVE_EDD_VALIDATION_CDO is used for this purpose.
   Note:
   QLIB_REMOVE_EDD_VALIDATION_xx is also (potentially) defined in 
   ICreditEventOverrideName.hpp
*/
#define QLIB_REMOVE_EDD_VALIDATION_CDO


DRLIB_BEGIN_NAMESPACE


/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IRebateCalculator);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);

class MARKET_DLL ITrancheCreditEventOverride: public virtual IObject,
                                              public virtual IGetMarket
{

public:
    static CClassConstSP const TYPE;
    virtual ~ITrancheCreditEventOverride();

    /** Returns the contingent leg loss cashflows corresponding to a defaulted 
     * name, assuming that the default settles in one date.
     * This object contains the credit event parameters */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const = 0;

    /** Returns the fee leg loss cashflows corresponding to a defaulted name.
     * This object contains the credit event parameters.
     * CAUTION: The first two parameters are really the outputs of this method.
     * The original content of the SPs will be discarded. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        const double notional,
        const double recoveryRate,
        const DateTime& lastTriggerDate,
        const DateTime& creditEventDate,
        IBadDayAdjusterConstSP bda) const = 0;

    /* Returns the overall recovery rate accross all settlements, such that
     * total loss = overall defaulted notional * (1-overall recovery rate) */
    virtual double getOverallRecoveryRate(double recoveryRate) const = 0;

    /* Returns the overall defaulted notional accross all settlements, such that
     * total loss = overall defaulted notional * overall recovery rate */
    virtual double getOverallDefaultedNotionalFraction() const = 0;


    //--------------------------------------------------------------------------
    //                      Some utility static methods
    // Counter-intuitively, they live in this inteface class because they are
    // related to tranche credit event overrides and have to be in a file in the
    // market folder to be accessible where required
    //--------------------------------------------------------------------------
    /** Static method that rolls the originalDate passed in and adjusts it as
     * required. See comments within the method for further description.
     * Why is this method here instead of being in TrancheCreditEventOverride,
     * where it would fit more naturally? Because TrancheCreditEventOverride.hpp
     * lives in the products folder, so it could not be included by market files
     * requiring access to this method. On the other hand, 
     * ITrancheCreditEventOverride.hpp lives in the market folder (for this very
     * reason), so it does the job. */
    static DateTime rollAndAdjustDate(const DateTime& originalDate,
                                      CIntConstSP offset,
                                      const DateTime& valueDate,
                                      const DateTime& minDate,
                                      const DateTime& maxDate,
                                      IBadDayAdjusterConstSP bda);

    /** Static method that rolls the originalDate passed in and adjusts it as 
     * required. */
    static DateTime rollAndAdjustDate(const DateTime& originalDate,
                                      CIntConstSP offset,
                                      const DateTime& valueDate,
                                      const DateTime& minDate,
                                      IBadDayAdjusterConstSP bda);

private:
    /** For reflection */
    static void load (CClassSP& clazz);
};

typedef smartConstPtr<ITrancheCreditEventOverride> ITrancheCreditEventOverrideConstSP;
typedef smartPtr<ITrancheCreditEventOverride> ITrancheCreditEventOverrideSP;
typedef array<ITrancheCreditEventOverrideSP, ITrancheCreditEventOverride> ITrancheCreditEventOverrideArray;
typedef smartPtr<ITrancheCreditEventOverrideArray> ITrancheCreditEventOverrideArraySP;
typedef smartConstPtr<ITrancheCreditEventOverrideArray> ITrancheCreditEventOverrideArrayConstSP;


DRLIB_END_NAMESPACE

#endif
