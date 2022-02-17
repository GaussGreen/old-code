//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : ICreditEventOverrideName.hpp
//
//   Description : Interface used by credit instruments in order to override
//                 (at instrument level) a name's default-related parameters.
//                 Due to design changes, this interface now derives from
//                 ITrancheCreditEventOverride and defines just one new 
//                 method, so conceptually it is almost identical to 
//                 ITrancheCreditEventOverride.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICREDITEVENTOVERRIDENAME_HPP
#define QLIB_ICREDITEVENTOVERRIDENAME_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/ITrancheCreditEventOverride.hpp"

/* 
   Pyramid is currently not using the event handling functionality as it was
   designed to work. Whenever there is a default they keep rolling the market
   credit event date until the calculation date is reached, and at that time
   they mark the trade as triggered and set the calculation date. 
   However, if any trades on that name still have not been triggered, they 
   keep rolling the market credit event date. This results in defaults being 
   triggered before the actual market default date and needs to be accomodated
   in QLib - QLIB_REMOVE_EDD_VALIDATION_CDS is used for this purpose.
   Note:
   QLIB_REMOVE_EDD_VALIDATION_xx is also (potentially) defined in 
   ITrancheCreditEventOverride.hpp
*/
#define QLIB_REMOVE_EDD_VALIDATION_CDS

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ICreditEventOverrideName : 
    public virtual ITrancheCreditEventOverride 
{
public:
    static CClassConstSP const TYPE;
    virtual ~ICreditEventOverrideName();

    /** Name of this market object */
    virtual string getName() const = 0;

private:
    /** For reflection */
    static void load (CClassSP& clazz);
};

typedef smartConstPtr<ICreditEventOverrideName> ICreditEventOverrideNameConstSP;
typedef smartPtr<ICreditEventOverrideName> ICreditEventOverrideNameSP;
typedef array<ICreditEventOverrideNameSP, ICreditEventOverrideName> ICreditEventOverrideNameArray;
typedef smartPtr<ICreditEventOverrideNameArray> ICreditEventOverrideNameArraySP;
typedef smartConstPtr<ICreditEventOverrideNameArray> ICreditEventOverrideNameArrayConstSP;


DRLIB_END_NAMESPACE

#endif
