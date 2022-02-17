//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ICreditEventOverride.hpp
//
//   Description : Interface used by credit instruments (e.g., CIS) in order 
//                 to override (at instrument level) the names' default-related 
//                 parameters.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICREDITEVENTOVERRIDE_HPP
#define QLIB_ICREDITEVENTOVERRIDE_HPP

#include "edginc/Object.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICreditEventOverrideName);


class MARKET_DLL ICreditEventOverride: public virtual IObject,
                                       public virtual IGetMarket
{

public:
    static CClassConstSP const TYPE;
    virtual ~ICreditEventOverride();

    // Returns the override information for the supplied name. If no override
    // information is available returns ICreditEventOverrideNameSP()
    virtual ICreditEventOverrideNameSP getOverrideForName(
        const string& name) const = 0;

private:
    // For reflection
    static void load (CClassSP& clazz);
};

DECLARE(ICreditEventOverride);

DRLIB_END_NAMESPACE

#endif
