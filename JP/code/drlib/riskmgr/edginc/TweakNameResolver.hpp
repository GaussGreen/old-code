//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Author      : Mark A Robson
//
//   Date        : 31 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_TWEAK_NAME_RESOLVER_HPP
#define EDR_TWEAK_NAME_RESOLVER_HPP
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface defining methods which allows a specific type of object to be 
    identified by name  */
class RISKMGR_DLL ITweakNameResolver {
public:
    ITweakNameResolver(); // in SensMgr.cpp

    virtual ~ITweakNameResolver(); // in SensMgr.cpp

    /** returns the name identifying the market data to be shifted. Note that
        typically the name is the name of a piece of *market* data but does 
        not have to be (eg could be an instrument piece of data) */
    virtual OutputNameConstSP getMarketDataName() const = 0;

    /** does the given name match the name identifying the market data
        to be shifted. Could consider being const, but has large downstream
        consequences. */
    virtual bool nameMatches(const OutputName& name,
                             IObjectConstSP    obj) = 0;
};

DRLIB_END_NAMESPACE

#endif

