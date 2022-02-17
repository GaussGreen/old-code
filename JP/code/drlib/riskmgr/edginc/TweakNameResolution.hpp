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

#ifndef EDR_TWEAK_NAME_RESOLUTION_HPP
#define EDR_TWEAK_NAME_RESOLUTION_HPP
#include "edginc/TweakTypeID.hpp"

DRLIB_BEGIN_NAMESPACE
class ITweakNameResolver;

/** Interface defining how objects can be identified by name (as well as
    their type via being derived from ITweakTypeID) */
class RISKMGR_DLL ITweakNameResolution: public virtual ITweakTypeID {
public:
    ITweakNameResolution(); // in SensMgr.cpp

    virtual ~ITweakNameResolution(); // in SensMgr.cpp

    /** Returns either null (=> all names match) or an object which controls
        when an object matches. Do not delete returned pointer. (Could
        be made into a smart pointer if required) */
    virtual ITweakNameResolver* nameResolver() = 0;
};

DRLIB_END_NAMESPACE

#endif

