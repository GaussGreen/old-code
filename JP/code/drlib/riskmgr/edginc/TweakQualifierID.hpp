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

#ifndef EDR_TWEAK_QUALIFIER_ID_HPP
#define EDR_TWEAK_QUALIFIER_ID_HPP
#include "edginc/TweakNameResolution.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface defining methods which allow a 'qualifier' to be returned for an
    object of a specific type. Here 'qualifier' can be, for example, benchmark
    dates.  */
class RISKMGR_DLL ITweakQualifierID: public virtual ITweakNameResolution {
public:
    ITweakQualifierID(); // in SensMgr.cpp

    virtual ~ITweakQualifierID(); // in SensMgr.cpp

    /** Returns an object needed for qualifying sensitivity eg benchmark
        dates. Could consider being const, but has large downstream
        consequences. */
    virtual IObjectConstSP qualifier(IObjectConstSP obj) = 0;
};

DRLIB_END_NAMESPACE

#endif

