//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Author      : Mark A Robson
//
//   Date        : 25 July 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VECTOR_TWEAK_HPP
#define EDR_VECTOR_TWEAK_HPP
#include "edginc/ScalarTweak.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for tweaks that have a well defined shift size which is a 
    double together with the idea that it should be applied at a particular
    expiry */
class RISKMGR_DLL IVectorTweak: public virtual IScalarTweak{ // possibly dubious inheritance
public:
    IVectorTweak();  // In SensControl.cpp
    virtual ~IVectorTweak();  // In SensControl.cpp

    /** Returns the expiry which is should be tweaked */
    virtual ExpiryConstSP getExpiry() const = 0;
};

DRLIB_END_NAMESPACE

#endif

