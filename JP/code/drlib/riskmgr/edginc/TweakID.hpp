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

#ifndef EDR_TWEAK_ID_HPP
#define EDR_TWEAK_ID_HPP
#include "edginc/TweakNameResolution.hpp"

DRLIB_BEGIN_NAMESPACE
class ITweakNameResolver;

/** Interface that specifies all the information needed to apply a tweak  */
class RISKMGR_DLL ITweakID: public virtual ITweakNameResolution {
public:

    ITweakID();  // in SensMgr.cpp

    virtual ~ITweakID(); // in SensMgr.cpp

    /** resets internally stored values associated with tweaking. This is 
        because this type of
        object (or more precisely derived types) can store information about
        what has been tweaked (eg initial values) */
    virtual void reset() = 0;

    /** Shifts the object (which supports being tweaked
        by this type) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    virtual bool shift(IObjectSP obj) = 0;
};

DRLIB_END_NAMESPACE

#endif

