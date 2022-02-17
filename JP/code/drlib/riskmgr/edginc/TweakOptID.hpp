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

#ifndef EDR_TWEAK_OPT_ID_HPP
#define EDR_TWEAK_OPT_ID_HPP

#include "edginc/TweakID.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface that specifies all the information needed to apply a tweak which
    can be undone in a potentially optimal manner */
class RISKMGR_DLL ITweakOptID: public virtual ITweakID {
public:
    ITweakOptID(); // in SensMgr.cpp

    virtual ~ITweakOptID(); // in SensMgr.cpp

    /** whether the object (which supports being tweaked
        by this type of sens control) implements a restorable shift */
    virtual bool restorableShift(IObjectConstSP obj) const = 0;

    /** Restores the object (which supports being tweaked
        by this type of ITweakOptID) to its original form */
    virtual void restore(IObjectSP obj) = 0;
};

FORWARD_DECLARE(ITweakOptID)

DRLIB_END_NAMESPACE

#endif

