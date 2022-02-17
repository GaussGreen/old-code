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

#ifndef EDR_TWEAK_NAME_LIST_ID_HPP
#define EDR_TWEAK_NAME_LIST_ID_HPP
#include "edginc/OutputName.hpp"
#include "edginc/TweakTypeID.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface defining methods which allows a list of names of objects that
    implement a specific interface to be returned. For example, if you want
    to calculate delta, then you want a list of underlyings to tweak - so
    your greek should implement this interface so that the infrastructure can
    find the objects to tweak */
class RISKMGR_DLL ITweakNameListID: public virtual ITweakTypeID {
public:
    ITweakNameListID(); // in SensMgr.cpp

    virtual ~ITweakNameListID(); // in SensMgr.cpp

    /** Appends the name(s) of the supplied object with respect to
        this tweak to the supplied list. Could consider being const,
        but has large downstream consequences. */
    virtual void appendName(OutputNameArray& namesList,
                            IObjectConstSP   obj) = 0;
};

DRLIB_END_NAMESPACE

#endif

