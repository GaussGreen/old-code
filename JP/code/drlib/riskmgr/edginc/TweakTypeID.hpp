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

#ifndef EDR_TWEAK_TYPE_ID_HPP
#define EDR_TWEAK_TYPE_ID_HPP
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface defining method which specifies the type of object to look
    for when doing tweak related operations */
class RISKMGR_DLL ITweakTypeID {
public:
    ITweakTypeID(); // in SensMgr.cpp
    virtual ~ITweakTypeID(); // in SensMgr.cpp

    /** returns the interface identifying type of objects to look for */
    virtual CClassConstSP shiftInterface() const = 0;

};

DRLIB_END_NAMESPACE

#endif

