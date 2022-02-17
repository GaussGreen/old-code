//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Additive.hpp
//
//   Description : Marker interface for sensitivities to show if they can
//                 be added together
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef ADDITIVE_HPP
#define ADDITIVE_HPP

#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

/** Marker interface for sensitivities to show if they can be added together */
class RISKMGR_DLL Additive {
public:
    static CClassConstSP const TYPE;

    virtual ~Additive();

protected:
    Additive();
};


DRLIB_END_NAMESPACE
#endif
