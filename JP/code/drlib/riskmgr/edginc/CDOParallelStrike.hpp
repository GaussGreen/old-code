//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CDOParallelStrike.hpp
//
//   Description : Tag class for parallel shift in CDO strikes
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#ifndef DRLIB_CDO_PARALLEL_SHIFT_H
#define DRLIB_CDO_PARALLEL_SHIFT_H

#include "edginc/Void.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct RISKMGR_DLL CDOParallelStrike: CObject {
    static CClassConstSP const TYPE;
    CDOParallelStrike(); ~CDOParallelStrike();
    typedef Void Qualifier;
    //Is continuous, i.e. tweaks to it can be made arbitrarily small
    enum { discrete = 0};
};

DRLIB_END_NAMESPACE

#endif
