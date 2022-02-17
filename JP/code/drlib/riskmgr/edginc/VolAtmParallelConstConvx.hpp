/**
 * @file VolAtmParallelConstConvx.hpp
 */

#ifndef DRLIB_VolAtmParallelConstConvx_H
#define DRLIB_VolAtmParallelConstConvx_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "general level of vol curve" as a property of
 * market names.
 */

struct RISKMGR_DLL VolAtmParallelConstConvx: CObject {
    static CClassConstSP const TYPE;
    VolAtmParallelConstConvx(); ~VolAtmParallelConstConvx();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
