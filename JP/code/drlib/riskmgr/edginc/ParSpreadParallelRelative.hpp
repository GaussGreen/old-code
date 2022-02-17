/**
 * @file ParSpreadParallelRelative.hpp
 */

#ifndef DRLIB_ParSpreadParallelRelative_H
#define DRLIB_ParSpreadParallelRelative_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "par spread level" as a property of market names.
 *
 * Similar to ParSpreadParallel, but the shift is relative.
 */

struct RISKMGR_DLL ParSpreadParallelRelative: CObject {
    static CClassConstSP const TYPE;
    ParSpreadParallelRelative(); ~ParSpreadParallelRelative();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

#endif
