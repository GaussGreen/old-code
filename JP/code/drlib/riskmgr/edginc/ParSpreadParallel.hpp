/**
 * @file ParSpreadParallel.hpp
 */

#ifndef DRLIB_ParSpreadParallel_H
#define DRLIB_ParSpreadParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "par spread level" as a property of market names.
 *
 * Cf. Spot, and see ParSpreadRhoParallel.cpp.
 */

struct RISKMGR_DLL ParSpreadParallel: CObject {
    static CClassConstSP const TYPE;
    ParSpreadParallel(); ~ParSpreadParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
