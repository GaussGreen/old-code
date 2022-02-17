/**
 * @file CRVolParallel.hpp
 */

#ifndef DRLIB_CRVolParallel_H
#define DRLIB_CRVolParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "credit vol level" as a property of market names.
 *
 * Cf. Spot, and see CRVegaParallel.cpp.
 */

struct RISKMGR_DLL CRVolParallel: CObject {
    static CClassConstSP const TYPE;
    CRVolParallel(); ~CRVolParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
