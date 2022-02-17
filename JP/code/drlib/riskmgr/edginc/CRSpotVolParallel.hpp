/**
 * @file CRSpotVolParallel.hpp
 */

#ifndef DRLIB_CRSpotVolParallel_H
#define DRLIB_CRSpotVolParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class.
 *
 * Cf. Spot, and see CRSpotVegaParallel.cpp.
 */

struct RISKMGR_DLL CRSpotVolParallel: CObject {
    static CClassConstSP const TYPE;
    CRSpotVolParallel(); ~CRSpotVolParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
