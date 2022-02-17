/**
 * @file CorrSwapBasisAdjAbsolute.hpp
 */

#ifndef DRLIB_CorrSwapBasisAdjAbsolute_H
#define DRLIB_CorrSwapBasisAdjAbsolute_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "corr swap basis adjustment" as a property of market names.
 */

struct RISKMGR_DLL CorrSwapBasisAdjAbsolute: CObject {
    static CClassConstSP const TYPE;
    CorrSwapBasisAdjAbsolute(); ~CorrSwapBasisAdjAbsolute();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
