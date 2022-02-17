/**
 * @file CorrSwapSamplingAdjAbsolute.hpp
 */

#ifndef DRLIB_CorrSwapSamplingAdjAbsolute_H
#define DRLIB_CorrSwapSamplingAdjAbsolute_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "corr swap sampling adjustment" as a property of market names.
 */

struct RISKMGR_DLL CorrSwapSamplingAdjAbsolute: CObject {
    static CClassConstSP const TYPE;
    CorrSwapSamplingAdjAbsolute(); ~CorrSwapSamplingAdjAbsolute();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
