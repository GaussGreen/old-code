/**
 * @file CRMeanReversionParallelProperty.hpp
 */

#ifndef DRLIB_CRMeanReversionParallelProperty_H
#define DRLIB_CRMeanReversionParallelProperty_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "credit mean reversion" as a property of market names.
 *
 * Cf. Spot, and see CRMeanReversionParallel.cpp.
 */

struct RISKMGR_DLL CRMeanReversionParallelProperty: CObject {
    static CClassConstSP const TYPE;
    CRMeanReversionParallelProperty(); ~CRMeanReversionParallelProperty();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
