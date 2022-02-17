/**
 * @file IRMeanReversionParallelProperty.hpp
 */

#ifndef DRLIB_IRMeanReversionParallelProperty_H
#define DRLIB_IRMeanReversionParallelProperty_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "interest rate mean reversion" as a property of market names.
 *
 * Cf. Spot, and see IRMeanReversionParallel.cpp.
 */

struct RISKMGR_DLL IRMeanReversionParallelProperty: CObject {
    static CClassConstSP const TYPE;
    IRMeanReversionParallelProperty(); ~IRMeanReversionParallelProperty();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
