/**
 * @file CreditIndexBasisParallel.hpp
 */

#ifndef DRLIB_CreditIndexBasisParallel_H
#define DRLIB_CreditIndexBasisParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "par spread level" as a property of market names.
 *
 * Cf. Spot, and see ParSpreadRhoParallel.cpp.
 */

struct RISKMGR_DLL CreditIndexBasisParallel: CObject {
    static CClassConstSP const TYPE;
    CreditIndexBasisParallel(); ~CreditIndexBasisParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
