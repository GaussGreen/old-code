/**
 * @file BasketSpot.hpp
 */

#ifndef DRLIB_BasketSpot_H
#define DRLIB_BasketSpot_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "basket spot" as a property of market names.
 *
 * Cf. Spot, and see BasketDelta.
 */

struct RISKMGR_DLL BasketSpot: CObject {
    static CClassConstSP const TYPE;
    BasketSpot(); ~BasketSpot();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
