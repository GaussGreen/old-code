/**
 * @file EnergyFuturesCurveParallel.hpp
 */

#ifndef DRLIB_EnergyFuturesCurveParallel_H
#define DRLIB_EnergyFuturesCurveParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "par spread level" as a property of market names.
 *
 * Similar to ParSpreadParallel, but the shift is relative.
 */

struct RISKMGR_DLL EnergyFuturesCurveParallel: CObject {
    static CClassConstSP const TYPE;
    EnergyFuturesCurveParallel(); ~EnergyFuturesCurveParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

#endif
