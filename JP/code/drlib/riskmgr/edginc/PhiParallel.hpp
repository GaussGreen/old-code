#ifndef EDG_PHI_PARALLEL_HPP
#define EDG_PHI_PARALLEL_HPP

#include "edginc/Sensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Default shift size for PhiParallel sensitivity
 */

const double PhiParallel_defaultShift = 0.01;

/**
 * A newly built PhiParallel Sensitivity
 */

SensitivitySP PhiParallel_SP(double shiftSize, bool onlyEqEq);

DRLIB_END_NAMESPACE

#endif
