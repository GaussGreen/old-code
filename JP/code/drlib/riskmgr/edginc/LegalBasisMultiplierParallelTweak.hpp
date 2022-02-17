/**
 *
 */

#ifndef EDG_LegalBasisMultiplierParallelTweak_H
#define EDG_LegalBasisMultiplierParallelTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE


/**
 * Parameterises a bulk shift to the multiplier coefficients of the Legal Basis
 * adjustment done to a CDSParSpreadCurve
 */

struct LegalBasisMultiplierParallelTweak: public virtual IScalarTweak {
    virtual ~LegalBasisMultiplierParallelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
