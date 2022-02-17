/**
 *
 */

#ifndef EDG_LegalBasisAdditiveParallelTweak_H
#define EDG_LegalBasisAdditiveParallelTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a bulk shift to the additive coefficients of the Legal Basis
 * adjustment done to a CDSParSpreadCurve
 */

struct LegalBasisAdditiveParallelTweak: public virtual IScalarTweak {
    virtual ~LegalBasisAdditiveParallelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
