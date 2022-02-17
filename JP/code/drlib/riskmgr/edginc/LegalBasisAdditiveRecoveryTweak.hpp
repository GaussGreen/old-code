/**
 *
 */

#ifndef EDG_LegalBasisAdditiveRecoveryTweak_H
#define EDG_LegalBasisAdditiveRecoveryTweak_H

#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a shift to the additive recovery coefficient of the Legal Basis
 * adjustment done to a CDSParSpreadCurve
 */

struct LegalBasisAdditiveRecoveryTweak: public virtual IScalarTweak {
    virtual ~LegalBasisAdditiveRecoveryTweak() {}
};

DRLIB_END_NAMESPACE

#endif
