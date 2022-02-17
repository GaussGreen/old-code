/**
 *
 */

#ifndef EDG_LegalBasisMultiplierRecoveryTweak_H
#define EDG_LegalBasisMultiplierRecoveryTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a shift to the multiplier recovery coefficient of the Legal Basis
 * adjustment done to a CDSParSpreadCurve
 */

struct LegalBasisMultiplierRecoveryTweak: public virtual IScalarTweak {
    virtual ~LegalBasisMultiplierRecoveryTweak() {}
};

DRLIB_END_NAMESPACE

#endif
