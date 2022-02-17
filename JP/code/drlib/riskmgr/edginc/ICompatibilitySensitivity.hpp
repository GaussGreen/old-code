/**
 * @file ICompatibilitySensitivity.hpp
 */

#ifndef QLIB_ICompatibilitySensitivity_H
#define QLIB_ICompatibilitySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICompatibilitySensitivity)
FORWARD_DECLARE(AtomicHypothesis)
FORWARD_DECLARE(AbstractPropertyTweakHypothesis)

/**
 * Interim support for model/instrument mechanisms based on
 * Control::getCurrentSensitivity() to work with "declarative" sensitivities.
 *
 * See IRiskQuantityFactory for an overview of the "declarative" sensitivities
 * framework, and the note on RiskQuantityEvaluator::storeResults() in
 * IRiskQuantityFactory.cpp for an account of the issue motivating the
 * existence of this interface.
 */

class RISKMGR_DLL ICompatibilitySensitivity: public virtual IObject {

    friend class HypothesisTree;

    virtual void setCurrentHypothesis(AtomicHypothesisArraySP history,
                                      AbstractPropertyTweakHypothesisConstSP current,
                                      double oldValue) = 0;

public:

    ICompatibilitySensitivity();
    ~ICompatibilitySensitivity();

    static CClassConstSP const TYPE;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
