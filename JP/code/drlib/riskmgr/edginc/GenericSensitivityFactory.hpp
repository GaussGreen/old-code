/**
 * \file GenericSensitivityFactory.hpp
 *
 *
 * $History$
 *
 */

#ifndef DRLIB_GenericSensitivityFactory_H
#define DRLIB_GenericSensitivityFactory_H

#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Standardised SensitivityFactory, parameterised by the class of Sensitivity to
 * be constructed
 *
 * This avoids the need to create a specially cut-and-pasted "Factory"
 * class to go with every Sensitivity.  You can just put a GenericSensitivityFactory
 * directly in your SensitivityFactory::addSens() call:
 *
 * <PRE>
 *     void MySensitivity::load(...) {
 *         ...
 *         SensitivityFactory::addSens(
 *             ...,
 *             new GenericSensitivityFactory<MySensitivity>(),
 *             ...);
 *     }
 * </PRE>
 */

template <class SENSITIVITY>
struct GenericSensitivityFactory: public SensitivityFactory::IDefault,
                                  public SensitivityFactory::IScalar {

    virtual Sensitivity* createDefault() {
        return new SENSITIVITY(SENSITIVITY::DEFAULT_SHIFT);
    }

    virtual Sensitivity* createScalar(double shiftSize) {
        return new SENSITIVITY(shiftSize);
    }
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
