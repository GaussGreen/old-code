/**
 *
 */

#ifndef EDG_CompositeVegaParallelTweak_H
#define EDG_CompositeVegaParallelTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a bulk shift to the composite volatilities
 */

struct CompositeVegaParallelTweak: public virtual IScalarTweak {
    virtual ~CompositeVegaParallelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
