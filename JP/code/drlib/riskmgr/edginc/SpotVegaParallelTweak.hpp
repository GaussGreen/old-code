/**
 *
 */

#ifndef EDG_SpotVegaParallelTweak_H
#define EDG_SpotVegaParallelTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a bulk shift to the spot volatility
 */

struct SpotVegaParallelTweak: public virtual IScalarTweak {
    virtual ~SpotVegaParallelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
