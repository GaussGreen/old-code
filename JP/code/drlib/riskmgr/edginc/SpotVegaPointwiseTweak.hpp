/**
 *
 */

#ifndef EDG_SpotVegaPointwiseTweak_H
#define EDG_SpotVegaPointwiseTweak_H

#include "edginc/VectorTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a pointwise shift to the spot volatilities
 */

struct SpotVegaPointwiseTweak: public virtual IVectorTweak {
    virtual ~SpotVegaPointwiseTweak() {}
};

DRLIB_END_NAMESPACE

#endif
