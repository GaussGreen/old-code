/**
 *
 */

#ifndef EDG_CompositeVegaPointwiseTweak_H
#define EDG_CompositeVegaPointwiseTweak_H

#include "edginc/VectorTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a pointwise shift to the composite volatilities
 */

struct CompositeVegaPointwiseTweak: public virtual IVectorTweak{
    virtual ~CompositeVegaPointwiseTweak() {}
};


DRLIB_END_NAMESPACE

#endif
