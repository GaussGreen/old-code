/**
 *
 *
 */

#ifndef EDG_CurrencyBasisRhoPointwiseTweak_H
#define EDG_CurrencyBasisRhoPointwiseTweak_H

#include "edginc/VectorTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a pointwise shift to the spread curve of a currency basis
 */
struct CurrencyBasisRhoPointwiseTweak: public virtual IVectorTweak {
    virtual ~CurrencyBasisRhoPointwiseTweak() {}
};

DRLIB_END_NAMESPACE

#endif
