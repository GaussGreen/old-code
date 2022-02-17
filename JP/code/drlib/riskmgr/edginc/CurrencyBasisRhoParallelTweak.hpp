/**
 *
 *
 */

#ifndef EDG_CurrencyBasisRhoParallelTweak_H
#define EDG_CurrencyBasisRhoParallelTweak_H

#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a bulk shift to a currency basis spread curve
 */
struct CurrencyBasisRhoParallelTweak: public virtual IScalarTweak  {
    virtual ~CurrencyBasisRhoParallelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
