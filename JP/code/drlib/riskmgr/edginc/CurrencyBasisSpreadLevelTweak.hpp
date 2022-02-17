/**
 *
 *
 */

#ifndef EDG_CurrencyBasisSpreadLevelTweak_H
#define EDG_CurrencyBasisSpreadLevelTweak_H
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises a scenario to set level to a Currency Basis
 */
struct CurrencyBasisSpreadLevelTweak: public virtual IScalarTweak {
    virtual ~CurrencyBasisSpreadLevelTweak() {}
};

DRLIB_END_NAMESPACE

#endif
