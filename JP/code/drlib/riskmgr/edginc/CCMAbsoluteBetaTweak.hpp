/**
 * @file CCMAbsoluteBetaTweak.hpp
 *
 *
 */

#ifndef EDG_CCMAbsoluteBetaTweak_H
#define EDG_CCMAbsoluteBetaTweak_H

#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Parameterises an absolute (not proportional) tweak to CCM "beta" parameter
 *
 * This defines the parameters (actually just getShiftSize()) required by
 * ParSpreadCurve when it applies the same scaling to all its par spreads
 *
 * See ParSpreadCurve::shiftSens(CCMAbsoluteBetaTweak *), which is
 * an implementation of TweakableWith<CCMAbsoluteBetaTweak>.
 *
 * The sensitivities defined with respect to this tweak are
 *
 *    -  CCMDBasketBetaDSpread.cpp
 */

struct CCMAbsoluteBetaTweak: public CreditTweak,
                             public virtual IScalarTweak {

    virtual ~CCMAbsoluteBetaTweak() {}
};

DRLIB_END_NAMESPACE

#endif
