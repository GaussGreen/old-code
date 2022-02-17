/**
 * \file CreditTweak.hpp
 *
 *
 */

#ifndef DRLIB_CreditTweak_H
#define DRLIB_CreditTweak_H

DRLIB_BEGIN_NAMESPACE

/**
 * "Tag" interface for tweaks with respect to credit
 *
 * CDO::endDate() looks at whether the SensControl "implements" this to decide
 * what date to return.
 *
 * Implementing classes are
 *
 *    -  ParSpreadRhoPointwiseTweak (and therefore ParSpreadRhoPointwise,
 *       ParSpreadRhoPointwiseTwoSided),
 *    -  ParSpreadRhoParallelTweak (and therefore ParSpreadRhoParallel,
 *       ParSpreadRhoParallelTwoSided)
 *    -  NOT LiquiditySpreadRhoPointwise: perhaps it should, but that's
 *       not how CDO::endDate() was set up originally
 */

struct CreditTweak {
    virtual ~CreditTweak() {}
};

DRLIB_END_NAMESPACE

#endif
