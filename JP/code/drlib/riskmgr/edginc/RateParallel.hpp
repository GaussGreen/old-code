/**
 * @file RateParallel.hpp
 */

#ifndef DRLIB_RateParallel_H
#define DRLIB_RateParallel_H

#include "edginc/Void.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "general level of vol curve" as a property of
 * market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 *    -  RiskProperty<RateParallel> represents vol, the property which a market
 *       name can have (see RiskProperty, and VegaParallel for where it's used
 *       to define "sensitivity to vol")
 *
 *    -  ITweakableWithRespectTo<RateParallel> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<RateParallel> parameterises 1D changes to spot of a
 *       particular market name
 *
 *    -  PropertyTweakHypothesis<RateParallel> represents a hypothesis that the
 *       spot of a particular market name takes a tweaked value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define RateParallel as a scalar property (has no term structure
 * or other qualifiers).  For an example of a term-structured property, see
 * RatePointwise.
 */

struct RISKMGR_DLL RateParallel: CObject {
    static CClassConstSP const TYPE;
    RateParallel(); ~RateParallel();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for parallel vol (none)
     *
     * Some properties exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties---so RatePointwise::Qualifier is
     * ExpiryWindow.
     */

    typedef Void Qualifier;

    enum {
        /**
         * RateParallel is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_RATEPARALLEL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<RateParallel>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<RateParallel>);
#endif

DRLIB_END_NAMESPACE

#endif
