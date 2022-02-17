/**
 * @file Smile2QElementwise.hpp
 */

#ifndef QLIB_Smile2QElementwise_HPP
#define QLIB_Smile2QElementwise_HPP

#include "edginc/Void.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/BoxedInt.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "general level of vol curve" as a property of
 * market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 *    -  RiskProperty<Smile2QElementwise> represents vol, the property which a market
 *       name can have (see RiskProperty, and VegaParallel for where it's used
 *       to define "sensitivity to vol")
 *
 *    -  ITweakableWithRespectTo<Smile2QElementwise> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<Smile2QElementwise> parameterises 1D changes to spot of a
 *       particular market name
 *
 *    -  PropertyTweakHypothesis<Smile2QElementwise> represents a hypothesis that the
 *       spot of a particular market name takes a tweaked value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define Smile2QElementwise as a scalar property (has no term structure
 * or other qualifiers).  For an example of a term-structured property, see
 * RatePointwise.
 */

struct RISKMGR_DLL Smile2QElementwise: CObject {
    static CClassConstSP const TYPE;
    Smile2QElementwise(); ~Smile2QElementwise();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for parallel vol (none)
     *
     * Some properties exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties---so RatePointwise::Qualifier is
     * ExpiryWindow.
     */

    typedef BoxedInt Qualifier;

    enum {
        /**
         * Smile2QElementwise is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_Smile2QElementwise_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Smile2QElementwise>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Smile2QElementwise>);
#endif

DRLIB_END_NAMESPACE

#endif
