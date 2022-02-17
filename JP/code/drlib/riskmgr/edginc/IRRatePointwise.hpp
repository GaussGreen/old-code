/**
 * @file IRRatePointwise.hpp
 */

#ifndef DRLIB_IRRatePointwise_H
#define DRLIB_IRRatePointwise_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "general level of vol curve" as a property of
 * market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 *    -  RiskProperty<IRRatePointwise> represents vol, the property which a market
 *       name can have (see RiskProperty, and VegaParallel for where it's used
 *       to define "sensitivity to vol")
 *
 *    -  ITweakableWithRespectTo<IRRatePointwise> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<IRRatePointwise> parameterises 1D changes to spot of a
 *       particular market name
 *
 *    -  PropertyTweakHypothesis<IRRatePointwise> represents a hypothesis that the
 *       spot of a particular market name takes a tweaked value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define IRRatePointwise as a scalar property (has no term structure
 * or other qualifiers).  For an example of a term-structured property, see
 * IRRatePointwise.
 */

struct RISKMGR_DLL IRRatePointwise: CObject {
    static CClassConstSP const TYPE;
    IRRatePointwise(); ~IRRatePointwise();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for parallel vol (none)
     *
     * Some properties exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties---so IRRatePointwise::Qualifier is
     * ExpiryWindow.
     */

    typedef ExpiryWindow Qualifier;

    enum {
        /**
         * IRRatePointwise is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_IRRATEPOINTWISE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<IRRatePointwise>);
#endif

DRLIB_END_NAMESPACE

#endif
