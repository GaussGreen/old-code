/**
 * @file VolParallel.hpp
 */

#ifndef DRLIB_VolParallel_H
#define DRLIB_VolParallel_H

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
 *    -  RiskProperty<VolParallel> represents vol, the property which a market
 *       name can have (see RiskProperty, and VegaParallel for where it's used
 *       to define "sensitivity to vol")
 *
 *    -  ITweakableWithRespectTo<VolParallel> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<VolParallel> parameterises 1D changes to spot of a
 *       particular market name
 *
 *    -  PropertyTweakHypothesis<VolParallel> represents a hypothesis that the
 *       spot of a particular market name takes a tweaked value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define VolParallel as a scalar property (has no term structure
 * or other qualifiers).  For an example of a term-structured property, see
 * VolPointwise.
 */

struct RISKMGR_DLL VolParallel: CObject {
    static CClassConstSP const TYPE;
    VolParallel(); ~VolParallel();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for parallel vol (none)
     *
     * Some properties exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties---so VolPointwise::Qualifier is
     * ExpiryWindow.
     */

    typedef Void Qualifier;

    enum {
        /**
         * VolParallel is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_VOLPARALLEL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<VolParallel>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<VolParallel>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
