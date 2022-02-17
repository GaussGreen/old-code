/**
 * @file VolPointwise.hpp
 */

#ifndef QLIB_VolPointwise_H
#define QLIB_VolPointwise_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "vol at an expiry" as a property of market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 *    -  RiskProperty<VolPointwise> represents spot, the property which a market
 *       name can have (see RiskProperty, and VegaPointwise for where it's used
 *       to define "sensitivity to spot")
 *
 *    -  ITweakableWithRespectTo<VolPointwise> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<VolPointwise> parameterises 1D changes to spot of a
 *       particular market name at a particular expiry (see PropertyRiskAxis)
 *
 *    -  PropertyTweakHypothesis<VolPointwise> represents a hypothesis that the
 *       spot of a particular market name at a particular expiry takes a tweaked
 *       value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define VolPointwise as an term-structured property.  For an
 * example of a property that doesn't require any qualifiers, see Spot.
 */

struct RISKMGR_DLL VolPointwise: CObject {
    static CClassConstSP const TYPE;
    VolPointwise(); ~VolPointwise();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for spot (an ExpiryWindow).
     *
     * VolPointwise exists in multiple instances for each name, and requires an
     * ExpiryWindow for full specification.
     */

    typedef ExpiryWindow Qualifier;

    enum {
        /**
         * VolPointwise is continuous, i.e. tweaks to it can be made
         * arbitrarily small.
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_VOLPOINTWISE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<VolPointwise>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<VolPointwise>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<VolPointwise>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
