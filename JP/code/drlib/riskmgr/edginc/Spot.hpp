/**
 * @file Spot.hpp
 */

#ifndef DRLIB_Spot_H
#define DRLIB_Spot_H

#include "edginc/Void.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "spot" as a property of market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 *    -  RiskProperty<Spot> represents spot, the property which a market
 *       name can have (see RiskProperty, and Delta for where it's used
 *       to define "sensitivity to spot")
 *
 *    -  ITweakableWithRespectTo<Spot> is the interface identifying
 *       such market objects
 *
 *    -  PropertyRiskAxis<Spot> parameterises 1D changes to spot of a
 *       particular market name (see PropertyRiskAxis)
 *
 *    -  PropertyTweakHypothesis<Spot> represents a hypothesis that the
 *       spot of a particular market name takes a tweaked value
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * are a part, see IRiskQuantityFactory (particularly <I>The new setup: property
 * tags</I>).
 *
 * The typedefs define Spot as a scalar property (has no term structure or
 * other qualifiers).  For an example of a term-structured property, see
 * VolPointwise.
 */

struct Spot: CObject {
    static CClassConstSP const TYPE;
    Spot(); ~Spot();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for spot (none)
     *
     * Some properties exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties---so VolPointwise::Qualifier is
     * ExpiryWindow.  Spot is a "scalar" property so we just specify Void here.
     */

    typedef Void Qualifier;

    enum {
        /**
         * Spot is continuous, i.e. tweaks to it can be made arbitrarily small
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_SPOT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Spot>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Spot>);
#endif
DRLIB_END_NAMESPACE

#endif
