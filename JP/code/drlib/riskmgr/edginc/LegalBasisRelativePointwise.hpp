/**
 * @file LegalBasisRelativePointwise.hpp
 */

#ifndef QLIB_LegalBasisRelativePointwise_H
#define QLIB_LegalBasisRelativePointwise_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "par spread legal basis at an expiry" as a property of
 * market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 * See VolPointwise for a fully documented analogous example.
 */

struct RISKMGR_DLL LegalBasisRelativePointwise: CObject {
    static CClassConstSP const TYPE;
    LegalBasisRelativePointwise(); ~LegalBasisRelativePointwise();

    /**
     * Info needed in addition to market data name to make a one-dimensional
     * "risk axis" for spot (an ExpiryWindow).
     *
     * LegalBasisRelativePointwise exists in multiple instances for each name, and requires an
     * ExpiryWindow for full specification.
     */

    typedef ExpiryWindow Qualifier;

    enum {
        /**
         * LegalBasisRelativePointwise is continuous, i.e. tweaks to it can be made
         * arbitrarily small.
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_LEGALBASISRELATIVEPOINTWISE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<LegalBasisRelativePointwise>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
