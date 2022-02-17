/**
 * @file LegalBasisAdditiveRelativeParallelTweak.hpp
 */

#ifndef QLIB_LegalBasisAdditiveRelativeParallelTweak_H
#define QLIB_LegalBasisAdditiveRelativeParallelTweak_H

#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "relative parallel par spread legal basis" as a property of
 * market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 * See VolPointwise for a fully documented analogous example.
 */

struct RISKMGR_DLL LegalBasisAdditiveRelativeParallelTweak: CObject {
    static CClassConstSP const TYPE;
    LegalBasisAdditiveRelativeParallelTweak(); ~LegalBasisAdditiveRelativeParallelTweak();

  
    typedef Void Qualifier;

    enum {
        /**
         * LegalBasisAdditiveRelativeParallelTweak is continuous, i.e. tweaks to it can be made
         * arbitrarily small.
         *
         * This is carried through to Sensitivity::discreteShift()
         * via RiskProperty<PROPERTY>::discrete().
         */

        discrete = 0
    };
};

#ifndef QLIB_LegalBasisAdditiveRelativeParallelTweak_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<LegalBasisAdditiveRelativeParallelTweak>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
