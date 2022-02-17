/**
 * @file CorrSwapBasisAdjTweak.hpp
 */

#ifndef QLIB_CorrSwapBasisAdjTweak_H
#define QLIB_CorrSwapBasisAdjTweak_H

#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CorrSwapBasisAdjTweak)

/**
 * A greek calculated by tweaking each market name's corr swap basis adjustment at all its
 * defined expiries simultaneously.
 */

class RISKMGR_DLL CorrSwapBasisAdjTweak: public ScalarRiskPropertySensitivity,
                             public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CorrSwapBasisAdjTweak(double shiftSize = DEFAULT_SHIFT);
    ~CorrSwapBasisAdjTweak();
};

typedef smartPtr<CorrSwapBasisAdjTweak> CorrSwapBasisAdjTweakSP;

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
