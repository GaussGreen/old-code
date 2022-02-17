/**
 * @file CorrSwapSamplingAdjTweak.hpp
 */

#ifndef QLIB_CorrSwapSamplingAdjTweak_H
#define QLIBCorrSwapSamplingAdjTweakl_H

#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CorrSwapSamplingAdjTweak)

/**
 * A greek calculated by tweaking each market name's corr swap sampling adjustment  
 */

class RISKMGR_DLL CorrSwapSamplingAdjTweak: public ScalarRiskPropertySensitivity,
                                public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CorrSwapSamplingAdjTweak(double shiftSize = DEFAULT_SHIFT);
    ~CorrSwapSamplingAdjTweak();
};

typedef smartPtr<CorrSwapSamplingAdjTweak> CorrSwapSamplingAdjTweakSP;

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
