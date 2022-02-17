/**
 * @file VegaAtmParallelConstConvx.hpp
 */

#ifndef QLIB_VegaAtmParallelConstConvx_H
#define QLIB_VegaAtmParallelConstConvx_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VegaAtmParallelConstConvx)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries simultaneously.
 *
 * VegaAtmParallel has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL VegaAtmParallelConstConvx: public ScalarRiskPropertySensitivity,
                    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    VegaAtmParallelConstConvx(double shiftSize = DEFAULT_SHIFT);
    ~VegaAtmParallelConstConvx();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
