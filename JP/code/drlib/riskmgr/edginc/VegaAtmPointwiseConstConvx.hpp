/**
 * @file VegaAtmPointwiseConstConvx.hpp
 */

#ifndef DRLIB_VegaAtmPointwiseConstConvx_H
#define DRLIB_VegaAtmPointwiseConstConvx_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VegaAtmPointwiseConstConvx)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries.
 */

class RISKMGR_DLL VegaAtmPointwiseConstConvx: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    VegaAtmPointwiseConstConvx(const string& name, double shiftSize);
    VegaAtmPointwiseConstConvx(double shiftSize = DEFAULT_SHIFT);
    ~VegaAtmPointwiseConstConvx();

};

DRLIB_END_NAMESPACE

#endif
