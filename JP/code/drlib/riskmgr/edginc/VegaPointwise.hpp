/**
 * @file VegaPointwise.hpp
 */

#ifndef DRLIB_VegaPointwise_H
#define DRLIB_VegaPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VegaPointwise)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries.
 *
 * VegaPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL VegaPointwise: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    VegaPointwise(const string& name, double shiftSize);
    VegaPointwise(double shiftSize = DEFAULT_SHIFT);
    ~VegaPointwise();

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
