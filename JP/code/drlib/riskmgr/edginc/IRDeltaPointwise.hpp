/**
 * @file IRDeltaPointwise.hpp
 */

#ifndef DRLIB_IRDeltaPointwise_H
#define DRLIB_IRDeltaPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IRDeltaPointwise)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries.
 *
 * IRDeltaPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL IRDeltaPointwise: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double SENSITIVITY_UNIT;
    static const double DEFAULT_SHIFT;

    IRDeltaPointwise(const string& name, double shiftSize);
    IRDeltaPointwise(double shiftSize = DEFAULT_SHIFT);
    ~IRDeltaPointwise();

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
