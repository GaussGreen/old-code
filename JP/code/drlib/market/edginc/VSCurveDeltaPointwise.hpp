#ifndef DRLIB_VSCurveDeltaPointwise_H
#define DRLIB_VSCurveDeltaPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VSCurveDeltaPointwise)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries
 */

class MARKET_DLL VSCurveDeltaPointwise:
        public PerNameRiskPropertySensitivity<ExpiryWindow>,
        public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    VSCurveDeltaPointwise(const string& name, double shiftSize);
    VSCurveDeltaPointwise(double shiftSize = DEFAULT_SHIFT);
    ~VSCurveDeltaPointwise();

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
