#ifndef DRLIB_VSCurveDeltaParallel_H
#define DRLIB_VSCurveDeltaParallel_H

#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VSCurveDeltaParallel)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries
 */

class MARKET_DLL VSCurveDeltaParallel:
        public PerNameRiskPropertySensitivity<Void>,
        public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    VSCurveDeltaParallel(const string& name, double shiftSize);
    VSCurveDeltaParallel(double shiftSize = DEFAULT_SHIFT);
    ~VSCurveDeltaParallel();

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
