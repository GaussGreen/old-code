/**
 * @file ParSpreadUpfrontPointwise.hpp
 */

#ifndef DRLIB_ParSpreadUpfrontPointwise_H
#define DRLIB_ParSpreadUpfrontPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ParSpreadUpfrontPointwise);

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries.
 */

class RISKMGR_DLL ParSpreadUpfrontPointwise: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                 public virtual Additive 
{
private:
    static void load(CClassSP& clazz);
    Deriv deriv() const;
    
public:
    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    ParSpreadUpfrontPointwise(const string& name, double shiftSize);
    ParSpreadUpfrontPointwise(double shiftSize = DEFAULT_SHIFT);
    ~ParSpreadUpfrontPointwise();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
