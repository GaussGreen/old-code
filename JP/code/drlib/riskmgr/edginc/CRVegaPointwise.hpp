/**
 * Group : Credit QR NY
 * 
 * File  : CRVegaPointwise.hpp
 * 
 * Description : the tweak tag that describes a pointwise tweak for spread volatility cubes. The greek is calculated by tweaking each market volatility smile at all its
 * defined expiries simultaneously.
 * 
 * Author : Adrian Bozdog
 */

#ifndef DRLIB_CRVegaPointwise_H
#define DRLIB_CRVegaPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CRVegaPointwise)

class RISKMGR_DLL CRVegaPointwise: public PerNameRiskPropertySensitivity<ExpiryPair>,
                       public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    CRVegaPointwise(const string& name, double shiftSize);
    CRVegaPointwise(double shiftSize = DEFAULT_SHIFT);
    ~CRVegaPointwise();

};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
