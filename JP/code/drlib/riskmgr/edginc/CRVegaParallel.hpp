/**
 * Group : Credit QR NY
 * 
 * File  : CRVegaParallel.hpp
 * 
 * Description : the tweak tag that describes a parallel tweak for spread volatility cubes. The greek is calculated by tweaking all market volatility smiles at all its
 * defined expiries simultaneously.
 * 
 * Author : Adrian Bozdog
 */


#ifndef QLIB_CRVegaParallel_H
#define QLIB_CRVegaParallel_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CRVegaParallel)

class RISKMGR_DLL CRVegaParallel: public ScalarRiskPropertySensitivity,
                    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CRVegaParallel(double shiftSize = DEFAULT_SHIFT);
    ~CRVegaParallel();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
