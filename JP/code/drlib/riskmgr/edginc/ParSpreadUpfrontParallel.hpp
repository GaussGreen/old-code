//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ParSpreadUpfrontParallel.hpp
//
//   Description : ParSpread curve Upfront parallel shift
//
//   Date        : Nov 2007
//
//----------------------------------------------------------------------------

#ifndef DRLIB_PARSPREADUPFRONTPARALLEL_H
#define DRLIB_PARSPREADUPFRONTPARALLEL_H

#include "edginc/Additive.hpp"
#include "edginc/ParSpreadUpfrontParallelTP.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A greek calculated by parallel tweaking parspread curve upfront fees.
 */

class RISKMGR_DLL ParSpreadUpfrontParallel: public ScalarRiskPropertySensitivity,
                                            public virtual Additive 
{
private:
    static void load(CClassSP& clazz);
    ScalarRiskPropertySensitivity::Deriv deriv() const;
    
public:
    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    ParSpreadUpfrontParallel(double shiftSize, const string& name);
    ParSpreadUpfrontParallel(double shiftSize = DEFAULT_SHIFT);
    ~ParSpreadUpfrontParallel();
};

DRLIB_END_NAMESPACE

#endif
