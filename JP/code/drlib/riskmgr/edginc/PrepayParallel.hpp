//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : PrepayParallel.hpp
//
//   Description : ABCDS prepay curve horizontal parallel shift
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef DRLIB_PREPAYPARALLEL_H
#define DRLIB_PREPAYPARALLEL_H

#include "edginc/Additive.hpp"
#include "edginc/PrepayParallelTP.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A greek calculated by tweaking each market name's prepay curve in datetime.
 */

class RISKMGR_DLL PrepayParallel: public ScalarRiskPropertySensitivity,
                      public virtual Additive 
{
private:
    static void load(CClassSP& clazz);
    ScalarRiskPropertySensitivity::Deriv deriv() const;
    
public:
    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    PrepayParallel(double shiftSize, const string& name);
    PrepayParallel(double shiftSize = DEFAULT_SHIFT);
    ~PrepayParallel();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
