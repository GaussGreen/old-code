/**
 * @file LocalCorrSensParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/LocalCorrVoid.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(LocalCorrSensParallel)

class RISKMGR_DLL LocalCorrSensParallel : public PerNameRiskPropertySensitivity<Void>,
                                          public virtual Additive {
public:
    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;


    LocalCorrSensParallel(CClassConstSP            type, 
                          double                   shiftSize, 
                          string                   name,
                          LocalCorrVoidConstSP     tag);
    
    LocalCorrVoidConstSP tag;

    Deriv deriv() const;

    static void load(CClassSP& clazz);
    
    void validatePop2Object();

};

DRLIB_END_NAMESPACE
