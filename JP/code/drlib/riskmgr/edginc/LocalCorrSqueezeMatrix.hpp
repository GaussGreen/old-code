/**
 * @file LocalCorrSqueezeMatrix.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(LocalCorrSqueezeMatrix)

class RISKMGR_DLL LocalCorrSqueezeMatrix: 
    public PerNameRiskPropertySensitivity<ExpiryAndStrike>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

    void validatePop2Object();

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    LocalCorrSqueezeMatrix(const string& name, double shiftSize);
    LocalCorrSqueezeMatrix(double shiftSize = DEFAULT_SHIFT);    
    ~LocalCorrSqueezeMatrix();

};

DRLIB_END_NAMESPACE
