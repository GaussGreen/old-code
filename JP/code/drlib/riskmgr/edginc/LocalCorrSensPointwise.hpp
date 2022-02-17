/**
 * @file LocalCorrSensPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/LocalCorrExpiry.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(LocalCorrSensPointwise)

class RISKMGR_DLL LocalCorrSensPointwise : public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                           public virtual Additive {
public:
    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    LocalCorrSensPointwise(CClassConstSP            type, 
                           double                   shiftSize, 
                           string                   name,
                           LocalCorrExpiryConstSP   tag);

    LocalCorrExpiryConstSP tag;

    Deriv deriv() const;

    static void load(CClassSP& clazz);
    
    void validatePop2Object();
};

DRLIB_END_NAMESPACE
