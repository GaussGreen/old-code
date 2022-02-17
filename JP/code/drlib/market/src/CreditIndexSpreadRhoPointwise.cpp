/**
 * @file CreditIndexSpreadRhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/CreditIndexSpreadPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/TRACE.hpp"

DRLIB_BEGIN_NAMESPACE

static const double ONE_BASIS_POINT = 0.0001;
static const double SENSITIVITY_UNIT = ONE_BASIS_POINT;

class MARKET_DLL CreditIndexSpreadRhoPointwise : public PerNameRiskPropertySensitivity<ExpiryWindow>,
                             public virtual Additive,
                             public CreditTweak {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CreditIndexSpreadRhoPointwise, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CreditIndexSpreadRhoPointwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CreditIndexSpreadRhoPointwise>(), 
                                    new CreditIndexSpreadRhoPointwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<CreditIndexSpreadPointwise>::TYPE);
    }

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<CreditIndexSpreadPointwise>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const char* NAME;

    CreditIndexSpreadRhoPointwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
    {}
};

const char* CreditIndexSpreadRhoPointwise::NAME = "CREDIT_INDEX_SPREAD_RHO_POINTWISE";
const double CreditIndexSpreadRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

CClassConstSP const CreditIndexSpreadRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "CreditIndexSpreadRhoPointwise", typeid(CreditIndexSpreadRhoPointwise), load);


bool CreditIndexSpreadRhoPointwiseLinkIn() {
    return CreditIndexSpreadRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE
