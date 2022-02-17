/**
 * @file CreditIndexSpreadRhoParallel.cpp
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
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/CreditIndexSpreadPointwise.hpp"
#include "edginc/CreditIndexSpreadParallel.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/TRACE.hpp"

DRLIB_BEGIN_NAMESPACE

static const double ONE_BASIS_POINT = 0.0001;
static const double SENSITIVITY_UNIT = ONE_BASIS_POINT;

class MARKET_DLL CreditIndexSpreadRhoParallel: public ScalarRiskPropertySensitivity,
                            public virtual Additive,
                            public virtual CreditTweak {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CreditIndexSpreadRhoParallel, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CreditIndexSpreadRhoParallel>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CreditIndexSpreadRhoParallel>(), 
                                    new CreditIndexSpreadRhoParallel(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<CreditIndexSpreadParallel>::TYPE);
    }

    ScalarRiskPropertySensitivity::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                    RiskProperty<CreditIndexSpreadParallel>::SP(),
                    IScalarDerivative::oneSided(),
                    ONE_BASIS_POINT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const char* NAME;

    CreditIndexSpreadRhoParallel(double shiftSize = DEFAULT_SHIFT,
                                 const string& packetName = NAME):
    ScalarRiskPropertySensitivity(TYPE, packetName, shiftSize) {}

    ~CreditIndexSpreadRhoParallel() {}

};

const char* CreditIndexSpreadRhoParallel::NAME = "CREDIT_INDEX_SPREAD_RHO_PARALLEL";
const double CreditIndexSpreadRhoParallel::DEFAULT_SHIFT = ONE_BASIS_POINT;

CClassConstSP const CreditIndexSpreadRhoParallel::TYPE = CClass::registerClassLoadMethod(
    "CreditIndexSpreadRhoParallel", typeid(CreditIndexSpreadRhoParallel), load);

bool CreditIndexSpreadRhoParallelLinkIn() {
    return CreditIndexSpreadRhoParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

