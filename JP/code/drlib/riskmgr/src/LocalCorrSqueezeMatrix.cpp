/**
 * @file LocalCorrSqueezeMatrix.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/LocalCorrExpiryAndStrikewise.hpp"
#include "edginc/LocalCorrSqueezeMatrix.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

const string LocalCorrSqueezeMatrix::NAME = "LOCAL_CORR_SQUEEZE_MATRIX";
const double LocalCorrSqueezeMatrix::DEFAULT_SHIFT = 0.01;
const double SENSITIVITY_UNIT = 0.01;


LocalCorrSqueezeMatrix::LocalCorrSqueezeMatrix(const string& name, double shiftSize) :
    PerNameRiskPropertySensitivity<ExpiryAndStrike>(TYPE, shiftSize, name) 
{}   

LocalCorrSqueezeMatrix::LocalCorrSqueezeMatrix(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryAndStrike>(TYPE, shiftSize, NAME) 
{}

PerNameRiskPropertySensitivity<ExpiryAndStrike>::Deriv LocalCorrSqueezeMatrix::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<LocalCorrExpiryAndStrikewise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

void LocalCorrSqueezeMatrix::validatePop2Object() {
    static const string method("LocalCorrSqueezeMatrix::validatePop2Object");
    try {
        if (Maths::isPositive(shiftSize - 1.0) || Maths::isNegative(shiftSize + 1.0) ) {
            throw ModelException(method, "Shiftsize has to be between -1.0 and +1.0, but equals" +
                Format::toString(shiftSize) + ".");
        }
    } catch (exception &e) {  
        throw ModelException(e, method);
    }
}

LocalCorrSqueezeMatrix::~LocalCorrSqueezeMatrix() {}

void LocalCorrSqueezeMatrix::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LocalCorrSqueezeMatrix, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryAndStrike>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<LocalCorrSqueezeMatrix>::iObject);
    SensitivityFactory::addSens(LocalCorrSqueezeMatrix::NAME, 
                                new GenericSensitivityFactory<LocalCorrSqueezeMatrix>(), 
                                new LocalCorrSqueezeMatrix(
                                        LocalCorrSqueezeMatrix::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<LocalCorrExpiryAndStrikewise>::TYPE);
}

CClassConstSP const LocalCorrSqueezeMatrix::TYPE = CClass::registerClassLoadMethod(
    "LocalCorrSqueezeMatrix", typeid(LocalCorrSqueezeMatrix), load);

bool LocalCorrSqueezeMatrixLinkIn() {
    return LocalCorrSqueezeMatrix::TYPE != NULL;
}

DRLIB_END_NAMESPACE
