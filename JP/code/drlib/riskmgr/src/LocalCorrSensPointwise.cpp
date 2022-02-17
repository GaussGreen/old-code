/**
 * @file LocalCorrSensPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/LocalCorrExpiry.hpp"
#include "edginc/LocalCorrSensPointwise.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const LocalCorrSensPointwise::TYPE = CClass::registerClassLoadMethod(
    "LocalCorrSensPointwise", typeid(LocalCorrSensPointwise), load);

const double LocalCorrSensPointwise::DEFAULT_SHIFT = 0.01;
const double LocalCorrSensPointwise::SENSITIVITY_UNIT = 0.01;

LocalCorrSensPointwise::LocalCorrSensPointwise(CClassConstSP             type, 
                                               double                    shiftSize, 
                                               string                    name,
                                               LocalCorrExpiryConstSP    tag) : 
    PerNameRiskPropertySensitivity<ExpiryWindow>(type, shiftSize, name), tag(tag) {}   

PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv LocalCorrSensPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<LocalCorrExpiry>::SP(tag),
                 IScalarDerivative::oneSided(),
                 LocalCorrSensPointwise::SENSITIVITY_UNIT);
}

void LocalCorrSensPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LocalCorrSensPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);    
    FIELD(tag, "tag");
    FIELD_MAKE_TRANSIENT(tag);
}

void LocalCorrSensPointwise::validatePop2Object() {
    static const string method("LocalCorrSensPointwise::validatePop2Object");
    try {
        if (Maths::isPositive(shiftSize - 1.0) || Maths::isNegative(shiftSize + 1.0) ) {
            throw ModelException(method, "Shiftsize has to be between -1.0 and +1.0, but equals" +
                Format::toString(shiftSize) + ".");
        }
    } catch (exception &e) {  
        throw ModelException(e, method);
    }
}

    
bool LocalCorrSensPointwiseLinkIn() {
    return LocalCorrSensPointwise::TYPE != NULL;
}

#define DECL(CLASS, OP, NAME)                                        \
                                                                                \
    struct CLASS: LocalCorrSensPointwise {                                                     \
                                                                                \
        static CClassConstSP const TYPE;                                        \
                                                                                \
        CLASS(double shiftSize = DEFAULT_SHIFT):                                \
            LocalCorrSensPointwise(TYPE, shiftSize, NAME,                                      \
                    LocalCorrExpiry::SP(OP))                                  \
        {}                                                                      \
                                                                                \
        static void load(CClassSP& clazz) {                                     \
            clazz->setPublic();                                                 \
            REGISTER(CLASS, clazz);                                             \
            SUPERCLASS(LocalCorrSensPointwise);                                                \
            EMPTY_SHELL_METHOD(&DefaultConstructor<CLASS>::iObject);            \
            SensitivityFactory::addSens(NAME,                                   \
                                        new GenericSensitivityFactory<CLASS>(), \
                                        new CLASS(),                            \
                                        ITweakableWithRespectTo<LocalCorrExpiry>::TYPE); \
        }                                                                                   \
    };                                                                                      \
                                                                                            \
    CClassConstSP const CLASS::TYPE = CClass::registerClassLoadMethod(#CLASS, typeid(CLASS), load);

/**
 * 
 */

DECL(LocalCorrSqueezePointwise, LocalCorrExpiry::PARALLEL_SHIFT, "LOCAL_CORR_SQUEEZE_POINTWISE")

/**
 * 
 */

DECL(LocalCorrSkewPointwise, LocalCorrExpiry::SKEWED_SHIFT, "LOCAL_CORR_SKEW_POINTWISE")

DRLIB_END_NAMESPACE
