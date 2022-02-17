/**
 * @file LocalCorrSensParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/LocalCorrVoid.hpp"
#include "edginc/LocalCorrSensParallel.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const LocalCorrSensParallel::TYPE = CClass::registerClassLoadMethod(
    "LocalCorrSensParallel", typeid(LocalCorrSensParallel), load);

const double LocalCorrSensParallel::DEFAULT_SHIFT = 0.01;
const double LocalCorrSensParallel::SENSITIVITY_UNIT = 0.01;

LocalCorrSensParallel::LocalCorrSensParallel(CClassConstSP           type, 
                                                   double                  shiftSize, 
                                                   string                  name,
                                                   LocalCorrVoidConstSP    tag) : 
    PerNameRiskPropertySensitivity<Void>(type, shiftSize, name), tag(tag) {}   

PerNameRiskPropertySensitivity<Void>::Deriv LocalCorrSensParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<LocalCorrVoid>::SP(tag),
                 IScalarDerivative::oneSided(),
                 LocalCorrSensParallel::SENSITIVITY_UNIT);
}

void LocalCorrSensParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(LocalCorrSensParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<Void>);    
    FIELD(tag, "tag");
    FIELD_MAKE_TRANSIENT(tag);
}

void LocalCorrSensParallel::validatePop2Object() {
    static const string method("LocalCorrSensParallel::validatePop2Object");
    try {
        if (Maths::isPositive(shiftSize - 1.0) || Maths::isNegative(shiftSize + 1.0) ) {
            throw ModelException(method, "Shiftsize has to be between -1.0 and +1.0, but equals" +
                Format::toString(shiftSize) + ".");
        }
    } catch (exception &e) {  
        throw ModelException(e, method);
    }
}

bool LocalCorrSensParallelLinkIn() {
    return LocalCorrSensParallel::TYPE != NULL;
}

#define DECL(CLASS, OP, NAME)                                        \
                                                                                \
    struct CLASS: LocalCorrSensParallel {                                                     \
                                                                                \
        static CClassConstSP const TYPE;                                        \
                                                                                \
        CLASS(double shiftSize = DEFAULT_SHIFT):                                \
            LocalCorrSensParallel(TYPE, shiftSize, NAME,                                      \
                    LocalCorrVoid::SP(OP))                                  \
        {}                                                                      \
                                                                                \
        static void load(CClassSP& clazz) {                                     \
            clazz->setPublic();                                                 \
            REGISTER(CLASS, clazz);                                             \
            SUPERCLASS(LocalCorrSensParallel);                                                \
            EMPTY_SHELL_METHOD(&DefaultConstructor<CLASS>::iObject);            \
            SensitivityFactory::addSens(NAME,                                   \
                                        new GenericSensitivityFactory<CLASS>(), \
                                        new CLASS(),                            \
                                        ITweakableWithRespectTo<LocalCorrVoid>::TYPE); \
        }                                                                                   \
    };                                                                                      \
                                                                                            \
    CClassConstSP const CLASS::TYPE = CClass::registerClassLoadMethod(#CLASS, typeid(CLASS), load);

/**
 * 
 */

DECL(LocalCorrSqueezeParallel, LocalCorrVoid::PARALLEL_SHIFT, "LOCAL_CORR_SQUEEZE_PARALLEL")

/**
 * 
 */

DECL(LocalCorrSkewParallel, LocalCorrVoid::SKEWED_SHIFT, "LOCAL_CORR_SKEW_PARALLEL")

DRLIB_END_NAMESPACE
