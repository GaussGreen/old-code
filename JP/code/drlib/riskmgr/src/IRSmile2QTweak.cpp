#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/Additive.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/IRSmile2QTweak.hpp"
#include "edginc/RiskProperty_TYPES.hpp"


DRLIB_BEGIN_NAMESPACE

const double ONE_BASIS_POINT = 0.0001;

/**
 * Calculates derivative of instrument price w.r.t. bulk shifts in
 * the mean reversion
 *
 * The output generated by this sensitivity is:
 * <DL>
 * <DT>CR_MEAN_REVERSION_PARALLEL 
 * <DD>Change in price when the mean reversion is shifted up by 1bp
 * </DL>
 *
 * The derivative is reported in 20%-units, i.e. the value reported is 
 * 5 times smaller than the actual partial derivative.
 */

void IRSmile2QTweak::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(IRSmile2QTweak, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(ScalarRiskPropertySensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<IRSmile2QTweak>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<IRSmile2QTweak>(), 
                                new IRSmile2QTweak(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<IRSmile2QTweak::Property>::TYPE);
}

ScalarRiskPropertySensitivity::Deriv IRSmile2QTweak::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<IRSmile2QTweak::Property>::SP(),
                 IScalarDerivative::oneSided(),
                 1.0,
                 IScalarDerivative::second(),
                 1.0);
}

IRSmile2QTweak::IRSmile2QTweak(double shiftSize):
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

const string IRSmile2QTweak::NAME = "IR_SMILE2Q_TWEAK";
const double IRSmile2QTweak::DEFAULT_SHIFT = 0.01; // Requirement: 1%

CClassConstSP const IRSmile2QTweak::TYPE = CClass::registerClassLoadMethod(
    "IRSmile2QTweak", typeid(IRSmile2QTweak), load);

void IRSmile2QTweak::Property::load(CClassSP& clazz) {
    REGISTER(IRSmile2QTweak::Property, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(&DefaultConstructor<IRSmile2QTweak::Property>::iObject);
}

CClassConstSP const IRSmile2QTweak::Property::TYPE = CClass::registerClassLoadMethod(
    "IRSmile2QTweak::Property", typeid(IRSmile2QTweak::Property), IRSmile2QTweak::Property::load);

RiskProperty_TYPES(IRSmile2QTweak::Property);

/**
   Included in RiskMgrLib::linkInClasses() to force linkage into the Windows
   exe.
 */
bool IRSmile2QTweakLinkIn() {
    return IRSmile2QTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE
