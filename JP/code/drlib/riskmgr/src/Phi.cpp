/**
 * @file Phi.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/Correl.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/Additive.hpp"
#include "edginc/Phi.hpp"

DRLIB_BEGIN_NAMESPACE

struct PhiBase: ScalarRiskPropertySensitivity,
                virtual Additive {

    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;

    PhiBase(CClassConstSP type, double shiftSize, string name,
            CorrelConstSP tag):
        ScalarRiskPropertySensitivity(type, name, shiftSize),
        tag(tag)
    {}

    CorrelConstSP tag;

    ScalarRiskPropertySensitivity::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<Correl>::SP(tag),
                     IScalarDerivative::oneSided(),
                     Phi_sensitivityUnit);
    }

    static void load(CClassSP& clazz) {
        REGISTER(PhiBase, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        FIELD(tag, "tag");
        FIELD_MAKE_TRANSIENT(tag);
    }

    void validatePop2Object() {
        if (tag->op == Correl::SQUEEZE &&
            (getShiftSize() < 0 || 1 < getShiftSize())) {
            throw ModelException(getClass()->getName() + "::validatePop2Object",
                getClass()->getName() + " shiftSize " +
                "(" + Format::toString(getShiftSize()) + ") "
                "must be between 0.0 and 1.0");
        }
    }
};

const double PhiBase::DEFAULT_SHIFT = 0.01;

CClassConstSP const PhiBase::TYPE = CClass::registerClassLoadMethod(
    "PhiBase", typeid(PhiBase), load);

bool PhiLinkIn() {
    return PhiBase::TYPE != NULL;
}

#define DECL(CLASS, ASSET, FX, OTHER, OP, NAME)                                 \
                                                                                \
    struct CLASS: PhiBase {                                                     \
                                                                                \
        static CClassConstSP const TYPE;                                        \
                                                                                \
        static const string name;                                               \
                                                                                \
        CLASS(double shiftSize = DEFAULT_SHIFT):                                \
            PhiBase(TYPE, shiftSize, NAME,                                      \
                    Correl::SP(ASSET, FX, OTHER, OP))                           \
        {}                                                                      \
                                                                                \
        static void load(CClassSP& clazz) {                                     \
            clazz->setPublic();                                                 \
            REGISTER(CLASS, clazz);                                             \
            SUPERCLASS(PhiBase);                                                \
            EMPTY_SHELL_METHOD(&DefaultConstructor<CLASS>::iObject);            \
            SensitivityFactory::addSens(NAME,                                   \
                                        new GenericSensitivityFactory<CLASS>(), \
                                        new CLASS(),                            \
                                        ITweakableWithRespectTo<Correl>::TYPE); \
        }                                                                       \
    };                                                                          \
                                                                                \
    CClassConstSP const CLASS::TYPE = CClass::registerClassLoadMethod(#CLASS, typeid(CLASS), load);

/**
 * First derivative of instrument price with respect to a change towards 0 in
 * each asset-asset correlation in the market in turn
 */

DECL(Phi, true, false, false, Correl::ABSOLUTE_INWARD, "PHI")

/**
 * First derivative of instrument price with respect to a change towards 0 in
 * each FX correlation in the market in turn
 */

DECL(FXPhi, false, true, false, Correl::ABSOLUTE_INWARD, "FX_PHI")

/**
 * First derivative of instrument price with respect to a change towards 0 in
 * each correlation in the market in turn (including asset, FX and other
 * correls)
 */

DECL(PhiPlus, true, true, true, Correl::ABSOLUTE_INWARD, "PHIPLUS")

/**
 * First derivative of instrument price with respect to a "squeeze" towards 1
 * in each correlation in the market in turn (including asset, FX and other
 * correls)
 */

DECL(CorrelationSqueeze, true, true, true, Correl::SQUEEZE, "CORRELATION_SQUEEZE")

DRLIB_END_NAMESPACE
