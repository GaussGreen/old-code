/**
 * @file PhiParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Phi.hpp"
#include "edginc/PhiParallel.hpp"
#include "edginc/Correl.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * First derivative of instrument price with respect to a simultaneous change
 * towards 0 in all the asset-asset correlations in the market.
 */

struct PhiParallel: AllNamesRiskPropertySensitivity,
                    virtual Additive {

    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    static const string NAME;

    bool onlyEqEq;

    PhiParallel(double shiftSize = DEFAULT_SHIFT,
                bool onlyEqEq = true):
        AllNamesRiskPropertySensitivity(TYPE, NAME, shiftSize),
        onlyEqEq(onlyEqEq)
    {}

    AllNamesRiskPropertySensitivity::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<Correl>::SP(Correl::SP(true, // asset
                                                         !onlyEqEq, // fx
                                                         !onlyEqEq, // other
                                                         Correl::ABSOLUTE)),
                     IScalarDerivative::oneSided(),
                     Phi_sensitivityUnit);
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(PhiParallel, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(AllNamesRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<PhiParallel>::iObject);
        FIELD(onlyEqEq, "Whether to shift Eq-Eq only or all correlations");
        FIELD_MAKE_OPTIONAL(onlyEqEq);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<PhiParallel>(), 
                                    new PhiParallel(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Correl>::TYPE);
    }
};

const string PhiParallel::NAME = "PHI_PARALLEL";
const double PhiParallel::DEFAULT_SHIFT = PhiParallel_defaultShift;

CClassConstSP const PhiParallel::TYPE = CClass::registerClassLoadMethod("PhiParallel", typeid(PhiParallel), load);

SensitivitySP PhiParallel_SP(double shiftSize, bool onlyEqEq) {
    return SensitivitySP(new PhiParallel(shiftSize, onlyEqEq));
}

bool PhiParallelLinkIn() {
    return PhiParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
