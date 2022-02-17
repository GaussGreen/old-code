/**
 * @file LiquiditySpreadRhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/ModeHypothesis.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/LiquiditySpreadPointwise.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/E2CModel.hpp"

DRLIB_BEGIN_NAMESPACE

static const char* NAME = "LIQUIDITY_SPREAD_RHO_POINTWISE";
static const double ONE_BASIS_POINT = 0.0001;

// 
// *******************************
//  ParSpreadPricingOffHypothesis
// *******************************
// 

class ParSpreadPricingOffHypothesis: public ModeHypothesis {

    static IObject* defaultOne() {
        return new ParSpreadPricingOffHypothesis();
    }

    static const CClassConstSP TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(ParSpreadPricingOffHypothesis, clazz);
        SUPERCLASS(ModeHypothesis);
        EMPTY_SHELL_METHOD(defaultOne);
    }

public:

    ParSpreadPricingOffHypothesis():
        ModeHypothesis(TYPE)
    {}

    double applyToWorld(MultiTweakGroupSP world, bool* changed) const {
        if (IE2CModel::TYPE->isInstance(world->getModelSP())) {
            dynamic_cast<IE2CModel&>(*world->getModel()).setParSpreadPricing(false);
            if (changed) *changed = true;
            return 1;
        }
        else {
            if (changed) *changed = false;
            return 0;
        }
    }

    void undo(MultiTweakGroupSP world) const {
        if (IE2CModel::TYPE->isInstance(world->getModelSP())) {
            dynamic_cast<IE2CModel&>(*world->getModel()).setParSpreadPricing(true);
        }
    }
};

CClassConstSP const ParSpreadPricingOffHypothesis::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadPricingOffHypothesis", typeid(ParSpreadPricingOffHypothesis), load);

// 
// *****************************
//  LiquiditySpreadRhoPointwise
// *****************************
// 

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * LiquiditySpreadRhoPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class LiquiditySpreadRhoPointwise: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                   public virtual Additive {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(LiquiditySpreadRhoPointwise, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<LiquiditySpreadRhoPointwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<LiquiditySpreadRhoPointwise>(), 
                                    new LiquiditySpreadRhoPointwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<LiquiditySpreadPointwise>::TYPE);
    }

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     ParSpreadMaxTweakSize::adapted(
                         RiskProperty<LiquiditySpreadPointwise>::SP()),
                     IScalarDerivative::oneSided()->underScenario(
                         IHypothesisSP(new ParSpreadPricingOffHypothesis())),
                     ONE_BASIS_POINT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;

    LiquiditySpreadRhoPointwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
    {}
};

CClassConstSP const LiquiditySpreadRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "LiquiditySpreadRhoPointwise", typeid(LiquiditySpreadRhoPointwise), load);

const double LiquiditySpreadRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

bool LiquiditySpreadRhoPointwiseLinkIn() {
    return LiquiditySpreadRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
