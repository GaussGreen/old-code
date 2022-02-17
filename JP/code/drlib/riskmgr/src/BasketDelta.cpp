/**
 * @file BasketDelta.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Delta.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/IQualifiedRiskAxis.hpp"

DRLIB_BEGIN_NAMESPACE

static const char* NAME = "DELTA_BASKET";
const double BasketDelta::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

BasketDelta::BasketDelta(double shiftSize):
    ScalarShiftTwoSided(TYPE, false, 1., NAME,
                        new RiskProperty<BasketSpot>(), shiftSize),
    useDeltaNames(false)
{}

BasketDelta::BasketDelta(double shiftSize,
                         IModel* model, CControl* control, bool useDeltaNames):
    ScalarShiftTwoSided(TYPE, false, 1.,
                        useDeltaNames ? Delta::NAME : NAME,
                        new RiskProperty<BasketSpot>(), shiftSize),
    useDeltaNames(useDeltaNames)
{
    this->algorithm = model;
    this->control = control;
}

BasketDelta::~BasketDelta() {}

double BasketDelta::divisor() const{
    static const char routine[] = "BasketDelta::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try{
            initSpot = getInitialValue();
        } catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "delta calculation");
            throw e;
        }
        // then just scale by the shift size
        shiftSize = getShiftSize();
        if (fabs(shiftSize) < DBL_EPSILON){
            throw ModelException(routine, "Shift size is zero");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return (shiftSize * initSpot);
}

const string& BasketDelta::getSecondOrderSensOutputName() const {
    static string SECOND_ORDER_NAME("GAMMA_BASKET");
    return useDeltaNames ? Delta::SECOND_ORDER_NAME : SECOND_ORDER_NAME;
}

static IObject* defaultBasketDelta() {
    return new BasketDelta(BasketDelta::DEFAULT_SHIFT);
}

void BasketDelta::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(BasketDelta, clazz);
    SUPERCLASS(ScalarShiftTwoSided);
    EMPTY_SHELL_METHOD(defaultBasketDelta);
    FIELD(useDeltaNames, "use Delta Names");
    FIELD_MAKE_OPTIONAL(useDeltaNames);
    SensitivityFactory::addSens(NAME, 
                                new GenericSensitivityFactory<BasketDelta>(), 
                                new BasketDelta(BasketDelta::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<BasketSpot>::TYPE);
}

CClassConstSP const BasketDelta::TYPE = CClass::registerClassLoadMethod(
    "BasketDelta", typeid(BasketDelta), load);

/**
 * Included in RiskMgrLib::linkInClasses() to force BasketDelta
 * to get linked into the Windows exe.
 */

bool BasketDeltaLinkIn() {
    return BasketDelta::TYPE != NULL;
}

DRLIB_END_NAMESPACE
