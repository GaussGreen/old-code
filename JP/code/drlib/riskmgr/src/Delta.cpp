/**
 * @file Delta.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/Spot.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/Delta.hpp"

DRLIB_BEGIN_NAMESPACE

const string Delta::NAME = "DELTA";
const string Delta::SECOND_ORDER_NAME = "GAMMA";
const double Delta::DEFAULT_SHIFT = 0.005;
const double Delta::MINIMUM_SHIFT = 0.0001;

Delta::Delta(double shiftSize):
    ScalarShiftTwoSided(TYPE, false, 1., NAME,
                        new RiskProperty<Spot>(), shiftSize)
{}

Delta::Delta(double shiftSize, IModel* model, CControl* control):
    ScalarShiftTwoSided(TYPE, false, 1., NAME,
                        new RiskProperty<Spot>(), shiftSize)
{
    this->algorithm = model;
    this->control = control;
}

Delta::~Delta() {}

const string& Delta::getSecondOrderSensOutputName() const {
    return SECOND_ORDER_NAME;
}

void Delta::calculate(TweakGroup*      tweakGroup,
                      CResults*        results){
    try{
        calculateTwoSidedDeriv(getSecondOrderSensOutputName(), tweakGroup,
                               results);
        BasketDeltaSP basketDelta(new BasketDelta(getShiftSize(),
                                                  getModel(),
                                                  getControl(),
                                                  true));
        basketDelta->calculate(tweakGroup, results);
    } catch (exception& e){
        throw ModelException(&e,  "Delta::calculate");
    }
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double Delta::divisor() const{
    static const char routine[] = "Delta::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try {
            initSpot = getInitialValue();            
        }
        catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "delta calculation");
            throw e;
        }

        if (Maths::isZero(initSpot)) {
            throw ModelException(routine, "spot price is zero");
        }
        // then just scale by the shift size
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(routine, "Shift size is zero");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return (shiftSize * initSpot);
}

NamedRiskQuantityArraySP Delta::riskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {

    NamedRiskQuantityArraySP rqs =
         ScalarShiftTwoSided::riskQuantities(world, riskMapping);

    BasketDeltaSP basketDelta(new BasketDelta(getShiftSize(), algorithm,
                                              control, true));

    NamedRiskQuantityArraySP brqs =
        basketDelta->riskQuantities(world, riskMapping);
    
    int nrqs = rqs->size();

    for (int b = 0; b < brqs->size(); ++b) {
        const IResultsIdentifier* bname = (*brqs)[b]->resultsName.get();
        if ((bname->packet() == NAME || bname->packet() == SECOND_ORDER_NAME) &&
                !bname->entry()) {
            continue;
        }
        int r;
        for (r = 0; r < nrqs; ++r) {
            if ((*rqs)[r]->resultsName->equalTo(bname)) break;
        }
        if (r == nrqs) {
            rqs->push_back((*brqs)[b]);
        }
    }

    return rqs;
}

static IObject* defaultDelta() {
    return new Delta(Delta::DEFAULT_SHIFT);
}

void Delta::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(Delta, clazz);
    SUPERCLASS(ScalarShiftTwoSided);
    EMPTY_SHELL_METHOD(defaultDelta);
    SensitivityFactory::addSens(Delta::NAME, 
                                new GenericSensitivityFactory<Delta>(), 
                                new Delta(Delta::DEFAULT_SHIFT),
                                ITweakableWithRespectTo<Spot>::TYPE);
}

CClassConstSP const Delta::TYPE = CClass::registerClassLoadMethod(
    "Delta", typeid(Delta), load);

bool DeltaLinkIn() {
    return Delta::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
