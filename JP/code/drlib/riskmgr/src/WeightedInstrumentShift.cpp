//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : WeightedInstrumentShift.cpp
//
//   Description : Sensitivity to a shift in instrument weight 
//                 (implemented by WeightedInstrumentCollection)
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/WeightedInstrumentShift.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
const string WeightedInstrumentShift::NAME = "WEIGHTED_INSTRUMENT_SHIFT";
const double WeightedInstrumentShift::DEFAULT_SHIFT = 0.001;
const double WeightedInstrumentShift::SENSITIVITY_UNIT = 1;

IHypothesis::AlternateWorldSP WeightedInstrumentShift::appliedTo(
        OutputNameConstSP name,
        double shiftSize,
        IObjectSP world) {
    try {
        return property()->axisFor(name, VoidSP())->hypothesis(shiftSize)->
                   appliedTo(world);
    }
    catch (exception& e) {
        throw ModelException(e, "ScalarRiskPropertySensitivity::appliedTo()");
    }
}

double WeightedInstrumentShift::getShiftSize() const {
    return shiftSize;
}

OutputNameArrayConstSP WeightedInstrumentShift::allNames(const IObject* object) const {
    return OutputNameArrayConstSP(new OutputNameArray());//Empty
}

WeightedInstrumentShift::WeightedInstrumentShift(double shiftSize):
    AllNamesRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

AllNamesRiskPropertySensitivity::Deriv WeightedInstrumentShift::deriv() const {
    smartPtr<WeightedInstrumentTweak> tag(new WeightedInstrumentTweak(
        IntArrayConstSP(new IntArray(instrumentList))));

    return Deriv(IResultsFunction::price(),
                 RiskProperty<WeightedInstrumentTweak>::SP(tag),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

WeightedInstrumentShift::~WeightedInstrumentShift() {}

void WeightedInstrumentShift::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(WeightedInstrumentShift, clazz);
    SUPERCLASS(AllNamesRiskPropertySensitivity);
    IMPLEMENTS(IScalarPerNameShift);
    EMPTY_SHELL_METHOD(&DefaultConstructor<WeightedInstrumentShift>::iObject);
    FIELD(instrumentList, "Which instruments in the InstrumentCollection to shift");
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<WeightedInstrumentShift>(), 
                                new WeightedInstrumentShift(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<WeightedInstrumentTweak>::TYPE);
}

CClassConstSP const WeightedInstrumentShift::TYPE = CClass::registerClassLoadMethod(
    "WeightedInstrumentShift", typeid(WeightedInstrumentShift), load);

bool WeightedInstrumentShiftLinkIn() {
    return WeightedInstrumentShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
