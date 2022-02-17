//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CDOParallelStrikeShift.cpp
//
//   Description : Sensitivity to a parallel shift in CDO strikes
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDOParallelStrikeShift.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CDOParallelStrikeShift::NAME = "CDO_PARALLEL_STRIKE_SHIFT";
const double CDOParallelStrikeShift::DEFAULT_SHIFT = 0.001;
const double CDOParallelStrikeShift::SENSITIVITY_UNIT = 1;

IHypothesis::AlternateWorldSP CDOParallelStrikeShift::appliedTo(
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

double CDOParallelStrikeShift::getShiftSize() const {
    return shiftSize;
}

OutputNameArrayConstSP CDOParallelStrikeShift::allNames(const IObject* object) const {
    return OutputNameArrayConstSP(new OutputNameArray());//Empty
}

CDOParallelStrikeShift::CDOParallelStrikeShift(double shiftSize):
    AllNamesRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

AllNamesRiskPropertySensitivity::Deriv CDOParallelStrikeShift::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CDOParallelStrike>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CDOParallelStrikeShift::~CDOParallelStrikeShift() {}

void CDOParallelStrikeShift::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDOParallelStrikeShift, clazz);
    SUPERCLASS(AllNamesRiskPropertySensitivity);
    IMPLEMENTS(IScalarPerNameShift);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CDOParallelStrikeShift>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CDOParallelStrikeShift>(), 
                                new CDOParallelStrikeShift(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CDOParallelStrike>::TYPE);
}

CClassConstSP const CDOParallelStrikeShift::TYPE = CClass::registerClassLoadMethod(
    "CDOParallelStrikeShift", typeid(CDOParallelStrikeShift), load);

bool CDOParallelStrikeShiftLinkIn() {
    return CDOParallelStrikeShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
