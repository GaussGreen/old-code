/**
 * @file ScalarRiskPropertySensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

ScalarRiskPropertySensitivity::ScalarRiskPropertySensitivity(
        CClassConstSP type,
        const string& outputName,
        double shiftSize):
    PerNameRiskPropertySensitivity<Void>(type, shiftSize, outputName, "")
{}

ScalarRiskPropertySensitivity::ScalarRiskPropertySensitivity(
        CClassConstSP type,
        const string& outputName,
        const string& outputName2,
        double shiftSize):
    PerNameRiskPropertySensitivity<Void>(type, shiftSize, outputName, outputName2)
{}

ScalarRiskPropertySensitivity::~ScalarRiskPropertySensitivity() {}

bool ScalarRiskPropertySensitivity::findAndShift(
        IObjectSP objectToShift, OutputNameConstSP name) {
    try {
        bool changed;
        property()->axisFor(name, VoidSP())->hypothesis(shiftSize)->
            applyTo(objectToShift, &changed);
        return changed;
    }
    catch (exception& e) {
        throw ModelException(e,
                             "ScalarRiskPropertySensitivity::findAndShift()");
    }
}

IHypothesis::AlternateWorldSP ScalarRiskPropertySensitivity::appliedTo(
        OutputNameConstSP name,
        double shiftSize,
        IObjectSP world) {
    try {
        return property()->axisFor(name, VoidSP())->hypothesis(shiftSize)->
                   appliedTo(world);
    }
    catch (exception& e) {
        throw ModelException(e,
                             "ScalarRiskPropertySensitivity::appliedTo()");
    }
}

double ScalarRiskPropertySensitivity::getShiftSize() const {
    return shiftSize;
}

OutputNameArrayConstSP ScalarRiskPropertySensitivity::allNames(
        const IObject* object) const {
    return PerNameRiskPropertySensitivity<Void>::allNames(object);
}

void ScalarRiskPropertySensitivity::load(CClassSP& clazz) {
    REGISTER(ScalarRiskPropertySensitivity, clazz);
    SUPERCLASS(PerNameRiskPropertySensitivity<Void>);
    IMPLEMENTS(IPerturbation);
    IMPLEMENTS(IScalarPerNameShift);
}

CClassConstSP const ScalarRiskPropertySensitivity::TYPE = CClass::registerClassLoadMethod(
    "ScalarRiskPropertySensitivity", typeid(ScalarRiskPropertySensitivity), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
