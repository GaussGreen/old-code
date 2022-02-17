/**
 * @file AllNamesRiskPropertySensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/Results.hpp"
#include "edginc/Additive.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

AllNamesRiskPropertySensitivity::AllNamesRiskPropertySensitivity(
        CClassConstSP type,
        const string& outputName,
        double shiftSize):
    RiskPropertySensitivity<Void>(type, shiftSize, outputName, "")
{}

AllNamesRiskPropertySensitivity::AllNamesRiskPropertySensitivity(
        CClassConstSP type,
        const string& outputName,
        const string& outputName2,
        double shiftSize):
    RiskPropertySensitivity<Void>(type, shiftSize, outputName, outputName2)
{}

AllNamesRiskPropertySensitivity::~AllNamesRiskPropertySensitivity() {}

NamedRiskQuantityArraySP AllNamesRiskPropertySensitivity::nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {

    NamedRiskQuantityArraySP rqs = NamedRiskQuantityArray::SP(derivatives()->size());
    IRiskAxisConstSP axis = property()->axisFor(OutputNameSP(), VoidSP());

    for (int d = 0; d < derivatives()->size(); ++d) {
        TRACE_BLOCK("Creating NamedRiskQuantity("
            << getPacketName() << ", "
            << (d == 0 ? getSensOutputName() : (*outputNames)[d])
            << ")");
        (*rqs)[d] = NamedRiskQuantity::SP(
            (*derivatives())[d]->riskQuantity(derivand(), axis, shiftSize),
            IResultsIdentifier::SP(
            getPacketName(),
            d == 0 ? getSensOutputName() : (*outputNames)[d]),
        pow((*sensitivityUnits())[d], (*derivatives())[d]->order()));
    }
    return rqs;
}

const string& AllNamesRiskPropertySensitivity::getPacketName() const {
    return Results::INSTRUMENT_PACKET;
}

void AllNamesRiskPropertySensitivity::load(CClassSP& clazz) {
    REGISTER(AllNamesRiskPropertySensitivity, clazz);
    SUPERCLASS(RiskPropertySensitivity<Void>);
}

CClassConstSP const AllNamesRiskPropertySensitivity::TYPE = CClass::registerClassLoadMethod(
    "AllNamesRiskPropertySensitivity", typeid(AllNamesRiskPropertySensitivity), load);

class GeneralAllNamesRiskPropertySensitivity:
        public AllNamesRiskPropertySensitivity,
        public virtual Additive {

    static CClassConstSP const TYPE;

private:

    GeneralAllNamesRiskPropertySensitivity():
        AllNamesRiskPropertySensitivity(TYPE, "", 0.) {}

    static IObject* emptyShell() {
        return new GeneralAllNamesRiskPropertySensitivity();
    }        

    static void load(CClassSP& clazz) {
        REGISTER(GeneralAllNamesRiskPropertySensitivity, clazz);
        SUPERCLASS(AllNamesRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(emptyShell);
        IMPLEMENTS(Additive);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunctionSP(),
                     smartPtr<IRiskProperty<Void> >(),
                     IScalarDerivativeSP(),
                     0.);
    }

public:

    GeneralAllNamesRiskPropertySensitivity(
            const string& outputName1,
            const string& outputName2,
            IResultsFunctionConstSP derivand,
            smartConstPtr<IRiskProperty<Void> > property,
            IScalarDerivativeConstSP derivative1,
            IScalarDerivativeConstSP derivative2,
            double sensitivityUnit1,
            double sensitivityUnit2,
            double shiftSize):
        AllNamesRiskPropertySensitivity(
            TYPE, outputName1, !derivative2 ? string() : outputName2, shiftSize)
    {
        _derivand = derivand;
        _property = property;
        _derivatives.reset(new IScalarDerivativeArray(
            1, IScalarDerivativeSP::constCast(derivative1)));
        _sensitivityUnits.reset(new DoubleArray(1, sensitivityUnit1));

        if (!!derivative2) {
            IScalarDerivativeArraySP::constCast(_derivatives)->push_back(
                IScalarDerivativeSP::constCast(derivative2));
            DoubleArraySP::constCast(_sensitivityUnits)->push_back(
                sensitivityUnit2);
        }
    }
};

CClassConstSP const GeneralAllNamesRiskPropertySensitivity::TYPE = CClass::registerClassLoadMethod(
    "GeneralAllNamesRiskPropertySensitivity", typeid(GeneralAllNamesRiskPropertySensitivity), load);

AllNamesRiskPropertySensitivitySP AllNamesRiskPropertySensitivity::SP(
        const string& outputName,
        IResultsFunctionConstSP derivand,
        IScalarRiskPropertyConstSP property,
        IScalarDerivativeConstSP derivative,
        double sensitivityUnit,
        double shiftSize) {
    return AllNamesRiskPropertySensitivitySP(
        new GeneralAllNamesRiskPropertySensitivity(
            outputName, "", derivand, property,
            derivative, IScalarDerivativeConstSP(), sensitivityUnit, 0, shiftSize));
}

AllNamesRiskPropertySensitivitySP AllNamesRiskPropertySensitivity::SP(
        const string& outputName1,
        const string& outputName2,
        IResultsFunctionConstSP derivand,
        IScalarRiskPropertyConstSP property,
        IScalarDerivativeConstSP derivative1,
        IScalarDerivativeConstSP derivative2,
        double sensitivityUnit1,
        double sensitivityUnit2,
        double shiftSize) {
    return AllNamesRiskPropertySensitivitySP(
        new GeneralAllNamesRiskPropertySensitivity(
            outputName1, outputName2, derivand, property,
            derivative1, derivative2,
            sensitivityUnit1, sensitivityUnit2, shiftSize));
}

DRLIB_END_NAMESPACE
