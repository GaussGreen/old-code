/**
 * @file FieldSensitivityDefinition.cpp
 */

#include <map>
#include "edginc/config.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ICrossDerivative.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include "edginc/RiskQuantityFactorySensitivity.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"
#include "edginc/CrossRiskPropertySensitivity.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

FieldSensitivityDefinition::FieldSensitivityDefinition():
    CObject(TYPE),
    unit(0.),
    defaultShiftSize(0.),
    defaultShiftSize2(0.),
    _useOverride(true),
    _requireOverride(false)
{}

FieldSensitivityDefinition::~FieldSensitivityDefinition() {}

static IScalarDerivativeConstSP deriv(const string& n) {
    if (n == "ONE_SIDED") {
        return IScalarDerivative::oneSided();
    }
    else if (n == "TWO_SIDED") {
        return IScalarDerivative::twoSided();
    }
    else if (n == "SECOND") {
        return IScalarDerivative::second();
    }
    else {
        throw ModelException(
            "\"" + n + "\" is not a valid derivative type: "
            "should be ONE_SIDED, TWO_SIDED or SECOND");
    }
}

static const char* resolution_ALLNAMES = "ALLNAMES",
                 * resolution_SCALAR = "SCALAR",
                 * resolution_POINTWISE = "POINTWISE",
                 * resolution_EXPIRYANDSTRIKE = "EXPIRYANDSTRIKE",
                 * resolution_ELEMENTWISE = "ELEMENTWISE",
                 * resolution_POINTWISExPOINTWISE = "POINTWISExPOINTWISE",
                 * resolution_SCALARxSCALAR = "SCALARxSCALAR";

void FieldSensitivityDefinition::validatePop2Object() {
    try {
        if (implementation == "DEFAULT") {
            _useOverride = true;
            _requireOverride = false;
        }
        else if (implementation == "INTERNAL") {
            _useOverride = true;
            _requireOverride = true;
        }
        else if (implementation == "FLEXIBLE") {
            _useOverride = false;
            _requireOverride = false;
        }
        else {
            throw ModelException(
                "'implementation' should be \"DEFAULT\", \"INTERNAL\" or "
                "\"FLEXIBLE\", but is \"" + implementation + "\"");
        }

        if (resolution != resolution_ALLNAMES &&
            resolution != resolution_SCALAR &&
            resolution != resolution_POINTWISE &&
            resolution != resolution_EXPIRYANDSTRIKE &&
            resolution != resolution_ELEMENTWISE &&
            resolution != resolution_POINTWISExPOINTWISE &&
            resolution != resolution_SCALARxSCALAR) {
            throw ModelException(
                "'resolution' should be \"ALLNAMES\", \"SCALAR\", "
                "\"SCALARxSCALAR\", \"POINTWISE\", \"POINTWISExPOINTWISE\" "
                "or \"ELEMENTWISE\" but is \"" + resolution + "\"");
        }

        _derivative1 = deriv(derivative1);

        if (packetName1 == "") {
            throw ModelException("'packetName1' may not be empty");
        }

        if (derivative2 == "") {
            if (packetName2 != "") {
                throw ModelException(
                    "If 'packetName2' is specified, 'derivative2' must be too");
            }
        }
        else {
            _derivative2 = deriv(derivative2);
            if (packetName2 == "") {
                throw ModelException(
                    "If 'derivative2' is specified, 'packetName2' must be too");
            }
        }

        if (derivand == "PRICE") {
            _derivand = IResultsFunction::price();
        }
        else {
            try {
                _derivand = IResultsFunction::outputRequest(
                                OutputRequest::SP(derivand));
            }
            catch (exception& e) {
                throw ModelException(e,
                    "\"" + derivand + "\" is not a valid derivand: "
                    "should be \"PRICE\" or the name of an output request");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(
            e, "FieldSensitivityDefinition::validatePop2Object()");
    }
}

static map<string,
           RiskQuantityFactorySensitivitySP (*)(const double* shiftSize,
                                                const string* packetName1,
                                                const string* packetName2)>
    builtins;

void FieldSensitivityDefinition::registerBuiltin(
        const string& name,
        RiskQuantityFactorySensitivitySP (*builtinConstructor)(
            const double* shiftSize,
            const string* packetName,
            const string* packetName2)) {
    if (builtins.find(name) != builtins.end()) {
        throw ModelException(
            "FieldSensitivityDefinition::registerBuiltin",
            "An override implementation has already been registered "
                "for the FieldSensitivityDefinition \"" + name + "\"");
    }

    builtins[name] = builtinConstructor;
}

RiskQuantityFactorySensitivitySP FieldSensitivityDefinition::sensitivity(
        OutputNameArrayConstSP toTweak,
        const double* overrideShiftSize,
        const string* overridePacketName1,
        const string* overridePacketName2) const {

    try {
        if (builtins.find(name) != builtins.end()) {
            if (_useOverride) {
                return builtins[name](overrideShiftSize,
                                      overridePacketName1, overridePacketName2);
            }
        }
        else if (_requireOverride) {
            throw ModelException(
                "'implementation' was specified as INTERNAL but no internal "
                "implementation for \"" + name + "\" exists");
        }

        double coeff = overrideShiftSize ? *overrideShiftSize
                                         : defaultShiftSize;
        string name1 = overridePacketName1 ? *overridePacketName1
                                           : packetName1;
        string name2 = overridePacketName2 ? *overridePacketName2
                                           : packetName2;

        RiskQuantityFactorySensitivitySP it;

        CDoubleConstSP one(CDouble::create(1.));

        if (resolution == resolution_SCALAR) {
            it = PerNameRiskPropertySensitivity<Void>::SP(
                 name1, name2, _derivand, property->scalarProperty(one),
                 _derivative1, _derivative2, unit, unit, coeff);
        }
        else if (resolution == resolution_POINTWISE) {
            it = PerNameRiskPropertySensitivity<ExpiryWindow>::SP(
                 name1, name2, _derivand, property->pointwiseProperty(one),
                 _derivative1, _derivative2, unit, unit, coeff);
        }
        else if (resolution == resolution_EXPIRYANDSTRIKE) {
            it = PerNameRiskPropertySensitivity<ExpiryAndStrike>::SP(
                 name1, name2, _derivand, property->expiryAndStrikewiseProperty(one),
                 _derivative1, _derivative2, unit, unit, coeff);
        }
        else if (resolution == resolution_ELEMENTWISE) {
            it = PerNameRiskPropertySensitivity<BoxedInt>::SP(
                 name1, name2, _derivand, property->elementwiseProperty(one),
                 _derivative1, _derivative2, unit, unit, coeff);
        }
        else if (resolution == resolution_ALLNAMES) {
            it = AllNamesRiskPropertySensitivity::SP(
                 name1, name2, _derivand, property->scalarProperty(one),
                 _derivative1, _derivative2, unit, unit, coeff);
        }
        else if (resolution == resolution_SCALARxSCALAR) {
            it = CrossRiskPropertySensitivity<Void, Void>::SP(
                 name1, _derivand,
                 property->scalarProperty(one),
                 (!property2 ? property : property2)->scalarProperty(one),
                 ICrossDerivative::cross(), unit,
                 coeff,
                 defaultShiftSize2 == 0. ? coeff : defaultShiftSize2);
        }
        else {
            ASSERT(resolution == resolution_POINTWISExPOINTWISE);
            it = CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>::SP(
                 name1, _derivand,
                 property->pointwiseProperty(one),
                 (!property2 ? property : property2)->pointwiseProperty(one),
                 ICrossDerivative::cross(), unit,
                 coeff,
                 defaultShiftSize2 == 0. ? coeff : defaultShiftSize2);
        }

        it->storeOverrideNames(
            OutputNameArraySP(!toTweak ? 0 : toTweak.clone()));

        return it;
    }
    catch (exception& e) {
        // If we throw 'e' out, the whole pricing run will fail ---
        // better to show the error in the results

        return RiskQuantityFactorySensitivity::singleton(
             name,
             NamedRiskQuantity::SP(
                 RiskQuantity::untweakable(ModelException(
                     e,
                     "FieldSensitivityDefinition::sensitivity()",
                     "Constructing sensitivity from "
                         "FieldSensitivityDefinition \"" + name + "\"")),
                 IResultsIdentifier::SP(name, OutputNameSP())));
    }
}

IObject* FieldSensitivityDefinition::emptyShell() {
    return new FieldSensitivityDefinition();
}

void FieldSensitivityDefinition::load(CClassSP& clazz) {
    REGISTER(FieldSensitivityDefinition, clazz);
    clazz->setPublic();
    SUPERCLASS(CObject);
    FIELD(name, "Unique name for the sensitivity");
    FIELD(implementation, "FLEXIBLE, INTERNAL, DEFAULT");
    FIELD(derivand, "What the greek is a derivative of, e.g. PRICE");
    FIELD(property, "What the greek is a derivative with respect to");
    FIELD(property2, "Second property with respect to which to take deriv, for cross-gammas");
    FIELD_MAKE_OPTIONAL(property2);
    FIELD(resolution, "SCALAR, POINTWISE, ELEMENTWISE, POINTWISExPOINTWISE, or SCALARxSCALAR");
    FIELD(unit, "Units in which to report greek, e.g. 0.0001 for \"per bp\"");
    FIELD(derivative1, "Deriv to compute: ONE_SIDED, TWO_SIDED or SECOND");
    FIELD(packetName1, "Packet name under which to report the greek");
    FIELD(derivative2, "Deriv to compute in addition to 'derivative1'");
    FIELD_MAKE_OPTIONAL(derivative2);
    FIELD(packetName2, "Packet name under which to report derivative2");
    FIELD_MAKE_OPTIONAL(packetName2);
    FIELD(defaultShiftSize, "Shift size for estimating derivative1");
    FIELD(defaultShiftSize2, "Shift size for estimating derivative2");
    FIELD_MAKE_OPTIONAL(defaultShiftSize2);
    FIELD(_derivative1, "_derivative1");
    FIELD_MAKE_TRANSIENT(_derivative1);
    FIELD(_derivative2, "_derivative2");
    FIELD_MAKE_TRANSIENT(_derivative2);
    FIELD(_derivand, "_derivand");
    FIELD_MAKE_TRANSIENT(_derivand);
    FIELD(_useOverride, "_useOverride");
    FIELD_MAKE_TRANSIENT(_useOverride);
    FIELD(_requireOverride, "_requireOverride");
    FIELD_MAKE_TRANSIENT(_requireOverride);
    EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const FieldSensitivityDefinition::TYPE = CClass::registerClassLoadMethod(
    "FieldSensitivityDefinition", typeid(FieldSensitivityDefinition), load);

DRLIB_END_NAMESPACE
