/**
 * @file IFieldRiskPropertyDefinition.cpp
 */

#include <map>
#include "edginc/config.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/**
 * @name Mechanism for overriding externally defined flexible properties with
 * internal hardwired ones
 */

//@{

/**
 * The type ("scalar", "parallel", "pointwise", "multi-element", "elementwise")
 * of each registered FieldRiskProperty::name
 */

static map<string, string> builtinTypes;

/**
 * Register against FieldRiskProperty::name of scalarProperty() overrides
 */

static map<string,
           IScalarRiskPropertyConstSP (*)(IObjectConstSP)>
    scalarBuiltins;

/**
 * Register against FieldRiskProperty::name of parallelProperty() overrides
 */

static map<string,
           IScalarRiskPropertyConstSP (*)(IArrayConstSP,
                                          IArrayConstSP)>
    parallelBuiltins;

/**
 * Register against FieldRiskProperty::name of pointwiseProperty() overrides
 */

static map<string,
           IExpiryRiskPropertyConstSP (*)(IObjectConstSP)>
    pointwiseBuiltins;

/**
 * Register against FieldRiskProperty::name of expiryAndStrikewiseProperty() overrides
 */

static map<string,
           IExpiryAndStrikeRiskPropertyConstSP (*)(IObjectConstSP)>
    expiryAndStrikewiseBuiltins;

/**
 * Register against FieldRiskProperty::name of elementsProperty() overrides
 */

static map<string,
           IScalarRiskPropertyConstSP (*)(IntArrayConstSP indices,
                                          IArrayConstSP)>
    elementsBuiltins;

/**
 * Register against FieldRiskProperty::name of elementwiseProperty() overrides
 */

static map<string,
           IIntRiskPropertyConstSP (*)(IObjectConstSP)>
    elementwiseBuiltins;

/**
 * Checks for common error cases
 */

static void noOverride(string type, string name, bool requireOverride) {
    if (builtinTypes.find(name) != builtinTypes.end()) {
        throw ModelException(
            "An internal implementation for \"" + name + "\" "
            "exists, but it is " + builtinTypes[name] +
            ", not " + type);
    }
    else if (requireOverride) {
        throw ModelException(
            "'implementation' was specified as INTERNAL but no "
            "internal implementation for \"" + name + "\" exists");
    }
}

template <class T>
static smartConstPtr<T> override(
        const map<string, smartConstPtr<T>(*)(IObjectConstSP) >& reg,
        const string& name,
        const string& type,
        IObjectConstSP arg,
        bool requireOverride) {
    try {
        if (reg.find(name) != reg.end()) {
            return reg.find(name)->second(arg);
        }

        if (name.size() > 1 && name[0] == '=') {
            try {
                CClassConstSP c = CClass::forName(name.substr(1));

                if (!T::TYPE->isAssignableFrom(c)) {
                    throw ModelException(
                        "Expected a " + T::TYPE->getName() + " because the "
                        "property is specified to be " + type + ", but got a " +
                        c->getName() + " which isn't");
                }

                IObjectConstSP o(c->newInstance());
                ASSERT(!!o);

                return smartConstPtr<T>::dynamicCast(o);
            }
            catch (exception& e) {
                throw ModelException(e,
                    "Constructing an instance of " + name.substr(1) +
                    " as internal FieldRiskProperty implementation");
            }
        }

        noOverride(type, name, requireOverride);

        return smartConstPtr<T>();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

IFieldRiskPropertyDefinition::IFieldRiskPropertyDefinition() {}
IFieldRiskPropertyDefinition::~IFieldRiskPropertyDefinition() {}

static void setBuiltinType(string name, string type) {
    if (builtinTypes.find(name) != builtinTypes.end()) {
        throw ModelException(
            string("An override implementation (") +
            (type == builtinTypes[name] ? "also " : "") +
            "of type " + builtinTypes[name] + ")"
            "has already been registered for the "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }

    builtinTypes[name] = type;
}

void IFieldRiskPropertyDefinition::registerScalarBuiltin(
        const string& name,
        IScalarRiskPropertyConstSP (*builtin)(IObjectConstSP)) {
    setBuiltinType(name, "scalar");
    scalarBuiltins[name] = builtin;
}

void IFieldRiskPropertyDefinition::registerParallelBuiltin(
        const string& name,
        IScalarRiskPropertyConstSP (*builtin)(IArrayConstSP, IArrayConstSP)) {
    setBuiltinType(name, "parallel");
    parallelBuiltins[name] = builtin;
}

void IFieldRiskPropertyDefinition::registerPointwiseBuiltin(
        const string& name,
        IExpiryRiskPropertyConstSP (*builtin)(IObjectConstSP)) {
    setBuiltinType(name, "pointwise");
    pointwiseBuiltins[name] = builtin;
}

void IFieldRiskPropertyDefinition::registerExpiryAndStrikewiseBuiltin(
        const string& name,
        IExpiryAndStrikeRiskPropertyConstSP (*builtin)(IObjectConstSP)) {
    setBuiltinType(name, "expiry-and-strikewise");
    expiryAndStrikewiseBuiltins[name] = builtin;
}

void IFieldRiskPropertyDefinition::registerElementsBuiltin(
        const string& name, 
        IScalarRiskPropertyConstSP (*builtin)(IntArrayConstSP, IArrayConstSP)) {
    setBuiltinType(name, "multi-element");
    elementsBuiltins[name] = builtin;
}

void IFieldRiskPropertyDefinition::registerElementwiseBuiltin(
        const string& name,
        IIntRiskPropertyConstSP (*builtin)(IObjectConstSP)) {
    setBuiltinType(name, "elementwise");
    elementwiseBuiltins[name] = builtin;
}

//@}

/**
 * What NumericFieldRiskPropertyDefinition and FieldRiskPropertyDefinition have
 * in common
 */

class BaseFieldRiskPropertyDefinition: public CObject,
                                       public virtual IFieldRiskPropertyDefinition {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

protected:

    BaseFieldRiskPropertyDefinition(CClassConstSP type):
        CObject(type),
        qualifierDimension(0),
        beforeGetMarket(false),
        _qualifiersAreIndices(false),
        _useOverride(true),
        _requireOverride(false)
    {}

    /**
     * A name for the property
     *
     * See IFieldRiskPropertyDefinition::registerScalarBuiltin() and friends.
     */

    string name; // $required(Unique name for the property)

    /**
     * Whether to use a hard-coded implementation of the property, or follow
     * this "flexible" definition
     *
     * FLEXIBLE, INTERNAL or DEFAULT.  See
     * IFieldRiskPropertyDefinition::registerScalarBuiltin() and friends.
     */

    string implementation; // $required(FLEXIBLE, INTERNAL or DEFAULT)

    /**
     * The type of object which carries the property
     */

    string type; // $required(Name of type to be tweaked)

    /**
     * The field on 'type' in which the property is stored, or (if more than
     * one entry) a chain of fields leading from 'type' to the field
     *
     * Maps to the @a field constructor argument to e.g.
     * PointwiseFieldRiskProperty::PointwiseFieldRiskProperty()
     */

    CStringArrayConstSP field; // $required(Path to field to tweak, starting from 'type')

    /**
     * If the property is array-valued, 'field' will lead us to an array rather
     * than an object; if the array itself has object values rather than just
     * doubles, we need to know which field on those objects to edit.
     */

    CStringArrayConstSP subField; // $optional(For tweaking object arrays, the field within each object that we should actually tweak)

    /**
     * Optionally: the field on 'type' in which the expiry labels for the
     * property's array value is stored, or (if more than one entry) a
     * chain of fields leading from 'type' to the field
     *
     * Maps to the @a expiriesField constructor argument to e.g.
     * PointwiseFieldRiskProperty::PointwiseFieldRiskProperty().  If it's
     * omitted, then when parallelProperty() or pointwiseProperty() are
     * called, we'll look for a @a getExpiriesMethod set for 'field' via
     * Calibrator::IAdjustable::registerBootstrappableField before
     * giving up and throwing an exception.
     */

    CStringArrayConstSP qualifierField; // $optional(Path to field in which expiry labels are stored, starting from 'type')

    /**
     * If qualifierField leads us to an object array rather than an ExpiryArray
     * we need to know which field holds the actual Expiry
     */

    CStringArrayConstSP qualifierSubField; // $optional(If the expiry labels are stored as a field on the objects in an object array rather than just as an ExpiryArray)

    /**
     * Not actually used yet
     */

    CStringArrayConstSP qualifier2Field; // $optional

    /**
     * Not actually used yet
     */

    CStringArrayConstSP qualifier2SubField; // $optional

    /**
     * Which dimension of the matrix-valued @a field is indexed by Expiry
     *
     * Maps to the @a expiryQualifier constructor argument to e.g.
     * PointwiseFieldRiskProperty::PointwiseFieldRiskProperty()
     */

    int qualifierDimension; // $optional(For matrix-valued fields: which dimension corresponds to 'qualifierField')

    /**
     * The scale in which changes to the property should be measured
     */

    string measure; // $required(Terms for reporting deriv: NOMINAL or ABSOLUTE)

    /**
     * When the property is used to generate scenarios: whether the scenario is
     * to be applied before or after the "getMarket" phase
     */

    bool beforeGetMarket; // $optional(Whether to apply scenarios before the getMarket phase)

    bool _qualifiersAreIndices;           // $transient
    bool _useOverride, _requireOverride;  // $transient
    bool _absoluteDistance;               // $transient

    void validatePop2Object();

    virtual IFieldTweak::IOperatorConstSP operat0r() const = 0;

public:

    ~BaseFieldRiskPropertyDefinition() {}

    IScalarRiskPropertyConstSP scalarProperty(IObjectConstSP argument) const;

    IScalarRiskPropertyConstSP parallelProperty(
        IArrayConstSP qualifiers,
        IArrayConstSP arguments) const;

    IExpiryRiskPropertyConstSP pointwiseProperty(IObjectConstSP argument) const;

    IExpiryAndStrikeRiskPropertyConstSP expiryAndStrikewiseProperty(
        IObjectConstSP argument) const;

    IScalarRiskPropertyConstSP elementsProperty(
        IntArrayConstSP indices,
        IArrayConstSP arguments) const;

    IIntRiskPropertyConstSP elementwiseProperty(IObjectConstSP argument) const;

    bool applyBeforeGetMarket() const { return beforeGetMarket; }
};

void BaseFieldRiskPropertyDefinition::validatePop2Object() {
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

        if (field->empty()) {
            throw ModelException("'field' must be a non-empty field path");
        }

        if (!!subField && subField->empty()) {
            throw ModelException("'subField' must be a non-empty field path");
        }

        if (!!qualifierField) {
            if (qualifierField->empty()) {
                throw ModelException(
                    "'qualifierField' must be a non-empty field path");
            }

            if ((*qualifierField)[0] == "#") {
                _qualifiersAreIndices = true;
            }
        }

        if (!!qualifierSubField && qualifierSubField->empty()) {
            throw ModelException(
                "'qualifierSubField' must be a non-empty field path");
        }

        if (measure == "NOMINAL") {
            _absoluteDistance = false;
        }
        else if (measure == "ABSOLUTE") {
            _absoluteDistance = true;
        }
        else {
            throw ModelException(
                "'measure' should be \"NOMINAL\" or \"ABSOLUTE\", but is \"" +
                measure + "\"");
        }

        // We can get some more validation by doing the following ---
        // in particular detecting invalid overrides --- but we deliberately
        // don't so as to defer errors till later and get them in the
        // results as Untweakable ...

        // scalarProperty();
    }
    catch (ModelException& e) {
        throw ModelException(e, __FUNCTION__,
            "Validating FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

static FieldPathSP pathOrNULL(CStringArrayConstSP fields) {
    return !fields ? FieldPathSP() : FieldPath::SP(fields);
}

static FieldPathSP pathOrEmpty(CStringArrayConstSP fields) {
    return !fields ? FieldPath::SP(CStringArray::SP()) :
                     FieldPath::SP(fields);
}

IScalarRiskPropertyConstSP BaseFieldRiskPropertyDefinition::scalarProperty(
        IObjectConstSP argument) const {
    try {
        if (_useOverride) {
            IScalarRiskPropertyConstSP o = override(
                scalarBuiltins, name, "scalar", argument, _requireOverride);
            if (!!o) return o;
        }

        return fieldRiskProperty::scalar(
            CClass::forName(type),
            IFieldTweak::bulk(FieldPath::SP(field), pathOrNULL(subField),
                              argument, operat0r()),
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing scalar risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

IScalarRiskPropertyConstSP BaseFieldRiskPropertyDefinition::parallelProperty(
        IArrayConstSP qualifiers,
        IArrayConstSP arguments) const {
    try {
        if (_qualifiersAreIndices) {
            IntArraySP inds = IntArray::SP(qualifiers->getLength());
            try {
                for (int q = 0; q < inds->size(); ++q) {
                    inds->set(q, IObjectSP::constCast(qualifiers->get(q)));
                }
            }
            catch (exception& e) {
                throw ModelException(e, "Translating qualifier list to "
                                        "integer index list");
            }

            return elementsProperty(inds, arguments);
        }

        if (_useOverride) {
            if (parallelBuiltins.find(name) != parallelBuiltins.end()) {
                return parallelBuiltins[name](qualifiers, arguments);
            }

            noOverride("parallel", name, _requireOverride);
        }

        return fieldRiskProperty::scalar(
            CClass::forName(type),
            IFieldTweak::elementwise(
                FieldPath::SP(field), pathOrNULL(subField),
                IFieldTweak::IIndices::lookupOrCalibratorExpiries(
                    FieldPath::SP(field),
                    pathOrNULL(qualifierField),
                    pathOrNULL(qualifierSubField),
                    qualifiers),
                qualifierDimension,
                arguments,
                operat0r()),
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing parallel risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

IExpiryRiskPropertyConstSP BaseFieldRiskPropertyDefinition::pointwiseProperty(
        IObjectConstSP argument) const {
    try {
        if (_useOverride) {
            IExpiryRiskPropertyConstSP o = override(
                pointwiseBuiltins, name, "pointwise", argument,
                _requireOverride);
            if (!!o) return o;
        }

        if (_qualifiersAreIndices) {
            throw ModelException("Can't be a pointwise property because "
                                 "'qualifier' is \"#\"");
        }

        return fieldRiskProperty::pointwise(
            CClass::forName(type),
            FieldPath::SP(field), pathOrNULL(subField),
            pathOrNULL(qualifierField), pathOrEmpty(qualifierSubField),
            qualifierDimension,
            operat0r(),
            argument,
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing pointwise risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

IExpiryAndStrikeRiskPropertyConstSP BaseFieldRiskPropertyDefinition::expiryAndStrikewiseProperty(
        IObjectConstSP argument) const {
    try {
        if (_useOverride) {
            IExpiryAndStrikeRiskPropertyConstSP o = override(
                expiryAndStrikewiseBuiltins, name, "expiry-and-strikewise",
                argument, _requireOverride);
            if (!!o) return o;
        }

        if (_qualifiersAreIndices) {
            throw ModelException("Can't be an expiry-and-strikewise property "
                                 "because 'qualifier' is \"#\"");
        }

        return fieldRiskProperty::expiryAndStrikewise(
            CClass::forName(type),
            FieldPath::SP(field), pathOrNULL(subField),
            pathOrNULL(qualifierField), pathOrEmpty(qualifierSubField),
            qualifierDimension,
            operat0r(),
            argument,
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing expiry-and-strikewise risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

IScalarRiskPropertyConstSP BaseFieldRiskPropertyDefinition::elementsProperty(
        IntArrayConstSP indices,
        IArrayConstSP arguments) const {
    try {
        if (_useOverride) {
            if (elementsBuiltins.find(name) != elementsBuiltins.end()) {
                return elementsBuiltins[name](indices, arguments);
            }

            noOverride("multi-element", name, _requireOverride);
        }

        return fieldRiskProperty::scalar(
            CClass::forName(type),
            IFieldTweak::elementwise(
                FieldPath::SP(field), pathOrNULL(subField),
                IFieldTweak::IIndices::direct(indices),
                qualifierDimension,
                arguments,
                operat0r()),
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing multi-element risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

IIntRiskPropertyConstSP BaseFieldRiskPropertyDefinition::elementwiseProperty(
        IObjectConstSP argument) const {
    try {
        if (_useOverride) {
            IIntRiskPropertyConstSP o = override(
                elementwiseBuiltins, name, "elementwise", argument,
                _requireOverride);
            if (!!o) return o;
        }

        return fieldRiskProperty::elementwise(
            CClass::forName(type),
            FieldPath::SP(field), pathOrNULL(subField),
            qualifierDimension,
            operat0r(),
            argument,
            _absoluteDistance);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Constructing elementwise risk property from "
            "FieldRiskPropertyDefinition \"" + name + "\"");
    }
}

/**
 * @name FieldRiskPropertyDefinition
 */

struct FieldRiskPropertyDefinition: BaseFieldRiskPropertyDefinition { // $public

    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;
    static IObject* emptyShell();

    string xml1;                // $optional
    string xml2;                // $optional
    string xml3;                // $optional
    string xml4;                // $optional

    FieldRiskPropertyDefinition():
        BaseFieldRiskPropertyDefinition(TYPE)
    {}

    IFieldTweak::IOperatorConstSP operat0r() const {
        string xml = xml1 + xml2 + xml3 + xml4;
        return xml.empty() ? IFieldTweak::IOperator::setter() :
                             IFieldTweak::IOperator::xmlSetter(xml);
    }
};

/**
 * @name NumericFieldRiskPropertyDefinition
 */

struct NumericFieldRiskPropertyDefinition:
        BaseFieldRiskPropertyDefinition { // $public

    static CClassConstSP const TYPE;

private:

    static IObject* emptyShell();
    static void load(CClassSP&);

    NumericFieldRiskPropertyDefinition():
        BaseFieldRiskPropertyDefinition(TYPE),
        clip(false)
    {}

    /**
     * The operation to be applied when the property is tweaked
     *
     * "+", "*", "AUTO" or "=".
     *
     * Together with floor/floorInclusive/cap/capInclusive,
     * maps to the @a tweak constructor argument to e.g.
     * PointwiseFieldRiskProperty::PointwiseFieldRiskProperty().
     */

    string operator_; // $required(How to tweak: 'AUTO', '+', '*', '=')

    /**
     * Low end of valid range for field (defaults to -inf)
     */

    CDoubleConstSP floor; // $optional(Level below which 'field' is invalid)

    /**
     * Whether 'floor' is included in the field's valid range
     */

    CBoolConstSP floorInclusive; // $optional(Is field = 'floor' valid?)

    /**
     * Upper end of valid range for field (defaults to +inf)
     */

    CDoubleConstSP cap; // $optional(Level above which 'field' is invalid)

    /**
     * Whether 'cap' is included in the field's valid range
     */

    CBoolConstSP capInclusive; // $optional(Is field = 'cap' valid?)

    /**
     * Whether tweaks should be clipped to the field's valid range
     *
     * If false: tweaks that go outside the range will fail.  If true: they
     * will be clipped to the endpoint if that's included in the range,
     * else reset to halfway between the endpoint and the initial value.
     */

    bool clip; // $required(Whether to clip out-of-range tweaks or just fail)

    // Unpacked object representation, set in validatePop2Object()

    TweakFunctionConstSP _operator;       // $transient
    double _floor, _cap;                  // $transient
    bool _floorInclusive, _capInclusive;  // $transient
    bool _useCalibratorRange;             // $transient

protected:

    void validatePop2Object();

    IFieldTweak::IOperatorConstSP operat0r() const {
        return IFieldTweak::IOperator::numeric(
            _operator,
            Range(_floor, _floorInclusive, _cap, _capInclusive),
            _useCalibratorRange,
            clip);
    }
};

void NumericFieldRiskPropertyDefinition::validatePop2Object() {
    BaseFieldRiskPropertyDefinition::validatePop2Object();
    try {
        _floor = !floor ? -HUGE_VAL : floor->doubleValue();
        _floorInclusive = !floorInclusive ? true : floorInclusive->boolValue();
        _cap = !cap ? HUGE_VAL : cap->doubleValue();
        _capInclusive = !capInclusive ? true : capInclusive->boolValue();

        if (!(_floor < _cap)) {
            throw ModelException("Range defined by cap and floor is degenerate");
        }

        _useCalibratorRange =
            !floor && !floorInclusive && !cap && !capInclusive;

        if (operator_ == "=") {
            _operator = TweakFunction::setter();
        }
        else if (operator_ == "+") {
            _operator = TweakFunction::additive();
        }
        else if (operator_ == "*") {
            _operator = TweakFunction::multiplicative();
        }
        else if (operator_ == "e" || operator_ == "E") {
            _operator = TweakFunction::exponential();
        }
        else if (operator_ == "AUTO") {
            _operator = TweakFunction::adaptive();
        }
        else {
            throw ModelException(
                "Operator code \"" + operator_ + "\" is invalid: "
                "should be \"=\", \"+\", \"*\" or \"AUTO\"");
        }
    }
    catch (ModelException& e) {
        throw ModelException(e, __FUNCTION__,
            "Validating NumericFieldRiskPropertyDefinition \"" + name + "\"");
    }
}












void BaseFieldRiskPropertyDefinition::load(CClassSP& clazz) {
  REGISTER(BaseFieldRiskPropertyDefinition, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldRiskPropertyDefinition);
  FIELD(name, "Unique name for the property");
  FIELD(implementation, "FLEXIBLE, INTERNAL or DEFAULT");
  FIELD(type, "Name of type to be tweaked");
  FIELD(field, "Path to field to tweak, starting from 'type'");
  FIELD(subField, "Path to field to tweak, starting from 'type'");
  FIELD_MAKE_OPTIONAL(subField);
  FIELD(qualifierField, "Path to field in which expiry labels are stored, starting from 'type'");
  FIELD_MAKE_OPTIONAL(qualifierField);
  FIELD(qualifierSubField, "Path to field in which expiry labels are stored, starting from 'type'");
  FIELD_MAKE_OPTIONAL(qualifierSubField);
  FIELD(qualifier2Field, "qualifier2Field");
  FIELD_MAKE_OPTIONAL(qualifier2Field);
  FIELD(qualifier2SubField, "qualifier2SubField");
  FIELD_MAKE_OPTIONAL(qualifier2SubField);
  FIELD(qualifierDimension, "For matrix-valued fields: which dimension corresponds to 'qualifierField'");
  FIELD_MAKE_OPTIONAL(qualifierDimension);
  FIELD(measure, "Terms for reporting deriv: NOMINAL or ABSOLUTE");
  FIELD(beforeGetMarket, "Whether to apply scenarios before the getMarket phase");
  FIELD_MAKE_OPTIONAL(beforeGetMarket);
  FIELD(_qualifiersAreIndices, "_qualifiersAreIndices");
  FIELD_MAKE_TRANSIENT(_qualifiersAreIndices);
  FIELD(_useOverride, "_useOverride");
  FIELD_MAKE_TRANSIENT(_useOverride);
  FIELD(_requireOverride, "_requireOverride");
  FIELD_MAKE_TRANSIENT(_requireOverride);
  FIELD(_absoluteDistance, "_absoluteDistance");
  FIELD_MAKE_TRANSIENT(_absoluteDistance);
}

CClassConstSP const BaseFieldRiskPropertyDefinition::TYPE = CClass::registerClassLoadMethod(
  "BaseFieldRiskPropertyDefinition", typeid(BaseFieldRiskPropertyDefinition), BaseFieldRiskPropertyDefinition::load);

IObject* FieldRiskPropertyDefinition::emptyShell() {
  return new FieldRiskPropertyDefinition();
}

void FieldRiskPropertyDefinition::load(CClassSP& clazz) {
  clazz->setPublic();
  REGISTER(FieldRiskPropertyDefinition, clazz);
  SUPERCLASS(BaseFieldRiskPropertyDefinition);
  FIELD(xml1, "xml1");
  FIELD_MAKE_OPTIONAL(xml1);
  FIELD(xml2, "xml2");
  FIELD_MAKE_OPTIONAL(xml2);
  FIELD(xml3, "xml3");
  FIELD_MAKE_OPTIONAL(xml3);
  FIELD(xml4, "xml4");
  FIELD_MAKE_OPTIONAL(xml4);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const FieldRiskPropertyDefinition::TYPE = CClass::registerClassLoadMethod(
  "FieldRiskPropertyDefinition", typeid(FieldRiskPropertyDefinition), FieldRiskPropertyDefinition::load);

IObject* NumericFieldRiskPropertyDefinition::emptyShell() {
  return new NumericFieldRiskPropertyDefinition();
}

void NumericFieldRiskPropertyDefinition::load(CClassSP& clazz) {
  clazz->setPublic();
  REGISTER(NumericFieldRiskPropertyDefinition, clazz);
  SUPERCLASS(BaseFieldRiskPropertyDefinition);
  FIELD_USING_ALIAS(operator_, operator,
                           "How to tweak: 'AUTO', '+', '*', '='");
  FIELD(floor, "Level below which 'field' is invalid");
  FIELD_MAKE_OPTIONAL(floor);
  FIELD(floorInclusive, "Is field = 'floor' valid?");
  FIELD_MAKE_OPTIONAL(floorInclusive);
  FIELD(cap, "Level above which 'field' is invalid");
  FIELD_MAKE_OPTIONAL(cap);
  FIELD(capInclusive, "Is field = 'cap' valid?");
  FIELD_MAKE_OPTIONAL(capInclusive);
  FIELD(clip, "Whether to clip out-of-range tweaks or just fail");
  FIELD(_operator, "_operator");
  FIELD_MAKE_TRANSIENT(_operator);
  FIELD(_floor, "_floor");
  FIELD_MAKE_TRANSIENT(_floor);
  FIELD(_cap, "_cap");
  FIELD_MAKE_TRANSIENT(_cap);
  FIELD(_floorInclusive, "_floorInclusive");
  FIELD_MAKE_TRANSIENT(_floorInclusive);
  FIELD(_capInclusive, "_capInclusive");
  FIELD_MAKE_TRANSIENT(_capInclusive);
  FIELD(_useCalibratorRange, "_useCalibratorRange");
  FIELD_MAKE_TRANSIENT(_useCalibratorRange);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const NumericFieldRiskPropertyDefinition::TYPE = CClass::registerClassLoadMethod(
  "NumericFieldRiskPropertyDefinition", typeid(NumericFieldRiskPropertyDefinition), NumericFieldRiskPropertyDefinition::load);

void IFieldRiskPropertyDefinition::load(CClassSP& clazz) {
  REGISTER_INTERFACE(IFieldRiskPropertyDefinition, clazz);
  clazz->setPublic();
  EXTENDS(IObject);
}

CClassConstSP const IFieldRiskPropertyDefinition::TYPE = CClass::registerInterfaceLoadMethod(
  "IFieldRiskPropertyDefinition", typeid(IFieldRiskPropertyDefinition), IFieldRiskPropertyDefinition::load);

DRLIB_END_NAMESPACE
