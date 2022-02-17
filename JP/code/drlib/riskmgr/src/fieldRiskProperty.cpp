/**
 * @file fieldRiskProperty.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Format.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/SimpleTweakNameListID.hpp"
#include "edginc/SimpleTweakNameResolver.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IFieldTweak.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/FieldRiskAxis.hpp"

DRLIB_BEGIN_NAMESPACE

// 
// *******************
//  FieldRiskProperty
// *******************
// 

struct FieldRiskProperty: CObject,
                          virtual IAbstractRiskProperty {

    static CClassConstSP const TYPE;

private:

    static void load(CClassSP& clazz);

    const string className;      // $required
    mutable CClassConstSP clazz; // $unregistered

protected:

    const bool absoluteDistance; // $required

    /**
     * Constructor.
     *
     * @param type    [Type of implementing subclass, as usual]
     *
     * @param clazz   The classs to which the property applies
     */

    FieldRiskProperty(CClassConstSP type,
                      CClassConstSP clazz,
                      bool absoluteDistance):
        CObject(type),
        className(clazz->getName()),
        clazz(clazz),
        absoluteDistance(absoluteDistance)
    {}

    FieldRiskProperty(CClassConstSP type):
        CObject(type),
        clazz(0),
        absoluteDistance(false)
    {}

    void validatePop2Object() {
        try {
            const_cast<CClassConstSP&>(clazz) = CClass::forName(className);
        }
        catch (ModelException& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

public:

    ~FieldRiskProperty() {}

    virtual bool discrete() const {
        return false;
    }

    virtual CClassConstSP subjectInterface() const {
        try {
            if (!clazz) {
                clazz = CClass::forName(className);
            }

            return clazz;
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

    virtual OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        try {
            SimpleTweakNameListID namesListID(subjectInterface());
            return OutputName::trim(SensMgrConst(world).allNames(&namesListID));
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }
};

// 
// ********
//  scalar
// ********
// 

struct ScalarFieldRiskProperty: FieldRiskProperty,
                                virtual IRiskProperty<Void> {

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    IFieldTweakConstSP tweak;             // $required

    ScalarFieldRiskProperty(CClassConstSP clazz,
                            IFieldTweakConstSP tweak,
                            bool absoluteDistance):
        FieldRiskProperty(TYPE, clazz, absoluteDistance),
        tweak(tweak)
    {}

    ScalarFieldRiskProperty():
        FieldRiskProperty(TYPE)
    {}

    VoidArrayConstSP subjectQualifiers(IObjectConstSP world,
                                       OutputNameConstSP name) const {
        return VoidArrayConstSP();
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP name, VoidConstSP) const {
        return IRiskAxisConstSP(new FieldRiskAxis(
            IAbstractRiskPropertyConstSP::attachToRef(this),
            name, tweak, absoluteDistance));
    }

    string toString() const {
        return subjectInterface()->getName() + " " + tweak->toString();
    }
};

IScalarRiskPropertyConstSP fieldRiskProperty::scalar(CClassConstSP clazz,
                                                     IFieldTweakConstSP tweak,
                                                     bool absoluteDistance) {
    return IScalarRiskPropertyConstSP(new ScalarFieldRiskProperty(
        clazz, tweak, absoluteDistance));
}

// 
// ************************
//  ArrayFieldRiskProperty
// ************************
// 

struct ArrayFieldRiskProperty: FieldRiskProperty {
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    FieldPathConstSP field;                  // $required
    FieldPathConstSP subField;               // $required
    const int indexDimension;                // $required
    const IFieldTweak::IOperatorConstSP op;  // $required
    const IObjectConstSP arg;                // $optional

    void validate() {
        if (!field || field->path()->empty()) {
            throw ModelException(__FUNCTION__,
                                 "'field' must be a nonempty array");
        }
    }
    
    /**
     * @param field             May not be NULL or empty
     *
     * @param subField          May be NULL to mean "empty"
     *
     * @param arg               May be NULL to denote NULL as a value
     */

    ArrayFieldRiskProperty(CClassConstSP type,
                           CClassConstSP clazz,
                           FieldPathConstSP field,
                           FieldPathConstSP subField,
                           int indexDimension,
                           IFieldTweak::IOperatorConstSP op,
                           IObjectConstSP arg,
                           bool absoluteDistance):
        FieldRiskProperty(type, clazz, absoluteDistance),
        field(field),
        subField(!subField ? FieldPath::SP() : subField),
        indexDimension(indexDimension),
        op(op),
        arg(arg)
    {
        validate();
    }

    void validatePop2Object() {
        validate();
    }

    ArrayFieldRiskProperty(CClassConstSP type):
        FieldRiskProperty(type),
        indexDimension(0)
    {}

    string toString() const {
        return subjectInterface()->getName() + field->toString() +
               subField->toString() + " (" + op->toString() +
               " " + IObject::stringOf(arg) + ")";
    }
};

// 
// *************
//  elementwise
// *************
// 

struct ElementwiseFieldRiskProperty: ArrayFieldRiskProperty,
                                     virtual IRiskProperty<BoxedInt> {

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    /**
     * @param field             May not be NULL or empty
     *
     * @param subField          May be NULL to mean "empty"
     *
     * @param arg               May be NULL to denote NULL as a value
     */

    ElementwiseFieldRiskProperty(CClassConstSP clazz,
                                 FieldPathConstSP field,
                                 FieldPathConstSP subField,
                                 int indexDimension,
                                 IFieldTweak::IOperatorConstSP op,
                                 IObjectConstSP arg,
                                 bool absoluteDistance):
        ArrayFieldRiskProperty(TYPE, clazz, field, subField,
                               indexDimension, op, arg, absoluteDistance)
    {}

    ElementwiseFieldRiskProperty():
        ArrayFieldRiskProperty(TYPE)
    {}

    BoxedIntArrayConstSP subjectQualifiers(IObjectConstSP world,
                                           OutputNameConstSP name) const {
        try {
            SimpleTweakNameResolver nameres(name);
            IObjectConstSP x = SensMgrConst(world).theFirst(subjectInterface(),
                                                            &nameres);
            try {
                IObjectConstSP v = field->get(x);
                if (!v) throw ModelException("Field's value is NULL");

                const CDoubleMatrix* m =
                    dynamic_cast<const CDoubleMatrix*>(v.get());
                int n = m ? indexDimension == 0 ? m->numCols() : m->numRows() :
                            IArrayConstSP::dynamicCast(v)->getLength();
                BoxedIntArraySP is(new BoxedIntArray(n));
                for (int i = 0; i < n; ++i) {
                    (*is)[i].reset(BoxedInt::create(i));
                }
                return is;
            }
            catch (exception& e) {
                throw ModelException(e, "Examining " + field->contextMessage(x));
            }
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                                 "Finding length of array to be tweaked");
        }
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP name,
                             BoxedIntConstSP index) const {
        return IRiskAxisConstSP(new FieldRiskAxis(
            IAbstractRiskPropertyConstSP::attachToRef(this),
            name,
            IFieldTweak::elementwise(
                field, subField,
                IFieldTweak::IIndices::direct(
                    IntArray::SP(1, index->intValue())),
                indexDimension,
                ObjectArray::SP(1, IObjectSP::constCast(arg)),
                op),
            absoluteDistance));
    }
};

IIntRiskPropertyConstSP fieldRiskProperty::elementwise(
        CClassConstSP clazz,
        FieldPathConstSP field,
        FieldPathConstSP subField,
        int indexDimension,
        IFieldTweak::IOperatorConstSP op,
        IObjectConstSP arg,
        bool absoluteDistance) {
    return IIntRiskPropertyConstSP(new ElementwiseFieldRiskProperty(
        clazz, field, subField, indexDimension, op, arg, absoluteDistance));
}

// 
// *************************
//  LookupFieldRiskProperty
// *************************
// 

struct AbstractLookupFieldRiskProperty: ArrayFieldRiskProperty {

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    const FieldPathConstSP keyListField;  // $optional
    const FieldPathConstSP keyField;      // $required

    /**
     * @param field             May not be NULL or empty
     *
     * @param subField          May be NULL to mean "empty"
     *
     * @param keyListField
     *
     *     May be NULL to mean "use field registered via
     *     Calibrator::IAdjustable::registerBootstrappableField
     *
     * @param keyField          May be NULL to mean "empty"
     *
     * @param arg               May be NULL to denote NULL as a value
     */

    AbstractLookupFieldRiskProperty(CClassConstSP type,
                                    CClassConstSP clazz,
                                    FieldPathConstSP field,
                                    FieldPathConstSP subField,
                                    FieldPathConstSP keyListField,
                                    FieldPathConstSP keyField,
                                    int indexDimension,
                                    IFieldTweak::IOperatorConstSP op,
                                    IObjectConstSP arg,
                                    bool absoluteDistance):
        ArrayFieldRiskProperty(type, clazz, field, subField, indexDimension,
                               op, arg, absoluteDistance),
        keyListField(keyListField),
        keyField(!keyField ? FieldPath::SP() : keyField)
    {}

    AbstractLookupFieldRiskProperty(CClassConstSP type):
        ArrayFieldRiskProperty(type)
    {}

    IArrayConstSP qualifierObjects(IObjectConstSP root,
                                   CClassConstSP qualifierType) const {
        try {
            if (!keyListField) {
                try {
                    FieldPath::ObjectChain chain =
                        field->objectChain(IObjectSP::constCast(root));
                    return Calibrator::IAdjustable::
                        getGetExpiriesMethod(chain.fields.back())(
                            (*chain.objects)[chain.fields.size() - 1].get());
                }
                catch (exception& e) {
                    throw ModelException(e,
                        "Looking up Calibrator-registered expiries",
                        "Did you forget to specify an explicit expiryQualifier "
                        "in your FieldRiskPropertyDefinition?");
                }
            }

            IArrayConstSP keyList;
            try {
                keyList = IArrayConstSP::dynamicCast(keyListField->get(root));
            }
            catch (exception& e) {
                throw ModelException(e, "Field is not an array");
            }

            IArraySP keys(
                qualifierType->newArrayInstanceByComponent(keyList->getLength()));

            for (int k = 0; k < keyList->getLength(); ++k) {
                IObjectConstSP keyObj = keyList->get(k);
                try {
                    IObjectConstSP key = keyField->get(keyObj);
                    if (!qualifierType->isInstance(key)) {
                        throw ModelException(
                            "Key " + (!key ? string("NULL") : key->toString()) +
                            " is not of type " +
                            qualifierType->getName());
                    }

                    keys->set(k, IObjectSP::constCast(key));
                }
                catch (exception& e) {
                    throw ModelException(e,
                        "Index #" + Format::toString(k) + " in array (" +
                        (!keyObj ? string("NULL") : keyObj->toString()) + ")");
                }
            }

            return keys;
        }        
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                                 "Looking up indices in " +
								 !keyListField ? "calibrator-registered field" :
			                                     keyListField->contextMessage(root, -1));
        }
    }

    IRiskAxisConstSP axisForObject(OutputNameConstSP name,
                                   IObjectConstSP qualifier) const {
        IObjectSP q = IObjectSP::constCast(qualifier);

        return IRiskAxisConstSP(new FieldRiskAxis(
            IAbstractRiskPropertyConstSP::attachToRef(this),
            name,
            IFieldTweak::elementwise(
                field, subField,
                IFieldTweak::IIndices::lookupOrCalibratorExpiries(
                    field, keyListField, keyField, ObjectArray::SP(1, q)),
                indexDimension,
                ObjectArray::SP(1, IObjectSP::constCast(arg)),
                op),
            absoluteDistance));
    }
};

template <class QUALIFIER>
struct LookupFieldRiskProperty: AbstractLookupFieldRiskProperty,
                                virtual IRiskProperty<QUALIFIER> {

    typedef QUALIFIER Qualifier;
    DECLARE(Qualifier)

    static CClassConstSP const TYPE;
    static IObject* emptyShell() {
        return new LookupFieldRiskProperty();
    }
    static void load(CClassSP& clazz) {
        REGISTER(LookupFieldRiskProperty, clazz);
        SUPERCLASS(AbstractLookupFieldRiskProperty);
        IMPLEMENTS(IRiskProperty<QUALIFIER>);
        EMPTY_SHELL_METHOD(emptyShell);
    }

    /**
     * @param field             May not be NULL or empty
     *
     * @param subField          May be NULL to mean "empty"
     *
     * @param keyListField
     *
     *     May be NULL to mean "use field registered via
     *     Calibrator::IAdjustable::registerBootstrappableField
     *
     * @param keyField          May be NULL to mean "empty"
     *
     * @param arg               May be NULL to denote NULL as a value
     */

    LookupFieldRiskProperty(CClassConstSP type,
                            CClassConstSP clazz,
                            FieldPathConstSP field,
                            FieldPathConstSP subField,
                            FieldPathConstSP keyListField,
                            FieldPathConstSP keyField,
                            int indexDimension,
                            IFieldTweak::IOperatorConstSP op,
                            IObjectConstSP arg,
                            bool absoluteDistance):
        AbstractLookupFieldRiskProperty(type, clazz, field, subField,
                                        keyListField, keyField, indexDimension,
                                        op, arg, absoluteDistance)
    {}

    LookupFieldRiskProperty(CClassConstSP clazz,
                            FieldPathConstSP field,
                            FieldPathConstSP subField,
                            FieldPathConstSP keyListField,
                            FieldPathConstSP keyField,
                            int indexDimension,
                            IFieldTweak::IOperatorConstSP op,
                            IObjectConstSP arg,
                            bool absoluteDistance):
        AbstractLookupFieldRiskProperty(TYPE, clazz, field, subField,
                                        keyListField, keyField, indexDimension,
                                        op, arg, absoluteDistance)
    {}

    LookupFieldRiskProperty(CClassConstSP type = TYPE):
        AbstractLookupFieldRiskProperty(type)
    {}

    virtual QualifierArrayConstSP qualifiers(IObjectConstSP root) const {
        return QualifierArrayConstSP::dynamicCast(
            qualifierObjects(root, Qualifier::TYPE));
    }

    QualifierArrayConstSP subjectQualifiers(IObjectConstSP world,
                                            OutputNameConstSP name) const {
        try {
            SimpleTweakNameResolver nameres(name);
            return qualifiers(SensMgrConst(world).theFirst(
                                  subjectInterface(), &nameres));
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP name,
                             QualifierConstSP qualifier) const {
        return axisForObject(name, qualifier);
    }
};

// 
// ***********
//  pointwise
// ***********
// 

template struct LookupFieldRiskProperty<ExpiryWindow>;

struct PointwiseFieldRiskProperty: LookupFieldRiskProperty<ExpiryWindow> {

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    /**
     * @param field             May not be NULL or empty
     *
     * @param subField          May be NULL to mean "empty"
     *
     * @param keyListField
     *
     *     May be NULL to mean "use field registered via
     *     Calibrator::IAdjustable::registerBootstrappableField
     *
     * @param keyField          May be NULL to mean "empty"
     *
     * @param arg               May be NULL to denote NULL as a value
     */

    PointwiseFieldRiskProperty(CClassConstSP clazz,
                               FieldPathConstSP field,
                               FieldPathConstSP subField,
                               FieldPathConstSP keyListField,
                               FieldPathConstSP keyField,
                               int indexDimension,
                               IFieldTweak::IOperatorConstSP op,
                               IObjectConstSP arg,
                               bool absoluteDistance):
        LookupFieldRiskProperty<ExpiryWindow>(
            TYPE, clazz, field, subField,
            keyListField, keyField, indexDimension,
            op, arg, absoluteDistance)
    {}

    PointwiseFieldRiskProperty():
        LookupFieldRiskProperty<ExpiryWindow>(TYPE)
    {}

    ExpiryWindowArrayConstSP qualifiers(IObjectConstSP root) const {
        TRACE_METHOD;
        return ExpiryWindow::series(ExpiryArrayConstSP::dynamicCast(
                   qualifierObjects(root, Expiry::TYPE)));
    }
};

IExpiryRiskPropertyConstSP fieldRiskProperty::pointwise(
        CClassConstSP clazz,
        FieldPathConstSP field,
        FieldPathConstSP subField,
        FieldPathConstSP keyListField,
        FieldPathConstSP keyField,
        int indexDimension,
        IFieldTweak::IOperatorConstSP op,
        IObjectConstSP arg,
        bool absoluteDistance) {
    return IExpiryRiskPropertyConstSP(new PointwiseFieldRiskProperty(
        clazz, field, subField, keyListField, keyField,
        indexDimension, op, arg, absoluteDistance));
}

// 
// ****************
//  expiryPairwise
// ****************
// 

template struct LookupFieldRiskProperty<ExpiryPair>;

template <>
CClassConstSP const LookupFieldRiskProperty<ExpiryPair>::TYPE =
    CClass::registerClassLoadMethod(
        "LookupFieldRiskProperty<ExpiryPair>",
        typeid(LookupFieldRiskProperty<ExpiryPair>), load);

IExpiryPairRiskPropertyConstSP fieldRiskProperty::expiryPairwise(
        CClassConstSP clazz,
        FieldPathConstSP field,
        FieldPathConstSP subField,
        FieldPathConstSP keyListField,
        FieldPathConstSP keyField,
        int indexDimension,
        IFieldTweak::IOperatorConstSP op,
        IObjectConstSP arg,
        bool absoluteDistance) {
    return IExpiryPairRiskPropertyConstSP(new LookupFieldRiskProperty<ExpiryPair>(
        clazz, field, subField, keyListField, keyField,
        indexDimension, op, arg, absoluteDistance));
}

// 
// *********************
//  expiryAndStrikewise
// *********************
// 

template struct LookupFieldRiskProperty<ExpiryAndStrike>;

template <>
CClassConstSP const LookupFieldRiskProperty<ExpiryAndStrike>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "LookupFieldRiskProperty<ExpiryAndStrike>",
        typeid(LookupFieldRiskProperty<ExpiryAndStrike>), load);

IExpiryAndStrikeRiskPropertyConstSP fieldRiskProperty::expiryAndStrikewise(
        CClassConstSP clazz,
        FieldPathConstSP field,
        FieldPathConstSP subField,
        FieldPathConstSP keyListField,
        FieldPathConstSP keyField,
        int indexDimension,
        IFieldTweak::IOperatorConstSP op,
        IObjectConstSP arg,
        bool absoluteDistance) {
    return IExpiryAndStrikeRiskPropertyConstSP(new LookupFieldRiskProperty<ExpiryAndStrike>(
        clazz, field, subField, keyListField, keyField,
        indexDimension, op, arg, absoluteDistance));
}



void AbstractLookupFieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(AbstractLookupFieldRiskProperty, clazz);
  SUPERCLASS(ArrayFieldRiskProperty);
  FIELD(keyListField, "keyListField");
  FIELD_MAKE_OPTIONAL(keyListField);
  FIELD(keyField, "keyField");
}

CClassConstSP const AbstractLookupFieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "AbstractLookupFieldRiskProperty", typeid(AbstractLookupFieldRiskProperty), AbstractLookupFieldRiskProperty::load);

IObject* ScalarFieldRiskProperty::emptyShell() {
  return new ScalarFieldRiskProperty();
}

void ScalarFieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(ScalarFieldRiskProperty, clazz);
  SUPERCLASS(FieldRiskProperty);
  IMPLEMENTS(IRiskProperty<Void>);
  FIELD(tweak, "tweak");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ScalarFieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "ScalarFieldRiskProperty", typeid(ScalarFieldRiskProperty), ScalarFieldRiskProperty::load);

IObject* ElementwiseFieldRiskProperty::emptyShell() {
  return new ElementwiseFieldRiskProperty();
}

void ElementwiseFieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(ElementwiseFieldRiskProperty, clazz);
  SUPERCLASS(ArrayFieldRiskProperty);
  IMPLEMENTS(IRiskProperty<BoxedInt>);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ElementwiseFieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "ElementwiseFieldRiskProperty", typeid(ElementwiseFieldRiskProperty), ElementwiseFieldRiskProperty::load);

IObject* PointwiseFieldRiskProperty::emptyShell() {
  return new PointwiseFieldRiskProperty();
}

void PointwiseFieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(PointwiseFieldRiskProperty, clazz);
  SUPERCLASS(LookupFieldRiskProperty<ExpiryWindow>);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const PointwiseFieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "PointwiseFieldRiskProperty", typeid(PointwiseFieldRiskProperty), PointwiseFieldRiskProperty::load);

void FieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(FieldRiskProperty, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IAbstractRiskProperty);
  FIELD(className, "className");
  FIELD(absoluteDistance, "absoluteDistance");
}

CClassConstSP const FieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "FieldRiskProperty", typeid(FieldRiskProperty), FieldRiskProperty::load);

template <>
CClassConstSP const LookupFieldRiskProperty<ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
  "LookupFieldRiskProperty<ExpiryWindow>", typeid(LookupFieldRiskProperty<ExpiryWindow>), LookupFieldRiskProperty<ExpiryWindow>::load);

void ArrayFieldRiskProperty::load(CClassSP& clazz) {
  REGISTER(ArrayFieldRiskProperty, clazz);
  SUPERCLASS(FieldRiskProperty);
  FIELD(field, "field");
  FIELD(subField, "subField");
  FIELD(indexDimension, "indexDimension");
  FIELD(op, "op");
  FIELD(arg, "arg");
  FIELD_MAKE_OPTIONAL(arg);
}

CClassConstSP const ArrayFieldRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "ArrayFieldRiskProperty", typeid(ArrayFieldRiskProperty), ArrayFieldRiskProperty::load);

DRLIB_END_NAMESPACE
