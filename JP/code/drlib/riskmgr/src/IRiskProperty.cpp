/**
 * @file IRiskProperty.cpp
 */

#include "edginc/config.hpp"
#define QLIB_IRISKPROPERTY_CPP
#include "edginc/TRACE.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

// 
// *************************
//  ConditionedRiskProperty
// *************************
// 

template <class QUALIFIER>
class ConditionedRiskProperty:
        public CObject,
        public virtual IRiskProperty<QUALIFIER> {

    ConditionedRiskProperty(): CObject(TYPE) {}

    static IObject* defaultOne() {
        return new ConditionedRiskProperty();
    }

    static void load(CClassSP& clazz) {
        REGISTER(ConditionedRiskProperty, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRiskProperty<QUALIFIER>);
        EMPTY_SHELL_METHOD(defaultOne);
        FIELD(condition, "condition");
        FIELD(property, "property");
    }

public:

    static CClassConstSP const TYPE;

    DECLARE(QUALIFIER)

    typedef IRiskProperty<QUALIFIER> Property;
    DECLARE(Property)

private:

    IHypothesisConstSP condition;
    PropertyConstSP property;

public:

    ConditionedRiskProperty(IHypothesisConstSP condition,
                            PropertyConstSP property):
        CObject(TYPE),
        condition(condition),
        property(property)
    {}

    /**
     * See RiskProperty<PROPERTY>::discrete()
     */

    bool discrete() const {
        return property->discrete();
    }

    /**
     * See RiskProperty<PROPERTY>::subjectInterface()
     */

    CClassConstSP subjectInterface() const {
        return property->subjectInterface();
    }

    /**
     * See RiskProperty<PROPERTY>::subjectNames()
     */

    OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        return property->subjectNames(world);
    }

    /**
     * See RiskProperty<PROPERTY>::axisFor()
     */

    IRiskAxisConstSP axisFor(OutputNameConstSP subjectName,
                             QUALIFIERConstSP qualifier) const {
        return IRiskAxis::conditioned(
            condition, property->axisFor(subjectName, qualifier));
    }

    /**
     * See RiskProperty<PROPERTY>::subjectQualifiers()
     */

    QUALIFIERArrayConstSP subjectQualifiers(IObjectConstSP world,
                                            OutputNameConstSP name) const {
        return property->subjectQualifiers(world, name);
    }

    /**
     * For error messages
     */

    string toString() const {
        return property->toString() + " assuming " + condition->toString();
    }
};

template class ConditionedRiskProperty<Void>;
template class ConditionedRiskProperty<ExpiryWindow>;
template class ConditionedRiskProperty<ExpiryPair>;
template class ConditionedRiskProperty<ExpiryAndStrike>;
template class ConditionedRiskProperty<BoxedInt>;

template <>
CClassConstSP const ConditionedRiskProperty<Void>::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskProperty<Void>", typeid(ConditionedRiskProperty<Void>), load);

template <>
CClassConstSP const ConditionedRiskProperty<ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskProperty<ExpiryWindow>", typeid(ConditionedRiskProperty<ExpiryWindow>), load);

template <>
CClassConstSP const ConditionedRiskProperty<ExpiryPair>::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskProperty<ExpiryPair>", typeid(ConditionedRiskProperty<ExpiryPair>), load);

template <>
CClassConstSP const ConditionedRiskProperty<ExpiryAndStrike>::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskProperty<ExpiryAndStrike>", typeid(ConditionedRiskProperty<ExpiryAndStrike>), load);

template <>
CClassConstSP const ConditionedRiskProperty<BoxedInt>::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskProperty<BoxedInt>", typeid(ConditionedRiskProperty<BoxedInt>), load);

// 
// ***************
//  IRiskProperty
// ***************
// 

template <class QUALIFIER>
IRiskProperty<QUALIFIER>::IRiskProperty() {}

template <class QUALIFIER>
IRiskProperty<QUALIFIER>::~IRiskProperty() {}

template <class QUALIFIER>
void IRiskProperty<QUALIFIER>::tweakAllSubjects(
        IObjectSP world, smartConstPtr<QUALIFIER> qualifier,
        double coefficient) const {
    try {
        OutputNameArrayConstSP names(subjectNames(world));
        for (int n = 0; n < names->size(); ++n) {
            axisFor((*names)[n], qualifier)->hypothesis(coefficient)->
                applyTo(world);
        }
    }
    catch (ModelException& e) {
        throw ModelException(e, "IRiskProperty::tweakAllSubjects");
    }
}

template <class QUALIFIER>
DoubleArrayConstSP IRiskProperty<QUALIFIER>::adaptiveCoefficients(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        QualifierArrayConstSP qualifiers,
        double targetCoefficient) const {
    return DoubleArray::SP(qualifiers->size(), targetCoefficient);
}

template <class QUALIFIER>
BoolArrayConstSP IRiskProperty<QUALIFIER>::mayHaveEffect(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        QualifierArrayConstSP qualifiers,
        const Sensitivity* sensitivity) const {
    return CBoolArraySP(new BoolArray(qualifiers->size(), true));
}

template <>
BoolArrayConstSP IRiskProperty<ExpiryWindow>::mayHaveEffect(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        ExpiryWindowArrayConstSP exps,
        const Sensitivity* sensitivity) const {
    TRACE_METHOD;

    DateTime valueDate = world->getInstruments()->getValueDate();
    DateTime endDate = world->getInstruments()->endDate(sensitivity);
    TRACE("Last sensitive date is " << endDate.toString());

    CBoolArraySP them = CBoolArraySP(new BoolArray(exps->size()));

    bool expired = false;

    for (int e = 0; e < exps->size(); ++e) {
        (*them)[e] = !expired;

        expired = (*exps)[e]->expiry->toDate(valueDate).isGreater(endDate);

        if (expired && e > 0 && !(*them)[e-1]) {
            TRACE("Last sens date passed at " <<
                  *(*exps)[e]->expiry << " mark; "
                  "subsequent expiries will be marked as insensitive");
        }
    }

    return them;
}

template <>
BoolArrayConstSP IRiskProperty<ExpiryPair>::mayHaveEffect(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        ExpiryPairArrayConstSP exps,
        const Sensitivity* sensitivity) const {
    TRACE_METHOD;

    DateTime valueDate = world->getInstruments()->getValueDate();
    DateTime endDate = world->getInstruments()->endDate(sensitivity);
    TRACE("Last sensitive date is " << endDate);

    CBoolArraySP them = CBoolArraySP(new BoolArray(exps->size()));

    bool expired = false;

    for (int e = 0; e < exps->size(); ++e) {
        (*them)[e] = !expired;

        expired = (*exps)[e]->optExpiry->toDate(valueDate).isGreater(endDate);

        if (expired && e > 0 && !(*them)[e-1]) {
            TRACE("Last sens date passed at " <<
                  *(*exps)[e]->optExpiry << " mark; "
                  "subsequent expiries will be marked as insensitive");
        }
    }

    return them;
}

template <>
BoolArrayConstSP IRiskProperty<ExpiryAndStrike>::mayHaveEffect(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        ExpiryAndStrikeArrayConstSP exps,
        const Sensitivity* sensitivity) const {
    TRACE_METHOD;

    DateTime valueDate = world->getInstruments()->getValueDate();
    DateTime endDate = world->getInstruments()->endDate(sensitivity);
    TRACE("Last sensitive date is " << endDate.toString());

    CBoolArraySP them = CBoolArraySP(new BoolArray(exps->size()));

    bool expired = false;

    for (int e = 0; e < exps->size(); ++e) {
        (*them)[e] = !expired;

        expired = (*exps)[e]->expiry->toDate(valueDate).isGreater(endDate);

        if (expired && e > 0 && !(*them)[e-1]) {
            TRACE("Last sens date passed at " <<
                  *(*exps)[e]->expiry << " mark; "
                  "subsequent expiries will be marked as insensitive");
        }
    }

    return them;
}

template <class QUALIFIER>
smartPtr<IRiskProperty<QUALIFIER> > IRiskProperty<QUALIFIER>::conditioned(
        IHypothesisConstSP condition,
        smartConstPtr<IRiskProperty<QUALIFIER> > property) {
    return smartPtr<IRiskProperty<QUALIFIER> >(
        new ConditionedRiskProperty<QUALIFIER>(condition, property));
}



template <class QUALIFIER>
void IRiskProperty<QUALIFIER>::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IRiskProperty<QUALIFIER>, clazz);
    EXTENDS(IAbstractRiskProperty);
}

template <>
CClassConstSP const IRiskProperty<Void>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IScalarRiskProperty", typeid(IRiskProperty<Void>), load);

template <>
CClassConstSP const IRiskProperty<ExpiryWindow>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IExpiryRiskProperty",
        typeid(IExpiryRiskProperty), load);

template <>
CClassConstSP const IRiskProperty<ExpiryPair>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IExpiryPairRiskProperty",
        typeid(IExpiryPairRiskProperty), load);

template <>
CClassConstSP const IRiskProperty<ExpiryAndStrike>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IExpiryAndStrikeRiskProperty",
        typeid(IExpiryAndStrikeRiskProperty), load);

template <>
CClassConstSP const IRiskProperty<BoxedInt>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IIntRiskProperty",
        typeid(IIntRiskProperty), load);

// Must use normal template instantiation since we want these instantiated
// for both debug and optimised
template class RISKMGR_DLL IRiskProperty<Void>;
template class RISKMGR_DLL IRiskProperty<ExpiryWindow>;
template class RISKMGR_DLL IRiskProperty<ExpiryPair>;
template class RISKMGR_DLL IRiskProperty<ExpiryAndStrike>;
template class RISKMGR_DLL IRiskProperty<BoxedInt>;

DEFINE_TEMPLATE_TYPE(IScalarRiskPropertyArray);
DEFINE_TEMPLATE_TYPE(IExpiryRiskPropertyArray);
DEFINE_TEMPLATE_TYPE(IExpiryPairRiskPropertyArray);
DEFINE_TEMPLATE_TYPE(IExpiryAndStrikeRiskPropertyArray);
DEFINE_TEMPLATE_TYPE(IIntRiskPropertyArray);

DRLIB_END_NAMESPACE
