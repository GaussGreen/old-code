//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiExpiryShift.cpp
//
//   Description : Specialised Perturbation where shifts for a range of 
//                 expiries are defined together in buckets. The interpretation
//                 of the buckets is not defined at this level - up to any
//                 derived classes. Really meant as a scenario not as a greek.
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Format.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/MultiExpiryShift.hpp"

DRLIB_BEGIN_NAMESPACE

MultiExpiryShift::~MultiExpiryShift(){}

MultiExpiryShift::MultiExpiryShift(const CClassConstSP& clazz):
    Perturbation(clazz){}

void MultiExpiryShift::validatePop2Object() {
    static const string method("MultiExpiryShift::validatePop2Object");
    try {
        if (expiries.empty()) {
            throw ModelException(method, "expiries are empty");
        }
        if (shifts.empty()) {
            throw ModelException(method, "shifts are empty");
        }
        if (shifts.size() != expiries.size()) {
            throw ModelException(method, 
                                 "shifts & expiries are different lengths");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns array of output names which need to be tweaked for this
    sensitivity. In particular, if there are toTweak names, then returns
    these otherwise generates a list based on object supplied (which
    is typically, but needn't be, a TweakGroup */
OutputNameArrayConstSP MultiExpiryShift::names(const IObjectConstSP tweakGroup) const{
    if (!toTweak || toTweak->size()==0) {
        return OutputName::trim(SensMgrConst(tweakGroup).allNames(
            const_cast<MultiExpiryShift*>(this)));
    }
    else {
        return OutputName::trim(toTweak);
    }
}

class MultiExpiryShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MultiExpiryShift, clazz);
        SUPERCLASS(Perturbation);
        FIELD(expiries, "benchmark regions");
        FIELD(shifts, "shift per region");
        FIELD(toTweak, "Ignored - do not use");
        FIELD_MAKE_OPTIONAL(toTweak);
    }
};

CClassConstSP const MultiExpiryShift::TYPE = CClass::registerClassLoadMethod(
    "MultiExpiryShift", typeid(MultiExpiryShift), MultiExpiryShiftHelper::load);

/**
 * @name Implementation of MultiExpiryShift::asRiskProperty()
 */

//@{

FORWARD_DECLARE(MultiExpiryShiftRiskAxis)

class MultiExpiryShiftRiskAxis: public CObject,
                                public virtual IRiskAxis,
                                public virtual RiskAxis {

    MultiExpiryShiftRiskAxis(): CObject(TYPE) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

public:

    MultiExpiryShiftConstSP perturbation;    // $required
    OutputNameConstSP name;                  // $required
    IAbstractRiskPropertyConstSP property;   // $required

    static CClassConstSP const TYPE;

    MultiExpiryShiftRiskAxis(MultiExpiryShiftConstSP perturbation,
                               OutputNameConstSP name,
                               IAbstractRiskPropertyConstSP property):
        CObject(TYPE),
        perturbation(perturbation),
        name(name),
        property(property)
    {}

    IHypothesisConstSP hypothesis(double coeff) const;

    IAbstractRiskPropertyConstSP abstractProperty() const {
        return property;
    }

    OutputNameConstSP marketDataName() const {
        return name;
    }

    string toString() const {
        return (!name ? string("all") : name->toString() + "'s") + " " + perturbation->toString();
    }

    RiskAxisConstSP frozen() const {
        return RiskAxisConstSP::attachToRef(this);
    }

    IRiskAxisConstSP thawed() const {
        return IRiskAxisConstSP::attachToRef(this);
    }
};

FORWARD_DECLARE(MultiExpiryShiftHypothesis)

class MultiExpiryShiftHypothesis: public AtomicHypothesis {

    MultiExpiryShiftHypothesis(): AtomicHypothesis(TYPE) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    MultiExpiryShiftRiskAxisConstSP _axis;   // $required

public:

    static CClassConstSP const TYPE;

    MultiExpiryShiftHypothesis(MultiExpiryShiftRiskAxisConstSP axis):
        AtomicHypothesis(TYPE),
        _axis(axis)
    {}

    AlternateWorldSP appliedTo(IObjectSP world) const {
        IObjectSP world_(copy(world.get()));
        bool changed;
        double dist = applyTo(world_, &changed);
        return AlternateWorldSP(new AlternateWorld(world_, dist, changed));
    }

    double applyTo(IObjectSP world, bool* changed) const {
        bool was = const_cast<MultiExpiryShift&>(*_axis->perturbation.get()).
                       findAndShift(world, _axis->name);
        if (changed) *changed = was;
        return was ? 1. : 0.;
    }

    IRiskAxisConstSP axis() const {
        return _axis;
    }

    double axisCoefficient() const {
        return 1.;
    }

    string toString() const {
        return _axis->toString();
    }
};

IHypothesisConstSP MultiExpiryShiftRiskAxis::hypothesis(
        double coeff) const {
    try {
        if (coeff == 0.) {
            return IHypothesis::noop();
        }
        else if (coeff == 1.) {
            return IHypothesisConstSP(new MultiExpiryShiftHypothesis(
                MultiExpiryShiftRiskAxisConstSP::attachToRef(this)));
        }
        else {
            throw ModelException(
                "Requested shift size is " + Format::toString(coeff) +
                " but if you're using a MultiExpiryShift (in this case a " +
                perturbation->getClass()->getName() + ") as a RiskProperty you "
                "can only use shift size == 1");
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

FORWARD_DECLARE(MultiExpiryShiftRiskProperty)

class MultiExpiryShiftRiskProperty: public CObject,
                                    public virtual IScalarRiskProperty {

    MultiExpiryShiftRiskProperty(): CObject(TYPE) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    MultiExpiryShiftConstSP perturbation;    // $required

public:

    static CClassConstSP const TYPE;

    MultiExpiryShiftRiskProperty(MultiExpiryShiftConstSP perturbation):
        CObject(TYPE),
        perturbation(perturbation)
    {}

    bool discrete() const {
        return false;
    }

    CClassConstSP subjectInterface() const {
        return perturbation->shiftInterface();
    }

    virtual OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        return perturbation->names(world);
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP name,
                             VoidConstSP) const {
        return IRiskAxisConstSP(new MultiExpiryShiftRiskAxis(
            perturbation, name,
            IAbstractRiskPropertyConstSP::attachToRef(this)));
    }

    VoidArrayConstSP subjectQualifiers(IObjectConstSP world,
                                       OutputNameConstSP name) const {
        return VoidArrayConstSP();
    }
};

IScalarRiskPropertyConstSP MultiExpiryShift::asRiskProperty() const {
    return IScalarRiskPropertyConstSP(new MultiExpiryShiftRiskProperty(
        MultiExpiryShiftConstSP::attachToRef(this)));
}

//@}

IObject* MultiExpiryShiftHypothesis::emptyShell() {
  return new MultiExpiryShiftHypothesis();
}

void MultiExpiryShiftHypothesis::load(CClassSP& clazz) {
  REGISTER(MultiExpiryShiftHypothesis, clazz);
  SUPERCLASS(AtomicHypothesis);
  FIELD(_axis, "_axis");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const MultiExpiryShiftHypothesis::TYPE = CClass::registerClassLoadMethod(
  "MultiExpiryShiftHypothesis", typeid(MultiExpiryShiftHypothesis), MultiExpiryShiftHypothesis::load);

IObject* MultiExpiryShiftRiskAxis::emptyShell() {
  return new MultiExpiryShiftRiskAxis();
}

void MultiExpiryShiftRiskAxis::load(CClassSP& clazz) {
  REGISTER(MultiExpiryShiftRiskAxis, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskAxis);
  IMPLEMENTS(RiskAxis);
  FIELD(perturbation, "perturbation");
  FIELD(name, "name");
  FIELD(property, "property");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const MultiExpiryShiftRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "MultiExpiryShiftRiskAxis", typeid(MultiExpiryShiftRiskAxis), MultiExpiryShiftRiskAxis::load);

IObject* MultiExpiryShiftRiskProperty::emptyShell() {
  return new MultiExpiryShiftRiskProperty();
}

void MultiExpiryShiftRiskProperty::load(CClassSP& clazz) {
  REGISTER(MultiExpiryShiftRiskProperty, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskProperty<Void>);
  FIELD(perturbation, "perturbation");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const MultiExpiryShiftRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "MultiExpiryShiftRiskProperty", typeid(MultiExpiryShiftRiskProperty), MultiExpiryShiftRiskProperty::load);

DRLIB_END_NAMESPACE

