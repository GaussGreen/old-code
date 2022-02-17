//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : ScalarPerturbation.cpp
//
//   Description : Defines the ability to shift a named double
//
//   Author      : Mark A Robson
//
//   Date        : 29 June 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/ScalarPerturbation.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE

ScalarPerturbation::~ScalarPerturbation(){}

/** Returns the scalar shift size */
double ScalarPerturbation::getShiftSize() const{
    return shiftSize;
}

/** Sets the scalar shift size */
void ScalarPerturbation::setShiftSize(double shiftSize){
    this->shiftSize = shiftSize;
}

DECLARE(ScalarPerturbation)

class ScalarPerturbation_AltWorld: public IHypothesis::AlternateWorld {

public:

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(ScalarPerturbation_AltWorld, clazz);
        SUPERCLASS(IHypothesis::AlternateWorld);
        // no EMPTY_SHELL_METHOD: can't copy sensMgr
    };

    SensMgrOpt sensMgr; // $unregistered
    ScalarPerturbationSP shift; // $unregistered
    double shiftSize; // $unregistered

    ScalarPerturbation_AltWorld(ScalarPerturbationSP shift,
                                double shiftSize,
                                OutputNameConstSP overrideName,
                                IObjectSP tweakGroup):
        IHypothesis::AlternateWorld(TYPE, tweakGroup, 0., true),
        sensMgr(tweakGroup.get()),
        shift(shift),
        shiftSize(shiftSize)
    {
        double oldShiftSize = shift->getShiftSize();
        shift->setShiftSize(shiftSize);
        try {
            world = sensMgr.shift(shift.get(), overrideName);
        }
        catch (...) {
            shift->setShiftSize(oldShiftSize);
            throw;
        }
        shift->setShiftSize(oldShiftSize);
    }

    void _undo() {
        double oldShiftSize = shift->getShiftSize();
        shift->setShiftSize(shiftSize);
        try {
            sensMgr.restore();
        }
        catch (...) {
            shift->setShiftSize(oldShiftSize);
            throw;
        }
        shift->setShiftSize(oldShiftSize);
    }
};

CClassConstSP const ScalarPerturbation_AltWorld::TYPE =
    CClass::registerClassLoadMethod(
        "ScalarPerturbation_AltWorld", typeid(ScalarPerturbation_AltWorld),
        load);

/** IScalarPerNameShift implementation, for ImpliedScalarPerturbation  */
IHypothesis::AlternateWorldSP ScalarPerturbation::appliedTo(
        OutputNameConstSP nameToShift,
        double shiftSize,
        IObjectSP tweakGroup) {
    try {
        return IHypothesis::AlternateWorldSP(
            new ScalarPerturbation_AltWorld(
                ScalarPerturbationSP::attachToRef(this),
                shiftSize, nameToShift, tweakGroup));
    }
    catch (exception& e) {
        throw ModelException(e, "ScalarPerturbation::appliedTo()");
    }
}

OutputNameArrayConstSP ScalarPerturbation::allNames(const IObject* tweakGroup) const {
    return OutputName::trim(SensMgrConst(tweakGroup).allNames(
        const_cast<ScalarPerturbation*>(this)));
}

/** Returns array of output names which need to be tweaked for this
    sensitivity. In particular, if there are toTweak names, then returns
    these otherwise generates a list based on object supplied (which
    is typically, but needn't be, a TweakGroup */
OutputNameArrayConstSP ScalarPerturbation::names(const IObject* tweakGroup) const{
    
    if (!toTweak || toTweak->size()==0) {
        return allNames(tweakGroup);
    }
    else {
        return OutputName::trim(toTweak);
    }
}

ScalarPerturbation::ScalarPerturbation(CClassConstSP clazz):
    Perturbation(clazz){}

ScalarPerturbation::ScalarPerturbation(CClassConstSP clazz, double shiftSize):
    Perturbation(clazz), shiftSize(shiftSize) {}

/**
 * @name Implementation of ScalarPerturbation::asRiskProperty()
 */

//@{

FORWARD_DECLARE(ScalarPerturbationRiskAxis)

class ScalarPerturbationRiskAxis: public CObject,
                                  public virtual IRiskAxis,
                                  public virtual RiskAxis {

    ScalarPerturbationRiskAxis(): CObject(TYPE) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

public:

    ScalarPerturbationConstSP perturbation;    // $required
    OutputNameConstSP name;                    // $required
    IAbstractRiskPropertyConstSP property;     // $required

    static CClassConstSP const TYPE;

    ScalarPerturbationRiskAxis(ScalarPerturbationConstSP perturbation,
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

FORWARD_DECLARE(ScalarPerturbationHypothesis)

class ScalarPerturbationHypothesis: public AtomicHypothesis {

    ScalarPerturbationHypothesis(): AtomicHypothesis(TYPE), coefficient(0.) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    ScalarPerturbationRiskAxisConstSP _axis;   // $required
    double coefficient;                        // $required

public:

    static CClassConstSP const TYPE;

    ScalarPerturbationHypothesis(ScalarPerturbationRiskAxisConstSP axis,
                                 double coefficient):
        AtomicHypothesis(TYPE),
        _axis(axis),
        coefficient(coefficient)
    {}

    AlternateWorldSP appliedTo(IObjectSP world) const {
        ScalarPerturbation& p = const_cast<ScalarPerturbation&>(
                                    *_axis->perturbation.get());

        return p.appliedTo(_axis->name, coefficient * p.getShiftSize(), world);
    }

    double applyTo(IObjectSP world, bool* changed) const {
        ScalarPerturbation& p = const_cast<ScalarPerturbation&>(
                                    *_axis->perturbation.get());

        double originalShiftSize = p.getShiftSize();
        p.setShiftSize(coefficient * originalShiftSize);
        bool was;
        try {
            was = p.findAndShift(world, _axis->name);
            p.setShiftSize(originalShiftSize);
        }
        catch (exception&) {
            p.setShiftSize(originalShiftSize);
            throw;
        }

        if (changed) *changed = was;
        return was ? coefficient : 0.;
    }

    IRiskAxisConstSP axis() const {
        return _axis;
    }

    double axisCoefficient() const {
        return coefficient;
    }

    string toString() const {
        return Format::toString(coefficient) + " shift of " + _axis->toString();
    }
};

IHypothesisConstSP ScalarPerturbationRiskAxis::hypothesis(
        double coeff) const {
    return IHypothesisConstSP(new ScalarPerturbationHypothesis(
        ScalarPerturbationRiskAxisConstSP::attachToRef(this), coeff));
}

FORWARD_DECLARE(ScalarPerturbationRiskProperty)

class ScalarPerturbationRiskProperty: public CObject,
                                      public virtual IScalarRiskProperty {

    ScalarPerturbationRiskProperty(): CObject(TYPE) {}

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    ScalarPerturbationConstSP perturbation;    // $required

public:

    static CClassConstSP const TYPE;

    ScalarPerturbationRiskProperty(ScalarPerturbationConstSP perturbation):
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
        return perturbation->names(world.get());
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP name,
                             VoidConstSP) const {
        return IRiskAxisConstSP(new ScalarPerturbationRiskAxis(
            perturbation, name,
            IAbstractRiskPropertyConstSP::attachToRef(this)));
    }

    VoidArrayConstSP subjectQualifiers(IObjectConstSP world,
                                       OutputNameConstSP name) const {
        return VoidArrayConstSP();
    }
};

IScalarRiskPropertyConstSP ScalarPerturbation::asRiskProperty() const {
    return IScalarRiskPropertyConstSP(new ScalarPerturbationRiskProperty(
        ScalarPerturbationConstSP::attachToRef(this)));
}

//@}



void ScalarPerturbation::load(CClassSP& clazz){
    REGISTER(ScalarPerturbation, clazz);
    SUPERCLASS(Perturbation);
    IMPLEMENTS(IScalarPerNameShift);
    FIELD(shiftSize, "How big to make the tweak");
    FIELD(toTweak, "Ignored - do not use");
    FIELD_MAKE_OPTIONAL(toTweak);
};

CClassConstSP const ScalarPerturbation::TYPE = CClass::registerClassLoadMethod(
    "ScalarPerturbation", typeid(ScalarPerturbation), load);

IObject* ScalarPerturbationHypothesis::emptyShell() {
  return new ScalarPerturbationHypothesis();
}

void ScalarPerturbationHypothesis::load(CClassSP& clazz) {
  REGISTER(ScalarPerturbationHypothesis, clazz);
  SUPERCLASS(AtomicHypothesis);
  FIELD(_axis, "_axis");
  FIELD(coefficient, "coefficient");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ScalarPerturbationHypothesis::TYPE = CClass::registerClassLoadMethod(
  "ScalarPerturbationHypothesis", typeid(ScalarPerturbationHypothesis), ScalarPerturbationHypothesis::load);

IObject* ScalarPerturbationRiskAxis::emptyShell() {
  return new ScalarPerturbationRiskAxis();
}

void ScalarPerturbationRiskAxis::load(CClassSP& clazz) {
  REGISTER(ScalarPerturbationRiskAxis, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskAxis);
  IMPLEMENTS(RiskAxis);
  FIELD(perturbation, "perturbation");
  FIELD(name, "name");
  FIELD(property, "property");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ScalarPerturbationRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "ScalarPerturbationRiskAxis", typeid(ScalarPerturbationRiskAxis), ScalarPerturbationRiskAxis::load);

IObject* ScalarPerturbationRiskProperty::emptyShell() {
  return new ScalarPerturbationRiskProperty();
}

void ScalarPerturbationRiskProperty::load(CClassSP& clazz) {
  REGISTER(ScalarPerturbationRiskProperty, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskProperty<Void>);
  FIELD(perturbation, "perturbation");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ScalarPerturbationRiskProperty::TYPE = CClass::registerClassLoadMethod(
  "ScalarPerturbationRiskProperty", typeid(ScalarPerturbationRiskProperty), ScalarPerturbationRiskProperty::load);


DRLIB_END_NAMESPACE
