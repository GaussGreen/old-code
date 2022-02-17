/**
 * @file IHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_IHYPOTHESIS_CPP
#include "edginc/CompoundHypothesis.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/MarketObject.hpp"

// for NamesOnObjectHypothesis
#include "edginc/IRiskProperty.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/SimpleTweakNameResolver.hpp"

DRLIB_BEGIN_NAMESPACE

IHypothesis::IHypothesis() {}

IHypothesisConstSP IHypothesis::noop() {
    return CompoundHypothesis::empty();
}

IHypothesis::~IHypothesis() {}

bool IHypothesis::applyScenario(IObjectSP object) {
    bool changed;
    applyTo(object, &changed);
    return changed;
}

bool IHypothesis::preapplyScenario(IObjectSP object) {
    return false;
}

void IHypothesis::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IHypothesis, clazz);
    EXTENDS(IObject);
    EXTENDS(IScenarioShift);
}

CClassConstSP const IHypothesis::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IHypothesis", typeid(IHypothesis), load);

IHypothesis::AlternateWorld::AlternateWorld(IObjectSP world, double distance,
                                            bool found):
    CObject(TYPE),
    undone(false),
    oldValue(-666.66),
    world(world),
    distance(distance),
    found(found)
{}

IHypothesis::AlternateWorld::AlternateWorld(CClassConstSP type,
                                            IObjectSP world, double distance,
                                            bool found):
    CObject(type),
    undone(false),
    oldValue(-666.66),
    world(world),
    distance(distance),
    found(found)
{}

IHypothesis::AlternateWorld::AlternateWorld(CClassConstSP type):
    CObject(type),
    undone(false),
    oldValue(-666.66),
    distance(0),
    found(false)
{}

IHypothesis::AlternateWorld::~AlternateWorld() {}

void IHypothesis::AlternateWorld::undo() {
    ASSERT(!undone);
    _undo();
    undone = true;
}

void IHypothesis::AlternateWorld::_undo() {
}

static IObject* defaultAlternateWorld() {
    return new IHypothesis::AlternateWorld(IObjectSP(), 0., false);
}

void IHypothesis::AlternateWorld::load(CClassSP& clazz) {
    REGISTER(IHypothesis::AlternateWorld, clazz);
    SUPERCLASS(CObject);
    FIELD(world, "world");
    FIELD(distance, "distance");
    FIELD(undone, "undone");
    FIELD(found, "found");
    EMPTY_SHELL_METHOD(defaultAlternateWorld);
}

typedef IHypothesis::AlternateWorld IHypothesis_AlternateWorld;

CClassConstSP const IHypothesis_AlternateWorld::TYPE =
    CClass::registerClassLoadMethod(
        "IHypothesis::AlternateWorld", typeid(IHypothesis_AlternateWorld),
        load);

typedef IHypothesis::AlternateWorldArray IHypothesis_AlternateWorldArray;

DEFINE_TEMPLATE_TYPE_WITH_NAME("IHypothesis::AlternateWorldArray", IHypothesis_AlternateWorldArray);

IHypothesisConstSP IHypothesis::then(IHypothesisConstSP other) const {
    int as = numAtomics(), otherAs = other->numAtomics();

    AtomicHypothesisArraySP atomics(new AtomicHypothesisArray(as + otherAs));

    for (int a = 0; a < as; ++a) {
        (*atomics)[a] = AtomicHypothesisSP::constCast(atomic(a));
    }

    {for (int a = 0; a < otherAs; ++a) {
        (*atomics)[as + a] = AtomicHypothesisSP::constCast(other->atomic(a));
    }}

    return CompoundHypothesis::SP(atomics);
}

IHypothesis::IDistanceMetric::IDistanceMetric() {}
IHypothesis::IDistanceMetric::~IDistanceMetric() {}

CClassConstSP const IHypothesis::IDistanceMetric::TYPE = CClass::registerInterfaceLoadMethod(
    "IHypothesis::IDistanceMetric", typeid(IHypothesis::IDistanceMetric), 0);

class LastDistanceMetric: public CObject,
                          public virtual IHypothesis::IDistanceMetric {

public:

    static CClassConstSP const TYPE;

    LastDistanceMetric():
        CObject(TYPE)
    {}

    double operator()(const DoubleArray& distances) const {
        return distances.empty() ? 0. : distances.back();
    }

    static void load(CClassSP& clazz) {
        REGISTER(LastDistanceMetric, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IHypothesis::IDistanceMetric);
        EMPTY_SHELL_METHOD(&DefaultConstructor<LastDistanceMetric>::iObject);
    }
};

CClassConstSP const LastDistanceMetric::TYPE = CClass::registerClassLoadMethod(
    "LastDistanceMetric", typeid(LastDistanceMetric), load);

IHypothesis::IDistanceMetricConstSP IHypothesis::IDistanceMetric::last() {
    static IDistanceMetricConstSP it(new LastDistanceMetric());
    return it;
}

class ConstantDistanceMetric: public CObject,
                              public virtual IHypothesis::IDistanceMetric {

public:

    static CClassConstSP const TYPE;

    double c;

    ConstantDistanceMetric(double c = 0):
        CObject(TYPE),
        c(c)
    {}

    double operator()(const DoubleArray&) const {
        return c;
    }

    static void load(CClassSP& clazz) {
        REGISTER(ConstantDistanceMetric, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IHypothesis::IDistanceMetric);
        EMPTY_SHELL_METHOD(&DefaultConstructor<ConstantDistanceMetric>::iObject);
        FIELD(c, "c");
    }
};

IHypothesis::IDistanceMetricConstSP IHypothesis::IDistanceMetric::constant(double c) {
    return IDistanceMetricConstSP(new ConstantDistanceMetric(c));
}

CClassConstSP const ConstantDistanceMetric::TYPE = CClass::registerClassLoadMethod(
    "ConstantDistanceMetric", typeid(ConstantDistanceMetric), load);

DEFINE_TEMPLATE_TYPE(IHypothesisArray);

DRLIB_END_NAMESPACE
