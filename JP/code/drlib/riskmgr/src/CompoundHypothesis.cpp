/**
 * @file CompoundHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/CompoundHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

// 
// ==============
//  Constructors
// ==============
// 

CompoundHypothesis::CompoundHypothesis():
    CObject(TYPE)
{}

CompoundHypothesis::CompoundHypothesis(
        AtomicHypothesisArrayConstSP differencesWrtBaseState,
        IHypothesis::IDistanceMetricConstSP distanceMetric):
    CObject(TYPE),
    differencesWrtBaseState(differencesWrtBaseState),
    _distanceMetric(distanceMetric)
{}

CompoundHypothesisSP CompoundHypothesis::SP(
        AtomicHypothesisArrayConstSP differencesWrtBaseState,
        IHypothesis::IDistanceMetricConstSP distanceMetric) {
    return CompoundHypothesisSP(new CompoundHypothesis(
               differencesWrtBaseState, distanceMetric));
}

CompoundHypothesis::CompoundHypothesis(
        AtomicHypothesisArrayConstSP differencesWrtBaseState):
    CObject(TYPE),
    differencesWrtBaseState(differencesWrtBaseState),
    _distanceMetric(IHypothesis::IDistanceMetric::last())
{}

CompoundHypothesisSP CompoundHypothesis::SP(
        AtomicHypothesisArrayConstSP differencesWrtBaseState) {
    return CompoundHypothesisSP(
        new CompoundHypothesis(differencesWrtBaseState));
}

CompoundHypothesisSP CompoundHypothesis::SP(
        IHypothesisArrayConstSP differencesWrtBaseState,
        IHypothesis::IDistanceMetricConstSP distanceMetric) {
    int numAtomics = 0;
    for (int h = 0; h < differencesWrtBaseState->size(); ++h) {
        numAtomics += (*differencesWrtBaseState)[h]->numAtomics();
    }

    AtomicHypothesisArraySP atomics = AtomicHypothesisArray::SP(numAtomics);

    int b = 0;
    for (int h = 0; h < differencesWrtBaseState->size(); ++h) {
        IHypothesisConstSP hyp = (*differencesWrtBaseState)[h];
        for (int a = 0; a < hyp->numAtomics(); ++a) {
            (*atomics)[b++] = AtomicHypothesisSP::constCast(hyp->atomic(a));
        }
    }

    return SP(atomics, distanceMetric);
}

CompoundHypothesis::~CompoundHypothesis() {}

int CompoundHypothesis::numAtomics() const {
    return differencesWrtBaseState->size();
}

AtomicHypothesisConstSP CompoundHypothesis::atomic(int i) const {
    return (*differencesWrtBaseState)[i];
}

// 
// ===========
//  appliedTo
// ===========
// 

struct CompoundAlternateWorld: IHypothesis::AlternateWorld {

    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(CompoundAlternateWorld, clazz);
        SUPERCLASS(IHypothesis::AlternateWorld);
        EMPTY_SHELL_METHOD(defaultCompoundAlternateWorld);
        FIELD(alts, "alts");
    }        

    static IObject *defaultCompoundAlternateWorld() {
        return new CompoundAlternateWorld(
            IHypothesis::AlternateWorldArraySP(), IObjectSP(), 0, false);
    }

    IHypothesis::AlternateWorldArraySP alts;

    CompoundAlternateWorld(IHypothesis::AlternateWorldArraySP alts,
                           IObjectSP world, double distance, bool found):
        IHypothesis::AlternateWorld(TYPE, world, distance, found),
        alts(alts)
    {}

    void _undo() {
        for (int h = 0; h < alts->size(); ++h)
            (*alts)[h]->undo();
    }
};

IHypothesis::AlternateWorldSP CompoundHypothesis::appliedTo(IObjectSP world) const {
    IHypothesis::AlternateWorldArraySP alts(new AlternateWorldArray(
        differencesWrtBaseState->size()));

    CDoubleArray distances(differencesWrtBaseState->size());
    bool found = false;

    int h;
    try {
        for (h = 0; h < differencesWrtBaseState->size(); ++h) {
            (*alts)[h] = (*differencesWrtBaseState)[h]->appliedTo(world);
            world = (*alts)[h]->world;
            distances[h] = (*alts)[h]->distance;
            found = found || (*alts)[h]->found;
        }
    }
    catch (exception&) {
        for (--h; h >= 0; --h) (*alts)[h]->undo();
        throw;
    }

    return AlternateWorldSP(new CompoundAlternateWorld(
        alts, world, (*distanceMetric())(distances), found));
}

// 
// =========
//  applyTo
// =========
// 

double CompoundHypothesis::applyTo(IObjectSP world, bool* changed) const {
    CDoubleArray distances(differencesWrtBaseState->size());

    if (changed) *changed = false;
    for (int h = 0; h < differencesWrtBaseState->size(); ++h) {
        bool ch;
        distances[h] = (*differencesWrtBaseState)[h]->applyTo(world, &ch);
        if (ch && changed) *changed = true;
    }

    return (*distanceMetric())(distances);
}

// 
// ==============
//  Housekeeping
// ==============
// 

CompoundHypothesisSP CompoundHypothesis::empty() {
    static CompoundHypothesisSP it(new CompoundHypothesis(
        AtomicHypothesisArray::SP(), IHypothesis::IDistanceMetric::last()));
    return it;
}

IRiskAxisConstSP CompoundHypothesis::axis() const {
    return IRiskAxisConstSP();
}

double CompoundHypothesis::axisCoefficient() const {
    return 0.;
}

IHypothesis::IDistanceMetricConstSP CompoundHypothesis::distanceMetric() const {
    return _distanceMetric;
}

string CompoundHypothesis::toString() const {
    if (numAtomics() == 0) return "no changes, i.e. base case";
    string it;
    for (int h = 0; h < numAtomics(); ++h) {
        if (h) it += " & ";
        it += atomic(h)->toString();
    }
    return it;
}

// 
// =========
//  grouped
// =========
// 

struct GroupedHypothesis: AtomicHypothesis {
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(GroupedHypothesis, clazz);
        SUPERCLASS(AtomicHypothesis);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(compound, "compound");
    }        

    static IObject *emptyShell() {
        return new GroupedHypothesis(CompoundHypothesisSP());
    }

    CompoundHypothesisConstSP compound;

    GroupedHypothesis(CompoundHypothesisConstSP compound):
        AtomicHypothesis(TYPE),
        compound(compound)
    {}

    AlternateWorldSP appliedTo(IObjectSP world) const {
        IObjectSP world2(world.clone());

        bool changed;
        double distance = compound->applyTo(world2, &changed);

        return AlternateWorldSP(new AlternateWorld(world2, distance, changed));
    }

    double applyTo(IObjectSP world, bool* changed = 0) const {
        return compound->applyTo(world, changed);
    }

    IRiskAxisConstSP axis() const { return IRiskAxisConstSP(); }

    double axisCoefficient() const { return 0.; }

    string toString() const {
        return compound->toString() + " (all in one go)";
    }
};

AtomicHypothesisSP CompoundHypothesis::grouped() const {
    return AtomicHypothesisSP(new GroupedHypothesis(
                                  CompoundHypothesisConstSP::attachToRef(this)));
}

void CompoundHypothesis::setApproxOrder(int order){

	AtomicHypothesisArraySP dWBS(copy(differencesWrtBaseState.get()));
	for (int i = 0; i<dWBS->size(); i++){
		(*dWBS)[i]->setApproxOrder(order);
	}
    differencesWrtBaseState = dWBS;
}

// 
// ******************
//  Reflection stuff
// ******************
// 

IObject *CompoundHypothesis::defaultOne() {
    return new CompoundHypothesis();
}

void CompoundHypothesis::load(CClassSP& clazz) {
    REGISTER(CompoundHypothesis, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IHypothesis);
    EMPTY_SHELL_METHOD(defaultOne);
    FIELD(differencesWrtBaseState, "differencesWrtBaseState");
    FIELD(_distanceMetric, "distanceMetric");
}

CClassConstSP const CompoundHypothesis::TYPE =
    CClass::registerClassLoadMethod(
        "CompoundHypothesis", typeid(CompoundHypothesis), load);

CClassConstSP const CompoundAlternateWorld::TYPE =
    CClass::registerClassLoadMethod(
        "CompoundAlternateWorld", typeid(CompoundAlternateWorld), load);

CClassConstSP const GroupedHypothesis::TYPE =
    CClass::registerClassLoadMethod(
        "GroupedHypothesis", typeid(GroupedHypothesis), load);

DRLIB_END_NAMESPACE
