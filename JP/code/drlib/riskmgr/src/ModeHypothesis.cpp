/**
 * @file ModeHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/ModeHypothesis.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

class ModeHypothesis_AltWorld:
        public IHypothesis::AlternateWorld {

    ModeHypothesis_AltWorld():
        IHypothesis::AlternateWorld(TYPE)
    {}

    static IObject* defaultOne() {
        return new ModeHypothesis_AltWorld();
    }

    static const CClassConstSP TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(ModeHypothesis_AltWorld, clazz);
        SUPERCLASS(IHypothesis::AlternateWorld);
        FIELD(hyp, "hyp");
        EMPTY_SHELL_METHOD(defaultOne);
    }

    ModeHypothesisConstSP hyp;

public:

    ModeHypothesis_AltWorld(ModeHypothesisConstSP hyp, IObjectSP world,
                            double distance):
        AlternateWorld(TYPE, world, distance, true),
        hyp(hyp)
    {}

    void _undo() {
        try {
            MultiTweakGroupSP tg = MultiTweakGroupSP::dynamicCast(world);
            try {
                hyp->undo(tg);
            }
            catch (exception& e) {
                throw ModelException(
                    e, hyp->getClass()->getName() + "::undo(MultiTweakGroupSP)");
            }
        }
        catch (exception& e) {
            throw ModelException(e, "ModeHypothesis_AltWorld::_undo()");
        }
    }
};

ModeHypothesis::ModeHypothesis(CClassConstSP type):
    AtomicHypothesis(type)
{}

ModeHypothesis::~ModeHypothesis() {}

IRiskAxisConstSP ModeHypothesis::axis() const {
    return IRiskAxisConstSP();
}

double ModeHypothesis::axisCoefficient() const {
    return 0.;
}

IHypothesis::AlternateWorldSP ModeHypothesis::appliedTo(IObjectSP world) const {
    try {
        double distance = applyTo(world);
        return AlternateWorldSP(new ModeHypothesis_AltWorld(
            ModeHypothesisConstSP::attachToRef(this), world, distance));
        }
    catch (exception& e) {
        throw ModelException(e, "E2CParSpreadPricingOffHypothesis::appliedTo()");
    }
}

double ModeHypothesis::applyTo(IObjectSP world, bool* changed) const {
    try {
        MultiTweakGroupSP tg = MultiTweakGroupSP::dynamicCast(world);

        try {
            return applyToWorld(tg, changed);
        }
        catch (exception& e) {
            throw ModelException(
                e, getClass()->getName() + "::applyToWorld()");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "ModeHypothesis::applyTo()");
    }
}

void ModeHypothesis::load(CClassSP& clazz) {
    REGISTER(ModeHypothesis, clazz);
    SUPERCLASS(AtomicHypothesis);
}

CClassConstSP const ModeHypothesis::TYPE = CClass::registerClassLoadMethod(
    "ModeHypothesis", typeid(ModeHypothesis), load);

CClassConstSP const ModeHypothesis_AltWorld::TYPE = CClass::registerClassLoadMethod(
    "ModeHypothesis_AltWorld", typeid(ModeHypothesis_AltWorld), load);

DRLIB_END_NAMESPACE
