/**
 * @file PropertyTweakHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Format.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakOptID.hpp"
#include "edginc/TweakNameResolver.hpp"
#include "edginc/ICompatibilitySensitivity.hpp"
#include "edginc/MarketObject.hpp" // just for TRACE

DRLIB_BEGIN_NAMESPACE

AbstractPropertyTweakHypothesis::AbstractPropertyTweakHypothesis(
        CClassConstSP type):
    AtomicHypothesis(type)
{}

AbstractPropertyTweakHypothesis::~AbstractPropertyTweakHypothesis() {}

/**
 * Provides the SensMgr interfaces for AbstractPropertyTweakHypothesis.
 * We could very nearly just have AbstractPropertyTweakHypothesis implement them itself,
 * but we want "shift" to intercept the "distance" returned as part of the
 * TweakOutcome in a reentrant way.  Also we're able to keep the
 * public interface of AbstractPropertyTweakHypothesis a little simpler.
 */

class AbstractPropertyTweakHypothesisTweakOptID:
        public CObject,
        public virtual ITweakOptID,
        private virtual ITweakNameResolver {

    static IObject* defaultAbstractPropertyTweakHypothesisTweakOptID() {
        return new AbstractPropertyTweakHypothesisTweakOptID(
            AbstractPropertyTweakHypothesisSP(), NULL);
    }

    static void load(CClassSP& clazz) {
        REGISTER(AbstractPropertyTweakHypothesisTweakOptID, clazz);
        SUPERCLASS(CObject);
        //IMPLEMENTS(ITweakOptID);
        EMPTY_SHELL_METHOD(defaultAbstractPropertyTweakHypothesisTweakOptID);
        FIELD(tweak, "tweak");
        FIELD(found, "found");
        // FIELD(outcome, "outcome"); -- OK since only used by AbstractPropertyTweakHypothesis_AltWorld
    }

public:

    static CClassConstSP const TYPE;

    AbstractPropertyTweakHypothesisConstSP tweak;
    TweakOutcome* outcome; // $unregistered
    int found;

    AbstractPropertyTweakHypothesisTweakOptID(AbstractPropertyTweakHypothesisConstSP tweak,
                                              TweakOutcome* outcome):
        CObject(TYPE),
        tweak(tweak),
        outcome(outcome),
        found(0)
    {}

    bool restorableShift(IObjectConstSP obj) const {
        return !!tweak->restorableInterface() &&
               tweak->restorableInterface()->isInstance(obj);
    }

    void restore(IObjectSP obj) {
        tweak->sensRestore(obj);
    }

    void reset() {
    }

    bool shift(IObjectSP obj) {
        TRACE_METHOD;

        // Work around an inconsistency in underlying SensMgr: when we
        // run a PerNameRiskPropertySensitivity, objects whose sensName comes
        // back empty are excluded (are considered not to carry the property
        // being tweaked even though they are of the relevant type) --- but
        // when we run an AllNamesRiskPropertySensitivity they aren't excluded.
        // See top of SensMgrOpt::OptimalShift::invoke().

        if (!tweak->getMarketDataName() && // null => "all names"
                tweak->sensName(obj) == "") {
            TRACE("Not tweaking " << *obj << " since it doesn't want to be");
            return true;
        }
        else {
            TRACE("Tweaking " << *obj);

            TweakOutcome oc = tweak->sensShift(obj);
            if (outcome) *outcome = oc;
            ++found;
            return oc.tweakMembers();
        }
    }

    ITweakNameResolver* nameResolver() {
        return this;
    }

    CClassConstSP shiftInterface() const {
        return tweak->shiftableInterface();
    }

    OutputNameConstSP getMarketDataName() const {
        return tweak->getMarketDataName();
    }

    bool nameMatches(const OutputName &name, IObjectConstSP obj) {
        return name.equals(tweak->sensName(obj));
    }
};

CClassConstSP const AbstractPropertyTweakHypothesisTweakOptID::TYPE = CClass::registerClassLoadMethod(
    "AbstractPropertyTweakHypothesisTweakOptID", typeid(AbstractPropertyTweakHypothesisTweakOptID), load);

class AbstractPropertyTweakHypothesis_AltWorld:
        public IHypothesis::AlternateWorld { // $unregistered

    AbstractPropertyTweakHypothesis_AltWorld(
        const AbstractPropertyTweakHypothesis_AltWorld &other);
    AbstractPropertyTweakHypothesis_AltWorld &operator =(
        const AbstractPropertyTweakHypothesis_AltWorld &other);

public:

    TweakOutcome outcome;
    AbstractPropertyTweakHypothesisTweakOptID tweakid;
    SensMgrOpt sensMgr;

    AbstractPropertyTweakHypothesis_AltWorld(IObjectSP world,
                                             AbstractPropertyTweakHypothesisConstSP tweak):
        IHypothesis::AlternateWorld(TYPE, world, 0., true),
        outcome(0, 0, false),
        tweakid(tweak, &outcome),
        sensMgr(world)
    {
        world = sensMgr.shift(&tweakid);
        distance = outcome.distance();
        found = tweakid.found != 0;
        // the sentinel value is read by HypothesisTree to support Sensitivity
        // compatibility, "temporarily"
        oldValue = outcome.hasOldValue() ? outcome.oldValue() : -666e66;
    }

    void _undo() {
        sensMgr.restore();
    }
};

FORWARD_DECLARE(AbstractPropertyTweakHypothesis_AltWorld)

string AbstractPropertyTweakHypothesis::toString() const {
    // bit of a hack but hey

    string classname = getClass()->getName();
    size_t b = classname.rfind('<');
    string propname = b == (size_t)-1 ?
        classname : classname.substr(b + 1, (classname.size() - 1) - (b + 1));

    IObjectConstSP q = getQualifier();

    return Format::toString(axisCoefficient()) + " tweak of " +
           shiftableInterface()->getName() + " " +
           (!getMarketDataName() ? string()
                                 : getMarketDataName()->toString() + "'s ") +
           (!q ? "" : q->toString() + " ") +
           propname;
}

IHypothesis::AlternateWorldSP AbstractPropertyTweakHypothesis::appliedTo(
        IObjectSP world) const {
    TRACE_METHOD;
    try {
        IHypothesis::AlternateWorldSP it(
            new AbstractPropertyTweakHypothesis_AltWorld(
                world, AbstractPropertyTweakHypothesisConstSP::attachToRef(this)));
        return it;
    }
    catch (exception& e) {
        throw ModelException(
            e, "AbstractPropertyTweakHypothesis::appliedTo",
            string() + "Applying " + toString() + " to " +
                (!getMarketDataName() ?
                     "all names" : getMarketDataName()->toString()));
    }
}

TweakOutcome AbstractPropertyTweakHypothesis::applyTo_TweakOutcome(
        IObjectSP world) const {
    TRACE_METHOD;
    try {
        TweakOutcome it(0., 0., false);
        AbstractPropertyTweakHypothesisTweakOptID tweakid(
            AbstractPropertyTweakHypothesisConstSP::attachToRef(this), &it);
        SensMgr(world).shift(&tweakid);
        if (tweakid.found != 1) {
            throw ModelException(tweakid.found == 0 ?
                                     "No matching objects found" :
                                     "More than one object found");
        }

        return it;
    }
    catch (exception& e) {
        throw ModelException(
            e, "AbstractPropertyTweakHypothesis::appliedTo",
            string() + "Applying " + toString() + " to " +
                (!getMarketDataName() ?
                     "all names" : getMarketDataName()->toString()));
    }
}

double AbstractPropertyTweakHypothesis::applyTo(IObjectSP world,
                                                bool* changed) const {
    TRACE_METHOD;
    TRACE_SHOW(this->toString());
    try {
        TweakOutcome it(0., 0., false);
        AbstractPropertyTweakHypothesisTweakOptID tweakid(
            AbstractPropertyTweakHypothesisConstSP::attachToRef(this), &it);
        SensMgr(world).shift(&tweakid);
        if (changed) *changed = tweakid.found != 0;
        return it.distance();
    }
    catch (exception& e) {
        throw ModelException(
            e, "AbstractPropertyTweakHypothesis::appliedTo",
            string() + "Applying " + toString() + " to " +
                (!getMarketDataName() ?
                     "all names" : getMarketDataName()->toString()));
    }
}

ITweakOptID* AbstractPropertyTweakHypothesis::asTweakOptID() {
    return new AbstractPropertyTweakHypothesisTweakOptID(
        AbstractPropertyTweakHypothesisConstSP::attachToRef(this), NULL);
}

void AbstractPropertyTweakHypothesis::load(CClassSP& clazz) {
    REGISTER(AbstractPropertyTweakHypothesis, clazz);
    SUPERCLASS(AtomicHypothesis);
}

CClassConstSP const AbstractPropertyTweakHypothesis::TYPE = CClass::registerClassLoadMethod(
    "AbstractPropertyTweakHypothesis", typeid(AbstractPropertyTweakHypothesis), load);

DRLIB_END_NAMESPACE
