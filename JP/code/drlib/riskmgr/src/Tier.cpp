/**
 * @file Tier.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/IEqVolNamePair.hpp"
#include "edginc/SimpleTweakNameResolver.hpp"
#include "edginc/CompoundHypothesis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TierProperty)

/**
 * Specialised IRiskProperty for supporting "Tier" sensitivity
 *
 * The reason this is nonstandard is that Tier relies on tying together each
 * asset in the world with its volatility: see TierHypothesis for an account
 * of how we do this (via IEqVolNamePair).
 */

class TierProperty: public CObject,
                    public virtual IRiskProperty<Void> {

public:

    static CClassConstSP const TYPE;

private:

    TierProperty(): CObject(TYPE) {}
    static IObject* emptyShell();
    static void load(CClassSP&);

public:

    IScalarRiskPropertyConstSP spotProperty; // $required
    IScalarRiskPropertyConstSP volProperty;  // $required

    TierProperty(IScalarRiskPropertyConstSP spotProperty,
                 IScalarRiskPropertyConstSP volProperty):
        CObject(TYPE),
        spotProperty(spotProperty),
        volProperty(volProperty)
    {}

    bool discrete() const {
        return spotProperty->discrete() || volProperty->discrete();
    }

    CClassConstSP subjectInterface() const {
        return IEqVolNamePair::TYPE;
    }

    OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        map<string, OutputNameSP> pairs = IEqVolNamePair::namePairs(world);

        OutputNameArraySP names(new OutputNameArray(pairs.size()));

        int n = 0;
        for (map<string, OutputNameSP>::iterator p = pairs.begin();
             p != pairs.end(); ++p, ++n) {
            (*names)[n] = OutputName::SP(p->first);
        }

        return names;
    }

    VoidArrayConstSP subjectQualifiers(IObjectConstSP, OutputNameConstSP) const {
        return VoidArrayConstSP();
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP subjectName,
                             VoidConstSP) const;
};

FORWARD_DECLARE(TierHypothesis)

/**
 * Specialised IHypothesis and IRiskAxis for supporting Tier sensitivity
 *
 * The reason this is nonstandard is that Tier relies on tying together each
 * asset in the world with its volatility, which can't be done in a uniform
 * way.  So we reuse a mechanism originally put together for use in DDeltaDVol:
 * we rely on the "asset holder" objects like SimpleEquity to implement
 * IEqVolNamePair so that they can tell us
 *
 *    -  the names of the ITweakableWithRespectTo<Spot> objects representing
 *       the "asset" objects like Equity
 *
 *    -  the names of the corresponding vol objects
 *
 * See hypFor() for the guts of the implementation.
 *
 * Note that because the only permissable axisCoefficient() (shift size) is
 * one, we can collapse the IRiskAxis and IHypothesis stages into one class,
 * with hypothesis() just returning <tt>this</tt>.
 */

class TierHypothesis: public AbstractPropertyTweakHypothesis,
                      public virtual IRiskAxis,
                      public virtual RiskAxis {
public:

    static CClassConstSP const TYPE;

private:

    TierHypothesis(): AbstractPropertyTweakHypothesis(TYPE) {}
    static IObject* emptyShell();
    static void load(CClassSP&);

    TierPropertyConstSP property;  // $required
    OutputNameConstSP name;        // $optional

    IHypothesisSP hypFor(IObjectConstSP world) const {
        TRACE_METHOD;

        map<string, OutputNameSP> pairs = IEqVolNamePair::namePairs(world);

        IHypothesisArraySP hyps(new IHypothesisArray());

        map<string, OutputNameSP>::iterator p = pairs.find(name->idGet(0));
        if (p != pairs.end()) {
            hyps->push_back(IHypothesisSP::constCast(property->spotProperty->
                                axisFor(name)->hypothesis(1)));
            TRACE("Will apply " << *hyps->back() << " ...");
            hyps->push_back(IHypothesisSP::constCast(property->volProperty->
                                axisFor(p->second)->hypothesis(1)));
            TRACE("... and corresponding " << *hyps->back());
        }

        return CompoundHypothesis::SP(hyps,
                                      IHypothesis::IDistanceMetric::constant(1));
    }

public:

    virtual TweakOutcome sensShift(IObjectSP object) const {
        throw ModelException("Not Implemented");
    }
    virtual CClassConstSP restorableInterface() const {
        throw ModelException("Not Implemented");
    }
    virtual void sensRestore(IObjectSP object) const {
        throw ModelException("Not Implemented");
    }
    virtual OutputNameConstSP getMarketDataName() const {
        return name;
    }
    virtual IObjectConstSP getQualifier() const {
        throw ModelException("Not Implemented");
    }

    virtual string sensName(IObjectConstSP object) const {
        throw ModelException("Not Implemented");
    }
 
    virtual void setMarketDataName(OutputNameConstSP) {
        throw ModelException("Not Implemented");
    }

    virtual void setAxisCoefficient(double) {
        throw ModelException("Not Implemented");
    }
    virtual CClassConstSP shiftableInterface () const{
        throw ModelException("Not Implemented");
    }

    TierHypothesis(TierPropertyConstSP property,
                   OutputNameConstSP name):
        AbstractPropertyTweakHypothesis(TYPE),
        property(property),
        name(name)
    {}

    IAbstractRiskPropertyConstSP abstractProperty() const {
        return property;
    }

    OutputNameConstSP marketDataName() const {
        return name;
    }

    IRiskAxisConstSP thawed() const {
        return IRiskAxisConstSP::attachToRef(this);
    }

    RiskAxisConstSP frozen() const {
        return RiskAxisConstSP::attachToRef(this);
    }

    IHypothesisConstSP hypothesis(double coeff) const {
        if (coeff == 0.) {
            return IHypothesis::noop();
        }
        else if (coeff == 1.) {
            return IHypothesisConstSP::attachToRef(this);
        }
        else {
            throw ModelException(__FUNCTION__,
                "The only permissable shift sizes for a Tier sensitivity are "
                "0, 1");
        }
    }

    AlternateWorldSP appliedTo(IObjectSP world) const {
        TRACE_METHOD;
        try {
            return hypFor(world)->appliedTo(world);
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

    double applyTo(IObjectSP world, bool* changed) const {
        TRACE_METHOD;
        try {
            return hypFor(world)->applyTo(world, changed);
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

    IRiskAxisConstSP axis() const {
        return IRiskAxisConstSP::attachToRef(this);
    }

    double axisCoefficient() const {
        return 1.;
    }
};

IRiskAxisConstSP TierProperty::axisFor(OutputNameConstSP subjectName,
                                       VoidConstSP) const {
  return TierHypothesisConstSP(new TierHypothesis(
      TierPropertyConstSP::attachToRef(this),
      subjectName));
}

/**
 * A sensitivity which tells you the changes in your instrument's price when
 * each asset in the world has its spot tweaked (by a given SpotShift) and
 * simultaneously its vol tweaked (by a given VolAbsoluteShift)
 */

class Tier: public ScalarRiskPropertySensitivity,
            public virtual Additive {                        // $public

    static IObject* emptyShell();
    static void load(CClassSP&);

    ScalarRiskPropertySensitivity::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     TierPropertySP(new TierProperty(spotShift->asRiskProperty(),
                                                     volShift->asRiskProperty())),
                     IScalarDerivative::oneSided(),
                     1.);
    }

    string packetName;                 // $optional(Name under which to file this sensitivity in the results)
    SpotShiftConstSP spotShift;        // $required
    VolAbsoluteShiftConstSP volShift;  // $required

    Tier():
        ScalarRiskPropertySensitivity(TYPE, NAME, 1.)
    {}

public:

    virtual const string& getSensOutputName() const {
        return packetName == "" ? NAME : packetName;
    }

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const string NAME;
};

const string Tier::NAME("TIER");
const double Tier::DEFAULT_SHIFT = 1;

bool TierLoad() {
    return Tier::TYPE != NULL;
}


IObject* TierProperty::emptyShell() {
  return new TierProperty();
}

void TierProperty::load(CClassSP& clazz) {
  REGISTER(TierProperty, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskProperty<Void>);
  FIELD(spotProperty, "spotProperty");
  FIELD(volProperty, "volProperty");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const TierProperty::TYPE = CClass::registerClassLoadMethod(
  "TierProperty", typeid(TierProperty), TierProperty::load);

IObject* Tier::emptyShell() {
  return new Tier();
}

void Tier::load(CClassSP& clazz) {
  clazz->setPublic();
  REGISTER(Tier, clazz);
  SUPERCLASS(ScalarRiskPropertySensitivity);
  FIELD(packetName, "Name under which to file this sensitivity in the results");
  FIELD_MAKE_OPTIONAL(packetName);
  FIELD(spotShift, "spotShift");
  FIELD(volShift, "volShift");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const Tier::TYPE = CClass::registerClassLoadMethod(
  "Tier", typeid(Tier), Tier::load);

IObject* TierHypothesis::emptyShell() {
  return new TierHypothesis();
}

void TierHypothesis::load(CClassSP& clazz) {
  REGISTER(TierHypothesis, clazz);
  SUPERCLASS(AtomicHypothesis);
  IMPLEMENTS(IRiskAxis);
  IMPLEMENTS(RiskAxis);
  FIELD(property, "property");
  FIELD(name, "name");
  FIELD_MAKE_OPTIONAL(name);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const TierHypothesis::TYPE = CClass::registerClassLoadMethod(
  "TierHypothesis", typeid(TierHypothesis), TierHypothesis::load);

DRLIB_END_NAMESPACE
