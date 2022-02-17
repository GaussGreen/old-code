/**
 * @file ITweakableWithRespectTo.hpp
 */

#ifndef DRLIB_ITweakableWithRespectTo_H
#define DRLIB_ITweakableWithRespectTo_H

#include "edginc/PropertyTweak.hpp"
#include "edginc/TweakOutcome.hpp"
#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface designating market objects which have a particular "risk property".
 *
 * When a market object implements this interface, it puts itself in the domain
 * of PropertyTweakHypothesis<TAG>, i.e. has a different state in
 * alternate worlds constructed under that hypothesis.
 *
 * The object itself is responsible for implementing the actual state change
 * via the sensShift() method; PropertyTweakHypothesis is responsible for
 * finding objects the objects with the right sensName().
 *
 * For scalar properties like Spot (whose Qualifier type is specified as Void)
 * you don't have to implement sensQualifiers(), since that's done for you
 * in the specialised ITweakableWithRespectTo<Void> below.
 *
 * For an overview of how these classes fit into the "declarative"
 * sensitivities framework, see the IRiskQuantityFactory class documentation
 * (particularly <I>The new setup: property tags</I>).
 */

template <class TAG, class QUALIFIER = typename TAG::Qualifier>
class ITweakableWithRespectTo: public virtual IObject {

    static void load(CClassSP &clazz) {
        REGISTER_INTERFACE(ITweakableWithRespectTo, clazz);
        EXTENDS(IObject);
    }
public:

    static CClassConstSP const TYPE;

    /**
     * The market data name of this object (as far as the property is
     * concerned).
     *
     * Called by AbstractPropertyTweakHypothesis::appliedTo(), via
     * PropertyTweakHypothesis<TAG>::sensName(), to identify the objects
     * it wants to sensShift().
     *
     * Most objects, like Equity, have a single market name with which they're
     * associated, and just return that; a few return the name of an
     * appropriate sub-object which depends on the TAG tag.
     */

    virtual string sensName(const TAG*) const = 0;

    typedef typename TAG::Qualifier Qualifier;
    DECLARE(Qualifier)

    /**
     * The "qualifiers" defining the instances of the property on this
     * object.
     *
     * Some properties, like Spot, are scalar for each name in the market;
     * others exist in multiple instances for each name, and need
     * qualification.  The primary use is for term structured ("vector")
     * properties like VolPointwise, but others can be imagined.  See
     * axisFor().
     *
     * For scalar properties (e.g. Spot) you don't have to implement this
     * method: you get a trivial implemention for free (see the specialised
     * ITweakableWithRespectTo<Void> below).
     *
     * For term structured properties (e.g. VolPointwise) you must
     * return an ExpiryArrayConstSP giving the expiries at which vol is defined
     * for @a name.
     *
     * Called by PerNameRiskPropertySensitivity::nameRiskQuantities(), via
     * RiskProperty<TAG>::subjectQualifiers(), to get the list of
     * qualifiers to tweak on each object.
     */

    virtual QualifierArrayConstSP sensQualifiers(const TAG*) const = 0;

    /**
     * Alter the property on this object.
     *
     * This is the key method --- it's where you implement what a change to the
     * property does to the object's internal state, i.e. what the property
     * actually means.
     *
     * See for instance Equity::sensShift(const PropertyTweak<Spot> &).
     *
     * You get three arguments effectively:
     *
     *    -  TAG tells you which property to adjust
     *
     *    -  tweak.coefficient tells you how much to change it
     *       (see PropertyTweak<TAG>::coefficient)
     *
     *    -  tweak.qualifier tells you any other info you need to
     *       resolve TAG to a one-dimensional "risk axis".
     *       For term-structured properties (like VolPointwise),
     *       it's an ExpiryWindowConstSP; for scalar properties
     *       like Spot, it's Void and you can ignore it.
     *
     *    -  [tweak.marketDataName is generally redundant because it's
     *       your name by definition; it's used by a few basket-like
     *       objects which deliberately return names of their sub-components
     *       from sensName()
     *
     * If you are implementing IRestorableWithRespectTo, you need to save
     * enough of the object's previous state to undo the tweak:
     * see IRestorableWithRespectTo<TAG>::sensRestore().
     *
     * You return basically two things:
     *
     *    -  How far you actually moved the property value
     *       (@a distance param to TweakOutcome::TweakOutcome()) ---
     *       notionally the "distance" between the original and tweaked
     *       states of the world).  It's the divisor to be used in
     *       dp/dx ~= (p(x') - p(x)) / divisor, before any rescaling to
     *       e.g. express it in terms of basis points.  See also
     *       IHypothesis::AlternateWorld::distance.
     *
     *       Note that @a distance doesn't have to be equal to the
     *       PropertyTweak<TAG>::coefficient parameter described above.
     *       For instance Equity::sensShift(const PropertyTweak<Spot> &)
     *       treats 'coefficient' as a proportional change and returns
     *       the absolute change.
     *
     *    -  Whether you want the framework to tweak your sub-components
     *       (@a tweakMembers param to TweakOutcome::TweakOutcome()).
     *       Generally you do, but see
     *       XCB::sensShift(const PropertyTweak<BasketSpot> &) for an
     *       example of a basket-like object which doesn't.
     *
     *    -  [The three-argument
     *       TweakOutcome::TweakOutcome(double, double, bool), in which
     *       you specify oldValue and newValue rather than just distance,
     *       is for interoperability---see TweakOutcome::oldValue().]
     *
     * This method is called by AbstractPropertyTweakHypothesis::appliedTo(),
     * via PropertyTweakHypothesis<TAG>::sensShift().
     */

    virtual TweakOutcome sensShift(const PropertyTweak<TAG>& tweak) = 0;

    virtual ~ITweakableWithRespectTo() {}
private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG, class QUALIFIER> CClassConstSP const 
ITweakableWithRespectTo<TAG, QUALIFIER>::TYPE =
CClass::templateRegisterClass(
    typeid(ITweakableWithRespectTo<TAG, QUALIFIER>));
#endif

/**
 * Interface designating market objects which have a particular "risk property"
 * not requiring a qualifier
 *
 * This specialisation of ITweakableWithRespectTo allows market objects to omit
 * ITweakableWithRespectTo::sensQualifiers() if the qualifier is just "void",
 * i.e. the property is "scalar".
 */

template <class TAG>
class ITweakableWithRespectTo<TAG, Void>: public virtual IObject {
    static void load(CClassSP &clazz) {
        REGISTER_INTERFACE(ITweakableWithRespectTo, clazz);
        EXTENDS(IObject);
    }

public:

    static CClassConstSP const TYPE;

    virtual string sensName(const TAG*) const = 0;

    typedef Void Qualifier;
    DECLARE(Qualifier)

    virtual VoidArrayConstSP sensQualifiers(const TAG*) const {
        return VoidArrayConstSP();
    }

    virtual TweakOutcome sensShift(const PropertyTweak<TAG>& tweak) = 0;

    virtual ~ITweakableWithRespectTo() {}
private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const 
ITweakableWithRespectTo<TAG, Void>::TYPE =
CClass::templateRegisterClass(typeid(ITweakableWithRespectTo<TAG, Void>));
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
