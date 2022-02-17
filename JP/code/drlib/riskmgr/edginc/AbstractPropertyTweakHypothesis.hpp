/**
 * @file AbstractPropertyTweakHypothesis.hpp
 */

#ifndef QLIB_AbstractPropertyTweakHypothesis_H
#define QLIB_AbstractPropertyTweakHypothesis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/TweakOutcome.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(AbstractPropertyTweakHypothesis)
FORWARD_DECLARE(ICompatibilitySensitivity)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(OutputName)
class ITweakOptID;

/**
 * A function which, when applied to a "world", produces an alternate world in
 * which some property of a named market object has been modified.
 *
 * Most IHypothesis implementations involve tweaking "properties" (like spot)
 * of market objects; they're handled nice and uniformly via the
 * PropertyTweakHypothesis<PROPERTY> template, but to minimise the repetition
 * of templated code, we pull out the bits that don't need to know about
 * PROPERTY and stick them in this abstract base class.
 *
 * It also serves as a base class for FieldTweakHypothesis, which is
 * a "dynamically typed" counterpart to PropertyTweakHypothesis.
 *
 * For an overview of the "declarative" sensitivities framework of which these
 * classes are a part, see the IRiskQuantityFactory class documentation.
 */

class RISKMGR_DLL AbstractPropertyTweakHypothesis: public AtomicHypothesis {

    static void load(CClassSP& clazz);
    AbstractPropertyTweakHypothesis(const AbstractPropertyTweakHypothesis& rhs);
    AbstractPropertyTweakHypothesis& operator=(
        const AbstractPropertyTweakHypothesis& rhs);

public:

    static CClassConstSP const TYPE;

    ~AbstractPropertyTweakHypothesis();

    // Used to support product/engine logic based on Control::getCurrentSens(),
    // until we get round to re-expressing that in a more precise way.
    // See RiskQuantityEvaluator::storeResults().

    mutable ICompatibilitySensitivitySP _originatingSensitivity; // $unregistered

protected:

    friend class AbstractPropertyTweakHypothesisTweakOptID;
    friend class RiskQuantityFactorySensitivity;
    friend class ScalarShiftTwoSided;
    friend class HypothesisTree;

    AbstractPropertyTweakHypothesis(CClassConstSP type);

    /**
     * The interface which designates market objects affected by this hypothesis.
     *
     * It's actually ITweakableWithRespectTo<PROPERTY> (see
     * PropertyTweakHypothesis::shiftableInterface()).
     */

    virtual CClassConstSP shiftableInterface() const = 0;

    /**
     * Make a change to PROPERTY of a particular object.
     *
     * Calls ITweakableWithRespectTo<PROPERTY>::sensShift() on @a object.
     */

    virtual TweakOutcome sensShift(IObjectSP object) const = 0;

    /**
     * The interface which designates market objects which can be
     * mutated in place to apply this hypothesis.
     *
     * It's actually IRestorableWithRespectTo<PROPERTY> (see
     * PropertyTweakHypothesis::restorableInterface()).
     */

    virtual CClassConstSP restorableInterface() const = 0;

    /**
     * Undo the effect of a sensShift().
     *
     * Calls IRestorableWithRespectTo<PROPERTY>::sensRestore() on @a object.
     */

    virtual void sensRestore(IObjectSP object) const = 0;

    /**
     * The name of the market data object whose property is changed under this
     * hypothesis.
     *
     * Originates as @a subjectName argument to RiskProperty::axisFor().
     */

    virtual OutputNameConstSP getMarketDataName() const = 0;

    /**
     * Info needed in addition to getMarketDataName() to fully specify the
     * hypothesis.
     *
     * E.g. for pointwise tweaks (term-specific hypotheses) this is an
     * ExpiryWindow.  It's really of type PROPERTY::Qualifier.  You shouldn't
     * need to know this or what type this is; it's only used in toString().
     *
     * Originates as @a qualifier argument to RiskProperty::axisFor().
     */

    virtual IObjectConstSP getQualifier() const = 0;

public:

    /**
     * Description of this hypothesis (for use in error messages).
     */

    virtual string toString() const;

    /**
     * A version of @a world in which some property of a named market object
     * has been modified.
     *
     * This is the implementation of the key IHypothesis::appliedTo() method.
     * It uses shiftableInterface() and getMarketDataName() to find the object
     * which we want to tweak, then uses sensShift() to tweak it.
     */

    IHypothesis::AlternateWorldSP appliedTo(IObjectSP world) const;

    /**
     * Modifies @a world in place, altering a property of a named market
     * object.
     *
     * @return  The distance between the original and the modified worlds;
     *          see IHypothesis::AlternateWorld::distance.
     */

    virtual double applyTo(IObjectSP world, bool* changed = 0) const;

    /**
     * Modifies @a world in place, altering a property of a named market
     * object.
     *
     * As applyTo(), but provides TweakOutcome::oldValue(): necessary (for the
     * moment) for interoperating with SensControl and friends.
     */

    virtual TweakOutcome applyTo_TweakOutcome(IObjectSP world) const;

    /**
     * Methods used for implementing compatibility with
     * Control::getCurrentSensitivity()
     *
     * See RiskQuantityEvaluator::storeResults().  Eventually these will go away.
     */

    //@{

    virtual string sensName(IObjectConstSP object) const = 0; // would be nice not to need this

    virtual void setMarketDataName(OutputNameConstSP) = 0;

    virtual void setAxisCoefficient(double) = 0;

    // you own the returned pointer: this is for compatibility with VolBaseParamSurface
    ITweakOptID* asTweakOptID();

    //@}
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
