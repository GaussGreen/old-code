/**
 * @file CompoundHypothesis.hpp
 */

#ifndef DRLIB_CompoundHypothesis_H
#define DRLIB_CompoundHypothesis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(AtomicHypothesis)
FORWARD_DECLARE(CompoundHypothesis)

/**
 * A sequence of functions which can be applied to a "world" to produce an
 * alternate world.
 *
 * CompoundHypothesis is essentially an array of IHypothesis's which are
 * applied one after another to form an overall IHypothesis.  It's used when
 * computing conditional and higher-order sensitivities, e.g. a "delta next
 * day" tweak is a time tweak followed by a spot tweak, and a "d-delta-v" tweak
 * is a spot tweak followed by a vol tweak.  By breaking a CompoundHypothesis
 * into its atomic components, RiskQuantityEvaluator::storeResults() can optimise
 * the amount of tweaking it does.
 *
 * The only subtleties are:
 *
 *    -  The components of a CompoundHypothesis are AtomicHypothesis's,
 *       not (base) IHypothesis's.  This ensures that a CompoundHypothesis is
 *       always a one-level list, never a multi-level tree of compounds within
 *       compounds.
 *
 *    -  The distance between a world generated by a CompoundHypothesis
 *       and the original is defined to be a function of the distances between
 *       the intermediate worlds produced by the component atomic hypotheses:
 *       it's specified by the distanceMetric given at construct time.
 *
 * See IRiskQuantityFactory for the big picture of the subsystem in which
 * this class is used.
 */

class RISKMGR_DLL CompoundHypothesis: public CObject,
                                      public virtual IHypothesis {

    CompoundHypothesis();
    static IObject *defaultOne();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    AtomicHypothesisArrayConstSP differencesWrtBaseState;
    IHypothesis::IDistanceMetricConstSP _distanceMetric;

public:

    /**
     * Constructor
     *
     * @param differencesWrtBaseState   The underlying atomic hypotheses
     *
     * @param distanceMetric            Function for combining the
     *                                  TweakOutcome::distance()'s by which
     *                                  each constituent hypothesis "shifts"
     *                                  the world into an overall distance: see
     *                                  appliedTo().
     */

    CompoundHypothesis(AtomicHypothesisArrayConstSP differencesWrtBaseState,
                       IHypothesis::IDistanceMetricConstSP distanceMetric);

    /**
     * Constructor, with default distance metric
     *
     * The TweakOutcome::distance() reported by the resulting CompoundHypothesis
     * is just that reported by the last AtomicHypothesis.
     */

    CompoundHypothesis(AtomicHypothesisArrayConstSP differencesWrtBaseState);

    /**
     * smartPtr constructor
     */

    static CompoundHypothesisSP SP(
        AtomicHypothesisArrayConstSP differencesWrtBaseState,
        IHypothesis::IDistanceMetricConstSP distanceMetric);

    /**
     * smartPtr constructor, with default distance metric
     */

    static CompoundHypothesisSP SP(
        AtomicHypothesisArrayConstSP differencesWrtBaseState);

    /**
     * smartPtr constructor, from non-atomic hypotheses
     */

    static CompoundHypothesisSP SP(
        IHypothesisArrayConstSP differencesWrtBaseState,
        IHypothesis::IDistanceMetricConstSP distanceMetric);

    /**
     * An empty CompoundHypothesis (with zero component hypothesis,
     * it's a no-op).
     */

    static CompoundHypothesisSP empty();

    ~CompoundHypothesis();

    /**
     * Number of atomic hypotheses into which this decomposes.
     *
     * That's the length of the array passed to the CompoundHypothesis()
     * constructor.
     */

    int numAtomics() const;

    /**
     * Breaks hypotheses into numbered atomic components.
     *
     * That's just the components of the array passed to the CompoundHypothesis()
     * constructor.
     */

    AtomicHypothesisConstSP atomic(int) const;

    /**
     * A tweaked (perhaps in-place modified) version of an object.
     *
     * See IHypothesis::appliedTo().  In this implementation it's
     * the result of applying in turn each of the component hypotheses
     * supplied to the CompoundHypothesis() constructor.
     */

    IHypothesis::AlternateWorldSP appliedTo(IObjectSP world) const;

    /**
     * Modify an object in place.
     *
     * Each of the component hypotheses supplied to the CompoundHypothesis()
     * constructor is applied in turn to @a world.
     */

    double applyTo(IObjectSP world, bool* changed = 0) const;

    /**
     * The "risk axis" along which this hypothesis shifts the world, if
     * applicable.
     *
     * The default implementation returns a null pointer, since a
     * CompoundHypothesis in general isn't generated from a particular
     * IRiskAxis.  See IHypothesis::axis() for more info.
     */

    IRiskAxisConstSP axis() const;

    /**
     * The coefficient governing the amount by which this hypothesis shifts the
     * world along axis().
     *
     * The default implementation returns 0.  See axis().
     */

    double axisCoefficient() const;

    IHypothesis::IDistanceMetricConstSP distanceMetric() const;

    AtomicHypothesisSP grouped() const;

    /**
     * Used in making error messages
     */

    string toString() const;

	/**
     * Second order approximation for risk mapping
     */

    void setApproxOrder(int order);

};

DRLIB_END_NAMESPACE

#endif
