/**
 * @file IHypothesis.hpp
 */

#ifndef DRLIB_IHypothesis_H
#define DRLIB_IHypothesis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/ScenarioShift.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(AtomicHypothesis)
#ifndef QLIB_IHYPOTHESIS_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<IHypothesis>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<IHypothesis>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<IHypothesis>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<IHypothesis>);
#endif


/**
 * A function which maps a "world" to an alternate world
 *
 * When we apply a tweak to the TweakGroup (instrument/model assembly), we are
 * conceptually taking a "world" and constructing an "alternate world".  This
 * interface captures that idea.  For an overview of the "declarative"
 * sensitivities framework of which it's a key part, see IRiskQuantityFactory.
 *
 * The only subtleties are:
 *
 *    -  Implementations are free to construct the tweaked object either by
 *       returning an independent, altered copy of the object provided, or
 *       by altering it in place and providing an "undo" method.  It's up to
 *       the client* to call AlternateWorld::undo() on the AlternateWorld
 *       returned from appliedTo().  [* In fact HypothesisTree::evaluate()
 *       in IRiskQuantityFactory.cpp.]
 *
 *    -  There is a notion of "distance" between worlds---essentially the
 *       divisor to use in computing finite difference derivative estimates,
 *       (i.e. epsilon in dp/dx (f(x') - f(x)) / epsilon). The distance
 *       between the original and alternate world is available as
 *       AlternateWorld::distance.
 *
 *    -  IHypothesis's can be grouped into a CompoundHypothesis: these arise
 *       when computing conditional and higher-order sensitivities, e.g.
 *       a "delta next day" tweak is a time tweak followed by a spot tweak,
 *       and a "d-delta-v" tweak is a spot tweak followed by a vol tweak.
 *       See also numAtomics() and atomic().
 *
 *
 * <H3>Uses</H3>
 *
 * IHypothesis is used to specify the condition part of a HypotheticalQuantity,
 * and thereby plays a key role in defining RiskQuantity's.  See
 * IRiskQuantityFactory for the big picture.
 *
 * The chief implementation of IHypothesis is
 * PropertyTweakHypothesis<PROPERTY>: all tweaks to RiskProperty's of market
 * data names (like those used to calculate Delta, VegaPointwise, ...)  are
 * instances of that.  See also FieldTweakHypothesis for a "dynamically-typed"
 * version.
 */

class RISKMGR_DLL IHypothesis: public virtual IObject,
                               public virtual IScenarioShift {

public:

    static CClassConstSP const TYPE;

private:

    static void load(CClassSP& clazz);

public:

    IHypothesis();

    virtual ~IHypothesis();

    static IHypothesisConstSP noop();

    FORWARD_DECLARE(AlternateWorld)

    /**
     * Altered version of object originally passed to appliedTo(), together
     * with distance between the two and method for undoing any in-place
     * modifications
     */

    class RISKMGR_DLL AlternateWorld: public CObject {

    public:

        static CClassConstSP const TYPE;

    private:

        AlternateWorld();
        static void load(CClassSP& clazz);

        bool undone;

        // temporary support for Sensitivity compatibility
        double oldValue; // $unregistered
        friend class HypothesisTree;
        friend class AbstractPropertyTweakHypothesis_AltWorld;

    protected:

        /**
         * Underlying implemention of undo().  Chief implementation is
         * AbstractPropertyTweakHypothesis_AltWorld::_undo() in
         * AbstractPropertyTweakHypothesis.cpp.  Default implementation does
         * nothing, so assumes that world was constructed by non-destructive
         * copying.
         */

        virtual void _undo();

        /**
         * Constructor with specific @a type, for subclasses.
         */

        AlternateWorld(CClassConstSP type, IObjectSP world, double distance,
                       bool found);

        /**
         * (For reflection pre-construction only)
         */

        AlternateWorld(CClassConstSP type);

    public:

        /**
         * Base type constructor.
         *
         * See below for world, distance.
         */

        AlternateWorld(IObjectSP world, double distance, bool found);
        ~AlternateWorld();

        /**
         * Alternate version of world originally passed to
         * IHypothesis::appliedTo()
         */

        IObjectSP world;

        /**
         * Distance between world and that originally passed to
         * IHypothesis::appliedTo()
         *
         * This is essentially the divisor to be used in (f(x') - f(x)) /
         * epsilon (before any rescaling to account to e.g. express it
         * it terms of basis points).  Most often it comes directly
         * from TweakOutcome::distance() after an invocation of
         * ITweakableWithRespectTo<PROPERTY>::sensShift().
         */

        double distance;

        /**
         * If distance is 0, that may be because no objects in the world fell
         * into the domain of the hypothesis, or because the "tweaks" it
         * applied made no difference for some reason: this flag indicates the
         * latter case.
         *
         * We need to know, in order to report Untweakable or NotApplicable ...
         */

        bool found;

        /**
         * Undo any in-place modification made to world
         *
         * Before the client calls this, the world originally passed
         * to IHypothesis::appliedTo() may be in an invalid state.
         * After the client calls this, AlternateWorld::world may be
         * in an invalid state.
         *
         * The underlying implementation, which subclasses should
         * override, is _undo().
         */

        void undo();
    };

    /**
     * A tweaked (perhaps in-place modified) version of an object
     *
     * The altered world is returned as AlternateWorld::world, along with its
     * distance from the supplied @a world, and a method
     * (AlternateWorld::undo()) to undo any in-place modification.  Before you
     * call that, @a world may be invalid.  After you call it, the
     * returned AlternateWorld::world may be invalid.
     */

    virtual AlternateWorldSP appliedTo(IObjectSP world) const = 0;

    /**
     * Modify an object in place
     *
     * @return  The distance between the new state of the world and the old
     */

    virtual double applyTo(IObjectSP world, bool* changed = 0) const = 0;

    /**
     * Number of atomic hypotheses into which this decomposes
     *
     * CompoundHypothesis::numAtomics() returns the number of hypotheses in the
     * underlying array; all other implementations derive from
     * AtomicHypothesis which just returns 1.
     */

    virtual int numAtomics() const = 0;

    /**
     * Breaks hypothesis into numbered atomic components
     *
     * CompoundHypothesis::atomic() returns the corresponding entry in
     * the underlying array; all other implementations derive from
     * AtomicHypothesis which just returns itself having checked that
     * i == 0.
     */

    virtual AtomicHypothesisConstSP atomic(int i) const = 0;

    /**
     * The "risk axis" along which this hypothesis shifts the world, if
     * applicable
     *
     * May return a null pointer.  If not,
     * axis()->hypothesis(axisCoefficient()) is supposed to be equivalent to
     * this hypothesis --- see IRiskAxis::hypothesis().
     */

    virtual IRiskAxisConstSP axis() const = 0;

    /**
     * The coefficient governing the amount by which this hypothesis shifts the
     * world along axis()
     *
     * This is not necessary equal to the AlternateWorld::distance reported
     * after the hypothesis is applied: axisCoefficient() may be relative where
     * AlternateWorld::distance is absolute, the actual shift performed may be
     * adaptively determined, etc.
     */

    virtual double axisCoefficient() const = 0;

    FORWARD_DECLARE(IDistanceMetric)

    class RISKMGR_DLL IDistanceMetric: public virtual IObject {
    public:

        static CClassConstSP const TYPE;

        IDistanceMetric();

        ~IDistanceMetric();

        static IDistanceMetricConstSP last();
        static IDistanceMetricConstSP constant(double distance);

        virtual double operator ()(const DoubleArray&) const = 0;
    };

    /**
     * How the AlternateWorld::distance's reported by the hypothesis's atomic
     * components are to be combined into an overall distance
     *
     * For AtomicHypothesis's this is a noop.  For CompoundHypothesis's you
     * get to specify it (see CompoundHypothesis::CompoundHypothesis()).
     */

    virtual IDistanceMetricConstSP distanceMetric() const = 0;

    /**
     * A compound hypothesis, implemented by applying this one first followed
     * by @a other
     *
     * The AlternateWorld::distance reported by the resulting
     * CompoundHypothesis is just that reported by @a other.
     */

    IHypothesisConstSP then(IHypothesisConstSP other) const;

    /**
     * IScenarioShift implementation
     *
     * Calls applyTo() and returns true if the object was "shifted" by a
     * non-zero distance
     */

    bool applyScenario(IObjectSP object);
    bool preapplyScenario(IObjectSP object);

	/**
     * Second order approximation for risk mapping
     */

    virtual void setApproxOrder(int order) = 0;
};

DRLIB_END_NAMESPACE

#endif
