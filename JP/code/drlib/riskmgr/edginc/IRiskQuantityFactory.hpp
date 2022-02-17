/**
 * @file IRiskQuantityFactory.hpp
 */

#ifndef DRLIB_IRiskQuantityFactory_H
#define DRLIB_IRiskQuantityFactory_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/Control_forward.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MultiTweakGroup)
FORWARD_DECLARE(NamedRiskQuantity)
FORWARD_DECLARE(LazyRiskQuantityFactory)
FORWARD_DECLARE(IRiskQuantityFactory)
FORWARD_DECLARE(RiskMapping)

/**
 * Top level interface for "declarative" risk management implementation
 *
 * There's a presentation under "Declarative greeks" in the QLib doc db:
 * <A HREF=Notes://PPUSMC006/8525703B0051D832/55B689E36191AD7285256DFD00470F6B/900A04B2259FB097852570BB006C9386>here</A>.
 * 
 * 
 *
 * Table of contents:
 *
 *    -  <A HREF=#introduction>Introduction</A>
 *    -  <A HREF=#status>Current status</A>
 *    -  <A HREF=#using>Using the framework</A>, for the impatient
 *    -  <A HREF=#design>Design overview</A>
 *    -  <A HREF=#properties>"Risk properties" for your typical Greeks</A>
 *
 *
 * <A NAME=introduction><H3>Introduction</H3></A>
 *
 * This is the entry point to a new mechanism for calculating sensitivities,
 * which was implemented in order to support features required for RiskMapping,
 * such as efficient calculation of greeks for multiple instruments
 * simultaneously, but also has some general advantages over the current
 * scheme:
 *
 *    -  Doesn't require cutting and pasting lots of "boilerplate"
 *
 *    -  Cross-cutting concerns---like avoiding duplicate re-pricings, and
 *       reporting NotApplicable and Untweakable correctly---are handled
 *       cleanly by the framework, not by the individual sensitivity
 *       implementations.  This makes cross-gammas and other "complex"
 *       cases far simpler.
 *
 *    -  The multiple roles played by the current Sensitivity class and its
 *       descendants are split off into three simpler and clearer interfaces
 *
 * The key point about this design is that it's "declarative": it's centred
 * around a representation of what states the world needs to be tweaked into
 * and how the prices in each state are to be combined into a derivative ---
 * each of the different sensitivities works by describing what it wants to
 * compute, rather than by actually being a C++ implementation of the
 * computation.  That's what makes it possible to implement cross-cutting
 * concerns generically, and handle new requirements like efficient
 * multiple-instrument greeks, risk mapping, ...
 *
 * Other benefits which may be of interest:
 *
 *    -  It's easy to compute derivatives of things other than overall
 *       instrument FV
 *
 *    -  Coping with different "versions" of the same quantity (like
 *       risky/non-risky FV) is straightforward
 *
 *    -  The declarative representation would be amenable to parallisation
 *       if that were useful (I'm guessing it probably isn't)
 *
 *
 * <A NAME=status><H3>Current status</H3></A>
 *
 * It's intended that the new framework will gradually take over most (in the
 * long run perhaps all) of the jobs of the old one.  At the moment:
 *
 *    -  Some core Greeks for which we need RiskMapping have been "ported"
 *
 *          -  VegaParallel, VegaPointwise are fully in the new framework
 *
 *          -  Delta and BasketDelta are still in the old class hierarchy
 *             so that fiddly mechanisms which rely on the ScalarShift
 *             interface can be addressed at leisure --- but they are
 *             implemented using the new framework
 *
 *          -  Theta will follow
 *
 *    -  An assortment of credit Greeks have been ported as well
 *
 *          -  CreditDeltaPointwiseWithShift, LegalBasisAdditivePointwise,
 *             LegalBasisMultiplierPointwise, LiquiditySpreadRhoPointwise,
 *             ParSpreadRhoParallel, ParSpreadRhoPointwise,
 *             ParSpreadRhoPointwiseTwoSided
 *
 *    -  Logic in products/models/market data for doing different things
 *       depending on "the current sensitivity" has been left as far as
 *       possible the same, and is supported using a compatibility mechanism.
 *       The next step will be to leverage the new scheme to make these
 *       things simpler and more general---but, one thing at a time.
 *
 *
 * <A NAME=using><H3>Using the framework</H3></A>
 *
 * As well as facilitating complex greeks and facilities like risk mapping, the
 * new framework makes implementing common greeks even easier.  To implement a
 * typical sensitivity, you ...
 *
 *    -  Define a tag class for the "risk property" with respect to which
 *       the sensitivity is measured: let's call it PROPERTY
 *
 *          -  Examples: ParSpreadParallel.hpp / ParSpreadParallel.cpp,
 *             VegaPointwise.hpp / VegaPointwise.cpp
 *
 *    -  Make your tweakable market objects inherit from
 *       ITweakableWithRespectTo<PROPERTY> or
 *       IRestorableWithRespectTo<PROPERTY>, and implement the usual
 *       sensName() / sensShift() methods
 *
 *          -  Examples: Equity.hpp / Equity.cpp
 *
 *    -  Write a factory class on top of
 *       AllNamesRiskPropertySensitivity ("the new SensControlAllNames"),
 *       ScalarRiskPropertySensitivity ("the new ScalarShift") or
 *       PerNameRiskPropertySensitivity<ExpiryWindow> ("the new VectorShift"),
 *       or one of the other PerNameRiskPropertySensitivity variants,
 *       initialising it with RiskProperty<PROPERTY>
 *
 *          -  Examples: ParSpreadRhoParallel.cpp,
 *             VegaPointwise.hpp / VegaParallel.cpp
 *
 * The only slight subtlety arises from the requirement that <I>all</I> risk
 * properties must have a qualifier, even "scalar" ones: see
 * <A HREF=#qualifiers>Qualifiers and pointwise properties</A>.
 * When writing a scalar greek you have to
 *
 *    -  put <TT>typedef Void Qualifier;</TT> in your tag class, and
 *
 *    -  inherit your market objects from
 *       ITweakableWithRespectTo<PROPERTY> rather than
 *       ITweakableWithRespectTo<PROPERTY> (to get a default, empty
 *       implementation of ITweakableWithRespectTo<PROPERTY>::sensQualifiers()).
 *
 * When we move to Visual Studio 7.1 the latter wrinkle can go away, as
 * originally envisaged.
 *
 *
 * <A NAME=design><H3>Design overview</H3></A>
 *
 * Going right back to basics, the canonical Greek calculation is conceptually
 *
 * -     dValue/dX = (Value of instrument -
 *                    Value of instrument in a world in which X is different )
 *                         / amount by which X is different
 *
 * In the jargon of the new framework, this is a "risk quantity".  It's a
 * number computed by performing simple arithmetic (i.e. (a - b) / e) on two
 * "hypothetical quantities", namely, the instrument value in the real world,
 * and the instrument value under the "hypothesis" that X differs from its real
 * value by some tweak.
 *
 * A top-level sensitivity is understood as corresponding to a number of risk
 * quantities---for instance, VegaPointwise yields one for each name and
 * expiry, and Delta gives you two for each name (delta and gamma).
 *
 * In terms of C++ classes:
 *
 *    -  The "world" is a MultiTweakGroup, comprising a model and some
 *       instruments which it can price (and by inclusion all relevant
 *       market data)
 *
 *    -  A top level sensitivity is an IRiskQuantityFactory: given the
 *       world, its riskQuantities() method returns a list of
 *       NamedRiskQuantity's to be computed
 *
 *    -  A NamedRiskQuantity is a RiskQuantity together with an
 *       IResultsIdentifier under which it will eventually appear in the
 *       Results (e.g. "VEGA_POINTWISE IBM 3M")
 *
 *    -  A RiskQuantity comprises a list of HypotheticalQuantity's,
 *       and a function for combining them once they've been calculated
 *       (e.g. onesided(a, b, e) = (a - b) / e)
 *
 *    -  A HypotheticalQuantity is a function, together with an
 *       IHypothesis defining an altered state of the world, in which
 *       the function is to be evaluated [in practice the function is
 *       ResultsFunction::price(), i.e. "value of instrument", but
 *       it could be anything]
 *
 *    -  An IHypothesis is essentially a function taking a world
 *       and yielding an adjusted world (e.g. tweaking spot for
 *       some name), although for efficiency it supports undoable
 *       in-place updates: see IHypothesis::AlternateWorld::undo()
 *
 * The code which actually calculates the greeks is
 * RiskQuantityEvaluator::storeResults(): it
 *
 *    -  enumerates the NamedRiskQuantity's requested by the greeks
 *       (IRiskQuantityFactory's)
 *
 *    -  collates the IHypothesis's which need to be applied to the
 *       base-case world in order to compute those quantities;
 *
 *    -  runs through each hypothesis in turn, tweaking the world and
 *       recording the instrument prices which are required from that
 *       world-state;
 *
 *    -  finally, passing the prices thus gathered to the NamedRiskQuantity's
 *       so that they can compute the derivatives, and put them in the
 *       Results.
 *
 * You can get an "educational execution trace" showing the sequence of
 * events by e.g.
 *
 * <pre>
 *    cd c:/qlib-trunk/testing/riskmappinginp
 *    set QLIB_TRACE=c:/temp/trace.txt    (  export QLIB_TRACE=c:/temp/trace.txt  )
 *    ../../code/drlib/regtest/lib-nt.vc71/models.exe RiskMgr-VolSVJ-cfcomparison.xml
 * </pre>
 *
 * and opening c:/temp/trace.txt in an editor.
 *
 * Each world state is visited only once, and the instruments are priced once
 * in each state, so there is no need to worry about cacheing (and avoiding
 * returning cached but unwanted results, or reusing results computed using
 * different shift sizes or assumptions).
 *
 *
 * <H4>Atomic and compound hypotheses</H4>
 *
 * Hypotheses like "IBM spot up 1%" and "FT 3M vol down 2%" can be combined
 * into a CompoundHypothesis: these arise in practice from cross-gammas, and
 * from greeks computed under a conditioning scenario, such as delta-next-day.
 *
 * The set of world states that the RiskQuantityEvaluator needs visit then
 * actually becomes a tree, with AtomicHypothesis's as links---traversing the
 * tree is an easy way of visiting all world states while expending a minimum
 * amount of effort on tweaking.
 *
 *
 * <A NAME=properties><H3>"Risk properties" for your typical Greeks</H3></A>
 *
 * Most of the time, you're interested in sensitivities to what could be called
 * "risk properties" of market names---spot, vol, etc.---and the hypotheses you
 * want to try out are tweaks a property of a specific name, qualified if
 * necessary by an expiry (denoting a point in the term structure of the
 * property).  For instance, IBM spot or Exxon 3M par spread.
 *
 * <H4>The current setup</H4>
 *
 * In the current setup, these common cases are covered by ScalarShift and
 * VectorShift.  You define methods for identifying objects which carry the
 * property and for adjusting it in each case---conventionally called
 * sensName() and sensShift().  You group the methods into an interface.  And
 * you write a class inheriting from ScalarShift or VectorShift containing that
 * interface definition and a number of boilerplate methods which essentially
 * map generic "shift(IObject)" calls onto invocations of specific sensShift()
 * methods, etc.
 *
 * <A NAME=tags><H4>The new setup: property tags</H4></A>
 *
 * In the new setup, you define a tiny "tag" class to name the "risk property"
 * in question---for instance, Spot and VolPointwise.  These don't generally
 * include any actual method or fields of their own---although they can, see
 * Correl for an example---but by the magic of templates, you then get the
 * following for free:
 *
 *    -  RiskProperty<PROPERTY> represents the property in general,
 *       e.g. RiskProperty<VolPointwise> represents "pointwise volatility"
 *       
 *    -  PropertyRiskAxis<PROPERTY> represents a property of a particular
 *       name, qualified so as to be one-dimensional: for instance,
 *       PropertyRiskAxis<VolPointwise>("IBM", 6M) represents IBM
 *       6M vol
 *
 *    -  PropertyTweakHypothesis<PROPERTY> represents a shift along a
 *       PropertyRiskAxis: for instance
 *       PropertyTweakHypothesis<VolPointwise>("IBM", 6M, 0.01)
 *       represents a 1% tweak to IBM 6M vol
 *       
 * (RiskProperty<PROPERTY>::axisFor() is a factory method for
 * PropertyRiskAxis<PROPERTY>; and PropertyRiskAxis<PROPERTY>::hypothesis()
 * is a factory method for PropertyTweakHypothesis<PROPERTY>.)
 *
 *    -  ITweakableWithRespectTo<PROPERTY> is the interface identifying
 *       objects which have the property in question---it defines the
 *       sensName() and sensShift() methods which
 *       PropertyTweakHypothesis<PROPERTY> uses to find and perturb
 *       its subjects
 *
 * <H4>Abstract interfaces</H4>
 *
 * Underneath these are abstract interfaces IRiskProperty, IRiskAxis, and
 * AbstractPropertyTweakHypothesis (a refinement of AtomicHypothesis,
 * discussed above).  Because the code is written to these interfaces rather
 * than using the templates directly, the templates only have to be '#included'
 * in a couple of <TT>.cpp</TT>s---for instance, <TT>VolPointwise.cpp</TT> and
 * <TT>VegaPointwise.cpp</TT>---and you don't have to use them if they're not
 * appropriate for your particular application.
 *
 * <A NAME=qualifiers><H4>Qualifiers and pointwise properties</H4></A>
 *
 * The "qualifier" of a risk property is the extra info you have to add in
 * order to make it one-dimensional for a given market object.  For
 * term-structured properties like VolPointwise---those we conventionally call
 * "pointwise" or "vector" properties---this is essentially an Expiry.
 *
 * However, in some contexts we mean to tweak not literally a single expiry of
 * a property, but rather all expiries within a window around the "current"
 * point: back to the previous one and forward to the next.  So strictly
 * speaking the qualifier for a pointwise property needs to be an ExpiryWindow.
 * 
 * To handle pointwise properties uniformly with properties that don't require
 * qualification (like Spot and ParSpreadParallel), and potentially also
 * properties qualified by types other than ExpiryWindow, we require that
 * <I>all</I> <A HREF=#tags>property tag</A> classes include a
 * <TT>typedef</TT>:
 *
 *    -  Qualifier is either ExpiryWindow for pointwise properties,
 *       or Void for "scalar" ones [capital V: it's an empty class
 *       defined in toolkit/edginc/Void.hpp]
 *
 * Market data objects implementing ITweakableWithRespectTo<PROPERTY> advertise
 * the qualifiers they have available for tweaking through the method
 * ITweakableWithRespectTo<PROPERTY>::sensQualifiers().  The qualifier being
 * tweaked by a particular call to
 * ITweakableWithRespectTo<PROPERTY>::sensShift(shift) is available to the
 * object as <TT>shift.qualifier</TT>.
 *
 * See Equity::sensQualifiers(Spot*) and Equity::sensShift(const
 * PropertyTweak<Spot>&) for "pointwise" examples.  For "scalar" properties,
 * where Qualifier is Void, you get a default sensQualifiers implementation
 * for free by the magic of partial template specialisation, so you only
 * have to implement sensShift and sensName.
 */

class RISKMGR_DLL IRiskQuantityFactory: public virtual IObject {

public:

    static CClassConstSP const TYPE;

    IRiskQuantityFactory();
    virtual ~IRiskQuantityFactory();

    /**
     * The RiskQuantity's which this IRiskQuantityFactory says should be
     * computed for a given base-state world.
     *
     * This is called from RiskQuantityEvaluator::storeResults().
     * See the <A HREF=#design>design overview</A> above.
     *
     * For implementations, see
     * RiskQuantityFactorySensitivity::riskQuantities() and behind that
     * PerNameRiskPropertySensitivity::nameRiskQuantities(),
     * AllNamesRiskPropertySensitivity::nameRiskQuantities().
     *
     * For example, VegaPointwise might look at @a world, find names it knows
     * how to tweak, interrogate them for their vol expiries, and return risk
     * quantities representing "first derivative of price with respect to IBM
     * 3M vol", "... IBM 6M vol", "... JPM 1M vol".
     */

    virtual NamedRiskQuantityArraySP riskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const = 0;

    /**
     * Mechanism for handling the case in which the set of desired
     * RiskQuantity's can't be determined from the base-state world.
     *
     * All current and nearly all future implementations of
     * IRiskQuantityFactory need only provide riskQuantities(), with lazies()
     * just being a stub returning NULL (like
     * RiskQuantityFactorySensitivity::lazies()).
     *
     * What's it for?  Conceivably you might want to put the world into a
     * different state (perhaps do a theta tweak) and only then look at it
     * to decide exactly which RiskQuantity's you want to evaluate.  This
     * method allows you to return a list of such tweaks to apply, paired
     * with IRiskQuantityFactory's to invoke once they have been.
     *
     * Called from RiskQuantityEvaluator::storeResults().
     */

    virtual LazyRiskQuantityFactoryArraySP lazies(
        MultiTweakGroupConstSP world) const = 0;

    /**
     * A degenerate IRiskQuantityFactory which just returns the RiskQuantity's
     * you constructed it from
     */

    static IRiskQuantityFactorySP fromRiskQuantities(
                                      NamedRiskQuantityArraySP rqs);
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
