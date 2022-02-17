/**
 * @file IRiskAxis.hpp
 */

#ifndef DRLIB_IRiskAxis_H
#define DRLIB_IRiskAxis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IHypothesis.hpp" // just for 'metric' actually

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(RiskAxis)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(IAbstractRiskProperty)

#ifndef QLIB_RISKAXIS_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<IRiskAxis>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<IRiskAxis>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<IRiskAxis>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<IRiskAxis>);
#endif

/**
 * A property of the world (mostly of a particular market object), with respect
 * to which risk can be estimated
 *
 * By IRiskAxis we mean a 1D manifold along which we can change the state of
 * the world.  The most common implementation is PropertyRiskAxis<PROPERTY>,
 * which parameterises changes to a "hard-wired" property of a particular
 * market name.  For instance PropertyRiskAxis<Spot>("IBM") denotes IBM spot.
 * CalibratorFieldRiskAxis is a "dynamically-typed" counterpart.
 *
 * The key method is hypothesis(), which gives you an IHypothesis representing
 * a move along the risk axis your chosen magnitude.
 *
 * Risk axes corresponding to a property of a market data object (i.e. those
 * which can be generated using the IRiskProperty::axisFor() method of some
 * IRiskProperty implementation) will tell you what that property was: see
 * abstractProperty() and marketDataName().  The following should hold for
 * any IRiskProperty:
 *
 * -     if axis = property->axisFor(subjectName, qualifier) then:
 * -     axis->abstractProperty()->equalTo(property) == true
 * -     axis->marketDataName() == (subjectName)
 * -     axis->qualifier() == qualifier [assuming axis is an IQualifiedRiskAxis]
 *
 * For an overview of how these classes fit into the "declarative"
 * sensitivities framework, see the IRiskQuantityFactory class documentation.
 */

class RISKMGR_DLL IRiskAxis: public virtual IObject {

public:

    static CClassConstSP const TYPE;

    IRiskAxis();
    virtual ~IRiskAxis();

    /**
     * The IRiskProperty which changes as we move along this risk axis (or
     * NULL)
     *
     * IRiskAxis's corresponding to an IRiskProperty, i.e. generatable using
     * the latter's IRiskProperty::axisFor() method, should return something
     * equalTo that IRiskProperty.  If there's no suitable property, just
     * return a NULL smartPtr --- but note that you won't then be compatible
     * with RiskMapping.
     */

    virtual IAbstractRiskPropertyConstSP abstractProperty() const = 0;

    /**
     * The market data name to which the risk axis relates (or NULL if there
     * isn't a single one in particular)
     */

    virtual OutputNameConstSP marketDataName() const = 0;

    /**
     * See PropertyRiskAxis<PROPERTY>::hypothesis().
     */

    virtual IHypothesisConstSP hypothesis(double coeff) const = 0;

    /**
     * See RiskAxis.
     */

    virtual RiskAxisConstSP frozen() const = 0;

    /**
     * An "affine" IRiskAxis "shifted" by a given scenario
     *
     * Takes an IRiskAxis @a axis, and returns one representing that same risk
     * axis "after" a given scenario has been applied to the world.  The
     * IHypothesis's arising from IRiskAxis::hypothesis() are prefixed by @a
     * hypothesis.
     *
     * For instance, prefixing a PropertyRiskAxis<Spot> with a time-shift
     * hypothesis would give you a "delta next day" axis.
     *
     * You probably want to use IRiskProperty::conditioned(), rather than
     * calling this directly.
     */

    static IRiskAxisSP conditioned(IHypothesisConstSP hypothesis,
                                   IRiskAxisConstSP axis);

    /**
     * A compound IRiskAxis, parameterising parallel changes to a set of
     * underlying axes
     *
     * If @a axes is an IRiskAxisArray, then
     *
     * -   IHypothesis hyp = IRiskAxis::compound(axes)->hypothesis(coeff)
     *
     * is a CompoundHypothesis comprising
     *
     * -   { a->hypothesis(coeff) for a in axes }
     *
     * Generally you'll want to supply a null @a metric (which is the default).
     * Then the IHypothesis::AlternateWorld::distance() reported by the
     * returned @a hyp will just be @a coeff.
     *
     * Note: there is of course a class CompoundRiskAxis to implement this
     * method: see IRiskAxis.cpp.
     */

    static IRiskAxisSP compound(
        IRiskAxisArrayConstSP axes,
        IHypothesis::IDistanceMetricConstSP metric =
            IHypothesis::IDistanceMetricConstSP());
};

DRLIB_END_NAMESPACE

#endif
