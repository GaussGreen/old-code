/**
 * @file IScalarDerivative.hpp
 */

#ifndef DRLIB_IScalarDerivative_H
#define DRLIB_IScalarDerivative_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IResultsFunction)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(RiskQuantity)
FORWARD_DECLARE(IScalarDerivative)

/**
 * A derivative operator
 *
 * This interface represents the idea of a derivative operator, or more
 * correctly a rule for estimating a derivative by finite-difference
 * tweaking.
 *
 * Three are provided: oneSided(), twoSided(), and second().  They're used as
 * the @a derivative argument to RiskPropertySensitivity::Deriv::Deriv().
 *
 * The methods underScenario() and averagedOverScenarios() allow you to convert
 * one of the three "plain" derivatives into one which is evaluated under some
 * hypothetical condition.
 *
 * The concept of the derivative as an operator on functions is expressed in
 * the riskQuantity() method; you can see how it's used in
 * e.g. PerNameRiskPropertySensitivity.cpp.
 *
 * For an overview of the "declarative" sensitivities framework of
 * which these classes form a part, see IRiskQuantityFactory.
 */

class RISKMGR_DLL IScalarDerivative: public virtual IObject {

public:

    static CClassConstSP const TYPE;

    IScalarDerivative();
    virtual ~IScalarDerivative();

    /**
     * One-sided first derivative
     *
     * Its riskQuantity() calculates (v+ - v) ÷ d+.
     */

    static IScalarDerivativeConstSP oneSided();

    /**
     * Two-sided first derivative
     *
     * Its riskQuantity() calculates (v+ - v-) ÷ (d+ - d-).
     */

    static IScalarDerivativeConstSP twoSided();

    /**
     * Second derivative.
     *
     * Its riskQuantity() calculates (v+ + v- - 2v) ÷ (d+ × d+).
     */

    static IScalarDerivativeConstSP second();

    /**
     * A derivative computed after a "scenario" hypothesis has been applied
     *
     * For instance,
     *
     * <PRE>
     *       IScalarDerivative::oneSided()->underScenario(
     *            PropertyTweakHypothesis<ParSpreadParallelRelative>(0.1),
     *            0.1)
     * </PRE>
     *
     * is a one-sided derivative, to be computed for each market name after
     * a 10% parallel increase has been applied to the par spreads of all
     * names in the market.
     */

    IScalarDerivativeConstSP underScenario(IHypothesisConstSP scenario) const;

    /**
     * A derivative computed after a "scenario" hypothesis has been applied to
     * each applicable name
     *
     * For instance,
     *
     * <PRE>
     *       IScalarDerivative::oneSided()->underScenario(
     *            RiskProperty<ParSpreadParallelRelative>(),
     *            0.1)
     * </PRE>
     *
     * is a one-sided derivative, to be computed for each market name after
     * a 10% parallel increase in that name's par spreads has been applied.
     * See CreditDeltaPointwiseWithShift for an example.
     */

    IScalarDerivativeConstSP underScenario(
        IScalarRiskPropertyConstSP scenarioProperty,
        double scenarioShift) const;

    /**
     * A derivative averaged over a number of "scenario" hypotheses, for each
     * applicable name
     *
     * For instance,
     *
     * -     IScalarDerivative::oneSided()->averagedOverScenarios(
     *           RiskProperty<ParSpreadParallelRelative>(),
     *           0.1, -0.1)
     *
     * is a one-sided derivative, computed for each market name after a 10%
     * parallel increase in that name's par spreads has been applied, computed
     * again for a 10% decrease, and averaged.  See
     * EnergyDeltaPointwiseWithShift for an example.
     */

    IScalarDerivativeConstSP averagedOverScenarios(
        IScalarRiskPropertyConstSP scenarioProperty,
        double scenario1Shift,
        double scenario2Shift) const;

    /**
     * A RiskQuantity which evaluates to the derivative of a given function
     * with respect to a given property of the world
     *
     * Called from PerNameRiskPropertySensitivity::nameRiskQuantities() and
     * AllNamesRiskPropertySensitivity::nameRiskQuantities().
     *
     * In terms of the parameters explained below, we can define
     *
     *    -  v = value of @a function in the base-case world
     *
     *    -  v+ = value of @a function under the hypothesis that the property
     *            parameterised by @a axis differs by @a +coefficient
     *
     *    -  v- = value of @a function under the hypothesis that the property
     *            parameterised by @a axis differs by @a -coefficient
     *
     *    -  d+, d- = actual distances along @a axis between real and alternate
     *                worlds
     *
     * See oneSided(), twoSided(), second() for the formulas they
     * use to combine these numbers.
     * 
     * (For an account of the difference between @a coefficient
     * and d+, see ITweakableWithRespectTo::sensShift().)
     *
     * @param function     The function whose derivative is to be estimated
     *
     * @param axis         The IRiskAxis along which the derivative is to be
     *                     taken: typically this is a
     *                     PropertyRiskAxis<PROPERTY>, for instance the
     *                     second() derivative of price with respect to
     *                     PropertyRiskAxis<Spot> is "gamma"
     *
     * @param coefficient  Controls the size of the perturbation along the
     *                     @a axis to be used in estimating the derivative
     */

    virtual RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                        IRiskAxisConstSP axis,
                                        double coefficient) const = 0;

    /**
     * Whether the derivative is first or second order
     *
     * oneSided() and twoSided() return 1, second() returns 2.  This is used in
     * PerNameRiskPropertySensitivity::nameRiskQuantities() and
     * AllNamesRiskPropertySensitivity::nameRiskQuantities() to obtain a
     * scaling factor for transforming the result into the right "sensitivity
     * unit".
     */

    virtual int order() const = 0;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
