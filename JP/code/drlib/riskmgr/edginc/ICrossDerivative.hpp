/**
 * @file ICrossDerivative.hpp
 */

#ifndef DRLIB_ICrossDerivative_H
#define DRLIB_ICrossDerivative_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IResultsFunction)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(RiskQuantity)
FORWARD_DECLARE(ICrossDerivative)

/**
 * A cross-derivative operator
 *
 * This interface represents the idea of a cross-derivative operator, or more
 * correctly a rule for estimating a cross-derivative by finite-difference
 * tweaking.
 *
 * Only one is currently provided: cross().  It's used as
 * the @a derivative argument to CrossRiskPropertySensitivity::Deriv::Deriv().
 *
 * The concept of the derivative as an operator on functions is expressed in
 * the riskQuantity() method; you can see how it's used in
 * e.g. CrossRiskPropertySensitivity.cpp.
 *
 * For an overview of the "declarative" sensitivities framework of
 * which these classes form a part, see IRiskQuantityFactory.
 */

class RISKMGR_DLL ICrossDerivative: public virtual IObject {

public:

    static CClassConstSP const TYPE;

    ICrossDerivative();
    virtual ~ICrossDerivative();

    /**
     * Second derivative.
     *
     * Its riskQuantity() calculates
     *
     * -     (v++ + v-- - 2v00) + (g0 d0- d0+ + g1 d-0 d+0)  /
     *                ½ (d0+ - d0-) (d+0 - d-0)
     *
     * where
     *
     * -     g0 = (v0+ + v0- - 2 v00) / -d0- d0+
     * -     g1 = (v+0 + v-0 - 2 v00) / -d-0 d+0
     *
     * which shares all but two prices (v++, v--) with those required for
     * computing the IScalarDerivative::second()'s g0 and g1.  Cf. the comment
     * on calculateCrossDerivative in ScalarShift.cpp
     */

    static ICrossDerivativeConstSP cross();

    /**
     * A RiskQuantity which evaluates to the cross-derivative of a given function
     * with respect to a given pair of properties of the world
     *
     * Called from CrossRiskPropertySensitivity::nameRiskQuantities().
     *
     * In terms of the parameters explained below, we can define
     *
     *    -  v00 = value of @a function in the base-case world
     *
     *    -  v++ = value of @a function under the hypothesis that the property
     *             parameterised by @a axis0 differs by @a +coefficient0 and that
     *             parameterised by @a axis1 differs by @a +coefficient1
     *
     *    -  v-+ = value of @a function under the hypothesis that the property
     *             parameterised by @a axis0 differs by @a -coefficient0 and that
     *             parameterised by @a axis1 differs by @a +coefficient1
     *
     *    - ...
     *
     *    -  d++, d-+, ... = actual distances along @a axis between real and
     *             alternate worlds
     *
     * See cross() for the formulas used to combine these numbers.
     * 
     * (For an account of the difference between @a coefficient
     * and d+, see ITweakableWithRespectTo::sensShift().)
     *
     * @param function     The function whose derivative is to be estimated
     *
     * @param axis0        The first IRiskAxis along which the derivative is to
     *                     be taken: typically this is a
     *                     PropertyRiskAxis<PROPERTY>
     *
     * @param coefficient0 Controls the size of the perturbation along
     *                     @a axis0 to be used in estimating the derivative
     *
     * @param axis1        The second IRiskAxis along which the derivative is to
     *                     be taken
     *
     * @param coefficient0 Controls the size of the perturbation along
     *                     @a axis1 to be used in estimating the derivative
     */

    virtual RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                        IRiskAxisConstSP axis0,
                                        double coefficient0,
                                        IRiskAxisConstSP axis1,
                                        double coefficient1) const = 0;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
