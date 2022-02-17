/**
 * @file RiskQuantity.hpp
 */

#ifndef DRLIB_RiskQuantity_H
#define DRLIB_RiskQuantity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(RiskQuantity)
FORWARD_DECLARE(HypotheticalQuantity)
FORWARD_DECLARE(IRiskAxis)

/**
 * A quantity estimated by computing prices (or other values) in several
 * hypothetical states of the world, and performing some arithmetic on them.
 *
 * See IRiskQuantityFactory (<I>Design overview</I>) for a description of
 * the "declarative" sensitivities framework of which this class is a part.
 *
 * The process of estimating (say) a Delta by tweaking can be summarised as:
 *
 *    -  compute FV in the true state of the "world"
 *       (= model/instrument/market assembly);
 *
 *    -  compute FV in a hypothetical state of the world where spot
 *       has been tweaked;
 *
 *    -  calculate (hypothetical FV - base FV) / distance
 *       between the two worlds.
 *
 * A RiskQuantity is a specification for a slightly more general form of this
 * process:
 *
 *    -  The quantity computed in each world state doesn't have to be price.
 *       (Could be NAKED_BOND_PRICE, or anything.)
 *
 *    -  It may be computed in more than two world states, i.e. against
 *       more than one hypotheses.  That's needed for gammas, etc.
 *
 *    -  The hypotheses don't have to involve spot; they can be any
 *       IHypothesis.
 *
 *    -  The rule for combining the hypothetical quantities doesn't
 *       have to be "(x' - x) / distance"; it can be any function of their
 *       values and the distances of the corresponding worlds from
 *       the base world.
 *
 * parameters() defines the values to be computed and the world states in which
 * to compute them, as HypotheticalQuantity's; value() defines the rule for
 * combining them.
 *
 * You typically don't create RiskQuantity's for yourself: you get them from
 * IScalarDerivative::riskQuantity() implementations.  For instance,
 * IScalarDerivative::oneSided() ->riskQuantity(...) gives you a RiskQuantity
 * which calculates a one-sided Greek.
 *
 * The work of evaluating RiskQuantity's is carried out by
 * RiskQuantityEvaluator::storeResults().  (The reason that's not a method on
 * RiskQuantity is that for efficiency we want to share the work across
 * all the RiskQuantity's we have to compute.)
 */

class RISKMGR_DLL RiskQuantity: public CObject {
    RiskQuantity(const RiskQuantity& rhs);
    RiskQuantity& operator=(const RiskQuantity& rhs);

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    FORWARD_DECLARE(IFunction)

protected:

    HypotheticalQuantityArrayConstSP _parameters;

    RiskQuantity(CClassConstSP type,
                 HypotheticalQuantityArrayConstSP parameters);

    RiskQuantity(CClassConstSP type);

public:

    /**
     * A very simple RiskQuantity which evaluates to a constant.
     */

    static RiskQuantityConstSP constant(double c);

    /**
     * A very simple RiskQuantity which "evaluates" to "not applicable".
     *
     * Its value() throws NotApplicableException.  See
     * RiskQuantityFactorySensitivity::riskQuantities() for an example usage.
     */

    static RiskQuantityConstSP notApplicable();

    /**
     * A very simple RiskQuantity which "evaluates" to "untweakable".
     *
     * Its value() throws an exception.  See
     * PerNameRiskPropertySensitivity::nameRiskQuantities() for an example
     * usage.
     */

    static RiskQuantityConstSP untweakable(const ModelException&);

	/**
	 * Slight pun but used to differentiate between special NotApplicable
	 * and Untweakable varieties for exposure recording.
	 */

	virtual bool isExceptional() const;

	/**
	 * Allows test for notApplicable() without catching NotApplicableException
	 * which turns out to be slow.
	 */

	virtual bool isNotApplicable() const;

    ~RiskQuantity();

    /**
     * The hypothetical quantities to be evaluated and passed to value()
     * for combination.
     */

    HypotheticalQuantityArrayConstSP parameters() const;

    /**
     * Combine values of HypotheticalQuantity's, and the distances of the
     * corresponding hypothetical worlds from the base case, into a derivative
     * (or whatever).
     *
     * If parameters() contains N HypotheticalQuantity's, then both @a values
     * and @a distances are of length N, and
     *
     *    -  distances[i] is the amount by which (typically) some
     *       property of some market name in the world was tweaked
     *       under parameters()[i].hypothesis()
     *
     *    -  values[i] is the value which parameters()[i].quantity()
     *       was found to take under that hypothesis
     */

    virtual double value(const CDoubleArray& values,
                         const CDoubleArray& distances) const = 0;

    /**
     * Allows us to test for exposure as non-constant RiskQuantities are
     * defined to indicate that there is exposure.
     */
    virtual bool isConstant() const { return false; }
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
