/**
 * @file TweakFunction.hpp
 */

#ifndef QLIB_TweakFunction_H
#define QLIB_TweakFunction_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Range.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TweakFunction)

/**
 * A mutation which can be applied to the value of a numeric field as part of a
 * flexible scenario or greek
 *
 * Used by fieldRiskProperty (via IFieldTweak::IOperator::numeric()) to define
 * the operation to be performed on a numeric field value.  The key method is operator
 * ()().  There are four static factory functions additive(), multiplicative(),
 * automatic() and setter() which support the operations we actually use.
 *
 * For more information about the "flexible scenarios and greeks" framework,
 * see FlexibleSensitivity.
 */

class RISKMGR_DLL TweakFunction: public CObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    virtual double value(double argument, double x, const Range&) const = 0;

protected:

    TweakFunction(CClassConstSP type, bool clip);
    TweakFunction(CClassConstSP type);

public:

    /**
     * The mutated version of an original field value, given a "shift size" (or
     * other) argument
     *
     * For additive(), multiplicative() and automatic(), @a argument is a shift
     * size (either absolute or proportional).  setter() just returns @a
     * argument (@a x being ignored in this case).
     *
     * [Implementations don't override this method but rather value().]
     */

    double operator ()(double argument, double x, const Range& range,
                       bool clip) const;

    /**
     * Whether an @a argument of 0 to operator ()() implies that @a x will be
     * returned unchanged
     *
     * True for additive(), multiplicative(), automatic(), false for setter().
     */

    virtual bool zeroIsNoop() const = 0;

    ~TweakFunction();

    /**
     * An absolute tweak
     *
     * Its operator ()() maps (argument, x) to x + argument.  If the
     * result would fall outside @a range, then
     *
     *    -  if clip is true and the endpoint breached (floor/cap) is
     *       closed (floorInclusive/capInclusive), the endpoint itself
     *       is returned
     *
     *   -   if clip is true and the endpoint is open, the result is
     *       halfway between x and the endpoint
     *
     *   -   if clip is false, an exception is thrown
     *
     * An exception is also thrown if @a x is already outside @a range.
     */

    static TweakFunctionConstSP additive();

    /**
     * A proportional tweak
     *
     * Its operator ()() maps (argument, x) to x × (1 + argument).  If the
     * result would fall outside @a range, then
     *
     *    -  if clip is true and the endpoint breached (floor/cap) is
     *       closed (floorInclusive/capInclusive), the endpoint itself
     *       is returned
     *
     *   -   if clip is true and the endpoint is open, the result is
     *       halfway between x and the endpoint
     *
     *   -   if clip is false, an exception is thrown
     *
     * An exception is also thrown if @a x is already outside @a range.
     */

    static TweakFunctionConstSP multiplicative();

    /**
     * An exponential tweak
     *
     * Its operator ()() maps (argument, x) to x × exp argument.  If the
     * result would fall outside @a range, then
     *
     *    -  if clip is true and the endpoint breached (floor/cap) is
     *       closed (floorInclusive/capInclusive), the endpoint itself
     *       is returned
     *
     *   -   if clip is true and the endpoint is open, the result is
     *       halfway between x and the endpoint
     *
     *   -   if clip is false, an exception is thrown
     *
     * An exception is also thrown if @a x is already outside @a range.
     */

    static TweakFunctionConstSP exponential();

    /**
     * Either an absolute tweak, or an exponential one relative to one or both
     * endpoints of a given range
     *
     * What operator ()() on the returned TweakFunction does depends on the @a
     * range specified.  If the range is unbounded,
     *
     *    -  (argument, x) => x + argument  [same as additive()]
     *
     * If the range is bounded from above or below, x is tweaked exponentially
     * relative to the bound:
     *
     *    -  (argument, x) => bound + (x - bound) × exp argument
     *       [cf. multiplicative()]
     *       
     * If the range is bounded from both above and below, x is tweaked using
     * a sigmoid (squashing) function:
     *
     *    -  (argument, x) => s (s¯¹(x) + argument)
     *
     * where s maps the real line to the given range (a, b):
     *
     *    -  s(x) = a + (b - a) / (1 + exp (- x))
     *
     * An exception is thrown if the result would fall outside @a range (which
     * should never happen) or if @a x is already outside it.
     *
     * In the exponential or sigmoid case, if @a x is exactly on the boundary
     * of @a range, then it will map to itself whatever @a argument is.
     * In the context of a finite difference calculation, this will cause the
     * greek to fail as Untweakable with a "zero divisor" error.
     */

    static TweakFunctionConstSP adaptive();

    /**
     * A "tweak" which just maps everything to a fixed value
     *
     * Its operator () () maps (argument, x) to argument, ignoring x.
     */

    static TweakFunctionConstSP setter();

    /**
     * String representation of the TweakFunction, suitable for feeding back to
     * fromRepr()
     *
     * We use this in FieldRiskAxis to get a flat representation for simple db
     * storage.
     */

    virtual string repr() const;

    /**
     * Reconstruct a TweakFunction from a string representation
     */

    static TweakFunctionConstSP fromRepr(string repr);
};

DRLIB_END_NAMESPACE

#endif
