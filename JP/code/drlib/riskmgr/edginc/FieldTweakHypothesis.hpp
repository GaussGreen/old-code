/**
 * @file FieldTweakHypothesis.hpp
 */

#ifndef QLIB_FieldTweakHypothesis_H
#define QLIB_FieldTweakHypothesis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(FieldRiskAxis)
FORWARD_DECLARE(IFieldTweak)
FORWARD_DECLARE(FieldTweakHypothesis)

/**
 * A function which maps the world to an alternate world, in which a
 * dynamically specified field on a named market object has been modified
 *
 * For general information about the "flexible scenarios and greeks" framework,
 * see FlexibleSensitivity.
 */

class RISKMGR_DLL FieldTweakHypothesis: public AbstractPropertyTweakHypothesis {

    static void load(CClassSP& clazz);
    static IObject* emptyShell();

    FieldRiskAxisConstSP _axis;            // $required
    double _coefficient;                   // $required
    mutable IFieldTweakConstSP _tweak;     // $transient

    // for supporting setMarketDataName() --- temporary
    // compatibility mechanism
    OutputNameConstSP _marketDataName;     // $transient

    IFieldTweakConstSP tweak() const;

public:

    static CClassConstSP const TYPE;

    /**
     * Constructor
     */

    FieldTweakHypothesis(FieldRiskAxisConstSP axis,
                         double coefficient);

    /**
     * (For reflection pre-construction)
     */

    FieldTweakHypothesis();

    ~FieldTweakHypothesis();

    /**
     * The "risk axis" along which this hypothesis shifts the world
     *
     * Returns a FieldRiskAxis.
     */

    virtual IRiskAxisConstSP axis() const;

    /**
     * The interface which designates market objects affected by this hypothesis
     *
     * It's the CField::getDeclaringClass() of the first entry in the
     * FieldPath::path() of the @a field argument to
     * FieldRiskProperty::FieldRiskProperty().
     */

    virtual CClassConstSP shiftableInterface() const;

    /**
     * Make a change to the value of our field on a particular object
     *
     * Guts of this is implemented in ScalarFieldTweakHypothesis::sensShift(),
     * ExpiryFieldTweakHypothesis::sensShift().
     */

    virtual TweakOutcome sensShift(IObjectSP object) const;

    /**
     * The interface which designates market objects which can be
     * mutated in place to apply this hypothesis
     *
     * There isn't one currently, so this returns a null pointer.
     */

    virtual CClassConstSP restorableInterface() const;

    /**
     * Undo the effect of the last sensShift()
     *
     * Not implemented; throws an exception to that effect.
     */

    virtual void sensRestore(IObjectSP object) const;

    /**
     * "Market data name" of the object whose field value changes under this
     * hypothesis
     *
     * Originates as @a subjectName argument to FieldRiskProperty::axisFor().
     *
     * See sensName().
     */

    virtual OutputNameConstSP getMarketDataName() const;

    virtual IObjectConstSP getQualifier() const;

    /**
     * Amount by which the hypothesis moves the world along its FieldRiskAxis
     */

    virtual double axisCoefficient() const;

    /**
     * "Market data name" of an object
     *
     * Unlike PropertyTweakHypothesis::sensName() this doesn't do anything
     * tricky---it just casts @a object to MarketObject and returns
     * MarketObject::getName(), or returns "" if it isn't one.
     */

    virtual string sensName(IObjectConstSP object) const;

    /**
     * Support for Control::getCurrentSensitivity() compatibility
     *
     * These methods are provided so that we can implement a compatibility
     * mechanism for model/instrument logic which relies on interrogating
     * Control::getCurrentSensitivity() --- see <I>Current status</I> in the
     * IRiskQuantityFactory class documentation.
     */

    //@{

    virtual void setMarketDataName(OutputNameConstSP);

    virtual void setAxisCoefficient(double);

    //@}

    virtual string toString() const;
};

DRLIB_END_NAMESPACE

#endif
