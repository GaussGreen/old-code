/**
 * @file TweakOutcome.hpp
 */

#ifndef DRLIB_TweakOutcome_H
#define DRLIB_TweakOutcome_H

#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * The effect of a tweak to a property of a market object.
 *
 * This is what's returned from ITweakableWithRespectTo<PROPERTY>::sensShift().
 *
 * The important info is distance(); tweakMembers() is used in implementing
 * basket-type objects.  oldValue() is for interoperability with SensControl
 * and friends, and may go away in the medium term.
 */

class RISKMGR_DLL TweakOutcome: public CObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    double _oldValue;
    double _newValue;
    bool _tweakMembers;

public:

    /**
     * Report the amount by which a property of an object was altered, and
     * whether its component objects should be recursively altered.
     *
     * See ITweakableWithRespectTo<PROPERTY>::sensShift(), distance(),
     * tweakMembers().
     */

    TweakOutcome(double distance, bool tweakMembers);

    /**
     * Report the original and altered values of a property of an object,
     * and whether its component objects should be recursively altered.
     *
     * You want the simpler two-argument TweakOutcome(); see oldValue().
     */

    TweakOutcome(double oldValue, double newValue, bool tweakMembers);
    ~TweakOutcome();

    /**
     * The amount by which the property of the tweaked object was altered.
     *
     * See ITweakableWithRespectTo<PROPERTY>::sensShift().
     */

    double distance() const;

    /**
     * Whether 
     */

    bool tweakMembers() const;

    /**
     * See tweakMembers().
     *
     * Please don't use this method if you can avoid it.
     */

    void setTweakMembers(bool);

    /**
     * Whether oldValue() is valid (i.e. whether the 3-argument constructor was
     * used).
     */

    bool hasOldValue() const;

    /**
     * The original value of the property of the object that was tweaked.
     *
     * Only available if the 3-argument constructor was used.  Otherwise will
     * throw an exception.  oldValue() is provided to simplify interoperability
     * with SensControl::divisor() implementations which use the previous value
     * to back out the absolute change after a proportional tweak.  Mostly
     * you don't need to worry about it.
     */

    double oldValue() const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
