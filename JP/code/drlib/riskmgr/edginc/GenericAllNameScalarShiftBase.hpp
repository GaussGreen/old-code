/**
 * @file GenericAllNameScalarShift.hpp
 *
 *
 */

#ifndef DRLIB_GenericAllNameScalarShiftBase_H
#define DRLIB_GenericAllNameScalarShiftBase_H

#include "edginc/Additive.hpp"
#include "edginc/SensControlAllNames.hpp"
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A tweak applied to all names simultaneously (rather than to each in
 * turn, as is more usual). Note that this is a base non-templated class. It
 * contains the common code needed by GenericAllNameScalarShift
 *
 * This is not a particularly common requirement.  See
 * e.g. CCMAbsoluteBetaSens for an example instantiation.
 */

class RISKMGR_DLL GenericAllNameScalarShiftBase: public SensControlAllNames,
                                     public virtual IScalarTweak,
                                     public virtual Additive {
public:

    /** Reflection type of this class */
    static CClassConstSP const TYPE;

    /** Is this sensitivity made using a discrete shift (ie a jump) or
     * an approximately continuous one (returns false) */
    virtual bool discreteShift() const;

    /** Name of packet for output ("Instrument") */
    virtual const string& getPacketName() const;
    
    /** Units in which sensitivity outputs are expressed. Default
        implementation is 1.0 */
    virtual double getSensitivityUnit() const;

    /** Returns size of the shift */
    virtual double getShiftSize() const;

protected:

    GenericAllNameScalarShiftBase(CClassConstSP clazz,
                                  const string& outputName,
                                  double        shiftSize);

    /**
     * Compute the greek
     *
     * Calculates a one-sided first derivative.
     */
    virtual void calculate(TweakGroup* tweakGroup, Results* results);

    /**
     * The e to use in calculating (f(a')-f(a)) / e
     *
     * Note that this incorporates getSensitivityUnit() so that the output
     * sensitivities end up scaled into the desired units (e.g. per
     * basis point).
     */
    virtual double divisor() const;

private:
    static void load(CClassSP& clazz);
    /// fields ////
    double shiftSize;
};

DRLIB_END_NAMESPACE

#endif
