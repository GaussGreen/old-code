/**
 * @file GenericAllNameScalarShift.hpp
 *
 *
 */

#ifndef DRLIB_GenericAllNameScalarShift_H
#define DRLIB_GenericAllNameScalarShift_H

#include "edginc/GenericAllNameScalarShiftBase.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/RestorableWith.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A tweak applied to all names simultaneously (rather than to each in turn, as
 * is more usual)
 *
 * This is not a particularly common requirement.  See e.g. CCMAbsoluteBetaSens for
 * an example instantiation.
 */

template <class TWEAK, bool constantDivisor = false> class GenericAllNameScalarShift: 
    public GenericAllNameScalarShiftBase,
    public virtual TWEAK{

public:

    /** Reflection type of this class */
    static CClassConstSP const TYPE;

protected:

    GenericAllNameScalarShift(CClassConstSP clazz,
                              const string& outputName,
                              double        shiftSize):
        GenericAllNameScalarShiftBase(clazz, outputName, shiftSize) {}

public:

    /**
     * @name Boilerplate methods implemented automatically
     */

    //@{

    /**
     * The interface which the objects to be tweaked must implement
     *
     * (I.e. TweakableWith<TWEAK>.)
     */

    typedef TweakableWith<TWEAK> Tweakable;

    /**
     * The interface which the objects to be tweaked must implement if
     * in-place tweak/restore is to be used
     *
     * (I.e. RestorableWith<TWEAK>.)
     */

    typedef RestorableWith<TWEAK> Restorable;

    /**
     * The interface which the objects to be tweaked must implement
     *
     * (I.e. TweakableWith<TWEAK>.)
     */

    CClassConstSP shiftInterface() const {
        return Tweakable::TYPE;
    }

    /**
     * The interface which the objects to be tweaked must implement if
     * in-place tweak/restore is to be used
     *
     * (I.e. RestorableWith<TWEAK>.)
     */

    CClassConstSP restorableShiftInterface() const {
        return Restorable::TYPE;
    }

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */

    bool shift(IObjectSP obj) {
        return dynamic_cast<Tweakable &>(*obj).sensShift(this);
    }

    /** Restores the object (which supports being tweaked
        by this type of sens control) to its original form */

    void restore(IObjectSP obj) {
        dynamic_cast<Restorable &>(*obj).sensRestore(this);
    }

    /**
     * Invoked when Class is 'loaded'
     *
     * When you instantiate this template, you get this for free but
     * you still need to ensure it gets called.  See for instance
     * CCMAbsoluteBetaSens.cpp .
     */

    static void load(CClassSP& clazz) {
        REGISTER(GenericAllNameScalarShift, clazz);
        SUPERCLASS(GenericAllNameScalarShiftBase);
    }
    
    /**
     * The e to use in calculating (f(a')-f(a)) / e
     *
     * Note that this incorporates getSensitivityUnit() so that the output
     * sensitivities end up scaled into the desired units (e.g. per
     * basis point).
     */
    virtual double divisor() const {
        if (constantDivisor){
            return 1.0/getSensitivityUnit();
        } else {
            return GenericAllNameScalarShiftBase::divisor();
        }
    }
    

    //@}
};

DRLIB_END_NAMESPACE

#endif
