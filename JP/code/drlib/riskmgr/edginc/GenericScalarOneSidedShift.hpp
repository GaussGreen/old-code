/**
 * @file GenericScalarOneSidedShift.hpp
 */

#ifndef EDG_GENERICSCALAR_ONE_SIDED_SHIFT_HPP
#define EDG_GENERICSCALAR_ONE_SIDED_SHIFT_HPP

#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/Maths.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Template to automate the boilerplate aspects of a ScalarShift implementation
 *
 *
 * <H3>Purpose</H3>
 *
 * The point of this class is to absorb the boilerplate which is currently
 * copy-and-pasted between every implementation of ScalarShift.  See
 * ParSpreadRhoParallel.cpp, ParSpreadRhoParallelTwoSided.cpp for
 * examples.
 *
 *
 * <H3>How to use it</H3>
 *
 * -  Make a class defining the tweak with respect to which the sensitivity is
 *    calculated.  Let's call it TWEAK.  Example: ParSpreadRhoParallelTweak.
 *
 * -  Ensure that the (market) objects which get tweaked implement the interface
 *    TweakableWith<TWEAK>.  Example: ParSpreadCurve implements
 *    TweakableWith<ParSpreadRhoParallelTweak>.  (If desired,
 *    implement RestorableWith<TWEAK> as well.)
 *
 * -  Instantiate GenericScalarOneSidedShift using a typedef, with template 
 *    arguments TWEAK. Example: ParSpreadRhoParallel in
 *    ParSpreadRhoParallel.cpp.
 *
 * -  Define the NAME, DEFAULT_SHIFT,
 *    and SENSITIVITY_UNIT for the sensitivity.
 *    Example: ParSpreadRhoParallel.cpp.
 *
 * -  Register the class with the reflection system in the usual way
 *    (the load() method is predefined for you).
 *
 * This takes about 17 lines in total; all the necessary methods are
 * implemented automatically.  Note that the usage scheme is almost the same as
 * the informal cut-and-paste usage pattern for ScalarShift, but less work and
 * more explicit.
 *
 * The only significant difference is that the TWEAK class, which defines the
 * actual tweak made to the data and nothing else, is brought out separately
 * rather than being implicitly bundled into the Sensitivity along with
 * other things like one- versus two-sidedness, etc.  This makes it much
 * easier and cleaner to define several sensitivities against the same
 * tweak.
 *
 *
 * <H3>Template parameters</H3>
 *
 * <DL>
 * <DT>TWEAK
 * <DD>The interface defining the tweak with respect to which the sensitivity
 *     is calculated.  See for instance ParSpreadRhoParallelTweak.  It must have
 *     the following signature:
 *       <PRE>
 *         struct TWEAK {
 *           double getShiftSize() const;
 *         };
 *       </PRE>
 *       Or alternatively derive virtually from IScalarTweak
 * <DT>constantDivisor
 * <DD>With the default value of false, sensitivities will be computed as
 *     (shifted price - original price) * (shift size/SENSITIVITY_UNIT)
 *     If set to true, then instead
 *     (shifted price - original price) * (1.0/SENSITIVITY_UNIT)
 *     is used
 * </DL>
 * <DT>isAbstract
 * <DD>If isAbstract is true then no defaultConstructor will be registered in
 *     the load method nor will any sensitivity factory be registered
 * </DL>
 *
 *
 * <H3>See also</H3>
 * 
 *    -  ScalarShift is the (parent) class that actually does the work.
 *    -  Sample instantiations are ParSpreadRhoParallel.cpp,
 *       ParSpreadRhoParallelTwoSided.cpp.
 */
template <class TWEAK, bool constantDivisor = false, bool isAbstract = false>
class GenericScalarOneSidedShift: public ScalarShift,
                                  virtual public Additive,
                                  virtual public TWEAK {

    // Note: maybe we should own a TWEAK rather than be one; but this is more
    // or less dictated by the way things are set up.  If you stop inheriting
    // from TWEAK, do make sure to grep for "dynamic_cast<ParSpreadRho" (in the
    // whole library!).  William.

public:

    /**
     * @name Things you have to provide
     * See for example ParSpreadRhoParallelTwoSided.cpp.
     */

    //@{

    /**
     * Reflection type of this class
     */
    static CClassConstSP const TYPE;

    /**
     * Output name for derivative
     *
     * You must define this yourself.
     * See for example ParSpreadRhoParallelTwoSided.cpp.
     */
    const static string NAME;

    /**
     * Shift size to use if none provided
     *
     * You must define this yourself.
     * See for example ParSpreadRhoParallelTwoSided.cpp.
     */
    const static double DEFAULT_SHIFT;

    /**
     * Units in which sensitivity outputs are expressed
     *
     * You must define this yourself.  Normally it should just be 1.0, but
     * people want to see some sensitivities reported in terms of basis
     * points (0.0001). See for example ParSpreadRhoParallelTwoSided.cpp.
     */
    const static double SENSITIVITY_UNIT;

    //@}

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

    /**
     * Whether @a obj's name is @a name
     */
    bool nameMatches(const OutputName& name, IObjectConstSP obj) {
        return name.equals(dynamic_cast<const Tweakable&>(*obj).sensName(this));
    }

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    void appendName(OutputNameArray& namesList, IObjectConstSP obj) {
        const Tweakable& tweak = dynamic_cast<const Tweakable&>(*obj);
        OutputNameSP outputName(new OutputName(tweak.sensName(this)));;
        namesList.push_back(outputName);
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
     * @name Boilerplate methods implemented automatically:
     * essentially reflection and divisor() methods
     */

    //@{

    /** A GenericScalarOneSidedShift with shift size = DEFAULT_SHIFT */
    GenericScalarOneSidedShift(): ScalarShift(TYPE, NAME, DEFAULT_SHIFT) {}

    /** A GenericScalarOneSidedShift with given shift size */
    GenericScalarOneSidedShift(double shift): ScalarShift(TYPE, NAME, shift) {}
    
    /**
     * The e to use in calculating (f(a')-f(a)) / e
     *
     * Note that this incorporates SENSITIVITY_UNIT so that the output
     * sensitivities end up scaled into the desired units (eg per basis point).
     */
    virtual double divisor() const {
        if (constantDivisor){
            return 1.0/SENSITIVITY_UNIT;
        }
        double shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)) {
            throw ModelException("GenericScalarOneSidedShift::divisor",
                                 "Shift size is zero");
        }
        // We use this as a hook to apply SENSITIVITY_UNIT
        return shiftSize / SENSITIVITY_UNIT;
    } 
protected:
    /** For derived classes */
    GenericScalarOneSidedShift(
        CClassConstSP clazz,
        const string& outputName,
        double        shiftSize): ScalarShift(clazz, outputName, shiftSize) {}


private:
    /**
     * Invoked when Class is 'loaded'
     *
     * When you instantiate this template, you get this for free but
     * you still need to ensure it gets called.  See for instance
     * ParSpreadRhoParallelTwoSided.cpp .
     */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericScalarOneSidedShift, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        if (!isAbstract){
            EMPTY_SHELL_METHOD(defaultConstructor);
            // no fields
            // register how to build our sensitivity
            SensitivityFactory::addSens(
                NAME,
                new GenericSensitivityFactory<GenericScalarOneSidedShift>(),
                new GenericScalarOneSidedShift(DEFAULT_SHIFT),
                Tweakable::TYPE);
        }
    }

    //@}

    static IObject* defaultConstructor(){
        return new GenericScalarOneSidedShift();
    }
};



DRLIB_END_NAMESPACE

#endif
