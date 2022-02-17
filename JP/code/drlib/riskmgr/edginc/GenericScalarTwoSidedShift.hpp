#ifndef EDG_GENERICSCALAR_TWO_SIDED_SHIFT_HPP
#define EDG_GENERICSCALAR_TWO_SIDED_SHIFT_HPP

#include "edginc/ScalarShiftTwoSided.hpp"
#include "edginc/Additive.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/RiskProperty.hpp"

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
 * -  Instantiate GenericScalarTwoSidedShift using a typedef, with template 
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
template <class TAG, bool constantDivisor = false, bool isAbstract = false>
class GenericScalarTwoSidedShift: public ScalarShiftTwoSided {

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
     * Output name for 1st order derivative
     *
     * You must define this yourself.
     * See for example ParSpreadRhoParallelTwoSided.cpp.
     */
    const static string NAME;

    /**
     * Output name for 2nd order derivative
     *
     * You must define this yourself.
     * See for example ParSpreadRhoParallelTwoSided.cpp.
     */
    const static string SECOND_ORDER_NAME;

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

public:

    /** identifies the name used for storing the second order derivative 
        in the output*/
    virtual const string& getSecondOrderSensOutputName() const{
        return SECOND_ORDER_NAME;
    }

    /** A GenericScalarTwoSidedShift with shift size = DEFAULT_SHIFT */
    GenericScalarTwoSidedShift(): 
        ScalarShiftTwoSided(TYPE, constantDivisor, SENSITIVITY_UNIT, NAME,
                            new RiskProperty<TAG>(), DEFAULT_SHIFT) {}

    /** A GenericScalarTwoSidedShift with given shift size */
    GenericScalarTwoSidedShift(double shift): 
        ScalarShiftTwoSided(TYPE, constantDivisor, SENSITIVITY_UNIT, NAME,
                            new RiskProperty<TAG>(), shift) {}
    
protected:
    /** For derived classes */
    GenericScalarTwoSidedShift(
        CClassConstSP clazz,
        const string& outputName,
        double        shiftSize): 
            ScalarShiftTwoSided(
                clazz, constantDivisor, SENSITIVITY_UNIT, outputName,
                new RiskProperty<TAG>(), shiftSize)
    {}

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
        REGISTER(GenericScalarTwoSidedShift, clazz);
        SUPERCLASS(ScalarShiftTwoSided);
        if (!isAbstract){
            EMPTY_SHELL_METHOD(defaultConstructor);
            // no fields
            // register how to build our sensitivity
            SensitivityFactory::addSens(
                NAME,
                new GenericSensitivityFactory<GenericScalarTwoSidedShift>(),
                new GenericScalarTwoSidedShift(DEFAULT_SHIFT),
                ITweakableWithRespectTo<TAG>::TYPE);
        }
    }

    //@}

    static IObject* defaultConstructor(){
        return new GenericScalarTwoSidedShift();
    }
};

DRLIB_END_NAMESPACE

#endif
