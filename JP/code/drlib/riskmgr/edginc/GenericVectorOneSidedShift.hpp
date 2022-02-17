/**
 * @file GenericVectorOneSidedShift.hpp
 *
 *
 */

#ifndef EDG_GENERIC_VECTOR_ONE_SIDED_SHIFT_HPP
#define EDG_GENERIC_VECTOR_ONE_SIDED_SHIFT_HPP

#include "edginc/VectorShift.hpp"
#include "edginc/PointwiseRestorableWith.hpp"
#include "edginc/Additive.hpp"
#include "edginc/Maths.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Template to automate the boilerplate aspects of a VectorShift implementation
 */
template <class TWEAK, bool isAbstract = false>
class GenericVectorOneSidedShift: public VectorShift,
                                  public virtual Additive,
                                  public virtual TWEAK {
public:
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

    //// constructor with specified shift size
    GenericVectorOneSidedShift(double shiftSize): 
        VectorShift(TYPE, NAME, shiftSize) {}

    //// constructor with default shift size
    GenericVectorOneSidedShift(): VectorShift(TYPE, NAME, DEFAULT_SHIFT) {}

protected:
    //// for derived types
    GenericVectorOneSidedShift(const CClassConstSP& clazz,
                               const string&        sensName,
                               const double&        shiftSize):
        VectorShift(clazz, sensName, shiftSize) {}

public:
    /**
     * @name Boilerplate methods implemented automatically
     */

    //@{

    /**
     * The interface which the objects to be tweaked must implement
     *
     * (I.e. PointwiseTweakableWith<TWEAK>.)
     */

    typedef PointwiseTweakableWith<TWEAK> Tweakable;

    /**
     * The interface which the objects to be tweaked must implement if
     * in-place tweak/restore is to be used
     *
     * (I.e. PointwiseRestorableWith<TWEAK>.)
     */

    typedef PointwiseRestorableWith<TWEAK> Restorable;

    /**
     * The interface which the objects to be tweaked must implement
     *
     * (I.e. PointwiseTweakableWith<TWEAK>.)
     */

    CClassConstSP shiftInterface() const {
        return Tweakable::TYPE;
    }

    /**
     * The interface which the objects to be tweaked must implement if
     * in-place tweak/restore is to be used
     *
     * (I.e. PointwiseRestorableWith<TWEAK>.)
     */

    CClassConstSP restorableShiftInterface() const{
        return Restorable::TYPE;
    }

    /**
     * Whether @a obj's name is @a name
     *
     * Will crash if @a obj is not an instance of Tweakable.
     */
    bool nameMatches(const OutputName&     name, 
                     IObjectConstSP      obj) {
        return name.equals(dynamic_cast<const Tweakable&>(*obj).
                           sensName(this));
    }

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    void appendName(OutputNameArray&      namesList, 
                    IObjectConstSP      obj) {
        namesList.push_back(OutputNameSP(new OutputName(
            dynamic_cast<const Tweakable&>(*obj).sensName(this))));
    }

    /** The supplied object is queried for the expiries array needed
        for doing this tweak and this array is returned. The supplied
        object must implement the RhoPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj) {
        return dynamic_cast<const Tweakable&>(*obj).sensExpiries(this);
    }

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    bool shift(IObjectSP obj) {
        return dynamic_cast<Tweakable&>(*obj).sensShift(this);
    }

    /** Restores the object (which supports being tweaked
        by this type of sens control) to its original form */

    void restore(IObjectSP obj) {
        dynamic_cast<Restorable &>(*obj).sensRestore(this);
    }

    /**
     * The e to use in calculating (f(a+e)-f(a)) / e
     *
     * Note that this incorporates SENSITIVITY_UNIT so that the output
     * sensitivities end up scaled into the desired units (eg per basis point).
     */
    virtual double divisor() const {
        double shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)) {
            throw ModelException("GenericScalarOneSidedShift::divisor",
                                 "Shift size is zero");
        }
        // We use this as a hook to apply SENSITIVITY_UNIT
        return (shiftSize / SENSITIVITY_UNIT);
    } 
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP &clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericVectorOneSidedShift, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        if (!isAbstract){
            EMPTY_SHELL_METHOD(defaultConstructor);
            // no fields
            // register how to build our sensitivity
            SensitivityFactory::addSens(
                NAME,
                new GenericSensitivityFactory<GenericVectorOneSidedShift>(),
                new GenericVectorOneSidedShift(DEFAULT_SHIFT),
                Tweakable::TYPE);
        }
    }
    static IObject* defaultConstructor() {
        return new GenericVectorOneSidedShift();
    }
};

DRLIB_END_NAMESPACE

#endif
