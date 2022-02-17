
#ifndef EDG_VECTOR_SHIFT_TWO_SIDED_HPP
#define EDG_VECTOR_SHIFT_TWO_SIDED_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/TwoSidedDeriv.hpp"

DRLIB_BEGIN_NAMESPACE

/** Specialised version of VectorShift which implements a two sided calcuation
    of the derivative */
class RISKMGR_DLL VectorShiftTwoSided: public VectorShift,
                           public virtual ITwoSidedDeriv{
public:    
    static CClassConstSP const TYPE;

    virtual ~VectorShiftTwoSided();

    /** identifies the name used for storing the second order derivative 
        in the output*/
    virtual const string& getSecondOrderSensOutputName() const = 0;

    /** Overridden to perform two sided tweak */
    virtual void calculate(TweakGroup*  tweakGroup,
                           Results*     results);

    /** Overridden to scale 1st and 2nd order derivatives */
    virtual void scaleResult(Results*     results,
                             double       scaleFactor) const;

    /** Overridden to add 1st and 2nd order derivatives */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

    /** Returns the shift(s) which have been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;
protected:
    /** Note VectorShiftTwoSided is abstract. Create a scalar shift of
        type clazz, which uses outputName (eg VEGA_PARALLEL) to
        identify results and with given shiftSize */
    VectorShiftTwoSided(CClassConstSP clazz,
                        const string& outputName,
                        double        shiftSize);

    /** for reflection */
    VectorShiftTwoSided(CClassConstSP clazz,
                        const string& sensName);

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
