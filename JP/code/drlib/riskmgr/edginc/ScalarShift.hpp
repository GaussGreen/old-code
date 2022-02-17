
#ifndef EDG_SCALAR_SHIFT_H
#define EDG_SCALAR_SHIFT_H
#include "edginc/SensControlPerName.hpp"
#include "edginc/ScalarPerNameShift.hpp"
#include "edginc/IHypothesis.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/** This is a work around for solaris 2.5.1 - the assembler breaks if the
    class RISKMGR_DLL name is too long. Basically map<string, string> breaks */
typedef map<string, OutputNameSP> MapStringToName;

/** Slightly more specialised SensControl where the shifting
    information is encapsulated by a single double and is associated to a
    particular name */
class RISKMGR_DLL ScalarShift: public SensControlPerName,
                   public virtual IScalarPerNameShift{
public:    
    friend class ScalarShiftHelper;
    static CClassConstSP const TYPE;

    ~ScalarShift();

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** Returns the scalar shift size */
    virtual double getShiftSize() const;

    /** Sets the scalar shift size */
    virtual void setShiftSize(double shiftSize);

    /** implements a one sided scalar derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup*  tweakGroup,
                   CResults*    results);

    /** IScalarPerNameShift implementation, for ImpliedScalarShift */
    virtual IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                                    double shiftSize,
                                                    IObjectSP world);

    /** IScalarPerNameShift implementation, for ImpliedScalarShift */
    virtual OutputNameArrayConstSP allNames(const IObject*) const;

    /** Calculates cross derivative using the two supplied shifts. 'This'
    shift is used for storing results. Note that the second order
    derivatives must have been calculated already for each of shift1
    and shift2. Additionally, each one must have stored the shift size
    used in that calculation.

    None of the inputs may be null.

    Algorithm used:
    ****************************************************************
    * cross gamma calculator - gamma with respect to multiple stocks 
    *
    *    gamma(i,j) = f(+,+) - f(-,+) - f(+,-) + f(-,-)  = gamma(j,i)
    *               ---------------------------------
    *                     (2*S(i)*dS) * (2*S(j)*dS)
    *                                                              
    * [where f(+/-,+/-) = f( S(i) +/- S(i)*dS, S(j) +/- S(j)*dS )]
    *
    * This requires four reprices, but is equivalent to the following
    * that uses the already computed stock gammas and requires only
    * two additional pricings
    *
    *gmma(i,j)=f(+,+)+f(-,-)-2*f(0,0)-gamma(i)*(S(i)*dS)^2-gamma(j)*(S(j)*dS)^2
    *            -------------------------------------------------------------
    *                                 2 * (S(i)*dS) * (S(j)*dS)
    *
    * [where f(0,0) is the value of the option]
    ****************************************************************/
    /** added optional namePair filtering */
    static void calculateCrossDerivative(
        Sensitivity*  shift,         // used for storing results
        ScalarShift*  shift1,        // how to make the first shift
        const string& secondDeriv1,  // what 'gamma' is called for shift1
        ScalarShift*  shift2,        // how to make the second shift
        const string& secondDeriv2,  // what 'gamma' is called for shift2
        TweakGroup*   tweakGroup,
        CResults*     results);

    static void calculateCrossDerivative(
        Sensitivity*  shift,         // used for storing results
        ScalarShift*  shift1,        // how to make the first shift
        const string& secondDeriv1,  // what 'gamma' is called for shift1
        ScalarShift*  shift2,        // how to make the second shift
        const string& secondDeriv2,  // what 'gamma' is called for shift2
        TweakGroup*   tweakGroup,
        CResults*     results,
        MapStringToName* namePair);

protected:
    /** Note ScalarShift is abstract. Create a scalar shift of type clazz,
        which uses outputName (eg VEGA_PARALLEL) to identify results and
        with given shiftSize */
    ScalarShift(const CClassConstSP& clazz,
                const string&        outputName,
                const double&        shiftSize);


    /** for reflection */
    ScalarShift(const CClassConstSP& clazz,
                const string&        sensName);


    /** calculates delta/gamma like derivatives via 2 sided.
        ie calculates a first order two-sided scalar derivative. The current
        ScalarShift (ie this) is used to calculate one price. The value of the
        shiftsize is then negated when shifting for the second price */
    void calculateTwoSidedDeriv(
        const string&    secondDerivName, // eg GAMMA
        TweakGroup*      tweakGroup,
        CResults*        results);

    pair<double, double> twoSidedDerivs(TweakGroup*    tweakGroup,
                                        CResults*      results);

private:
    /* calculates single cross gamma - no checking for NULLs. Assumes
       that this ScalarShift has the required market data name set in
       it for the output and that shift1 and shift2 have the required
       market data names needed for the tweaking */
    static void singleCrossGamma(
        Sensitivity*   shift,             // used for storing results
        ScalarShift*   shifts[2],
        double         origPrice,         /* (I) original price */
        double         gammas[2],         /* (I) single gammas for the two */
        double         gammaShifts[2],    /* (I) shift size used for gammas */
        TweakGroup*    tweakGroup,
        OutputNameSP   outputNames[2],    // (I) result stored under both names
        CResults*      results);

    double shiftSize;
    ScalarShift(const ScalarShift& rhs);
    ScalarShift& operator=(const ScalarShift& rhs);
    void calcSingleTwoSidedFirstDeriv(
        const string&  secondDerivName, // eg GAMMA
        TweakGroup*    tweakGroup,
        CResults*      results);
};

typedef smartPtr<ScalarShift>      ScalarShiftSP;
typedef smartConstPtr<ScalarShift> ScalarShiftConstSP;
typedef array<ScalarShiftSP, ScalarShift> ScalarShiftArray;
typedef smartPtr<ScalarShiftArray>        ScalarShiftArraySP;
#ifndef QLIB_SCALARSHIFT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<ScalarShift>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<ScalarShift>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<ScalarShiftSP _COMMA_ ScalarShift>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<ScalarShiftArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<ScalarShift>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<ScalarShift>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<ScalarShiftSP _COMMA_ ScalarShift>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<ScalarShiftArray>);
#endif
DRLIB_END_NAMESPACE

#endif
