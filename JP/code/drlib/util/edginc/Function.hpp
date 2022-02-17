//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Function.hpp
//
//   Description :
//
//   Date        : 22 May 02
//
//
//----------------------------------------------------------------------------

#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Complex.hpp"
#include "edginc/Range.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IIntegrand.hpp"

#ifndef EDR_FUNCTION_HPP
#define EDR_FUNCTION_HPP

DRLIB_BEGIN_NAMESPACE

/** Function f: R^N --> R^M */
class UTIL_DLL MFunctionND: public virtual IIntegrand{
public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const = 0;

    MFunctionND(int N,      // nb of vars
                int M);     // nb of funcs

    MFunctionND(int N,      // nb of vars
                int M,      // nb of funcs
                const RangeArray& intervals);   // intervals of definition

    int getNbVars() const;
    int getNbFuncs() const;
    const RangeArray& getIntervals() const;

	void setIntervals(RangeArray& _ranges);
	void setNbVars(int nbVars);
	void setNbFuncs(int nbFuncs);
    /**
     * Split this integrand into a vector of "1D" integrands
     * (eg: if integrand is a function f:R^N->R^M, to1DIntegrands() will
     * return M integrands that will be functions R^N->R).
     * 
     * WARNING: This is a generic implementation where each R^N->R integrand
     * still makes a call to original f:R^N->R^M (and then projects result
     * R^M to R), meaning that performance is generally very poor.
     * This method should be overriden on a case by case basis.
     * */
    virtual IIntegrandArrayConstSP to1DIntegrands() const;

protected:
    virtual ~MFunctionND(){}

private:
    int N;  // nb of vars
    int M;  // nb of funcs
    RangeArray ranges;
};
DECLARE_REF_COUNT(MFunctionND);

/** Real-valued N dimensional function */
class UTIL_DLL FunctionNDDouble: public MFunctionND{
public:
    FunctionNDDouble(int N); // N = nb of vars
    FunctionNDDouble(int N, const RangeArray& intervals);

    virtual double operator()(const CDoubleArray&  x) const = 0;

    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const;
};
DECLARE_REF_COUNT(FunctionNDDouble);

/** Real-valued 1 dimensional function */
class UTIL_DLL Function1DDouble: public FunctionNDDouble{
public:
    Function1DDouble();
    Function1DDouble(const Range& interval);

    const Range& getInterval() const;

    virtual double operator()(double  x) const = 0;

    virtual double operator()(const CDoubleArray&  x) const;
};
DECLARE_REF_COUNT(Function1DDouble);

/** Complex-valued 1 dimensional function */
class UTIL_DLL Function1DComplex: public MFunctionND{
public:
    Function1DComplex();
    Function1DComplex(const Range& interval);

    const Range& getInterval() const;

    virtual Complex operator()(double  x) const = 0;

    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const;
};

DECLARE_REF_COUNT(Function1DComplex);

/** Function f: R^N --> R^N */
class UTIL_DLL FunctionND: public MFunctionND{
public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const = 0;

    FunctionND(int N);     // nb of vars == nb of funcs

    FunctionND(int N,      // nb of vars == nb of funcs
               const RangeArray& intervals);   // intervals of definition

    int getN() const;
};

/** Function f: R^N --> R^N that supports computation of Jacobian */
class UTIL_DLL FunctionNDWithJacobian: public FunctionND{
public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const = 0;

    /** Given x, returns f(x) and Jacobian fjac(x)
        Recall fjac_i,j = df_i / dx_j */
    virtual void operator()(const DoubleArray&  x,
                            DoubleArray&        f,
                            DoubleMatrix&       fjac) const;

    /** Given x and f(x), returns Jacobian fjac(x).
        By default  fjac is calculated using a finite difference
        scheme using numerical recipes' fdjac routine */
    virtual void jac(const DoubleArray&  x,
                     const DoubleArray&  f,
                     DoubleMatrix&       fjac) const;

    FunctionNDWithJacobian(int N);     // nb of vars == nb of funcs

    FunctionNDWithJacobian(int N,      // nb of vars == nb of funcs
                           const RangeArray& intervals);   // intervals of definition

private:
    static void vecfunc(int, double [], double []);

    // transient
    mutable DoubleArraySP x;
    mutable DoubleArraySP fvec;
    static const FunctionNDWithJacobian* me;
};

DRLIB_END_NAMESPACE

#endif
