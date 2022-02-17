//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Integrator.hpp
//
//   Description :
//
//   Date        : 23 Nov 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Function.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"

#include <math.h>

#ifndef EDR_INTEGRATOR_HPP
#define EDR_INTEGRATOR_HPP

DRLIB_BEGIN_NAMESPACE


class Integrator1D;
typedef smartPtr<Integrator1D> Integrator1DSP;
typedef smartConstPtr<Integrator1D> Integrator1DConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<Integrator1D>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<Integrator1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<Integrator1D>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<Integrator1D>);
#endif

/** Base Integrator class */
class UTIL_DLL Integrator1D: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    virtual double integrate(const Function1DDouble& func) const = 0;

    static Integrator1DSP createIntegrator(const string& name,
                                           double        absPrecision,
                                           double        relPrecision);

private:
    static void load(CClassSP& clazz);
};


/** Trapezium rule with n-points */
class UTIL_DLL Trapez1DSimple:public CObject, public Integrator1D{
public:
    static CClassConstSP const TYPE;

    Trapez1DSimple(int nbSteps);

    virtual double integrate(const Function1DDouble& func) const;

protected:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultTrapez1DSimple();

    Trapez1DSimple();

    Trapez1DSimple(CClassConstSP clazz);

    Trapez1DSimple(CClassConstSP clazz, int nbSteps);

    int nbSteps;
};



/** Trapezium rule in closed interval
    Numerical Recipes in C page 137  */
class UTIL_DLL Trapez1D:public CObject, public Integrator1D{
public:
    static CClassConstSP const TYPE;
    friend class Trapez1DHelper;

    virtual double integrate(const Function1DDouble& func) const;

    static const double default_accuracy;
    static const int    default_nbIterMax;

protected:
    double accuracy;
    int nbIterMax;

	Trapez1D(double eps = default_accuracy, int jmax = default_nbIterMax);
    Trapez1D(int notUsed);
    double trapzd(const Function1DDouble* func, const double a, const double b, const int n) const;
};

typedef smartPtr<Trapez1D> Trapez1DSP;
typedef smartConstPtr<Trapez1D> Trapez1DConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<Trapez1D>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<Trapez1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<Trapez1D>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<Trapez1D>);
#endif


/** Romberg rule in closed interval
    Numerical Recipes in C page 140 */
class UTIL_DLL ClosedRomberg1D: public Trapez1D{
public:
    static CClassConstSP const TYPE;
    friend class ClosedRomberg1DHelper;

    virtual double integrate(const Function1DDouble& func) const;
    virtual void validatePop2Object(); /* For initializing nbIterMaxP */

    static const int    default_RichardsonPolyOrder;

	ClosedRomberg1D(double accuracy = default_accuracy,
                    int nbIterMaxP = default_nbIterMax,
                    int RichardsonPolyOrder = default_RichardsonPolyOrder);

protected:
    ClosedRomberg1D(int notUsed);

    int RichardsonPolyOrder;
    int nbIterMaxP;
};

typedef smartPtr<ClosedRomberg1D> ClosedRomberg1DSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<ClosedRomberg1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<ClosedRomberg1D>);
#endif


/** Romberg rule in an open interval
    Numerical Recipes in C page 143 */
class UTIL_DLL OpenRomberg1D:public CObject, public Integrator1D{
public:
    static CClassConstSP const TYPE;
    friend class OpenRomberg1DHelper;

    virtual double integrate(const Function1DDouble& func) const;
    virtual void validatePop2Object(); /* For initializing nbIterMaxP */

    OpenRomberg1D(double eps = default_accuracy,
                  int jmax = default_nbIterMax,
                  int k = default_RichardsonPolyOrder);

    class UTIL_DLL Choose{
    public:
        virtual ~Choose(){}
        virtual double operator()(int n) const = 0;
    };

    DECLARE_REF_COUNT(Choose);

    /* Bounded intervals */
    class UTIL_DLL MidpntRule: public Choose{
    public:
        MidpntRule(const Function1DDouble* infunc,
                   double                  a,
                   double                  b);
        double operator()(int n) const;
    private:
        const Function1DDouble* func;
        double a;
        double b;
    };

    /* Positive halfline */
    class UTIL_DLL ExpRule: public Choose{
    public:
        ExpRule(const Function1DDouble* infunc,
                double                  a);
        double operator()(int n) const;
    private:
        const Function1DDouble* func;
        double a;
        double b;
    };

    static const double default_accuracy;
    static const int    default_nbIterMax;
    static const int    default_RichardsonPolyOrder;

protected:
    double accuracy;
    int    nbIterMax;
    int    RichardsonPolyOrder;
    int    nbIterMaxP;

    OpenRomberg1D(int notUsed);
};

typedef smartPtr<OpenRomberg1D> OpenRomberg1DSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<OpenRomberg1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<OpenRomberg1D>);
#endif

/** Integrates the fourier transform
        1 / 2PI Integral_{-infty}^{+infty} exp(-i * u * w) f(u) du
    of a complex function u |--> f(u), u in R, using an FFT algorithm.
    It is assumed that f(-u) is == to the transpose of f(u) */
class UTIL_DLL FFTIntegrator1D: public CObject{
public:
    static CClassConstSP const TYPE;
    friend class FFTIntegrator1DHelper;

    /** Type of object returned by integrate method.
        In order to evaluate the integral at some "w" in the interval [-range / 2, +range / 2],
        w must be passed to getValue. getValue will then interpolate the integral along an
        nbPoint-array of integral values. */
    class UTIL_DLL Integral: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        virtual double getValue(double weight) const = 0;

        virtual void getValue(const CDoubleArray& weights,
                              CDoubleArray&       values) const = 0;

    private:
        static void load(CClassSP& clazz);
    };

    typedef smartPtr<Integral> IntegralSP;
    typedef smartConstPtr<Integral> IntegralConstSP;

    virtual IntegralConstSP integrate(const Function1DComplex& func) const;

    virtual void validatePop2Object();

	FFTIntegrator1D(double range = default_range, int nbPoint = default_nbPoint);

    static FFTIntegrator1D* create(double range = default_range, int nbPoint = default_nbPoint);

private:
    double range;
    int nbPoint;

    static const double default_range;
    static const int default_nbPoint;
};

typedef smartPtr<FFTIntegrator1D> FFTIntegrator1DSP;
typedef smartConstPtr<FFTIntegrator1D> FFTIntegrator1DConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<FFTIntegrator1D>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<FFTIntegrator1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<FFTIntegrator1D>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<FFTIntegrator1D>);
#endif


/** Wrapper around imsl's int_func_inf  */
class UTIL_DLL IntFuncInf:public CObject, public Integrator1D{
public:
    static CClassConstSP const TYPE;
    friend class IntFuncInfHelper;

    virtual double integrate(const Function1DDouble& func) const;

    static const double default_err_abs;
    static const double default_err_rel;
    static const int    default_max_subinter;

	IntFuncInf(double err_abs = default_err_abs,
               double err_rel = default_err_rel,
               int max_subinter = default_max_subinter);

    static IntFuncInf* create(double err_abs = default_err_abs,
                              double err_rel = default_err_rel,
                              int max_subinter = default_max_subinter);

private:
    IntFuncInf(int notUsed);

    double err_abs;
    double err_rel;
    int max_subinter;
};

typedef smartPtr<IntFuncInf> IntFuncInfSP;
typedef smartConstPtr<IntFuncInf> IntFuncInfConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<IntFuncInf>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IntFuncInf>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<IntFuncInf>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IntFuncInf>);
#endif


/** Wrapper around imsl's int_fcn  */
class UTIL_DLL IntFuncGeneral:public CObject, public Integrator1D {
public:
    static CClassConstSP const TYPE;
    friend class IntFuncGeneralHelper;

    virtual double integrate(const Function1DDouble& func) const;

    static const double default_err_abs;
    static const double default_err_rel;
    static const int    default_max_subinter;

	IntFuncGeneral(double err_abs = default_err_abs,
                   double err_rel = default_err_rel,
                   int max_subinter = default_max_subinter);

    static IntFuncGeneral* create(double err_abs = default_err_abs,
                                  double err_rel = default_err_rel,
                                  int max_subinter = default_max_subinter);

private:
    IntFuncGeneral(int notUsed);

    double err_abs;
    double err_rel;
    int max_subinter;
};

typedef smartPtr<IntFuncGeneral> IntFuncGeneralSP;
typedef smartConstPtr<IntFuncGeneral> IntFuncGeneralConstSP;
#ifndef QLIB_INTEGRATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<IntFuncGeneral>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IntFuncGeneral>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<IntFuncGeneral>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IntFuncGeneral>);
#endif

DRLIB_END_NAMESPACE

#endif
