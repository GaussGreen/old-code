//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootFinder.hpp
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
#include "edginc/smartPtr.hpp"
#include <math.h>

#ifndef EDR_ROOT_FINDER_HPP
#define EDR_ROOT_FINDER_HPP

DRLIB_BEGIN_NAMESPACE
#if defined(_MSC_VER)
// disable warning about inheriting method from base class rather than virtual
// class
#pragma warning(disable : 4250)
#endif

class UTIL_DLL Func1D{
public:
    class UTIL_DLL NoDeriv{
    public:
        virtual double operator()(double  x) const = 0;
    protected:
        virtual ~NoDeriv(){}
    };

    class UTIL_DLL WtDeriv{
    public:
        virtual void operator()(double  x,
                                double& f,
                                double& df) const = 0;
    protected:
        virtual ~WtDeriv(){}
    };
};

class UTIL_DLL Bracketer1D: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    virtual void bracket(const Func1D::NoDeriv& func,
                         double&                x1,
                         double&                x2) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<Bracketer1D> Bracketer1DSP;
typedef smartConstPtr<Bracketer1D> Bracketer1DConstSP;
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<Bracketer1D>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<Bracketer1D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<Bracketer1D>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<Bracketer1D>);
#endif

class UTIL_DLL RootFinder1D{
public:
    template<class T>
    class TwoInitVal: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        virtual double solve(const T&    func,
                             double      x1,
                             double      x2) const = 0;

    private:
        static void load(CClassSP& clazz);
    };

    typedef TwoInitVal<Func1D::NoDeriv> TwoInitValNoDeriv;
    typedef smartPtr<TwoInitValNoDeriv> TwoInitValNoDerivSP;
    typedef smartConstPtr<TwoInitValNoDeriv> TwoInitValNoDerivConstSP;
    typedef TwoInitVal<Func1D::WtDeriv> TwoInitValWtDeriv;
    typedef smartPtr<TwoInitValWtDeriv> TwoInitValWtDerivSP;
    typedef smartConstPtr<TwoInitValWtDeriv> TwoInitValWtDerivConstSP;
};

#if defined(QLIB_BUILD_DLL) && defined(_MSC_VER) && (_MSC_VER >=1300)
// haven't had a chance to get to the bottom of this. But gcc does not like
// this at all.
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL RootFinder1D::TwoInitVal<Func1D::NoDeriv>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL RootFinder1D::TwoInitVal<Func1D::NoDeriv>);
#endif
#endif

template<> void RootFinder1D::TwoInitVal<
    Func1D::WtDeriv>::load(CClassSP& clazz);
template<> void RootFinder1D::TwoInitVal<
    Func1D::NoDeriv>::load(CClassSP& clazz);

/** Numerical recipes p. 352..
    Given a function func and an initial guessed range x1 to x2, the routine expands the range
    geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
    returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0). */
class UTIL_DLL ZBrac: public CObject,
             public Bracketer1D{
public:
    static CClassConstSP const TYPE;
    friend class ZBracHelper;

    ZBrac(double FACTOR = default_FACTOR,
          int    NUM_TRIES = default_NUM_TRIES);

    virtual void bracket(const Func1D::NoDeriv& func,
                         double&                x1,                 // not used on input
                         double&                x2) const;

    static const double default_FACTOR;
    static const int    default_NUM_TRIES;

private:

    double FACTOR;
    int    NUM_TRIES;
};

template<class Func1DNoDeriv>
void ZBrac_bracket(const Func1DNoDeriv& func,
                   double&              x1,
                   double&              x2,
                   double               FACTOR = ZBrac::default_FACTOR,
                   int                  NUM_TRIES = ZBrac::default_NUM_TRIES){
    static string method = "ZBrac_bracket";

    if (Maths::equals(x1, x2)){
        throw ModelException(method, 
                             Format::toString("Bad initial range. x1 (%f) must be different from x2 (%f).",
                                              x1,
                                              x2));
    }

    double f1 = func(x1);
    double f2 = func(x2);

    int j;
    for (j = 1; j <= NUM_TRIES; j++){
        if (f1*f2 < 0.0){
            return;
        }
        if (fabs(f1) < fabs(f2)){
            f1 = func(x1 += FACTOR * (x1 - x2));
        }
        else{
            f2 = func(x2 += FACTOR * (x2 - x1));
        }
    }
    
    /* Failure */
    throw ModelException(method,
                         "Number of iterations exceeded " + 
                         Format::toString(NUM_TRIES) + ".");
}

typedef smartPtr<ZBrac> ZBracSP;
typedef smartConstPtr<ZBrac> ZBracConstSP;
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<ZBrac>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<ZBrac>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<ZBrac>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<ZBrac>);
#endif

/** Adapted from Numerical recipes p. 352.
    Works for function defined in (0, +infty), as opposed to (-infty, +infty). 
    Do that by making a change of variable x->x0/(2.0-x/x0) in interval (-infty, x0], x0 > 0. */
class UTIL_DLL ZBracPositive: public CObject,
                     public Bracketer1D{
public:
    static CClassConstSP const TYPE;
    friend class ZBracPositiveHelper;

    ZBracPositive(double FACTOR = default_FACTOR,
                  int    NUM_TRIES = default_NUM_TRIES);

    virtual void bracket(const Func1D::NoDeriv& func,
                         double&                x1,                 // not used on input
                         double&                x2) const;

    static const double default_FACTOR;
    static const int    default_NUM_TRIES;

private:

    double FACTOR;
    int    NUM_TRIES;
};

/** Finds x1, x2 to bracket supplied func. If throwOnError is false then
    on failure, false is returned and no exception is thrown */
template<class Func1DNoDeriv>
bool ZBracPositive_bracket(
    const Func1DNoDeriv& func,
    double&              x1,    // not used on input
    double&              x2,
    bool                 throwOnError,
    double               FACTOR = ZBracPositive::default_FACTOR,
    int                  NUM_TRIES = ZBracPositive::default_NUM_TRIES){

    static string method = "ZBracPositive_bracket";

    if (!Maths::isPositive(x2)){
        if (!throwOnError){
            return false; // failure
        }
        throw ModelException(method, "Upper bound must be strictly positive.");
    }

    x1 = 0.0;
    double x0  = 0.5 * x2;

    try{
#define FUNC2(x) (func((x)))
#define FUNC1(x) (func(x0 / (2.0 - (x) / x0)))
        double f1 = FUNC1(x1);
        double f2 = FUNC2(x2);

        int j;
        for (j = 1; j <= NUM_TRIES; j++)
        {
            if (Maths::isNegative(f1 * f2)) {
                /* don't forget to make change of variable before exiting */
                x1 = x0 / (2.0 - x1 / x0);
                return true; // success
            }
            if (Maths::isPositive(fabs(f2) - fabs(f1))) {
                x1 += FACTOR * (x1 - x2);
                f1 = FUNC1(x1);
            }
            else {
                x2 += FACTOR * (x2 - x1);
                f2 = FUNC2(x2);
            }
        }
#undef FUNC2
#undef FUNC1
    } catch (exception&){
        if (!throwOnError){
            return false; // failure
        }
        throw;
    }
        
    /* Failure */
    if (!throwOnError){
        return false; // failure
    }
    throw ModelException(method,
                         "Number of iterations exceeded " + 
                         Format::toString(NUM_TRIES) + ".");
}

typedef smartPtr<ZBracPositive> ZBracPositiveSP;
typedef smartConstPtr<ZBracPositive> ZBracPositiveConstSP;
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<ZBracPositive>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<ZBracPositive>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<ZBracPositive>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<ZBracPositive>);
#endif

/** From Numerical recipes p. 366.
    Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
    between x1 and x2. The root, returned as the function value rtsafe, will be refined until
    its accuracy is known within xacc. funcd is a user-supplied routine that returns both the
    function value and the first derivative of the function. */
class UTIL_DLL RtSafe: public CObject,
              public RootFinder1D::TwoInitValWtDeriv{
public:
    static CClassConstSP const TYPE;
    friend class RtSafeHelper;

    RtSafe(double xacc,
           int    MAXIT = default_MAXIT);

    virtual double solve(const Func1D::WtDeriv&  func,
                         double                  x1,
                         double                  x2) const;

    static const int default_MAXIT;     // Maximum allowed number of iterations.

protected:
    RtSafe(const CClassConstSP& clazz);

private:
    RtSafe();

    double xacc;
    int MAXIT;     // Maximum allowed number of iterations.
};

/* If throwOnError is false, returns status (true = success) rather than
   throwing an exception if there is a 'numerical' error */
template<class Func1DWtDeriv>
bool RtSafe_solve(const Func1DWtDeriv&  func,
                    double                x1,
                    double                x2,
                    double                xacc,
                    bool                  throwOnError,
                    double&               result,
                    int                   MAXIT = RtSafe::default_MAXIT){
    static string method = "RtSafe_solve";

    double df,fh,fl;
    func(x1, fl, df);
    func(x2, fh, df);

    if ((Maths::isPositive(fl) && Maths::isPositive(fh)) 
        || (Maths::isNegative(fl) && Maths::isNegative(fh))) {
        if (!throwOnError){
            return false; // failure
        }
        throw ModelException(method, "Root must be bracketed.");
    }

    if (Maths::isZero(fl)) {
        result = x1;
        return true; // success
    }
    if (Maths::isZero(fh)) {
        result = x2;
        return true; // success
    }

    double xh,xl;
    if (Maths::isNegative(fl)) { // Orient the search so that f(xl) < 0.
        xl = x1;
        xh = x2;
    } 
    else {
        xh = x1;
        xl = x2;
    }

    double rts = 0.5 * (x1 + x2);         // Initialize the guess for root,
    double dxold = fabs(x2 - x1);         // the "stepsize before last,"
    double dx = dxold;                    // and the last step.
    double f;

    try{
        func(rts, f, df);
        
        int j;
        for (j = 1; j <= MAXIT; j++) { // Loop over allowed iterations.
            if ((Maths::isPositive(((rts - xh) * df - f) * ((rts - xl) * df - f)))     // Bisect if Newton out of range,
                || Maths::isPositive((fabs(2.0 * f) - fabs(dxold * df)))) {            // or not decreasing fast enough.
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (Maths::equals(xl, rts)) {
                    result = rts;
                    return true;        // Change in root is negligible.
                }
            } 
            else { // Newton step acceptable. Take it.
                dxold = dx;
                dx = f / df;
                double temp = rts;
                rts -= dx;
                if (Maths::equals(temp, rts)) {
                    result = rts;
                    return true;
                }
            }
            if (Maths::isNegative(fabs(dx) - xacc)) {
                result = rts;
                return true; // Convergence criterion.
            }
            func(rts, f, df);
            // The one new function evaluation per iteration.
            if (Maths::isNegative(f)) { // Maintain the bracket on the root.
                xl = rts;
            }
            else {
                xh = rts;
            }
        }
    } catch (exception&){
        if (!throwOnError){
            return false;
        }
        throw;
    }
       
    if (!throwOnError){
        return false;
    }
    // failure
    throw ModelException(method,
                         "Maximum number of iterations exceeded.");
}

typedef smartPtr<RtSafe> RtSafeSP;
typedef smartConstPtr<RtSafe> RtSafeConstSP;
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<RtSafe>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<RtSafe>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<RtSafe>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<RtSafe>);
#endif

/* Straight from NR.
   Using Brent's method, find the root of a function func known to lie between x1 and x2. The
   root, returned as zbrent, will be refined until its accuracy is tol. */
class UTIL_DLL ZBrent: public CObject,
              public RootFinder1D::TwoInitValNoDeriv{
public:
    static CClassConstSP const TYPE;
    friend class ZBrentHelper;

    ZBrent(double tol,
           int    ITMAX = default_ITMAX);

    virtual void validatePop2Object();

    virtual double solve(const Func1D::NoDeriv&  func,
                         double                  x1,
                         double                  x2) const;

    static const int default_ITMAX;
    static const double EPS;

protected:
    ZBrent(const CClassConstSP& clazz);

private:
    ZBrent();

    double tol;
    int ITMAX;      // Maximum allowed number of iterations.
};

template<class Func1DNoDeriv>
double ZBrent_solve(const Func1DNoDeriv&  func,
                    double                x1,
                    double                x2,
                    double                tol,
                    int                   ITMAX = ZBrent::default_ITMAX){
    static string method = "ZBrent_solve";

    int iter;
    double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
    double fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        throw ModelException(method, "Root must be bracketed.");
    }
    fc=fb;
    for (iter=1;iter<=ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a; // Rename a, b, c and adjust bounding interval d. 
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*ZBrent::EPS*fabs(b)+0.5*tol; // Convergence check.
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa; // Attempt inverse quadratic interpolation.
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q; // Check whether in bounds.
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d; // Accept interpolation.
                d=p/q;
            } else {
                d=xm; // Interpolation failed, use bisection.
                e=d;
            }
        } else { // Bounds decreasing too slowly, use bisection.
            d=xm;
            e=d;
        }
        a=b; // Move last best guess to a.
        fa=fb;
        if (fabs(d) > tol1) // Evaluate new trial root.
            b += d;
        else
            b += Maths::sign(tol1,xm);

        fb=func(b);
    }
    throw ModelException(method, "Maximum number of iterations exceeded.");
    return 0.0; // Never get here.
}

typedef smartPtr<ZBrent> ZBrentSP;
typedef smartConstPtr<ZBrent> ZBrentConstSP;
#ifndef QLIB_ROOTFINDER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<ZBrent>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<ZBrent>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<ZBrent>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<ZBrent>);
#endif

DRLIB_END_NAMESPACE

#endif




