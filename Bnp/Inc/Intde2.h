#ifndef __INTDE2_H__
#define __INTDE2_H__

/*
DE-Quadrature
Numerical Automatic Integrator for Improper Integral
    method    : Double Exponential (DE) Transformation
    dimension : one
    table     : use
functions
    intde  : integrator of f(x) over (a,b).
    intdei : integrator of f(x) over (a,infinity),
                 f(x) is non oscillatory function.
    intdeo : integrator of f(x) over (a,infinity),
                 f(x) is oscillatory function.
*/

/*
intde
    [description]
        I = integral of f(x) over (a,b)
    [declaration]							*/
void intdeini(int lenaw, double tiny, double eps, double* aw);
void intde(
    double (*f)(double, void*), double a, double b, double* aw, double* i, double* err, void* comm);
/*  [usage]
        intdeini(lenaw, tiny, eps, aw);  // initialization of aw
        ...
        intde(f, a, b, aw, &i, &err, &comm);
    [parameters]
        lenaw     : length of aw (int)
        tiny      : minimum value that 1/tiny does not
                    overflow (double)
        eps       : relative error requested (double)
        aw        : points and weights of the quadrature
                    formula, aw[0...lenaw-1] (double *)
        f         : integrand f(x) (double (*f)(double, void *))
        a         : lower limit of integration (double)
        b         : upper limit of integration (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
                comm      : pointer to a user-defined structure - passed to the function f (void *)
    [remarks]
        initial parameters
            lenaw > 1000,
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
        function
            f(x) needs to be analytic over (a,b).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^b |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination.
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,b).
                              you must divide the interval
                              (a,b) at this points.
                           2. relative error of f(x) is
                              greater than eps.
                           3. f(x) has oscillatory factor
                              and frequency of the oscillation
                              is very high.
*/

/*
intdei
    [description]
        I = integral of f(x) over (a,infinity),
            f(x) has not oscillatory factor.
    [declaration]							*/
void intdeiini(int lenaw, double tiny, double eps, double* aw);
void intdei(double (*f)(double, void*), double a, double* aw, double* i, double* err, void* comm);
/*  [usage]
        intdeiini(lenaw, tiny, eps, aw);  // initialization of aw
        ...
        intdei(f, a, aw, &i, &err, &comm);
    [parameters]
        lenaw     : length of aw (int)
        tiny      : minimum value that 1/tiny does not
                    overflow (double)
        eps       : relative error requested (double)
        aw        : points and weights of the quadrature
                    formula, aw[0...lenaw-1] (double *)
        f         : integrand f(x) (double (*f)(double, void *))
        a         : lower limit of integration (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
                comm      : pointer to a user-defined structure - passed to the function f (void *)
    [remarks]
        initial parameters
            lenaw > 1000,
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^infinity |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination.
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,infinity).
                              you must divide the interval
                              (a,infinity) at this points.
                           2. relative error of f(x) is
                              greater than eps.
                           3. f(x) has oscillatory factor
                              and decay of f(x) is very slow
                              as x -> infinity.
*/

/*
intdeo
    [description]
        I = integral of f(x) over (a,infinity),
            f(x) has oscillatory factor :
            f(x) = g(x) * sin(omega * x + theta) as x -> infinity.
    [declaration]							*/
void intdeoini(int lenaw, double tiny, double eps, double* aw);
void intdeo(
    double (*f)(double, void*),
    double  a,
    double  omega,
    double* aw,
    double* i,
    double* err,
    void*   comm);
/*  [usage]
        intdeoini(lenaw, tiny, eps, aw);  // initialization of aw
        ...
        intdeo(f, a, omega, aw, &i, &err, &comm);
    [parameters]
        lenaw     : length of aw (int)
        tiny      : minimum value that 1/tiny does not
                    overflow (double)
        eps       : relative error requested (double)
        aw        : points and weights of the quadrature
                    formula, aw[0...lenaw-1] (double *)
        f         : integrand f(x) (double (*f)(double, void *))
        a         : lower limit of integration (double)
        omega     : frequency of oscillation (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
                comm      : pointer to a user-defined structure - passed to the function f (void *)
    [remarks]
        initial parameters
            lenaw > 1000,
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^R |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination.
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,infinity).
                              you must divide the interval
                              (a,infinity) at this points.
                           2. relative error of f(x) is
                              greater than eps.
*/

#endif  // #ifndef __INTDE2_H__