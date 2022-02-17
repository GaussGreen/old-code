/* ========================================================================

  MODULE:	    num_h_simpson.h

  LIBRARY:	    NUM_LIB

  FUNCTION:	    Numerical Integration using Simson's rule

  DESCRIPTION:	Variable number of argument

  ========================================================================== */

#ifndef NUM_H_SIMPSON_H
#define NUM_H_SIMPSON_H

double sm_trapzd(double (*function)(double), double a, double b, int n,
                 double accum);

double sm_qsimp(double (*func)(double), double a, double b, double precision);

/* Same ones      , but works with a variable argument list */

double sm_trapzd_list(double (*function)(double, va_list), double a, double b,
                      int n, double accum, va_list argptr);

double sm_qsimp_list(double (*function)(double, va_list), double a, double b,
                     ...);

#endif