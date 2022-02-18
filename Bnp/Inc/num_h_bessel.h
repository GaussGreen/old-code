/* ==============================================================================
   FILENAME:  num_h_bessel.h

   PURPOSE: compute the non-integer Bessel functions
   ============================================================================== */

#ifndef NUM_H_BESSEL_H
#define NUM_H_BESSEL_H

double bessi0(double x);

/* Modified Bessel function K0(x) for any real x */

double bessk0(double x);

/* Modified Bessel function I1(x) for any real x */

double bessi1(double x);

/* Modified Bessel function K1(x) for any real x */

double bessk1(double x);

/* Modified Bessel function In(x) for any real x and n > 2 */

double bessk(int n, double x);

/* Modified Bessel function In(x) for any real x and n >=0 */

double bessi(int n, double x);

/* ========================================================================= */

/* -------------------------------------------------------------------------
   Modified Bessel functions of fractionnal order
   ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Evaluates Gamma1 and Gamma2 by Chebyshev expansion for fabs(x) < 0.5.
   Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
   -------------------------------------------------------------------------- */

Err beschb(double x, double* gam1, double* gam2, double* gampl, double* gammi);

/* -------------------------------------------------------------------------
   Returns the modified Bessel function ri=I(nu) , rk=K(nu) and their
   derivatives rip and rkp for positive x and for xnu >0.
   The relative accuracy is within one or two significant digits of BESSEL_EPS.
   FPMIN is a number close to the machine smallest floating point number.
   -------------------------------------------------------------------------- */

Err bessik(double x, double xnu, double* ri, double* rk, double* rip, double* rkp);

/* --------------------------------------------
         Wrapper function for bessik(NumRec) so we
         return only I_nu(x).
   -------------------------------------------- */
Err I_nu(double xnu, double x, double* i_nu);

/* --------------------------------------------
         wrapper function for bessik(NumRec) so we
         return only K_nu(x).
   -------------------------------------------- */
Err K_nu(double xnu, double x, double* k_nu);

#endif
