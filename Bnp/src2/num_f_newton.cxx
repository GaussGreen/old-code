/* ==========================================================
   FILENAME :   num_f_newton.cxx

   PURPOSE:     Newton Method to find the root of a function

   DESCRIPTION:
                this subroutine solves f(x)=yans by newton's method
                it receives initial x vector (x1        ,x2        ,x3)
                and y vector (y1        ,y2        ,y3) where y(i)=f(i)
                if niter<= 2        , then we only use first derivatives
                                and x1        ,y1 are  dummies.
                                at the end x1 is assigned x2        , x2 is
   assigned x3 and x3 is assigned x3+dx2        , y1 is assigned y2 and y2 is
   assigned y3 if nstop=1        , we do not proceed as x3-x2 is small

   See NUMERICAL RECIPES IN C page 365
   ========================================================== */

#include "utallhdr.h"
#include <math.h"

void newton(double yans, double niter, double x[3], double y[3],
            double *nstop) {
  double d1, derv1, error, dx2, dfudge, sderv;
  double precision = 1.0e-07;
  *nstop = 0;

  /* Check on the x values to see if convergence has been reached (RELATIVE
   * precision)*/
  if (fabs(x[1]) != 0) {
    if (fabs(x[2] / x[1] - 1) < precision) {
      *nstop = 1;
      return;
    }
  } else {
    if (fabs(x[2] - x[1]) < precision) {
      *nstop = 1;
      return;
    }
  }

  /* Check on the y values to see if function is not locally flat (numerically
   * null gradient)*/
  if (fabs(y[1]) != 0) {
    if (fabs(y[2] / y[1] - 1.0) < precision) {
      *nstop = 1;
      return;
    }
  } else {
    if (fabs(y[2] - y[1]) < precision) {
      *nstop = 1;
      return;
    }
  }

  /* Start the computations */
  derv1 = (y[2] - y[1]) / (x[2] - x[1]);
  error = yans - y[2];
  dx2 = error / derv1;

  if (error == 0) {
    *nstop = 1;
    return;
  }

  if (niter >= 3) {
    d1 = (y[1] - y[0]) / (x[1] - x[0]);
    sderv = (derv1 - d1) / (x[2] - x[1]);
    dfudge = 0.5 * sderv * dx2 * dx2;
    dx2 = dx2 * error / (error + dfudge);
  }

  x[0] = x[1];
  x[1] = x[2];
  x[2] = x[2] + dx2;
  y[0] = y[1];
  y[1] = y[2];
}
