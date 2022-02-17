/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optscdout(14)
 *
 * PURPOSE      	: DOUBLE BARRIER SCUD Option
 *
 * DESCRIPTION  	: Barrier on a secondary underlying
 *
 *
 *******************************************************************************/

static double minmax_dbscud(double x, double y, double a, double b, double mu1,
                            double mu2, double sig1, double sig2, double rho,
                            double mat, int maxmax, int minmin);

double srt_f_optdbscud(double fwdx, double fwdy, double spotx, double spoty,
                       double strike, double barrierx, double barriery,
                       double sigx, double sigy, double rho, double mat,
                       double disc, SrtCallPutType call_put, int scud_type1,
                       int scud_type2) {
  double mux;
  double muy;

  double l_x;
  double l_y;
  double l_k;

  double premium;

  int maxmax, minmin;
  int cp;
  int du1;
  int du2;

  if (mat < 0)
    return 0;
  if (mat == 0) {
    if (call_put == SRT_CALL) {
      if (spotx > strike)
        return (spotx - strike);
      else
        return 0;
    } else if (call_put == SRT_PUT) {
      if (strike > spotx)
        return (strike - spotx);
      else
        return 0;
    }
  }
  mux = log(fwdx / spotx) / mat - (sigx * sigx / 2);
  muy = log(fwdy / spoty) / mat - (sigy * sigy / 2);

  l_x = log(barrierx / spotx);
  l_y = log(barriery / spoty);

  l_k = log(strike / spotx);

  cp = (call_put == SRT_CALL) ? 1 : -1;
  du1 = 1 - (2 * scud_type1);
  du2 = 1 - (2 * scud_type2);

  if (l_x * du1 > 0)
    return 0;
  if (l_y * du2 > 0)
    return 0;

  /*premium = minmax_dbscud ( x         , y        , a        , b        , mu1
   * , mu2        , sig1        , sig2        , rho        , mat)*/
  if (du1 == 1 && du2 == -1) {

    maxmax = 1;
    minmin = 1;

    if (cp == 1) /*** call x-down y-up ***/
    {
      premium = fwdx * minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                     muy + rho * sigx * sigy, sigx, sigy, rho,
                                     mat, maxmax, minmin);

      premium -= strike * minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                        sigy, rho, mat, maxmax, minmin);

    } else /*** put x-down y-up ***/
    {

      premium = fwdx * (

                           minmax_dbscud(l_x, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin) -
                           minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin));
      premium -= strike * (

                              minmax_dbscud(l_x, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin) -
                              minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin));
    }

  }

  else if (du1 == -1 && du2 == 1) {

    maxmax = 1;
    minmin = 1;

    if (cp == -1) /*** put x-up y-down ***/
    {
      premium = fwdx * minmax_dbscud(l_y, l_k, l_y, l_x,
                                     muy + rho * sigx * sigy, mux + sigx * sigx,
                                     sigy, sigx, rho, mat, maxmax, minmin);

      premium -= strike * minmax_dbscud(l_y, l_k, l_y, l_x, muy, mux, sigy,
                                        sigx, rho, mat, maxmax, minmin);

    } else /*** call x-up y-down ***/
    {

      premium =
          fwdx * (

                     minmax_dbscud(l_y, l_x, l_y, l_x, muy + rho * sigx * sigy,
                                   mux + sigx * sigx, sigy, sigx, rho, mat,
                                   maxmax, minmin) -
                     minmax_dbscud(l_y, l_k, l_y, l_x, muy + rho * sigx * sigy,
                                   mux + sigx * sigx, sigy, sigx, rho, mat,
                                   maxmax, minmin));
      premium -= strike * (

                              minmax_dbscud(l_y, l_x, l_y, l_x, muy, mux, sigy,
                                            sigx, rho, mat, maxmax, minmin) -
                              minmax_dbscud(l_y, l_k, l_y, l_x, muy, mux, sigy,
                                            sigx, rho, mat, maxmax, minmin));
    }
  }

  if (du1 == 1 && du2 == 1) {
    maxmax = 1;
    minmin = -1;

    if (cp == 1) /*** put x-down y-down ***/
    {
      premium = fwdx * minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                     muy + rho * sigx * sigy, sigx, sigy, rho,
                                     mat, maxmax, minmin);

      premium -= strike * minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                        sigy, rho, mat, maxmax, minmin);

    } else /*** call x-down y-down ***/
    {

      premium = fwdx * (

                           minmax_dbscud(l_x, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin) -
                           minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin));
      premium -= strike * (

                              minmax_dbscud(l_x, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin) -
                              minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin));
    }

  }

  else if (du1 == -1 && du2 == -1) {

    maxmax = -1;
    minmin = 1;

    if (cp == -1) /*** put x-up y-up ***/
    {
      premium = fwdx * minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                     muy + rho * sigx * sigy, sigx, sigy, rho,
                                     mat, maxmax, minmin);

      premium -= strike * minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                        sigy, rho, mat, maxmax, minmin);

    } else /*** call x-up y-up ***/
    {

      premium = fwdx * (

                           minmax_dbscud(l_x, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin) -
                           minmax_dbscud(l_k, l_y, l_x, l_y, mux + sigx * sigx,
                                         muy + rho * sigx * sigy, sigx, sigy,
                                         rho, mat, maxmax, minmin));
      premium -= strike * (

                              minmax_dbscud(l_x, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin) -
                              minmax_dbscud(l_k, l_y, l_x, l_y, mux, muy, sigx,
                                            sigy, rho, mat, maxmax, minmin));
    }
  }

  premium *= disc * cp;

  return premium;

} /* END srt_f_optdbscud() */

/******************************************************************************/

static double minmax_dbscud(double x, double y, double a, double b, double mu1,
                            double mu2, double sig1, double sig2, double rho,
                            double mat, int maxmax, int minmin) {

  double result = 0.0;
  double sig1rt, sig2rt;

  sig1rt = sig1 * sqrt(mat);
  sig2rt = sig2 * sqrt(mat);

  result = bivar(maxmax * (-x + mu1 * mat) / sig1rt,
                 minmin * (y - mu2 * mat) / sig2rt, -minmin * maxmax * rho);

  result -=
      exp(2.0 * mu2 * b / (sig2 * sig2)) *
      bivar(maxmax * (-x + (2.0 * rho * sig1 * b / sig2) + mu1 * mat) / sig1rt,
            minmin * (y - 2.0 * b - mu2 * mat) / sig2rt,
            -minmin * maxmax * rho);

  result -=
      exp(2.0 * mu1 * a / (sig1 * sig1)) *
      bivar(maxmax * (-x + 2.0 * a + mu1 * mat) / sig1rt,
            minmin * (y - (2.0 * rho * sig2 * a / sig1) - (mu2 * mat)) / sig2rt,
            -minmin * maxmax * rho);

  result +=
      exp(2.0 * mu1 * a / (sig1 * sig1) + 2.0 * mu2 * b / (sig2 * sig2) +
          4.0 * a * b * rho / (sig1 * sig2 * mat)) *
      bivar(maxmax * (-x + 2.0 * rho * sig1 * b / sig2 + 2.0 * a + mu1 * mat) /
                sig1rt,
            minmin * (y - 2.0 * rho * sig2 * a / sig1 - 2.0 * b - mu2 * mat) /
                sig2rt,
            -minmin * maxmax * rho);

  return result;
}
