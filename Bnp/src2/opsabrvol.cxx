

/******************************************************************************/
/*******************************************************************************
*
* FUNCTION     	: srt_f_optsarbvol(...)
*
*
* PURPOSE      	: A Quick vol transformation to go from a BS vol to a BS beta
vol and vice versa
*
* DESCRIPTION  	: The volatility transformation is based on a small noise
expansion of the density for both BS model and its equivalent BS beta (see Pat
Hagan's paper on the subject)
*
*
* PARAMETERS
*	INPUT	    : fwd_price	    - forward underlying price
*              	: strike      	- strike price
*              	: Maturity      - initial time        , in years
*              	: Volinput      - volatility
*              	: alpha			-
*              	: beta         	-
*				: rho
*				: [Falg1]		- type of vol in input
*				: [flag2]		- type of vol in ouput
*				: the model
*								- dF =
a*F^{beta}dW *							      da = v a
dZ
*                                 <dZ        ,dW> = rho dt
* RETURNS      	: vol			- equivalent [flag2] vol
*
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

Err srt_f_cutstrike(double forward, double strike, double mat, double volinput,
                    double alpha, double beta, double rho, double lambda,
                    double *numericalparams, int nparams,
                    SrtDiffusionType input,
                    Err (*bsvol)(double F, double K, double T, double sigma,
                                 double alpha, double betaordelay, double rho,
                                 double lambda, double *numericalparams,
                                 int nparams, SrtDiffusionType typeinput,
                                 SrtDiffusionType typeoutput,
                                 double *voloutput),
                    int nstd, double *modstrike) {
  Err err = NULL;
  double ATMNormVol, ATMLogVol;

  if (nstd == 0) {
    *modstrike = strike;
  } else {
    if (strike < forward) {
      err = bsvol(forward, forward, mat, volinput, alpha, beta, rho, lambda,
                  numericalparams, nparams, input, SRT_LOGNORMAL, &ATMLogVol);

      *modstrike = DMAX(strike, forward * exp(-nstd * sqrt(mat) * ATMLogVol));
    } else {
      err = bsvol(forward, forward, mat, volinput, alpha, beta, rho, lambda,
                  numericalparams, nparams, input, SRT_NORMAL, &ATMNormVol);

      *modstrike = DMIN(strike, forward + nstd * sqrt(mat) * ATMNormVol);
    }
  }

  return err;
}

/* ------------------------------------------------------------------------------
 */

Err srt_f_cutvol(double forward, double strike, double mat, double volinput,
                 double alpha, double betaordelay, double rho, double lambda,
                 double *numericalparams, int nparams, SrtDiffusionType input,
                 SrtDiffusionType output, int nstd, double lowbetacut,
                 double highbetacut,
                 Err (*bsvol)(double F, double K, double T, double sigma,
                              double alpha, double betaordelay, double rho,
                              double lambda, double *numericalparams,
                              int nparams, SrtDiffusionType typeinput,
                              SrtDiffusionType typeoutput, double *voloutput),
                 double *vol) {
  double price;
  double modstrike;
  Err err = NULL;

  /* A few security checks */
  if ((mat <= 0.0) || (volinput == 0.0)) {
    *vol = 0.0;
    return NULL;
  }

  err = srt_f_cutstrike(forward, strike, mat, volinput, alpha, betaordelay, rho,
                        lambda, numericalparams, nparams, input, bsvol, nstd,
                        &modstrike);
  if (err) {
    goto FREE_RETURN;
  }

  if (modstrike == strike) {
    err = bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                lambda, numericalparams, nparams, input, output, vol);
  } else {
    if (strike < forward) {
      if (lowbetacut == 1) {
        if (output == SRT_LOGNORMAL) {
          err = bsvol(forward, modstrike, mat, volinput, alpha, betaordelay,
                      rho, lambda, numericalparams, nparams, input,
                      SRT_LOGNORMAL, vol);

          goto FREE_RETURN;
        } else if (output == SRT_NORMAL) {
          err = bsvol(forward, modstrike, mat, volinput, alpha, betaordelay,
                      rho, lambda, numericalparams, nparams, input,
                      SRT_LOGNORMAL, vol);

          price = srt_f_optblksch(forward, strike, *vol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_NORMAL, vol);
        } else if (output == SRT_BETAVOL) {
          err =
              bsvol(forward, forward, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_BETAVOL, vol);
        }
      } else if (lowbetacut == 0) {
        if (output == SRT_NORMAL) {
          err =
              bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_NORMAL, vol);

          goto FREE_RETURN;
        } else if (output == SRT_LOGNORMAL) {
          err =
              bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_NORMAL, vol);

          price = srt_f_optblknrm(forward, strike, *vol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_LOGNORMAL, vol);
        } else if (output == SRT_BETAVOL) {
          err =
              bsvol(forward, forward, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_BETAVOL, vol);
        }
      } else {
        err =
            bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                  lambda, numericalparams, nparams, input, SRT_LOGNORMAL, vol);

        err = vol_conv(*vol, SABR_STR_LOG, vol, SABR_ATM_BETA, forward,
                       modstrike, mat, 0.0, lowbetacut, 0.0);
        err = srt_f_optsarbvol(forward, strike, mat, *vol, 0.0, lowbetacut, 0.0,
                               SRT_BETAVOL, output, vol);
      }

    } else if (strike > forward) {
      if (highbetacut == 1) {
        if (output == SRT_LOGNORMAL) {
          err = bsvol(forward, modstrike, mat, volinput, alpha, betaordelay,
                      rho, lambda, numericalparams, nparams, input,
                      SRT_LOGNORMAL, vol);

          goto FREE_RETURN;
        } else if (output == SRT_NORMAL) {
          err = bsvol(forward, modstrike, mat, volinput, alpha, betaordelay,
                      rho, lambda, numericalparams, nparams, input,
                      SRT_LOGNORMAL, vol);

          price = srt_f_optblksch(forward, strike, *vol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_NORMAL, vol);
        } else if (output == SRT_BETAVOL) {
          err =
              bsvol(forward, forward, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_BETAVOL, vol);
        }
      } else if (highbetacut == 0) {
        if (output == SRT_NORMAL) {
          err =
              bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_NORMAL, vol);

          goto FREE_RETURN;
        } else if (output == SRT_LOGNORMAL) {
          err =
              bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_NORMAL, vol);

          price = srt_f_optblknrm(forward, strike, *vol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_LOGNORMAL, vol);
        } else if (output == SRT_BETAVOL) {
          err =
              bsvol(forward, forward, mat, volinput, alpha, betaordelay, rho,
                    lambda, numericalparams, nparams, input, SRT_BETAVOL, vol);
        }
      } else {
        err =
            bsvol(forward, modstrike, mat, volinput, alpha, betaordelay, rho,
                  lambda, numericalparams, nparams, input, SRT_LOGNORMAL, vol);

        err = vol_conv(*vol, SABR_STR_LOG, vol, SABR_ATM_BETA, forward,
                       modstrike, mat, 0.0, highbetacut, 0.0);
        err = srt_f_optsarbvol(forward, strike, mat, *vol, 0.0, highbetacut,
                               0.0, SRT_BETAVOL, output, vol);
      }

    } else {
      err = bsvol(forward, forward, mat, volinput, alpha, betaordelay, rho,
                  lambda, numericalparams, nparams, input, output, vol);
    }
  }

FREE_RETURN:

  /* Return a success message */
  return err;
}

/* -------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------
 */
Err srt_f_optsarbvol(double forward, double strike, double mat, double volinput,
                     double alpha, double beta, double rho,
                     SrtDiffusionType input, SrtDiffusionType output,
                     double *vol) {
  double tempvol;
  double tempvol2;
  double price;
  Err err = NULL;

  /* A few security checks */
  if ((mat <= 0.0) || (volinput == 0.0)) {
    *vol = 0.0;
    return NULL;
  }

  if ((beta > 1.0) || (beta < 0.0))
    return serror("Beta has to be between 0.0 and 1.0");

  if ((input == SRT_BETAVOL) && (input == output)) {
    *vol = volinput;
  } else if (alpha == 0.0 && beta < 1.0E-08) {
    /* we are in the normal case */

    /* Calculate the BETA vol */
    if (input == SRT_BETAVOL || input == SRT_NORMAL) {
      tempvol = volinput;
    } else if (input == SRT_LOGNORMAL) {
      /* the inputed vol is the ATM LOGNORMAL vol ! */
      price = srt_f_optblksch(forward, forward, volinput, mat, 1.0, SRT_CALL,
                              PREMIUM);

      err = srt_f_optimpvol(price, forward, forward, mat, 1.0, SRT_CALL,
                            SRT_NORMAL, &tempvol);
    } else {
      err = "SABR: invalid input";
      return err;
    }

    /* Then price the option */
    if (output == SRT_BETAVOL || output == SRT_NORMAL) {
      *vol = tempvol;
    } else if (output == SRT_LOGNORMAL) {
      price = srt_f_optblknrm(forward, strike, tempvol, mat, 1.0, SRT_CALL,
                              PREMIUM);

      err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                            SRT_LOGNORMAL, vol);
    } else {
      err = "SABR: invalid output";
      return err;
    }
  } else if (alpha == 0) {
    if (input == SRT_BETAVOL) {
      err = srt_f_optbetavoltoblkvol(forward, strike, volinput, mat, beta,
                                     &tempvol);

      if (output == SRT_NORMAL) {
        /* 	err = srt_f_optblkvoltobetavol(
                                                                        forward
                 , strike        , tempvol        , mat        , 0.0        ,
           vol); */
        price = srt_f_optblksch(forward, strike, tempvol, mat, 1.0, SRT_CALL,
                                PREMIUM);

        err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                              SRT_NORMAL, vol);
      } else
        *vol = tempvol;
    } else if ((input == SRT_LOGNORMAL)) {
      err = srt_f_optblkvoltobetavol(forward, forward, volinput, mat, beta,
                                     &tempvol);
      if (output == SRT_BETAVOL) {
        *vol = tempvol;
      } else {
        err = srt_f_optbetavoltoblkvol(forward, strike, tempvol, mat, beta,
                                       &tempvol2);
        if (output == SRT_NORMAL) {
          err = srt_f_optblkvoltobetavol(forward, strike, tempvol2, mat, 0.0,
                                         vol);
        } else
          *vol = tempvol2;
      }

    } else if (input == SRT_NORMAL) {

      /* err = srt_f_optbetavoltoblkvol(
                                                                      forward ,
                                                                      forward ,
                                                                      volinput ,
                                                                      mat , 0.0
         , &tempvol);
       */
      price = srt_f_optblknrm(forward, forward, volinput, mat, 1.0, SRT_CALL,
                              PREMIUM);

      err = srt_f_optimpvol(price, forward, forward, mat, 1.0, SRT_CALL,
                            SRT_LOGNORMAL, &tempvol);

      err = srt_f_optblkvoltobetavol(forward, forward, tempvol, mat, beta,
                                     &tempvol2);

      if (output == SRT_BETAVOL) {
        *vol = tempvol2;
      } else {
        err = srt_f_optbetavoltoblkvol(forward, strike, tempvol2, mat, beta,
                                       &tempvol);

        if (output == SRT_NORMAL) {
          /* err = srt_f_optblkvoltobetavol(
                                                          forward        ,
                                                          strike        ,
                                                          tempvol        ,
                                                          mat        ,
                                                          0.0        ,
                                                          vol);  */
          price = srt_f_optblksch(forward, strike, tempvol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_NORMAL, vol);

        } else
          *vol = tempvol;
      }

    } else {
      err = "option not availlable";
      return err;
    }

  } else {
    if (input == SRT_BETAVOL) {
      err = srt_f_optbetastochvoltoblkvol(forward, strike, volinput, alpha, rho,
                                          mat, beta, &tempvol);

      if (output == SRT_NORMAL) {
        /* 	err = srt_f_optblkvoltobetavol(
                        forward        ,
                        strike        ,
                        tempvol        ,
                        mat        ,
                        0.0        ,
                        vol); */

        price = srt_f_optblksch(forward, strike, tempvol, mat, 1.0, SRT_CALL,
                                PREMIUM);

        err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                              SRT_NORMAL, vol);
      } else
        *vol = tempvol;

    } else if ((input == SRT_LOGNORMAL)) {
      err = srt_f_optblkvolATMtobetavolStochVol(forward, forward, volinput, mat,
                                                alpha, beta, rho, &tempvol);

      if (output == SRT_BETAVOL) {
        *vol = tempvol;
      } else {
        err = srt_f_optbetastochvoltoblkvol(forward, strike, tempvol, alpha,
                                            rho, mat, beta, &tempvol2);
        if (output == SRT_NORMAL) {
          /*	err = srt_f_optblkvoltobetavol(
                                                                  forward ,
                                                                  strike ,
                                                                  tempvol2 , mat
             , 0.0        , vol); */
          price = srt_f_optblksch(forward, strike, tempvol2, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_NORMAL, vol);
        } else
          *vol = tempvol2;
      }

    } else if ((input == SRT_NORMAL)) {
      /*
      err = srt_f_optbetavoltoblkvol(
                                                                      forward ,
                                                                      forward ,
                                                                      volinput ,
                                                                      mat , 0.0
      , &tempvol);
    */

      price = srt_f_optblknrm(forward, forward, volinput, mat, 1.0, SRT_CALL,
                              PREMIUM);

      err = srt_f_optimpvol(price, forward, forward, mat, 1.0, SRT_CALL,
                            SRT_LOGNORMAL, &tempvol);

      err = srt_f_optblkvolATMtobetavolStochVol(forward, forward, tempvol, mat,
                                                alpha, beta, rho, &tempvol2);

      if ((output == SRT_LOGNORMAL) || (output == SRT_NORMAL)) {
        err = srt_f_optbetastochvoltoblkvol(forward, strike, tempvol2, alpha,
                                            rho, mat, beta, &tempvol);
        if (output == SRT_NORMAL) {

          /*	err = srt_f_optblkvoltobetavol(
                                                                                  forward        ,
                                                                                  strike        ,
                                                                                  tempvol        ,
                                                                                  mat        ,
                                                                                  0.0        ,
                                                                                  vol); */
          price = srt_f_optblksch(forward, strike, tempvol, mat, 1.0, SRT_CALL,
                                  PREMIUM);

          err = srt_f_optimpvol(price, forward, strike, mat, 1.0, SRT_CALL,
                                SRT_NORMAL, vol);
        } else
          *vol = tempvol;

      } else
        *vol = tempvol2;
    }
  }

  /* Return a success message */
  return NULL;
}

Err srt_f_optsabr_mr_vol(double forward, double strike, double mat,
                         double volinput, double alpha, double beta, double rho,
                         double lambda, double *numericalparams, int nparams,
                         SrtDiffusionType input, SrtDiffusionType output,
                         double *vol) {
  Err err = NULL;
  double modalpha;

  if (lambda < 0.0) {
    return serror("lambda has to be positive or null");
  } else if (lambda == 0) {
    modalpha = alpha;
  } else {
    modalpha = alpha * sqrt((1 - exp(-2 * lambda * mat)) / (2 * lambda * mat));
  }

  err = srt_f_optsarbvol(forward, strike, mat, volinput, modalpha, beta, rho,
                         input, output, vol);

  return err;
}

/* -------------------------------------------------------------------------------
 */
/* Formulas for quantoed SABR */
/* -------------------------------------------------------------------------------
 */

Err srt_f_optsabrvolq(
    double forward, // In the writeup        , this is f_init
    double strike,
    double maturity, // T
    double volFwd,   // sigma
    double alpha,    // alpha
    double
        beta, // beta TO DO: note that this is converted to shifted-log params
    double corrFwdVol, // rhoSigmaF
    double volFx,      // nu
    double corrFwdFx,  // rhoFX
    double corrVolFx,  // rhoSigmaX
    SrtDiffusionType output,
    double *forward_adj, // Quanto-adjusted forward
    double *vol_adj,     // Quanto-adjusted volatility ("sigma_beta")
    double *alpha_adj,   // Quanto-adjusted alpha
    double *rho_adj // Quanto-adjusted rho = corr(forward        , forward_vol)
) {
  // TO DO Generalize implementation
  // Currently only handles shifted log

  Err err = NULL;
  char *memErrMsg = "FATAL ERROR: Insufficient memory to create array. In "
                    "srt_f_optsabrvolq()";
  double t_init = 0.0;
  double y_init = 1.0;
  double Y_init;
  double F_init;

  // Shifted log parameters
  double a, b;

  // For numerical integration
  int i, n;
  double tol = 1.0e-6;
  double dt;
  double *t = NULL;
  double *g = NULL;
  double *h = NULL;
  double *sigma = NULL;
  double *rho = NULL;
  double *y1 = NULL;
  double *y2 = NULL;
  double *y3 = NULL;
  double *y4 = NULL;
  double *y5 = NULL;
  double *y6 = NULL;
  double *y7 = NULL;

  double c0, c1, c2, c3, c4, psi;
  double eta_hat, upsilon_hat, sigma_hat, rho_hat, alpha_hat, implied_vol;

  // Initialize
  *forward_adj = 0.0;
  *vol_adj = 0.0;
  *rho_adj = 0.0;
  *alpha_adj = 0.0;

  // Calculate a large number of grid points        , just to be safe
  n = (int)(20.0 * maturity) + 1;
  if (n < 25)
    n = 25;
  else if (n > 401)
    n = 401;

  // Convert beta to shifted-log parameters a        , b
  a = pow(forward, beta);
  b = beta * a / forward;
  a *= (1.0 - beta);

  // We need to do some numerical integration        , so let's define a uniform
  // grid of points between t_init and maturity
  t = (double *)malloc(n * sizeof(double));
  if (t == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  t[n - 1] = maturity;
  t[0] = t_init;
  for (i = n - 2; i != 0; --i) {
    dt = (t[i + 1] - t[0]) / (double)(i + 1);
    t[i] = t[i + 1] - dt;
  }

  g = (double *)malloc(n * sizeof(double));
  if (g == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  for (i = n - 1; i != -1; --i) {
    g[i] = exp(corrVolFx * volFx * alpha * (maturity - t[i]));
  }

  // Quanto adjusted volatility factor
  Y_init = y_init * g[0];

  h = (double *)malloc(n * sizeof(double));
  if (h == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  h[n - 1] = 0.0;
  for (i = n - 2; i != -1; --i) {
    c0 = Y_init * b * alpha / g[i + 1];
    c1 = corrFwdFx * volFx * volFwd;
    c2 = corrVolFx * volFx * alpha;
    c3 = corrFwdVol * volFwd * c0;
    c4 = c1 + h[i + 1] * (c2 + c3 + 0.5 * c0 * alpha * h[i + 1]);
    h[i] = h[i + 1] + (t[i + 1] - t[i]) * c4;
  }

  // Quanto adjusted forward
  // TO DO: This is only OK for the shifted log model
  if (fabs(b) > tol) {
    F_init = ((a + b * forward) * exp(b * h[0]) - a) / b;
  } else // Use expansion
  {
    c0 = b * h[0];
    c1 = 1.0 + c0 / 2.0 * (1.0 + c0 / 3.0);
    c2 = 1.0 + c1 * c0;
    F_init = (forward * c2 + a * h[0] * c1);
  }

  // Quanto adjusted term structures of sigma (the forward's vol) and rho
  // (correlation between the forward and its vol)
  sigma = (double *)malloc(n * sizeof(double));
  if (sigma == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  rho = (double *)malloc(n * sizeof(double));
  if (rho == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  for (i = n - 1; i != -1; --i) {
    c0 = alpha * h[i];
    psi = sqrt(volFwd * (volFwd + 2.0 * corrFwdVol * c0) + c0 * c0);
    sigma[i] = (Y_init / y_init) * psi / g[i];
    rho[i] = (corrFwdVol * volFwd + alpha * h[i]) / psi;
  }

  // Now solve ODEs to get the equivalent constant values
  y1 = (double *)malloc(n * sizeof(double));
  if (y1 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y2 = (double *)malloc(n * sizeof(double));
  if (y2 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y3 = (double *)malloc(n * sizeof(double));
  if (y3 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y4 = (double *)malloc(n * sizeof(double));
  if (y4 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y5 = (double *)malloc(n * sizeof(double));
  if (y5 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y6 = (double *)malloc(n * sizeof(double));
  if (y6 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y7 = (double *)malloc(n * sizeof(double));
  if (y7 == NULL) {
    err = serror(memErrMsg);
    goto free_return;
  }
  y1[n - 1] = y2[n - 1] = y3[n - 1] = y4[n - 1] = y5[n - 1] = y6[n - 1] =
      y7[n - 1] = 0.0;
  for (i = n - 2; i != -1; --i) {
    dt = t[i + 1] - t[i];
    y1[i] = y1[i + 1] + dt * sigma[i + 1] * sigma[i + 1];
    y2[i] = y2[i + 1] + dt * rho[i + 1] * sigma[i + 1] * alpha * y1[i + 1];
    y3[i] = y3[i + 1] + dt * alpha * alpha * y1[i + 1];
    y4[i] = y4[i + 1] + dt * 4.0 * sigma[i + 1] * sigma[i + 1] * y1[i + 1];
    y5[i] =
        y5[i + 1] + dt * (alpha * alpha * y1[i + 1] * y1[i + 1] +
                          6.0 * alpha * rho[i + 1] * sigma[i + 1] * y2[i + 1]);
    y6[i] = y6[i + 1] + dt * (2.0 * sigma[i + 1] * sigma[i + 1] * y2[i + 1] +
                              alpha * rho[i + 1] * sigma[i + 1] * y4[i + 1]);
    y7[i] = y7[i + 1] + dt * sigma[i + 1] * sigma[i + 1] * y4[i + 1];
  }

  // And finally...
  eta_hat = 2.0 * y2[0] / (y1[0] * y1[0]);
  upsilon_hat = 3.0 * (y5[0] - 3.0 * eta_hat * (y6[0] - eta_hat * y7[0])) /
                (y1[0] * y1[0] * y1[0]);
  sigma_hat =
      (y1[0] + y3[0] - 0.5 * upsilon_hat * y1[0] * y1[0]) / (maturity - t_init);
  upsilon_hat = sqrt(upsilon_hat);
  sigma_hat = sqrt(sigma_hat);
  alpha_hat = upsilon_hat * sigma_hat;
  rho_hat = eta_hat / upsilon_hat;

  // Now we can look up the implied vol at the strike using the above SABR
  // parameters
  err = srt_f_optsarbvol(F_init, strike, maturity, sigma_hat, alpha_hat, beta,
                         rho_hat, SRT_BETAVOL, output, &implied_vol);

  *forward_adj = F_init;
  *vol_adj = implied_vol;
  *alpha_adj = alpha_hat;
  *rho_adj = rho_hat;

free_return:

  if (y7 != NULL)
    free(y7);
  if (y6 != NULL)
    free(y6);
  if (y5 != NULL)
    free(y5);
  if (y4 != NULL)
    free(y4);
  if (y3 != NULL)
    free(y3);
  if (y2 != NULL)
    free(y2);
  if (y1 != NULL)
    free(y1);
  if (rho != NULL)
    free(rho);
  if (sigma != NULL)
    free(sigma);
  if (h != NULL)
    free(h);
  if (g != NULL)
    free(g);

  return err;
}