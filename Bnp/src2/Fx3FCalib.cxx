/* ==========================================================================
   FILE_NAME:	Fx3FCalib.cxx

   PURPOSE:		Modelling of the spot FX vol by taking in consideration
   a LGM 1 factor on the domestic and foreign market and a lognormal model on
   the Spot Fx

   DATE:		05/25/00

   AUTHOR:		L.C.
   ========================================================================== */

#include "LGMSVUtil.h"
#include "math.h"
#include "old_cpd_calib_wrapper.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

/*	A set of static functions designed to calculate exponential integrals */

double Phi_Func(double x, double T, double s, double t) {
  double result;

  result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;

  return result;
}

double Etha_Func(double x, double T, double s, double t) {
  double result;

  result = (t - s - Phi_Func(x, T, s, t)) / x;

  return result;
}

double Psi_Func(double x, double y, double T, double s, double t) {
  double result;

  result = 1.0 / (x * y) *
           (t - s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) +
            Phi_Func(x + y, T, s, t));

  return result;
}

double Gamma_Func(double x, double y, double T, double s, double t) {
  double result;

  result = 1.0 / x * (Phi_Func(y, T, s, t) - Phi_Func(x + y, T, s, t));

  return result;
}

double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t) {
  double result;

  result = 1.0 / (x * y) *
           (t - s - Phi_Func(x, Tx, s, t) - Phi_Func(y, Ty, s, t) +
            exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

  return result;
}

double Zeta_Func(double x, double y, double T, double s, double t) {
  double result;

  result = 1.0 / y * (Phi_Func(x, T, s, t) - Phi_Func(x + y, T, s, t));

  return result;
}

double Gamma2_Func(double x, double y, double Tx, double Ty, double s,
                   double t) {
  double result;

  result =
      1.0 / x *
      (Phi_Func(y, Ty, s, t) - exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

  return result;
}

static double H_func(double lam, double t1, double t2) {
  return exp(lam * t2) - exp(lam * t1);
}

static double Z_Func(double x, double t1, double t2) {
  return ((exp(x * t2) - exp(x * t1)) / x);
}

/*	A set of useful functions to calculate variances */

static Err Coefs_Partial_Var(double T, double T1, double T2, double sig_dom,
                             double lda_dom, double sig_for, double lda_for,
                             double correl_dom_for, double correl_dom_fx,
                             double correl_for_fx, double *a, double *b,
                             double *c) {
  /* coefficient of sig_fx 2 */
  (*a) = T2 - T1;

  /* coefficient of sig_fx   */
  (*b) = -2 * correl_for_fx * sig_for * Etha_Func(lda_for, T, T1, T2) +
         2 * correl_dom_fx * sig_dom * Etha_Func(lda_dom, T, T1, T2);

  /* constant coefficient    */
  (*c) = sig_for * sig_for * Psi_Func(lda_for, lda_for, T, T1, T2) +
         sig_dom * sig_dom * Psi_Func(lda_dom, lda_dom, T, T1, T2) -
         2 * correl_dom_for * sig_dom * sig_for *
             Psi_Func(lda_for, lda_dom, T, T1, T2);

  return NULL;
}

static Err Partial_Var(double T, double T1, double T2, double sig_dom,
                       double lda_dom, double sig_for, double lda_for,
                       double sig_fx, double correl_dom_for,
                       double correl_dom_fx, double correl_for_fx,
                       double *var) {
  double a;
  double b;
  double c;
  Err err = NULL;

  /* Get Coefs */
  err = Coefs_Partial_Var(T, T1, T2, sig_dom, lda_dom, sig_for, lda_for,
                          correl_dom_for, correl_dom_fx, correl_for_fx, &a, &b,
                          &c);

  if (err) {
    return err;
  }

  (*var) = a * sig_fx * sig_fx + b * sig_fx + c;

  return NULL;
}

/*	Check that there is no term structure of tau and get lambda */
Err get_unique_lambda(double *tau, int ntau, double *lda) {
  double the_tau;
  int i;

  if (tau == NULL || ntau < 1) {
    return "Uninitialised tau";
  }

  the_tau = tau[0];
  if (ntau > 1) {
    for (i = 1; i < ntau; i++) {
      if (tau[i] != the_tau) {
        return "No tau term structure is allowed";
      }
    }
  }

  if (the_tau < 1.0e-16) {
    return "No negative tau is allowed";
  }

  *lda = 1.0 / the_tau;

  return NULL;
}

/*	Merge rates term structures */
Err merge_rates_ts(double *sig_dom_mat, double *sig_dom, int sig_n_dom,
                   double *sig_for_mat, double *sig_for, int sig_n_for,
                   double **sig_merge_mat, double **sig_merge_dom,
                   double **sig_merge_for, int *sig_merge_n) {
  int i;
  Err err = NULL;

  *sig_merge_mat = NULL;
  *sig_merge_dom = NULL;
  *sig_merge_for = NULL;

  *sig_merge_mat = (double *)calloc(sig_n_dom, sizeof(double));
  if (!(*sig_merge_mat)) {
    err = "Memory allocation error in merge_rates_ts";
    goto FREE_RETURN;
  }

  memcpy(*sig_merge_mat, sig_dom_mat, sig_n_dom * sizeof(double));
  *sig_merge_n = sig_n_dom;
  num_f_concat_vector(sig_merge_n, sig_merge_mat, sig_n_for, sig_for_mat);
  num_f_sort_vector(*sig_merge_n, *sig_merge_mat);
  num_f_unique_vector(sig_merge_n, *sig_merge_mat);

  *sig_merge_dom = (double *)calloc(*sig_merge_n, sizeof(double));
  *sig_merge_for = (double *)calloc(*sig_merge_n, sizeof(double));

  if (!(*sig_merge_dom) || !(*sig_merge_for)) {
    err = "Memory allocation error (2) in merge_rates_ts";
    goto FREE_RETURN;
  }

  for (i = *sig_merge_n - 1; i >= 0; i--) {
    (*sig_merge_dom)[i] =
        sig_dom[Get_Index((*sig_merge_mat)[i], sig_dom_mat, sig_n_dom)];
    (*sig_merge_for)[i] =
        sig_for[Get_Index((*sig_merge_mat)[i], sig_for_mat, sig_n_for)];
  }

FREE_RETURN:

  if (err) {
    if (*sig_merge_mat)
      free(*sig_merge_mat);
    if (*sig_merge_dom)
      free(*sig_merge_dom);
    if (*sig_merge_for)
      free(*sig_merge_for);
  }

  return err;
}

/*	Merge rates        , fx and corr term structures */
Err merge_rates_fx_corr_ts(double *sig_dom_mat, double *sig_dom, int sig_n_dom,
                           double *sig_for_mat, double *sig_for, int sig_n_for,
                           double *sig_fx_mat, double *sig_fx, int sig_n_fx,
                           double *corr_mat, double *corr_dom_for,
                           double *corr_dom_fx, double *corr_for_fx,
                           int corr_n_mat, double **sig_merge_mat,
                           double **sig_merge_dom, double **sig_merge_for,
                           double **sig_merge_fx, double **corr_merge_dom_for,
                           double **corr_merge_dom_fx,
                           double **corr_merge_for_fx, int *sig_merge_n) {
  int i, index;
  Err err = NULL;

  *sig_merge_mat = NULL;
  *sig_merge_dom = NULL;
  *sig_merge_for = NULL;
  if (sig_fx_mat) {
    *sig_merge_fx = NULL;
  }
  *corr_merge_dom_for = NULL;
  *corr_merge_dom_fx = NULL;
  *corr_merge_for_fx = NULL;

  *sig_merge_mat = (double *)calloc(sig_n_dom, sizeof(double));
  if (!(*sig_merge_mat)) {
    err = "Memory allocation error in merge_rates_ts";
    goto FREE_RETURN;
  }

  memcpy(*sig_merge_mat, sig_dom_mat, sig_n_dom * sizeof(double));
  *sig_merge_n = sig_n_dom;
  num_f_concat_vector(sig_merge_n, sig_merge_mat, sig_n_for, sig_for_mat);

  if (sig_fx_mat) {
    num_f_concat_vector(sig_merge_n, sig_merge_mat, sig_n_fx, sig_fx_mat);
  }

  num_f_concat_vector(sig_merge_n, sig_merge_mat, corr_n_mat, corr_mat);
  num_f_sort_vector(*sig_merge_n, *sig_merge_mat);
  num_f_unique_vector(sig_merge_n, *sig_merge_mat);

  *sig_merge_dom = (double *)calloc(*sig_merge_n, sizeof(double));
  *sig_merge_for = (double *)calloc(*sig_merge_n, sizeof(double));

  if (sig_fx_mat) {
    *sig_merge_fx = (double *)calloc(*sig_merge_n, sizeof(double));
  }

  *corr_merge_dom_for = (double *)calloc(*sig_merge_n, sizeof(double));
  *corr_merge_dom_fx = (double *)calloc(*sig_merge_n, sizeof(double));
  *corr_merge_for_fx = (double *)calloc(*sig_merge_n, sizeof(double));

  if (!(*sig_merge_dom) || !(*sig_merge_for) ||
      (sig_fx_mat && !(*sig_merge_fx)) || !(*corr_merge_dom_for) ||
      !(*corr_merge_dom_fx) || !(*corr_merge_for_fx)) {
    err = "Memory allocation error (2) in merge_rates_fx_corr_ts";
    goto FREE_RETURN;
  }

  for (i = *sig_merge_n - 1; i >= 0; i--) {
    (*sig_merge_dom)[i] =
        sig_dom[Get_Index((*sig_merge_mat)[i], sig_dom_mat, sig_n_dom)];
    (*sig_merge_for)[i] =
        sig_for[Get_Index((*sig_merge_mat)[i], sig_for_mat, sig_n_for)];

    if (sig_fx_mat) {
      (*sig_merge_fx)[i] =
          sig_fx[Get_Index((*sig_merge_mat)[i], sig_fx_mat, sig_n_fx)];
    }

    index = Get_Index((*sig_merge_mat)[i], corr_mat, corr_n_mat);
    (*corr_merge_dom_for)[i] = corr_dom_for[index];
    (*corr_merge_dom_fx)[i] = corr_dom_fx[index];
    (*corr_merge_for_fx)[i] = corr_for_fx[index];
  }

FREE_RETURN:

  if (err) {
    if (*sig_merge_mat)
      free(*sig_merge_mat);
    if (*sig_merge_dom)
      free(*sig_merge_dom);
    if (*sig_merge_fx)
      free(*sig_merge_fx);
    if (*corr_merge_dom_for)
      free(*corr_merge_dom_for);
    if (*corr_merge_dom_fx)
      free(*corr_merge_dom_fx);
    if (*corr_merge_for_fx)
      free(*corr_merge_for_fx);
  }

  return err;
}

/*	Compute the variance of the forward factor */
static Err Fx3DtsFwdFactVar(double T0,    /*	Forward start date */
                            double Tval1, /*	Value date of the 1st forward */
                            double Tval2, /*	Value date of the 2nd forward */
                            double Tfix,  /*	Fix date of both forwards */
                            double *maturity_rates, long nbMat,
                            double *sig_curve_dom, double lda_dom,
                            double *sig_curve_for, double lda_for,
                            double correl_dom_for, double *varff) {
  long StartIndex, EndIndex, i;
  double t1, t2;
  double phidom, phifor, phidomfor;
  double lamdom, lamfor;
  Err err = NULL;

  phidom = phifor = phidomfor = 0.0;

  StartIndex = Get_Index(T0, maturity_rates, nbMat);
  EndIndex = Get_Index(Tfix, maturity_rates, nbMat);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      t1 = maturity_rates[i - 1];
    } else {
      /* First part */
      t1 = T0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      t2 = Tfix;
    } else {
      t2 = maturity_rates[i];
    }

    phidom += sig_curve_dom[i] * sig_curve_dom[i] *
              Phi_Func(2 * lda_dom, Tfix, t1, t2);
    phifor += sig_curve_for[i] * sig_curve_for[i] *
              Phi_Func(2 * lda_for, Tfix, t1, t2);
    phidomfor += sig_curve_dom[i] * sig_curve_for[i] *
                 Phi_Func(lda_dom + lda_for, Tfix, t1, t2);
  }

  phidomfor *= correl_dom_for;

  lamdom = (exp(-lda_dom * Tval1) - exp(-lda_dom * Tval2)) / lda_dom;
  lamfor = (exp(-lda_for * Tval1) - exp(-lda_for * Tval2)) / lda_for;

  *varff = phidom * lamdom * lamdom + phifor * lamfor * lamfor -
           2 * phidomfor * lamdom * lamfor;

  return err;
}

/*	IR Phi */
Err Fx3DtsPhi(double maturity, double *maturity_sig, long nbSig,
              double *sig_curve, double lda, double *phi) {
  double sig;
  double T1, T2;
  double var;
  double var_partial;
  int i;
  long StartIndex, EndIndex;

  if (maturity == 0) {
    *phi = 0.0;
    return NULL;
  }

  StartIndex = Get_Index(0.0, maturity_sig, nbSig);
  EndIndex = Get_Index(maturity, maturity_sig, nbSig);

  var = 0.0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_sig[i - 1];
    } else {
      /* First part */
      T1 = 0.0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = maturity;
    } else {
      T2 = maturity_sig[i];
    }

    sig = sig_curve[i];
    var_partial = sig * sig * Phi_Func(2 * lda, maturity, T1, T2);

    var += var_partial;
  }

  *phi = var;

  return NULL;
}

/*	Fx Implied volatility */
/*  The rates maturity dates have to be merged ! */
Err Fx3DtsImpliedVol(double opt_maturity, double start_date, double end_date,
                     double *maturity_rates, long nbMat, double *sig_curve_dom,
                     double lda_dom, double *sig_curve_for, double lda_for,
                     double *maturity_fx, double *sig_curve_fx, long nbrMat_fx,
                     double correl_dom_for, double correl_dom_fx,
                     double correl_for_fx, double *fx_vol) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2;
  double var;
  double var_partial;
  int i, j;
  long StartIndex, EndIndex, StartIndex2, EndIndex2;
  Err err = NULL;

  if (start_date > end_date) {
    err = "end_date before start_date in Fx3DtsImpliedVol";
    return err;
  }

  if (end_date == 0) {
    (*fx_vol) = 0;
    return err;
  }

  StartIndex = Get_Index(start_date, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(end_date, maturity_fx, nbrMat_fx);

  var = 0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = end_date;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }
      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_fx, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  }

  if (fabs(end_date - start_date) > 1.0e-08) {
    *fx_vol = sqrt(var / (end_date - start_date));
  } else {
    *fx_vol = 0.0;
  }

  return err;
}

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DtsCalibration(double *exercise_opt, double *maturity_opt,
                      double *vol_opt, long nbrOpt, double *maturity_rates,
                      long nbrMat, double *sig_curve_dom, double lda_dom,
                      double *sig_curve_for, double lda_for,
                      double correl_dom_for, double correl_dom_fx,
                      double correl_for_fx, double **fx_vol_curve)

{
  double a, b, c, a_part, b_part, c_part, c2, delta;
  double sig_dom, sig_for;
  double cumvar, vol;
  double T1, T2, t1, t2;
  double VolMin;
  double val_opt;
  int idxopt, i;
  long StartIndex, EndIndex;

  Err err = NULL;

  /* loop on the number of options */
  (*fx_vol_curve) = NULL;
  (*fx_vol_curve) = (double *)calloc(nbrOpt, sizeof(double));
  if (!(*fx_vol_curve)) {
    err = "Memory allocation error in Fx3DtsCalibration";
    goto FREE_RETURN;
  }

  for (idxopt = 0; idxopt <= nbrOpt - 1; idxopt++) {
    T2 = exercise_opt[idxopt];
    val_opt = maturity_opt[idxopt];

    /* Get the cumulated variance till the previous maturity T1*/
    if (idxopt > 0) {
      T1 = exercise_opt[idxopt - 1];

      err = Fx3DtsImpliedVol(
          val_opt, 0, T1, maturity_rates, nbrMat, sig_curve_dom, lda_dom,
          sig_curve_for, lda_for, exercise_opt, *fx_vol_curve, nbrOpt,
          correl_dom_for, correl_dom_fx, correl_for_fx, &vol);

      if (err) {
        goto FREE_RETURN;
      }

      cumvar = vol * vol * T1;
    } else {
      T1 = 0;
      cumvar = 0;
    }

    /* Calculate the coefficient of the last cumulative of the variance between
     * T1 and T2 */

    a = b = c2 = 0;
    StartIndex = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex = Get_Index(T2, maturity_rates, nbrMat);

    for (i = StartIndex; i < EndIndex + 1; i++) {
      if (i > StartIndex) {
        t1 = maturity_rates[i - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (i == EndIndex || StartIndex == EndIndex) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[i];
      }

      sig_dom = sig_curve_dom[i];
      sig_for = sig_curve_for[i];

      err = Coefs_Partial_Var(val_opt, t1, t2, sig_dom, lda_dom, sig_for,
                              lda_for, correl_dom_for, correl_dom_fx,
                              correl_for_fx, &a_part, &b_part, &c_part);

      if (err) {
        goto FREE_RETURN;
      }

      a += a_part;
      b += b_part;
      c2 += c_part;
    }

    c = c2 + cumvar;
    /* substract value to match */
    c -= vol_opt[idxopt] * vol_opt[idxopt] * T2;

    /* just solve the second order equation */
    delta = b * b - 4 * a * c;
    /* delta < 0 or solutions are negatives */
    if ((delta < 0) || ((c > 0) && (b > 0))) {
      VolMin = sqrt((c2 + cumvar - b * b / (4 * a)) / val_opt) * 100;

      err = serror("Cannot find solution to match the option %d (its vol "
                   "should be > %.2f %%)",
                   idxopt + 1, VolMin);

      goto FREE_RETURN;
    }

    /* if there is two positive solutions we are taking the smallest one */
    /* it would be easier to match the next one */
    if ((c > 0) && (b < 0)) {
      (*fx_vol_curve)[idxopt] = (-b - sqrt(delta)) / (2. * a);
    } else {
      (*fx_vol_curve)[idxopt] = (-b + sqrt(delta)) / (2. * a);
    }
  }

FREE_RETURN:

  if (err) {
    if (*fx_vol_curve) {
      free(*fx_vol_curve);
      *fx_vol_curve = NULL;
    }
  }

  return err;
}

/*	Forward fx */
Err Fx3DtsFwdFx(SrtUndPtr dom_und, SrtUndPtr for_und, SrtUndPtr fx_und,
                long maturity_date, double *fwd_fx) {
  char *dom_yc, *for_yc;
  double dom_dc, for_dc, spot;
  long spot_date;
  Err err = NULL;

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    return err;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    return err;
  }

  spot_date = get_spotdate_from_underlying(fx_und);

  if (maturity_date >= spot_date) {
    dom_dc = swp_f_df(spot_date, maturity_date, dom_yc);
    for_dc = swp_f_df(spot_date, maturity_date, for_yc);
  } else {
    dom_dc = 1.0 / swp_f_df(maturity_date, spot_date, dom_yc);
    for_dc = 1.0 / swp_f_df(maturity_date, spot_date, for_yc);
  }

  spot = get_spot_from_fxund(fx_und);

  (*fwd_fx) = spot * for_dc / dom_dc;

  return err;
}

/*	Implied vol direct from underlying */
Err Fx3DImpliedVol(char *fx_underlying, double val_time, double start_time,
                   double end_time, double *vol) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx,
      nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL;
  double *tau_date_dom = NULL, *tau_dom = NULL;
  double *sigma_date_for = NULL, *sigma_for = NULL;
  double *tau_date_for = NULL, *tau_for = NULL;
  double *sigma_date_fx = NULL, *sigma_fx = NULL;
  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double lda_dom, lda_for;
  double *sig_dom = NULL, *sig_for = NULL, *merge_dates = NULL;

  Err err = NULL;

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  err = Fx3DtsImpliedVol(val_time, start_time, end_time, merge_dates,
                         nb_merge_dates, sig_dom, lda_dom, sig_for, lda_for,
                         sigma_date_fx, sigma_fx, sigma_n_fx, correl_dom_for,
                         correl_dom_fx, correl_for_fx, vol);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

/*	Forward fx direct from underlying */
Err Fx3DFwdFx(char *fx_underlying, long maturity_date, double *fwd_fx) {
  SrtUndPtr und, dom_und, for_und;
  char *domname, *forname;
  Err err = NULL;

  und = lookup_und(fx_underlying);
  if (!und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  domname = get_domname_from_fxund(und);
  forname = get_forname_from_fxund(und);

  dom_und = lookup_und(domname);
  for_und = lookup_und(forname);

  err = Fx3DtsFwdFx(dom_und, for_und, und, maturity_date, fwd_fx);

  return err;
}

/*	Fx calibration direct from underlying */
Err Fx3DCalibration(char *dom_underlying, char *for_underlying,
                    double correl_dom_for, double correl_dom_fx,
                    double correl_for_fx, double *exercise_opt,
                    double *maturity_opt, double *vol_opt, long nbropt,
                    double **fx_vol_curve) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;
  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  err = Fx3DtsCalibration(exercise_opt, maturity_opt, vol_opt, nbropt,
                          merge_dates, nb_merge_dates, sig_dom, lda_dom,
                          sig_for, lda_for, correl_dom_for, correl_dom_fx,
                          correl_for_fx, fx_vol_curve);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

/*	European option price */
Err Fx3DFxOption(char *fx_underlying, double strike, long fix_date,
                 long val_date, long pay_date, char *pay_reveceive,
                 double *price) {
  double fx_vol, df, fx_fwd;
  double fix_time, val_time, pay_time;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  char *dom_yc;
  long today;
  Err err = NULL;

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    return err;
  }

  today = get_today_from_underlying(dom_und);
  fix_time = (fix_date - today) / 365.0;
  val_time = (val_date - today) / 365.0;
  pay_time = (pay_date - today) / 365.0;

  df = swp_f_df(today, pay_date, dom_yc);

  err = Fx3DImpliedVol_corr(fx_underlying, val_time, 0, fix_time, &fx_vol);
  if (err) {
    return err;
  }

  err =
      Fx3DFwdTpay_corr(fx_underlying, fix_time, val_time, pay_time, 1, &fx_fwd);
  if (err) {
    return err;
  }

  err = OptBlkSch(fx_fwd, strike, fx_vol, fix_time, df, pay_reveceive,
                  "PREMIUM", price);

  return err;
}

/*	Forward volatility of spot fx */
Err Fx3DSpotVol(char *fx_underlying, double maturity, double *fx_vol)

{

  SrtUndPtr fx_und;
  Err err = NULL;
  TermStruct *ts;

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  ts = get_ts_from_irund(fx_und);

  (*fx_vol) = find_fx_sig(maturity, ts);

  return err;
}

/*	Useful functions designed to get term structures */

Err Get_LGM_TermStructure(char *underlying, double **sigma_date, double **sigma,
                          long *sigma_n, double **tau_date, double **tau,
                          long *tau_n)

{
  SrtUndPtr und;
  TermStruct *ts;
  double *beta = NULL, *vovol = NULL, *rho = NULL, *meanvol = NULL;
  long today;
  int i;
  Err err;

  und = lookup_und(underlying);
  if (!und) {
    return serror("Couldn't fin underlying named %s", underlying);
  }

  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", underlying);
  }

  if (get_mdltype_from_irund(und) != LGM ||
      get_mdldim_from_irund(und) != ONE_FAC) {
    return serror("Underlying %s is not of type LGM 1F", underlying);
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = srt_f_display_IRM_OneFac_TermStruct(ts, sigma_date, sigma, &beta,
                                            &vovol, &rho, &meanvol, sigma_n,
                                            tau_date, tau, tau_n);

  for (i = 0; i < *sigma_n; i++) {
    (*sigma_date)[i] = ((*sigma_date)[i] - today) / 365.0;
  }
  for (i = 0; i < *tau_n; i++) {
    (*tau_date)[i] = ((*tau_date)[i] - today) / 365.0;
  }

  if (beta)
    free(beta);
  if (vovol)
    free(vovol);
  if (rho)
    free(rho);
  if (meanvol)
    free(meanvol);

  return err;
}

/* Get the LGM term structure and check it has a constant tau */

Err Get_LGM_TermStructure2(char *underlying, double **sigma_time,
                           double **sigma, long *sigma_n, double *fixed_tau)

{
  SrtUndPtr und;
  TermStruct *ts;
  long today;
  int i;

  double *beta = NULL, *vovol = NULL, *rho = NULL, *meanvol = NULL,
         *TauDates = NULL, *Tau = NULL;

  long nbTau;
  double tau;

  Err err = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't fin underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) != LGM ||
      get_mdldim_from_irund(und) != ONE_FAC) {
    err = serror("Underlying %s is not of type LGM", underlying);
    goto FREE_RETURN;
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = srt_f_display_IRM_OneFac_TermStruct(ts, sigma_time, sigma, &beta,
                                            &vovol, &rho, &meanvol, sigma_n,
                                            &TauDates, &Tau, &nbTau);
  if (err) {
    goto FREE_RETURN;
  }

  /* Now transform dates into maturities and check tau	*/
  for (i = 0; i < *sigma_n; i++) {
    (*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
  }

  tau = Tau[0];

  for (i = 1; i < nbTau; i++) {
    if (fabs(Tau[i] - tau) > 1.0E-08) {
      err = "Tau must be constant in the LGM term struct";
      goto FREE_RETURN;
    }
  }

  (*fixed_tau) = tau;

FREE_RETURN:

  if (beta)
    free(beta);
  if (vovol)
    free(vovol);
  if (rho)
    free(rho);
  if (meanvol)
    free(meanvol);
  if (TauDates)
    free(TauDates);
  if (Tau)
    free(Tau);

  return err;
}

/* Get the LGM 2F term structure and check it has:
   - a constant tau
   - a constant alpha        , gamma and rho */

Err Get_LGM2F_TermStructure(char *underlying, double **sigma_time,
                            double **sigma, long *sigma_n, double *fixed_tau,
                            double *fixed_alpha, double *fixed_gamma,
                            double *fixed_rho)

{
  SrtUndPtr und;
  TermStruct *ts;
  long today;
  int i;

  double *beta1 = NULL, *sigma2 = NULL, *beta2 = NULL, *Rho = NULL,
         *TauDates = NULL, *Tau1 = NULL, *Tau2 = NULL;

  long nbTau;
  double alpha, rho, tau1, tau2;

  Err err = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't fin underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) != LGM ||
      get_mdldim_from_irund(und) != TWO_FAC) {
    err = serror("Underlying %s is not of type LGM 2F", underlying);
    goto FREE_RETURN;
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = srt_f_display_IRM_TwoFac_TermStruct(ts, sigma_time, sigma, &beta1,
                                            &sigma2, &beta2, &Rho, sigma_n,
                                            &TauDates, &Tau1, &Tau2, &nbTau);
  if (err) {
    goto FREE_RETURN;
  }

  /* Now transform dates into maturities and check alpha        , gamma        ,
   * rho and tau
   */

  alpha = sigma2[0] / (*sigma)[0];
  rho = Rho[0];
  (*sigma_time)[0] = ((*sigma_time)[0] - today) / 365.0;

  for (i = 1; i < *sigma_n; i++) {
    (*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
    if (fabs(sigma2[i] / (*sigma)[i] - alpha) > 1.0E-08) {
      err = "Alpha must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
    if (fabs(Rho[i] - rho) > 1.0E-08) {
      err = "Rho must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
  }

  tau1 = Tau1[0];
  tau2 = Tau2[0];

  for (i = 1; i < nbTau; i++) {
    if (fabs(Tau1[i] - tau1) > 1.0E-08 || fabs(Tau2[i] - tau2) > 1.0E-08) {
      err = "Tau must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
  }

  (*fixed_tau) = tau1;
  (*fixed_alpha) = alpha;
  (*fixed_gamma) = 1.0 / tau2 - 1.0 / tau1;
  (*fixed_rho) = rho;

FREE_RETURN:

  if (beta1)
    free(beta1);
  if (sigma2)
    free(sigma2);
  if (beta2)
    free(beta2);
  if (Rho)
    free(Rho);
  if (TauDates)
    free(TauDates);
  if (Tau1)
    free(Tau1);
  if (Tau2)
    free(Tau2);

  return err;
}

Err Get_LGM2F_TermStructure_ts(char *underlying, double **sigma_time,
                               double **sigma, long *sigma_n, double *fixed_tau,
                               double *fixed_alpha, double *fixed_gamma,
                               double *fixed_rho, int *tau_ts)

{
  SrtUndPtr und;
  TermStruct *ts;
  long today;
  int i;

  double *beta1 = NULL, *sigma2 = NULL, *beta2 = NULL, *Rho = NULL,
         *TauDates = NULL, *Tau1 = NULL, *Tau2 = NULL;

  long nbTau;
  double alpha, rho, tau1, tau2;

  Err err = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't fin underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) != LGM ||
      get_mdldim_from_irund(und) != TWO_FAC) {
    err = serror("Underlying %s is not of type LGM 2F", underlying);
    goto FREE_RETURN;
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = srt_f_display_IRM_TwoFac_TermStruct(ts, sigma_time, sigma, &beta1,
                                            &sigma2, &beta2, &Rho, sigma_n,
                                            &TauDates, &Tau1, &Tau2, &nbTau);
  if (err) {
    goto FREE_RETURN;
  }

  /* Now transform dates into maturities and check alpha        , gamma        ,
   * rho and tau
   */

  alpha = sigma2[0] / (*sigma)[0];
  rho = Rho[0];
  (*sigma_time)[0] = ((*sigma_time)[0] - today) / 365.0;

  for (i = 1; i < *sigma_n; i++) {
    (*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
    if (fabs(sigma2[i] / (*sigma)[i] - alpha) > 1.0E-08) {
      err = "Alpha must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
    if (fabs(Rho[i] - rho) > 1.0E-08) {
      err = "Rho must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
  }

  tau1 = Tau1[0];
  tau2 = Tau2[0];

  *tau_ts = 0;
  for (i = 1; i < nbTau; i++) {
    if (fabs(Tau1[i] - tau1) > 1.0E-08 || fabs(Tau2[i] - tau2) > 1.0E-08) {

      // err = "Tau must be constant in the LGM2F term struct";
      // goto FREE_RETURN;
      *tau_ts = 1;
      break;
    }
  }

  (*fixed_tau) = tau1;
  (*fixed_alpha) = alpha;
  (*fixed_gamma) = 1.0 / tau2 - 1.0 / tau1;
  (*fixed_rho) = rho;

FREE_RETURN:

  if (beta1)
    free(beta1);
  if (sigma2)
    free(sigma2);
  if (beta2)
    free(beta2);
  if (Rho)
    free(Rho);
  if (TauDates)
    free(TauDates);
  if (Tau1)
    free(Tau1);
  if (Tau2)
    free(Tau2);

  return err;
}

/* Get the LGM 2F term structure and check it has:
   - a constant alpha        , gamma and rho
   Then merge Tau and Vol term structure		*/

Err Get_LGM2F_TermStructure2(char *underlying, double **sigma,
                             double **sigma_time, long *nb_sigma,
                             double **lambda, double **lambda_time,
                             long *nb_lambda, double *fixed_alpha,
                             double *fixed_gamma, double *fixed_rho) {
  SrtUndPtr und;
  TermStruct *ts;
  long today;
  int i, j, k;

  double *beta1 = NULL, *Sigma2 = NULL, *beta2 = NULL, *Rho = NULL,
         *Tau2 = NULL;

  double alpha, gamma, rho;

  Err err = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't fin underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) != LGM ||
      get_mdldim_from_irund(und) != TWO_FAC) {
    err = serror("Underlying %s is not of type LGM 2F", underlying);
    goto FREE_RETURN;
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = srt_f_display_IRM_TwoFac_TermStruct(
      ts, sigma_time, sigma, &beta1, &Sigma2, &beta2, &Rho, nb_sigma,
      lambda_time, lambda, &Tau2, nb_lambda);
  if (err) {
    goto FREE_RETURN;
  }

  /* Now transform dates into maturities and check alpha        , gamma        ,
   * and rho */

  alpha = Sigma2[0] / (*sigma)[0];
  rho = Rho[0];

  for (i = 0; i < *nb_sigma; i++) {
    (*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
    if (fabs(Sigma2[i] / (*sigma)[i] - alpha) > 1.0E-08) {
      err = "Alpha must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
    if (fabs(Rho[i] - rho) > 1.0E-08) {
      err = "Rho must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }

    /* Remove redondant points */
    j = 1;

    while (i + j < *nb_sigma && fabs((*sigma)[i + j] - (*sigma)[i]) < 1.0E-10) {
      j++;
    }
    j--;

    *nb_sigma -= j;

    if (j > 0) {
      (*sigma_time)[i] = ((*sigma_time)[i + j] - today) / 365.0;

      for (k = 1; k < *nb_sigma - i; k++) {
        (*sigma_time)[i + k] = (*sigma_time)[i + j + k];
        (*sigma)[i + k] = (*sigma)[i + j + k];
        Sigma2[i + k] = Sigma2[i + j + k];
        Rho[i + k] = Rho[i + j + k];
      }
    }
  }

  gamma = 1.0 / Tau2[0] - 1.0 / (*lambda)[0];

  for (i = 0; i < *nb_lambda; i++) {
    (*lambda_time)[i] = ((*lambda_time)[i] - today) / 365.0;
    if (fabs(1.0 / Tau2[i] - 1.0 / (*lambda)[i] - gamma) > 1.0E-08) {
      err = "Gamma must be constant in the LGM2F term struct";
      goto FREE_RETURN;
    }
    (*lambda)[i] = 1.0 / (*lambda)[i];

    /* Remove redondant points */
    j = 1;

    while (i + j < *nb_lambda &&
           fabs((*lambda)[i + j] - 1.0 / (*lambda)[i]) < 1.0E-10) {
      j++;
    }
    j--;

    *nb_lambda -= j;

    if (j > 0) {
      (*lambda_time)[i] = ((*lambda_time)[i + j] - today) / 365.0;

      for (k = 1; k < *nb_lambda - i; k++) {
        (*lambda_time)[i + k] = (*lambda_time)[i + j + k];
        (*lambda)[i + k] = (*lambda)[i + j + k];
        Tau2[i + k] = Tau2[i + j + k];
      }
    }
  }

  (*fixed_alpha) = alpha;
  (*fixed_gamma) = gamma;
  (*fixed_rho) = rho;

FREE_RETURN:

  if (beta1)
    free(beta1);
  if (Sigma2)
    free(Sigma2);
  if (beta2)
    free(beta2);
  if (Rho)
    free(Rho);
  if (Tau2)
    free(Tau2);

  return err;
}

/*	IR Phi */
Err LGM2FDtsPhi(double maturity, double *maturity_sig, long nbSig,
                double *sig_curve, double *maturity_lam, long nbLam,
                double *lam_curve, double alpha, double gamma, double rho,
                double *phi1, double *phi2, double *phi12) {
  double T1, T2, T3, T4;
  int i, j;
  long StartIndex, EndIndex;
  long StartIndex2, EndIndex2;
  double alpha2, gamma2, dt;

  LGMSVSolFunc Phi1, Phi2, Phi12;
  double vector[2];

  vector[0] = 0.0;
  vector[1] = 0.0;

  if (maturity == 0) {
    *phi1 = 0.0;
    *phi2 = 0.0;
    *phi12 = 0.0;
    return NULL;
  }

  alpha2 = alpha * alpha;
  gamma2 = 2.0 * gamma;

  Phi1.bIsft1 = 1;
  Phi1.b = 0.0;
  Phi1.bIsgt1 = 1;
  Phi1.cxx = 0.0;
  Phi1.bIsht1 = 1;
  Phi1.dXt1 = 0.0;

  Phi2.bIsft1 = 1;
  Phi2.b = 0.0;
  Phi2.bIsgt1 = 1;
  Phi2.cxx = 0.0;
  Phi2.bIsht1 = 1;
  Phi2.dXt1 = 0.0;

  Phi12.bIsft1 = 1;
  Phi12.b = 0.0;
  Phi12.bIsgt1 = 1;
  Phi12.cxx = 0.0;
  Phi12.bIsht1 = 1;
  Phi12.dXt1 = 0.0;

  StartIndex = Get_Index(0.0, maturity_sig, nbSig);
  EndIndex = Get_Index(maturity, maturity_sig, nbSig);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_sig[i - 1];
    } else {
      /* First part */
      T1 = 0.0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = maturity;
    } else {
      T2 = maturity_sig[i];
    }

    Phi1.a = sig_curve[i] * sig_curve[i];
    Phi2.a = alpha2 * Phi1.a;
    Phi12.a = alpha * Phi1.a;

    StartIndex2 = Get_Index(T1, maturity_lam, nbLam);
    EndIndex2 = Get_Index(T2, maturity_lam, nbLam);

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        T3 = maturity_lam[j - 1];
      } else {
        /* First part */
        T3 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        T4 = T2;
      } else {
        T4 = maturity_lam[j];
      }

      dt = (T4 - T3);

      Phi1.dLambda = 2.0 * lam_curve[j];
      Phi2.dLambda = Phi1.dLambda + gamma2;
      Phi12.dLambda = Phi1.dLambda + gamma;

      Phi1.dXt1 = LGMSVFuncValue(Phi1, dt, vector, 0);
      Phi2.dXt1 = LGMSVFuncValue(Phi2, dt, vector, 0);
      Phi12.dXt1 = LGMSVFuncValue(Phi12, dt, vector, 0);
    }
  }

  *phi1 = Phi1.dXt1;
  *phi2 = Phi2.dXt1;
  *phi12 = Phi12.dXt1;

  return NULL;
}

Err Get_FX_StochRate_TermStructures(
    char *underlying, double **sigma_date_dom, double **sigma_dom,
    long *sigma_n_dom, double **tau_date_dom, double **tau_dom, long *tau_n_dom,
    double **sigma_date_for, double **sigma_for, long *sigma_n_for,
    double **tau_date_for, double **tau_for, long *tau_n_for,
    double **sigma_date_fx, double **sigma_fx, long *sigma_n_fx,
    double *correl_dom_for, double *correl_dom_fx, double *correl_for_fx) {
  SrtUndPtr und;
  char *domname, *forname;
  long corr_date_n;
  double *corr_date = NULL, *corr = NULL;
  long today;
  int i;

  Err err = NULL;

  *sigma_date_dom = NULL;
  *sigma_dom = NULL;
  *tau_date_dom = NULL;
  *tau_dom = NULL;
  *sigma_date_for = NULL;
  *sigma_for = NULL;
  *tau_date_for = NULL;
  *tau_for = NULL;
  *sigma_date_fx = NULL;
  *sigma_fx = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't fin underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", underlying);
    goto FREE_RETURN;
  }

  domname = get_domname_from_fxund(und);
  err = Get_LGM_TermStructure(domname, sigma_date_dom, sigma_dom, sigma_n_dom,
                              tau_date_dom, tau_dom, tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  forname = get_forname_from_fxund(und);
  err = Get_LGM_TermStructure(forname, sigma_date_for, sigma_for, sigma_n_for,
                              tau_date_for, tau_for, tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = srt_f_display_FX_TermStruct(underlying, sigma_n_fx, sigma_date_fx,
                                    sigma_fx, &corr_date_n, &corr_date, &corr);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(und);
  for (i = 0; i < *sigma_n_fx; i++) {
    (*sigma_date_fx)[i] = ((*sigma_date_fx)[i] - today) / 365.0;
  }

  if (corr_date_n != 1) {
    err = "no correlation term structure allowed";
    goto FREE_RETURN;
  }

  *correl_dom_for = corr[0];
  *correl_dom_fx = corr[1];
  *correl_for_fx = corr[2];

FREE_RETURN:

  if (err) {
    if (*sigma_date_dom)
      free(*sigma_date_dom);
    *sigma_date_dom = NULL;
    if (*sigma_dom)
      free(*sigma_dom);
    *sigma_dom = NULL;

    if (*sigma_date_for)
      free(*sigma_date_for);
    *sigma_date_for = NULL;
    if (*sigma_for)
      free(*sigma_for);
    *sigma_for = NULL;

    if (*tau_date_dom)
      free(*tau_date_dom);
    *tau_date_dom = NULL;
    if (*tau_dom)
      free(*tau_dom);
    *tau_dom = NULL;

    if (*tau_date_for)
      free(*tau_date_for);
    *tau_date_for = NULL;
    if (*tau_for)
      free(*tau_for);
    *tau_for = NULL;

    if (*sigma_date_fx)
      free(*sigma_date_fx);
    *sigma_date_fx = NULL;
    if (*sigma_fx)
      free(*sigma_fx);
    *sigma_fx = NULL;
  }

  if (corr_date)
    free(corr_date);
  if (corr)
    free(corr);

  return err;
}

/*	Compute the expectation of log ( S (T) / S (0) )
                under a variety of measures */

/*	Calculate adjustment between Q-Tpay and Q-beta        , such that
        Q-beta expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay expect log ( S (
   Tfix ) / S ( 0 ) ) + adj */
Err Fx3DtsFwdBetaAdjustment(double T0,   /*	Forward start date */
                            double Tval, /*	Value date of the forward */
                            double Tpay, /*	Pay date of the forward */
                            double Tfix, /*	Fix date of the forward */
                            /*	Model data */
                            double *maturity_rates, long nbrMat,
                            double *sig_curve_dom, double lda_dom,
                            double *sig_curve_for, double lda_for,
                            double *maturity_fx, double *sig_curve_fx,
                            long nbrMat_fx, double correl_dom_for,
                            double correl_dom_fx, double correl_for_fx,
                            /*	Result */
                            double *adjust) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2;
  double adj;
  double adjust_partial;
  int i, j;
  long StartIndex, EndIndex, StartIndex2, EndIndex2;
  Err err = NULL;

  if (T0 > Tfix) {
    err = "end_date before start_date in Fx3DtsFwdBetaAdjustment";
    return err;
  }

  if (Tfix == 0) {
    (*adjust) = 0;
    return err;
  }

  StartIndex = Get_Index(T0, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(Tfix, maturity_fx, nbrMat_fx);

  adj = 0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = T0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = Tfix;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbrMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }

      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      adjust_partial =
          correl_dom_fx * sig_dom * Etha_Func(lda_dom, Tpay, t1, t2) * sig_fx -
          correl_dom_for * sig_dom * sig_for *
              Psi2_Func(lda_for, lda_dom, Tval, Tpay, t1, t2) +
          sig_dom * sig_dom * Psi2_Func(lda_dom, lda_dom, Tval, Tpay, t1, t2);

      adj += adjust_partial;
    }
  }

  *adjust = adj;

  return err;
}

/*	Calculate adjustment between Q-Tpay1 and Q-Tpay2        , such that
        Q-Tpay2 expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay1 expect log ( S (
   Tfix ) / S ( 0 ) ) + adj */
Err Fx3DtsFwdPayAdjustment(double T0,    /*	Forward start date */
                           double Tval,  /*	Value date of the forward */
                           double Tpay1, /*	Original pay date */
                           double Tpay2, /*	Adjusted pay date */
                           double Tfix,  /*	Fix date of the forward */
                           /*	Model data */
                           double *maturity_rates, long nbrMat,
                           double *sig_curve_dom, double lda_dom,
                           double *sig_curve_for, double lda_for,
                           double *maturity_fx, double *sig_curve_fx,
                           long nbrMat_fx, double correl_dom_for,
                           double correl_dom_fx, double correl_for_fx,
                           /*	Result */
                           double *adjust) {
  double adj1, adj2;
  Err err = NULL;

  err = Fx3DtsFwdBetaAdjustment(
      T0, Tval, Tpay1, Tfix, maturity_rates, nbrMat, sig_curve_dom, lda_dom,
      sig_curve_for, lda_for, maturity_fx, sig_curve_fx, nbrMat_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, &adj1);

  if (err) {
    return err;
  }

  err = Fx3DtsFwdBetaAdjustment(
      T0, Tval, Tpay2, Tfix, maturity_rates, nbrMat, sig_curve_dom, lda_dom,
      sig_curve_for, lda_for, maturity_fx, sig_curve_fx, nbrMat_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, &adj2);

  *adjust = adj1 - adj2;

  return err;
}

/*	Calculate cumulative covariance of F ( mat1 ) and F ( mat2 ) between 0
 * and T */
Err Fx3DtsFwdCumCovar(double T0,    /*	Forward start date */
                      double Tval1, /*	Value date of the 1st forward */
                      double Tval2, /*	Value date of the 2nd forward */
                      double Tfix,  /*	Fix date of both forwards */
                      /*	Model data */
                      double *maturity_rates, long nbrMat,
                      double *sig_curve_dom, double lda_dom,
                      double *sig_curve_for, double lda_for,
                      double *maturity_fx, double *sig_curve_fx, long nbrMat_fx,
                      double correl_dom_for, double correl_dom_fx,
                      double correl_for_fx,
                      /*	Result */
                      double *adjust) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2;
  double covar;
  double covar_partial;
  int i, j;
  long StartIndex, EndIndex, StartIndex2, EndIndex2;
  Err err = NULL;

  if (T0 > Tfix) {
    err = "end_date before start_date in Fx3DtsFwdCumCovar";
    return err;
  }

  if (Tfix == 0) {
    (*adjust) = 0;
    return err;
  }

  StartIndex = Get_Index(T0, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(Tfix, maturity_fx, nbrMat_fx);

  covar = 0.0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = T0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = Tfix;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbrMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }

      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      covar_partial = sig_fx * sig_fx * (t2 - t1) +
                      sig_fx * (correl_dom_fx * sig_dom *
                                    (Etha_Func(lda_dom, Tval1, t1, t2) +
                                     Etha_Func(lda_dom, Tval2, t1, t2))

                                - correl_for_fx * sig_for *
                                      (Etha_Func(lda_for, Tval1, t1, t2) +
                                       Etha_Func(lda_for, Tval2, t1, t2))) +
                      sig_dom * sig_dom *
                          Psi2_Func(lda_dom, lda_dom, Tval1, Tval2, t1, t2) +
                      sig_for * sig_for *
                          Psi2_Func(lda_for, lda_for, Tval1, Tval2, t1, t2) -
                      correl_dom_for * sig_dom * sig_for *
                          (Psi2_Func(lda_dom, lda_for, Tval1, Tval2, t1, t2) +
                           Psi2_Func(lda_for, lda_dom, Tval1, Tval2, t1, t2));

      covar += covar_partial;
    }
  }

  *adjust = covar;

  return err;
}

/*	Q beta */
Err Fx3DFwdBeta(char *fx_underlying,
                /*	Forward fixing time */
                double Tfix,
                /*	Forward value time */
                double Tval,
                /*	0: expect [log Fx]        , 1: expect [Fx] */
                int log_flag,
                /*	Answer */
                double *expect) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL;

  double lda_dom, lda_for;

  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double adjustment, fx_vol, fx_fwd;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  long today, fix_date, val_date;
  double spot_fx;

  Err err = NULL;

  /*	Read model data */

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  today = get_today_from_underlying(dom_und);
  fix_date = (long)(today + Tfix * DAYS_IN_YEAR + 1.0e-08);
  val_date = (long)(today + Tval * DAYS_IN_YEAR + 1.0e-08);

  /*	Calculate fwd        , spot and implied vol */

  /*	Forward */
  err = Fx3DFwdFx(fx_underlying, val_date, &fx_fwd);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Spot */
  err = Fx3DFwdFx(fx_underlying, today, &spot_fx);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Implied vol */
  err = Fx3DImpliedVol(fx_underlying, Tval, 0.0, Tfix, &fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculate expectation of the log under Q-Tfix */

  *expect = log(fx_fwd / spot_fx) - 0.5 * fx_vol * fx_vol * Tfix;

  /*	Adjust for Q-beta */

  err = Fx3DtsFwdBetaAdjustment(
      0.0, Tval, Tval, Tfix, merge_dates, nb_merge_dates, sig_dom, lda_dom,
      sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx, correl_dom_for,
      correl_dom_fx, correl_for_fx, &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  if (log_flag) {
    *expect = spot_fx * exp(*expect + 0.5 * fx_vol * fx_vol * Tfix);
  }

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  return err;
}

/*	Q Tpay */
Err Fx3DFwdTpay(char *fx_underlying,
                /*	Forward fixing time */
                double Tfix,
                /*	Forward value time */
                double Tval,
                /*	Payment time */
                double Tpay,
                /*	0: expect [log Fx]        , 1: expect [Fx] */
                int log_flag,
                /*	Answer */
                double *expect) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL;

  double lda_dom, lda_for;

  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double adjustment, fx_vol, fx_fwd;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  long today, fix_date, val_date, pay_date;
  double spot_fx;

  Err err = NULL;

  /*	Read model data */

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  today = get_today_from_underlying(dom_und);
  fix_date = (long)(today + Tfix * DAYS_IN_YEAR + 1.0e-08);
  val_date = (long)(today + Tval * DAYS_IN_YEAR + 1.0e-08);
  pay_date = (long)(today + Tpay * DAYS_IN_YEAR + 1.0e-08);

  /*	Calculate fwd        , spot and implied vol */

  /*	Forward */
  err = Fx3DFwdFx(fx_underlying, val_date, &fx_fwd);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Spot */
  err = Fx3DFwdFx(fx_underlying, today, &spot_fx);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Implied vol */
  err = Fx3DImpliedVol(fx_underlying, Tval, 0.0, Tfix, &fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculate expectation of the log under Q-Tfix */

  *expect = log(fx_fwd / spot_fx) - 0.5 * fx_vol * fx_vol * Tfix;

  /*	Adjust for Q-Tpay */

  err = Fx3DtsFwdPayAdjustment(
      0.0, Tval, Tval, Tpay, Tfix, merge_dates, nb_merge_dates, sig_dom,
      lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  if (log_flag) {
    *expect = spot_fx * exp(*expect + 0.5 * fx_vol * fx_vol * Tfix);
  }

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  return err;
}

/*	Q F */
Err Fx3DFwdQf(char *fx_underlying,
              /*	F1 fixing time */
              double Tbarfix,
              /*	F1 value time */
              double Tbarval,
              /*	F2 fixing time */
              double Tfix,
              /*	F2 value time */
              double Tval,
              /*	F2 Payment time */
              double Tpay,
              /*	0: expect [log Fx]        , 1: expect [Fx] */
              int log_flag,
              /*	Answer */
              double *expect) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL;

  double lda_dom, lda_for;

  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double adjustment, fx_vol, fx_fwd;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  long today, bar_fix_date, bar_val_date, fix_date, val_date, pay_date;
  double spot_fx;

  Err err = NULL;

  /*	Read model data */

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't find underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  today = get_today_from_underlying(dom_und);
  bar_fix_date = (long)(today + Tbarfix * DAYS_IN_YEAR + 1.0e-08);
  bar_val_date = (long)(today + Tbarval * DAYS_IN_YEAR + 1.0e-08);
  fix_date = (long)(today + Tfix * DAYS_IN_YEAR + 1.0e-08);
  val_date = (long)(today + Tval * DAYS_IN_YEAR + 1.0e-08);
  pay_date = (long)(today + Tpay * DAYS_IN_YEAR + 1.0e-08);

  /*	Calculate fwd        , spot and implied vol */

  /*	Forward */
  err = Fx3DFwdFx(fx_underlying, bar_val_date, &fx_fwd);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Spot */
  err = Fx3DFwdFx(fx_underlying, today, &spot_fx);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Implied vol */
  err = Fx3DImpliedVol(fx_underlying, Tbarval, 0.0, Tbarfix, &fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculate expectation of the log under Q-Tbar */

  *expect = log(fx_fwd / spot_fx) - 0.5 * fx_vol * fx_vol * Tbarfix;

  /*	Adjust for Q-Tpay */

  err = Fx3DtsFwdPayAdjustment(
      0.0, Tbarval, Tbarval, Tpay, Tbarfix, merge_dates, nb_merge_dates,
      sig_dom, lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  /*	Adjust for Q-F */

  err = Fx3DtsFwdCumCovar(0.0, Tval, Tbarval, Tbarfix, merge_dates,
                          nb_merge_dates, sig_dom, lda_dom, sig_for, lda_for,
                          sigma_date_fx, sigma_fx, sigma_n_fx, correl_dom_for,
                          correl_dom_fx, correl_for_fx, &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  if (log_flag) {
    *expect = spot_fx * exp(*expect + 0.5 * fx_vol * fx_vol * Tbarfix);
  }

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  return err;
}

/*	Compute covar ( log ( S ( T1 ) )         , log ( S ( T2 ) ) ) */
Err Fx3DCovar(char *fx_underlying, double Tfix1, double Tval1, double Tfix2,
              double Tval2, double *covar) {
  SrtUndPtr fx_und;
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL;

  double lda_dom, lda_for;

  double correl_dom_for, correl_dom_fx, correl_for_fx;

  double vol01, vol02, vol12, var01, var02, var12, varff;

  double temp;

  Err err = NULL;

  /*	Read model data */

  if (Tfix1 > Tfix2) {
    temp = Tfix1;
    Tfix1 = Tfix2;
    Tfix2 = temp;

    temp = Tval1;
    Tval1 = Tval2;
    Tval2 = temp;
  }

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  err = Fx3DtsImpliedVol(Tval1, 0, Tfix1, merge_dates, nb_merge_dates, sig_dom,
                         lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx,
                         sigma_n_fx, correl_dom_for, correl_dom_fx,
                         correl_for_fx, &vol01);

  var01 = vol01 * vol01 * Tfix1;

  err = Fx3DtsImpliedVol(Tval2, 0, Tfix2, merge_dates, nb_merge_dates, sig_dom,
                         lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx,
                         sigma_n_fx, correl_dom_for, correl_dom_fx,
                         correl_for_fx, &vol02);

  var02 = vol02 * vol02 * Tfix2;

  err = Fx3DtsImpliedVol(Tval2, Tfix1, Tfix2, merge_dates, nb_merge_dates,
                         sig_dom, lda_dom, sig_for, lda_for, sigma_date_fx,
                         sigma_fx, sigma_n_fx, correl_dom_for, correl_dom_fx,
                         correl_for_fx, &vol12);

  var12 = vol12 * vol12 * (Tfix2 - Tfix1);

  err = Fx3DtsFwdFactVar(0.0, Tval1, Tval2, Tfix1, merge_dates, nb_merge_dates,
                         sig_dom, lda_dom, sig_for, lda_for, correl_dom_for,
                         &varff);

  var12 += varff;

  *covar = var01 + var02 - var12;
  *covar *= 0.5;

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  return err;
}

/* x = lambda_rates y = lambda_vol */
static double Psi3_Func(double x, double y, double t0, double T, double t1,
                        double t2) {
  double res;

  if (x == y) {
    res = 1.0 / y * (exp(-y * (t1 - t0)) - exp(-y * (t2 - t0)));
    res -= exp(-y * (T - t0)) * (t2 - t1);
    return res / x;
  } else {
    res = 1.0 / y * (exp(-y * (t1 - t0)) - exp(-y * (t2 - t0)));
    res -= 1.0 / (y - x) *
           (exp(-x * T + y * t0 - (y - x) * t1) -
            exp(-x * T + y * t0 - (y - x) * t2));
    return res / x;
  }
}

Err Fx3DtsImpliedVolExtrapol(double opt_maturity, double start_date,
                             double end_date, double *maturity_rates,
                             long nbMat, double *sig_curve_dom, double lda_dom,
                             double *sig_curve_for, double lda_for,
                             double *maturity_fx, double *sig_curve_fx,
                             long nbrMat_fx, double sig_inf, double lda_vol,
                             double correl_dom_for, double correl_dom_fx,
                             double correl_for_fx, double *fx_vol) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2;
  double var;
  double var_partial;
  int i, j;
  long StartIndex, EndIndex, StartIndex2, EndIndex2;
  double t0, sig0;
  Err err = NULL;

  t0 = maturity_fx[nbrMat_fx - 1];
  T2 = start_date;

  if (start_date > end_date) {
    err = "end_date before start_date in Fx3DtsImpliedVol";
    return err;
  }

  if (end_date == 0) {
    (*fx_vol) = 0;
    return err;
  }

  StartIndex = Get_Index(start_date, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(end_date, maturity_fx, nbrMat_fx);

  if (t0 < end_date) {
    EndIndex += 1;
  }

  var = 0;

  for (i = StartIndex; i < EndIndex; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    T2 = maturity_fx[i];

    StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }
      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_fx, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  }

  /* Last part */
  T1 = T2;
  T2 = end_date;

  if (T2 <= t0) {
    /* classical calculation */
    StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }
      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_fx, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  } else {
    /* include the extrapolation of vol fx */
    StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbMat);

    sig0 = sig_fx;

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }
      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_inf, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      var_partial +=
          2.0 * (sig0 - sig_inf) *
          (correl_dom_fx * sig_dom *
               Psi3_Func(lda_dom, lda_vol, t0, opt_maturity, t1, t2) -
           correl_for_fx * sig_for *
               Psi3_Func(lda_for, lda_vol, t0, opt_maturity, t1, t2));

      var_partial +=
          2.0 * sig_inf * (sig0 - sig_inf) * Phi_Func(-lda_vol, t0, t1, t2);
      var_partial += (sig0 - sig_inf) * (sig0 - sig_inf) *
                     Phi_Func(-2.0 * lda_vol, t0, t1, t2);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  }

  if (fabs(end_date - start_date) > 1.0e-08) {
    *fx_vol = sqrt(var / (end_date - start_date));
  } else {
    *fx_vol = 0.0;
  }

  return err;
}

Err Fx3DtsImpliedVolExtrapol_corr(
    double opt_maturity, double start_date, double end_date,
    double *maturity_rates, long nbMat, double *sig_curve_dom, double lda_dom,
    double *sig_curve_for, double lda_for, double *maturity_fx,
    double *sig_curve_fx, long nbrMat_fx, double sig_inf, double lda_vol,
    double *correl_times, double *correl_dom_for_ts, double *correl_dom_fx_ts,
    double *correl_for_fx_ts, long nb_correl, double *fx_vol) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2;
  double var;
  double var_partial;
  int i, j;
  long Index, StartIndex, EndIndex, StartIndex2, EndIndex2;
  double t0, sig0;
  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double *mat_rates_corr = NULL, *new_sig_curve_dom = NULL,
         *new_sig_curve_for = NULL, *new_correl_dom_for_ts = NULL,
         *new_correl_dom_fx_ts = NULL, *new_correl_for_fx_ts = NULL;
  long nb_merged;

  Err err = NULL;

  /* first merge rates and correl TS */
  nb_merged = nbMat;
  mat_rates_corr = calloc(nb_merged, sizeof(double));

  if (!mat_rates_corr) {
    err = "Memory allocation faillure (1) in Fx3DtsImpliedVolExtrapol_corr";
    goto FREE_RETURN;
  }

  memcpy(mat_rates_corr, maturity_rates, nb_merged * sizeof(double));

  num_f_concat_vector(&nb_merged, &mat_rates_corr, nb_correl, correl_times);
  num_f_unique_vector(&nb_merged, mat_rates_corr);
  num_f_sort_vector(nb_merged, mat_rates_corr);

  new_sig_curve_dom = dvector(0, nb_merged - 1);
  new_sig_curve_for = dvector(0, nb_merged - 1);
  new_correl_dom_for_ts = dvector(0, nb_merged - 1);
  new_correl_dom_fx_ts = dvector(0, nb_merged - 1);
  new_correl_for_fx_ts = dvector(0, nb_merged - 1);

  if (!new_sig_curve_dom || !new_sig_curve_for || !new_correl_dom_for_ts ||
      !new_correl_dom_fx_ts || !new_correl_for_fx_ts) {
    err = "Memory allocation faillure (2) in Fx3DtsImpliedVolExtrapol_corr";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_merged; i++) {
    Index = Get_Index(mat_rates_corr[i], maturity_rates, nbMat);
    new_sig_curve_dom[i] = sig_curve_dom[Index];
    new_sig_curve_for[i] = sig_curve_for[Index];

    Index = Get_Index(mat_rates_corr[i], correl_times, nb_correl);
    new_correl_dom_for_ts[i] = correl_dom_for_ts[Index];
    new_correl_dom_fx_ts[i] = correl_dom_fx_ts[Index];
    new_correl_for_fx_ts[i] = correl_for_fx_ts[Index];
  }

  t0 = maturity_fx[nbrMat_fx - 1];
  T2 = start_date;

  if (start_date > end_date) {
    err = "end_date before start_date in Fx3DtsImpliedVol";
    return err;
  }

  if (end_date == 0) {
    (*fx_vol) = 0;
    return err;
  }

  StartIndex = Get_Index(start_date, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(end_date, maturity_fx, nbrMat_fx);

  if (t0 < end_date) {
    EndIndex += 1;
  }

  var = 0;

  for (i = StartIndex; i < EndIndex; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    T2 = maturity_fx[i];

    StartIndex2 = Get_Index(T1, mat_rates_corr, nb_merged);
    EndIndex2 = Get_Index(T2, mat_rates_corr, nb_merged);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = mat_rates_corr[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = mat_rates_corr[j];
      }

      sig_dom = new_sig_curve_dom[j];
      sig_for = new_sig_curve_for[j];
      correl_dom_for = new_correl_dom_for_ts[j];
      correl_dom_fx = new_correl_dom_fx_ts[j];
      correl_for_fx = new_correl_for_fx_ts[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_fx, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  }

  /* Last part */
  T1 = T2;
  T2 = end_date;

  if (T2 <= t0) {
    /* classical calculation */
    StartIndex2 = Get_Index(T1, mat_rates_corr, nb_merged);
    EndIndex2 = Get_Index(T2, mat_rates_corr, nb_merged);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = mat_rates_corr[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = mat_rates_corr[j];
      }

      sig_dom = new_sig_curve_dom[j];
      sig_for = new_sig_curve_for[j];
      correl_dom_for = new_correl_dom_for_ts[j];
      correl_dom_fx = new_correl_dom_fx_ts[j];
      correl_for_fx = new_correl_for_fx_ts[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_fx, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  } else {
    /* include the extrapolation of vol fx */
    StartIndex2 = Get_Index(T1, mat_rates_corr, nb_merged);
    EndIndex2 = Get_Index(T2, mat_rates_corr, nb_merged);

    sig0 = sig_fx;

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = mat_rates_corr[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = mat_rates_corr[j];
      }

      sig_dom = new_sig_curve_dom[j];
      sig_for = new_sig_curve_for[j];
      correl_dom_for = new_correl_dom_for_ts[j];
      correl_dom_fx = new_correl_dom_fx_ts[j];
      correl_for_fx = new_correl_for_fx_ts[j];

      err = Partial_Var(opt_maturity, t1, t2, sig_dom, lda_dom, sig_for,
                        lda_for, sig_inf, correl_dom_for, correl_dom_fx,
                        correl_for_fx, &var_partial);

      var_partial +=
          2.0 * (sig0 - sig_inf) *
          (correl_dom_fx * sig_dom *
               Psi3_Func(lda_dom, lda_vol, t0, opt_maturity, t1, t2) -
           correl_for_fx * sig_for *
               Psi3_Func(lda_for, lda_vol, t0, opt_maturity, t1, t2));

      var_partial +=
          2.0 * sig_inf * (sig0 - sig_inf) * Phi_Func(-lda_vol, t0, t1, t2);
      var_partial += (sig0 - sig_inf) * (sig0 - sig_inf) *
                     Phi_Func(-2.0 * lda_vol, t0, t1, t2);

      if (err) {
        return err;
      }

      var += var_partial;
    }
  }

  if (fabs(end_date - start_date) > 1.0e-08) {
    *fx_vol = sqrt(var / (end_date - start_date));
  } else {
    *fx_vol = 0.0;
  }

FREE_RETURN:

  if (mat_rates_corr)
    free(mat_rates_corr);
  if (new_sig_curve_dom)
    free_dvector(new_sig_curve_dom, 0, nb_merged - 1);
  if (new_sig_curve_for)
    free_dvector(new_sig_curve_for, 0, nb_merged - 1);
  if (new_correl_dom_for_ts)
    free_dvector(new_correl_dom_for_ts, 0, nb_merged - 1);
  if (new_correl_dom_fx_ts)
    free_dvector(new_correl_dom_fx_ts, 0, nb_merged - 1);
  if (new_correl_for_fx_ts)
    free_dvector(new_correl_for_fx_ts, 0, nb_merged - 1);

  return err;
}

/* all the term structures have to be merged... */
Err Fx3DtsFxFwdCov(double option_maturity, double start_date, double end_date,
                   double *maturity, long nbMat, double *sig_curve_dom1,
                   double lda_dom1, double *sig_curve_for1, double lda_for1,
                   double *sig_curve_fx1, double *sig_curve_dom2,
                   double lda_dom2, double *sig_curve_for2, double lda_for2,
                   double *sig_curve_fx2, double **correlations,
                   double *covar) {
  double *correl_fxfx = NULL;
  int StartIndex, EndIndex, i;
  double T1, T2;
  double sig_1_dom, sig_1_for, sig_1;
  double sig_2_dom, sig_2_for, sig_2;
  Err err = NULL;

  correl_fxfx = (double *)calloc(9, sizeof(double));

  if (!correl_fxfx) {
    err = "Memory allocation failure in Fx3DtsFxFwdCov";
    return err;
  }

  StartIndex = Get_Index(start_date, maturity, nbMat);
  EndIndex = Get_Index(end_date, maturity, nbMat);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = end_date;
    } else {
      T2 = maturity[i];
    }

    sig_1_dom = sig_curve_dom1[i];
    sig_2_dom = sig_curve_dom2[i];

    sig_1_for = sig_curve_for1[i];
    sig_2_for = sig_curve_for2[i];

    sig_1 = sig_curve_fx1[i];
    sig_2 = sig_curve_fx2[i];

    /* Z_1 with Z_2 */
    correl_fxfx[0] += sig_1 * sig_2 * (T2 - T1);

    /* Z_1 with Dom_2 */
    correl_fxfx[1] +=
        sig_1 * sig_2_dom * (Etha_Func(lda_dom2, option_maturity, T1, T2));

    /* Z_1 with For_2 */
    correl_fxfx[2] +=
        sig_1 * sig_2_for * (Etha_Func(lda_for2, option_maturity, T1, T2));

    /* Dom_1 with Z_2 */
    correl_fxfx[3] +=
        sig_1_dom * sig_2 * (Etha_Func(lda_dom1, option_maturity, T1, T2));

    /* Dom_1 with Dom_2 */
    correl_fxfx[4] += sig_1_dom * sig_2_dom *
                      (Psi_Func(lda_dom1, lda_dom2, option_maturity, T1, T2));

    /* Dom_1 with For_2 */
    correl_fxfx[5] += sig_1_dom * sig_2_for *
                      (Psi_Func(lda_dom1, lda_for2, option_maturity, T1, T2));

    /* For_1 with Z_2 */
    correl_fxfx[6] +=
        sig_1_for * sig_2 * (Etha_Func(lda_for1, option_maturity, T1, T2));

    /* For_1 with Dom_2 */
    correl_fxfx[7] += sig_1_for * sig_2_dom *
                      (Psi_Func(lda_for1, lda_dom2, option_maturity, T1, T2));

    /* For_1 with For_2 */
    correl_fxfx[8] += sig_1_for * sig_2_for *
                      (Psi_Func(lda_for1, lda_for2, option_maturity, T1, T2));
  }

  *covar = correlations[2][5] * correl_fxfx[0] +
           correlations[2][3] * correl_fxfx[1] -
           correlations[2][4] * correl_fxfx[2] +
           correlations[0][5] * correl_fxfx[3] +
           correlations[0][3] * correl_fxfx[4] -
           correlations[0][4] * correl_fxfx[5] -
           correlations[1][5] * correl_fxfx[6] -
           correlations[1][3] * correl_fxfx[7] +
           correlations[1][4] * correl_fxfx[8];

  if (correl_fxfx) {
    free(correl_fxfx);
  }

  return err;
}

/* calculates any correlation within the 3F framework */
/* can be fx/fx fx/lgm lgm/lgm */
/* put NULL vector when not needed */
/* all the term structures have to be merged... */

Err Fx3DtsFxFwdCov_corr(double option_maturity1, double option_maturity2,
                        double start_date, double end_date, double *maturity,
                        long nbMat,
                        /* FX1 or LGM1 */
                        double *sig_curve_dom1, double lda_dom1,
                        double *sig_curve_for1, double lda_for1,
                        double *sig_curve_fx1,
                        /* FX2 or LGM2 */
                        double *sig_curve_dom2, double lda_dom2,
                        double *sig_curve_for2, double lda_for2,
                        double *sig_curve_fx2, double ***correlations,
                        double *var1, double *var2, double *covar) {
  double correl_fxfx[9], var_fx1[9], var_fx2[9];

  int StartIndex, EndIndex, i;
  double T1, T2;
  double sig_1_dom, sig_1_for, sig_1;
  double sig_2_dom, sig_2_for, sig_2;
  int is_fx1, is_fx2;
  Err err = NULL;

  memset(correl_fxfx, 0, 9 * sizeof(double));
  memset(var_fx1, 0, 9 * sizeof(double));
  memset(var_fx2, 0, 9 * sizeof(double));

  /* Check if FX or LGM */
  if (!sig_curve_for1 || !sig_curve_fx1) {
    is_fx1 = 0;
  } else {
    is_fx1 = 1;
  }

  if (!sig_curve_for2 || !sig_curve_fx2) {
    is_fx2 = 0;
  } else {
    is_fx2 = 1;
  }

  StartIndex = Get_Index(start_date, maturity, nbMat);
  EndIndex = Get_Index(end_date, maturity, nbMat);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = end_date;
    } else {
      T2 = maturity[i];
    }

    sig_1_dom = sig_curve_dom1[i];
    sig_2_dom = sig_curve_dom2[i];

    if (is_fx1) {
      sig_1_for = sig_curve_for1[i];
      sig_1 = sig_curve_fx1[i];
    }

    if (is_fx2) {
      sig_2_for = sig_curve_for2[i];
      sig_2 = sig_curve_fx2[i];
    }

    if (is_fx1 && is_fx2) {
      /* For COVAR */
      /* ********* */

      /* Z_1 with Z_2 */
      correl_fxfx[0] += sig_1 * sig_2 * (T2 - T1) * correlations[2][5][i];

      /* Z_1 with For_2 */
      correl_fxfx[2] += sig_1 * sig_2_for *
                        (Etha_Func(lda_for2, option_maturity2, T1, T2)) *
                        correlations[2][4][i];

      /* For_1 with Z_2 */
      correl_fxfx[6] += sig_1_for * sig_2 *
                        (Etha_Func(lda_for1, option_maturity1, T1, T2)) *
                        correlations[1][5][i];

      /* For_1 with For_2 */
      correl_fxfx[8] += sig_1_for * sig_2_for *
                        (Psi2_Func(lda_for1, lda_for2, option_maturity1,
                                   option_maturity2, T1, T2)) *
                        correlations[1][4][i];
    }

    if (is_fx1) {
      /* For VAR1 */
      /* ******** */

      /* Z_1 with Z_1 */
      var_fx1[0] += sig_1 * sig_1 * (T2 - T1);

      /* Z_1 with Dom_1 */
      var_fx1[1] += sig_1 * sig_1_dom *
                    (Etha_Func(lda_dom1, option_maturity1, T1, T2)) *
                    correlations[0][2][i];

      /* Z_1 with For_1 */
      var_fx1[2] += sig_1 * sig_1_for *
                    (Etha_Func(lda_for1, option_maturity1, T1, T2)) *
                    correlations[1][2][i];

      /* For_1 with Dom_1 */
      var_fx1[7] += sig_1_for * sig_1_dom *
                    (Psi_Func(lda_for1, lda_dom1, option_maturity1, T1, T2)) *
                    correlations[0][1][i];

      /* For_1 with For_1 */
      var_fx1[8] += sig_1_for * sig_1_for *
                    (Psi_Func(lda_for1, lda_for1, option_maturity1, T1, T2));

      /* For COVAR */
      /* ********* */

      if (is_fx2) {
        /* Z_1 with Dom_2 */
        correl_fxfx[1] += sig_1 * sig_2_dom *
                          (Etha_Func(lda_dom2, option_maturity2, T1, T2)) *
                          correlations[2][3][i];

        /* For_1 with Dom_2 */
        correl_fxfx[7] += sig_1_for * sig_2_dom *
                          (Psi2_Func(lda_for1, lda_dom2, option_maturity1,
                                     option_maturity2, T1, T2)) *
                          correlations[1][3][i];

      } else {
        /* Z_1 with Dom_2 */
        correl_fxfx[1] += sig_1 * sig_2_dom *
                          (Phi_Func(lda_dom2, option_maturity2, T1, T2)) *
                          correlations[2][3][i];

        /* For_1 with Dom_2 */
        // correl_fxfx[7] += sig_1_for * sig_2_dom * (Zeta_Func(lda_dom2 ,
        // lda_for1        , option_maturity        , T1        , T2))
        //	* correlations[1][3][i];
        correl_fxfx[7] += sig_1_for * sig_2_dom *
                          (Gamma2_Func(lda_for1, lda_dom2, option_maturity1,
                                       option_maturity2, T1, T2)) *
                          correlations[1][3][i];
      }
    }

    if (is_fx2) {
      /* For VAR2 */
      /* ******** */

      /* Z_2 with Z_2 */
      var_fx2[0] += sig_2 * sig_2 * (T2 - T1);

      /* Z_2 with For_2 */
      var_fx2[2] += sig_2 * sig_2_for *
                    (Etha_Func(lda_for2, option_maturity2, T1, T2)) *
                    correlations[4][5][i];

      /* Z_2 with Dom_2 */
      var_fx2[3] += sig_2 * sig_2_dom *
                    (Etha_Func(lda_dom2, option_maturity2, T1, T2)) *
                    correlations[3][5][i];

      /* Dom_2 with For_2 */
      var_fx2[5] += sig_2_dom * sig_2_for *
                    (Psi_Func(lda_dom2, lda_for2, option_maturity2, T1, T2)) *
                    correlations[3][4][i];

      /* For_2 with For_2 */
      var_fx2[8] += sig_2_for * sig_2_for *
                    (Psi_Func(lda_for2, lda_for2, option_maturity2, T1, T2));

      /* For COVAR */
      /* ********* */

      if (is_fx1) {
        /* Dom_1 with Z_2 */
        correl_fxfx[3] += sig_1_dom * sig_2 *
                          (Etha_Func(lda_dom1, option_maturity1, T1, T2)) *
                          correlations[0][5][i];

        /* Dom_1 with For_2 */
        correl_fxfx[5] += sig_1_dom * sig_2_for *
                          (Psi2_Func(lda_dom1, lda_for2, option_maturity1,
                                     option_maturity2, T1, T2)) *
                          correlations[0][4][i];
      } else {
        /* Dom_1 with Z_2 */
        correl_fxfx[3] += sig_1_dom * sig_2 *
                          (Phi_Func(lda_dom1, option_maturity1, T1, T2)) *
                          correlations[0][5][i];

        /* Dom_1 with For_2 */
        // correl_fxfx[5] += sig_1_dom * sig_2_for * (Zeta_Func(lda_dom1 ,
        // lda_for2        , option_maturity        , T1        , T2))
        //	* correlations[0][4][i];
        correl_fxfx[5] += sig_1_dom * sig_2_for *
                          (Gamma2_Func(lda_for2, lda_dom1, option_maturity2,
                                       option_maturity1, T1, T2)) *
                          correlations[0][4][i];
      }
    }

    /* For COVAR */
    /* ********* */

    /* Dom_1 with Dom_2 */
    if (is_fx1 && is_fx2) {
      correl_fxfx[4] += sig_1_dom * sig_2_dom *
                        (Psi2_Func(lda_dom1, lda_dom2, option_maturity1,
                                   option_maturity2, T1, T2)) *
                        correlations[0][3][i];
    } else if (is_fx1) {
      // correl_fxfx[4] += sig_1_dom * sig_2_dom * (Zeta_Func(lda_dom2        ,
      // lda_dom1        , option_maturity        , T1        , T2))
      //	* correlations[0][3][i];
      correl_fxfx[4] += sig_1_dom * sig_2_dom *
                        (Gamma2_Func(lda_dom1, lda_dom2, option_maturity1,
                                     option_maturity2, T1, T2)) *
                        correlations[0][3][i];
    } else if (is_fx2) {
      // correl_fxfx[4] += sig_1_dom * sig_2_dom * (Zeta_Func(lda_dom1        ,
      // lda_dom2        , option_maturity        , T1        , T2))
      //	* correlations[0][3][i];
      correl_fxfx[4] += sig_1_dom * sig_2_dom *
                        (Gamma2_Func(lda_dom2, lda_dom1, option_maturity2,
                                     option_maturity1, T1, T2)) *
                        correlations[0][3][i];
    } else {
      // correl_fxfx[4] += sig_1_dom * sig_2_dom * (Phi_Func(lda_dom1 + lda_dom2
      //       , option_maturity        , T1        , T2))
      //	* correlations[0][3][i];
      correl_fxfx[4] += sig_1_dom * sig_2_dom *
                        (Phi_Func(lda_dom1 + lda_dom2,
                                  (lda_dom1 * option_maturity1 +
                                   lda_dom2 * option_maturity2) /
                                      (lda_dom1 + lda_dom2),
                                  T1, T2)) *
                        correlations[0][3][i];
    }

    /* For VAR1 */
    /* ******** */

    /* Dom_1 with Dom_1 */
    if (is_fx1) {
      var_fx1[4] += sig_1_dom * sig_1_dom *
                    (Psi_Func(lda_dom1, lda_dom1, option_maturity1, T1, T2));
    } else {
      var_fx1[4] += sig_1_dom * sig_1_dom *
                    (Phi_Func(2.0 * lda_dom1, option_maturity1, T1, T2));
    }

    /* For VAR2 */
    /* ******** */

    /* Dom_2 with Dom_2 */
    if (is_fx2) {
      var_fx2[4] += sig_2_dom * sig_2_dom *
                    (Psi_Func(lda_dom2, lda_dom2, option_maturity2, T1, T2));
    } else {
      var_fx2[4] += sig_2_dom * sig_2_dom *
                    (Phi_Func(2.0 * lda_dom2, option_maturity2, T1, T2));
    }
  }

  *covar = correl_fxfx[0] + correl_fxfx[1] - correl_fxfx[2] + correl_fxfx[3] +
           correl_fxfx[4] - correl_fxfx[5] - correl_fxfx[6] - correl_fxfx[7] +
           correl_fxfx[8];

  *var1 = var_fx1[0] + 2.0 * var_fx1[1] - 2.0 * var_fx1[2] + var_fx1[4] -
          2.0 * var_fx1[7] + var_fx1[8];

  *var2 = var_fx2[0] - 2.0 * var_fx2[2] + 2.0 * var_fx2[3] + var_fx2[4] -
          2.0 * var_fx2[5] + var_fx2[8];

  return err;
}

Err display_FX1F_TermStruct(char *szFxUndName, long *sigma_n,
                            double **sigma_time, double **sigma)

{

  SrtLst *l;
  TermStruct *ts;
  Err err = NULL;
  SrtUndPtr sFxUnd;

  double today;

  /*-------------- Get the local volatility term
   * structure---------------------------*/

  sFxUnd = lookup_und(szFxUndName);
  if (!sFxUnd)
    return serror("Can not get the FX underlying name ");

  today = get_today_from_underlying(sFxUnd);
  err = get_underlying_ts(sFxUnd, &ts);
  if (err)
    serror("Can not get the underlying term struct ");

  l = ts->head;
  *sigma_n = 0;

  /* Allocate memory for the initial pointers */
  *sigma_time = (double *)malloc(sizeof(double));
  *sigma = (double *)malloc(sizeof(double));

  /* Loop on all the elements of the Term Sturcture */
  while (l != NULL) {

    (*sigma_n)++;
    *sigma_time = realloc(*sigma_time, (*sigma_n) * sizeof(double));
    *sigma = realloc(*sigma, (*sigma_n) * sizeof(double));
    (*sigma_time)[*sigma_n - 1] =
        ((double)(((FxTermStructVal *)l->element->val.pval)->date) - today) *
        YEARS_IN_DAY;
    (*sigma)[*sigma_n - 1] = ((FxTermStructVal *)l->element->val.pval)->sigx;

    l = l->next;
  }

  return err;
}

/*	Term correlation between a domestic IR index and the forward Fx */
Err quanto_correl(
    long fwd_start_date, /*	Fwd start date */
    long fix_date,       /*	Fixing date */
    long start_date,     /*	Index start date */
    long end_date,       /*	Index end date */
    long pay_date,       /*	Payment time */
    double fx_ivol_pay,  /*	Implied vol of Fx for payment time */
    char *dom_yc_name,   /*	Name of the domestic yield curve */
    char *dom_vc_name,   /*	Name of the domestic market vol curve */
    char *dom_ref_name,  /*	Name of the domestic reference rate */
    char *dom_swap_freq, /*	Frequency and basis of underlying swaptions */
    char *dom_swap_basis,
    char *for_yc_name,   /*	Name of the foreign yield curve */
    char *for_vc_name,   /*	Name of the foreign market vol curve */
    char *for_ref_name,  /*	Name of the foreign reference rate */
    char *for_swap_freq, /*	Frequency and basis of underlying swaptions */
    char *for_swap_basis,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Lambdas */
    double dom_lam, double for_lam,
    /*	The 3 standard correlation inputs */
    double corr_dom_for, double corr_spot_fx_dom, double corr_spot_fx_for,
    /*	The "exotic" correl term structure        , between domestic yield and
       domestic index */
    int num_corr, double *corr_mats, double *corr_dom_yld_idx,
    /*	The output */
    double *the_corr, double *quanto_adj) {
  long today, fx_fix_date;
  double fix_time, pay_time, fx_fix_time, fwd_start_time, start_time, end_time;
  int ndomsig, nforsig, nmergesig;
  double *domsigtime = NULL, *forsigtime = NULL, *domsig = NULL, *forsig = NULL,
         *mergesigtime = NULL, *mergedomsig = NULL, *mergeforsig = NULL,
         *mergecorr = NULL, *temp_fx_vol_curve = NULL;
  double sigma_s;
  int StartIndex, EndIndex, i;
  double t1, t2;
  double I1, I2, I3, I4, I5;
  double denomivol;

  SrtCurvePtr yc_ptr;

  Err err = NULL;

  /*	Get Today */
  yc_ptr = lookup_curve(dom_yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);
  yc_ptr = lookup_curve(for_yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  if (today != get_today_from_curve(yc_ptr)) {
    err = "Inconsistent reference dates";
    goto FREE_RETURN;
  }

  if (fwd_start_date == 0) {
    fwd_start_date = today;
  }

  /*	Calculate times from dates */
  if (pay_date < fix_date || start_date < fix_date || end_date <= start_date ||
      fwd_start_date > fix_date) {
    err = "Inconsistent dates";
    goto FREE_RETURN;
  }
  fix_time = (fix_date - today) * YEARS_IN_DAY;
  pay_time = (pay_date - today) * YEARS_IN_DAY;
  start_time = (start_date - today) * YEARS_IN_DAY;
  end_time = (end_date - today) * YEARS_IN_DAY;

  fwd_start_time = (fwd_start_date - today) * YEARS_IN_DAY;

  fx_fix_date = add_unit(pay_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
  fx_fix_time = (fx_fix_date - today) * YEARS_IN_DAY;

  /*	If fixing <= 1.5y just return correl between spot and index */
  if (fix_time < 1.5) {
    *the_corr = corr_spot_fx_dom;
    return NULL;
  }

  /*	Calibrate IR models */
  /*	Domestic */
  err = old_cpd_calib_diagonal_wrapper(
      dom_yc_name, dom_vc_name, dom_ref_name, get_cash_vol, 0.0, 0, 0, NULL,
      pay_date, NULL, NULL, 0, 0.0, 0.0, dom_swap_freq, dom_swap_basis, 1, 0, 0,
      0, 0, 0, 0, NULL, &dom_lam, 1, 0.0, 0.0, 0.0, &ndomsig, &domsigtime,
      &domsig, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Foreign */
  err = old_cpd_calib_diagonal_wrapper(
      for_yc_name, for_vc_name, for_ref_name, get_cash_vol, 0.0, 0, 0, NULL,
      pay_date, NULL, NULL, 0, 0.0, 0.0, for_swap_freq, for_swap_basis, 1, 0, 0,
      0, 0, 0, 0, NULL, &for_lam, 1, 0.0, 0.0, 0.0, &nforsig, &forsigtime,
      &forsig, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Merge rates ts */

  err = merge_rates_ts(domsigtime, domsig, ndomsig, forsigtime, forsig, nforsig,
                       &mergesigtime, &mergedomsig, &mergeforsig, &nmergesig);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Update correlations */
  mergecorr = (double *)calloc(nmergesig, sizeof(double));
  if (!mergecorr) {
    err = "Memory allocation failure in quanto_correl";
    goto FREE_RETURN;
  }
  for (i = 0; i < nmergesig; i++) {
    mergecorr[i] =
        interp(corr_mats, corr_dom_yld_idx, num_corr,
               (fix_time > mergesigtime[i] ? fix_time - mergesigtime[i] : 0.0),
               0, mergecorr + i);
  }

  /*	Calibrate Fx */
  err = Fx3DtsCalibration(&fx_fix_time, &pay_time, &fx_ivol_pay, 1,
                          mergesigtime, nmergesig, mergedomsig, dom_lam,
                          mergeforsig, for_lam, corr_dom_for, corr_spot_fx_dom,
                          corr_spot_fx_for, &temp_fx_vol_curve);
  if (err) {
    goto FREE_RETURN;
  }
  sigma_s = temp_fx_vol_curve[0];
  free(temp_fx_vol_curve);
  temp_fx_vol_curve = NULL;

  /*	Calculate the integrals	*/
  I1 = I2 = I3 = I4 = I5 = 0.0;
  StartIndex = Get_Index(fwd_start_time, mergesigtime, nmergesig);
  EndIndex = Get_Index(fix_time, mergesigtime, nmergesig);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      t1 = mergesigtime[i - 1];
    } else {
      /* First part */
      t1 = fwd_start_time;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      t2 = fix_time;
    } else {
      t2 = mergesigtime[i];
    }

    I1 += mergedomsig[i] * (exp(dom_lam * t2) - exp(dom_lam * t1));
    I2 += mergecorr[i] * mergedomsig[i] * mergedomsig[i] *
          (exp(dom_lam * t2) - exp(dom_lam * t1) -
           exp(-dom_lam * pay_time) *
               (exp(2.0 * dom_lam * t2) - exp(2.0 * dom_lam * t1)) / 2.0);
    I3 += mergedomsig[i] * mergeforsig[i] *
          ((exp(dom_lam * t2) - exp(dom_lam * t1)) / dom_lam -
           exp(-for_lam * pay_time) *
               (exp((dom_lam + for_lam) * t2) - exp((dom_lam + for_lam) * t1)) /
               (dom_lam + for_lam));
    I4 += mergedomsig[i] * mergedomsig[i] *
          (exp(2.0 * dom_lam * t2) - exp(2.0 * dom_lam * t1));
  }

  I1 *= corr_spot_fx_dom * sigma_s *
        (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
        dom_lam / (end_time - start_time);
  I2 *= (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
        dom_lam / dom_lam / (end_time - start_time);
  I3 *= corr_dom_for * (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) /
        dom_lam / for_lam / (end_time - start_time);
  I4 *= ((exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
         (end_time - start_time)) *
        ((exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
         (end_time - start_time)) /
        2.0 / dom_lam;

  /*	Implied vol of Fx */
  err = Fx3DtsImpliedVol(pay_time, fwd_start_time, fix_time, mergesigtime,
                         nmergesig, mergedomsig, dom_lam, mergeforsig, for_lam,
                         &fx_fix_time, &sigma_s, 1, corr_dom_for,
                         corr_spot_fx_dom, corr_spot_fx_for, &denomivol);
  if (err) {
    goto FREE_RETURN;
  }
  I5 = denomivol * denomivol * fix_time;

  /*	Return result */

  *the_corr = (I1 + I2 - I3) / sqrt(I4 * I5);
  *quanto_adj = (I1 + I2 - I3);

FREE_RETURN:

  if (domsigtime)
    free(domsigtime);
  if (domsig)
    free(domsig);
  if (forsigtime)
    free(forsigtime);
  if (forsig)
    free(forsig);
  if (mergesigtime)
    free(mergesigtime);
  if (mergedomsig)
    free(mergedomsig);
  if (mergeforsig)
    free(mergeforsig);
  if (mergecorr)
    free(mergecorr);
  if (temp_fx_vol_curve)
    free(temp_fx_vol_curve);

  return err;
}

/*	Term correlation between a domestic IR index and the forward Fx */
Err quanto_correl_corr(
    long fwd_start_date, /*	Fwd start date */
    long fix_date,       /*	Fixing time */
    long start_date,     /*	Index start date */
    long end_date,       /*	Index end date */
    long pay_date,       /*	Payment time */
    double fx_ivol_pay,  /*	Implied vol of Fx for payment time */
    char *dom_yc_name,   /*	Name of the domestic yield curve */
    char *dom_vc_name,   /*	Name of the domestic market vol curve */
    char *dom_ref_name,  /*	Name of the domestic reference rate */
    char *dom_swap_freq, /*	Frequency and basis of underlying swaptions */
    char *dom_swap_basis,
    char *for_yc_name,   /*	Name of the foreign yield curve */
    char *for_vc_name,   /*	Name of the foreign market vol curve */
    char *for_ref_name,  /*	Name of the foreign reference rate */
    char *for_swap_freq, /*	Frequency and basis of underlying swaptions */
    char *for_swap_basis,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Lambdas */
    double dom_lam, double for_lam,
    /*	The 3 standard correlation inputs */
    double *corr_times, double *corr_dom_for_ts, double *corr_dom_fx_ts,
    double *corr_for_fx_ts, long corr_n_times,
    /*	The "exotic" correl term structure        , between domestic yield and
       domestic index */
    int num_corr, double *corr_mats, double *corr_dom_yld_idx,
    /*	The output */
    double *the_corr, double *quanto_adj) {
  long today, fx_fix_date;
  double fix_time, pay_time, fx_fix_time, fwd_start_time, start_time, end_time;
  int ndomsig, nforsig, nmergesig;
  double *domsigtime = NULL, *forsigtime = NULL, *domsig = NULL, *forsig = NULL,
         *mergesigtime = NULL, *mergedomsig = NULL, *mergeforsig = NULL,
         *mergecorr = NULL, *mergecorrdomfor = NULL, *mergecorrdomfx = NULL,
         *mergecorrforfx = NULL, *temp_fx_vol_curve = NULL;
  double sigma_s;
  int StartIndex, EndIndex, i, index;
  double t1, t2;
  double I1, I2, I3, I4, I5;
  double denomivol;

  SrtCurvePtr yc_ptr;

  Err err = NULL;

  /*	Get Today */
  yc_ptr = lookup_curve(dom_yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);
  yc_ptr = lookup_curve(for_yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  if (today != get_today_from_curve(yc_ptr)) {
    err = "Inconsistent reference dates";
    goto FREE_RETURN;
  }

  if (fwd_start_date == 0) {
    fwd_start_date = today;
  }

  /*	Calculate times from dates */
  if (pay_date < fix_date || start_date < fix_date || end_date <= start_date ||
      fwd_start_date > fix_date) {
    err = "Inconsistent dates";
    goto FREE_RETURN;
  }
  fix_time = (fix_date - today) * YEARS_IN_DAY;
  pay_time = (pay_date - today) * YEARS_IN_DAY;
  start_time = (start_date - today) * YEARS_IN_DAY;
  end_time = (end_date - today) * YEARS_IN_DAY;

  fwd_start_time = (fwd_start_date - today) * YEARS_IN_DAY;

  fx_fix_date = add_unit(pay_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
  fx_fix_time = (fx_fix_date - today) * YEARS_IN_DAY;

  /*	If fixing <= 1.5y just return correl between spot and index */
  if (fix_time < 1.5) {
    I1 = 0.0;
    StartIndex = Get_Index(fwd_start_time, corr_times, corr_n_times);
    EndIndex = Get_Index(fix_time, corr_times, corr_n_times);

    for (i = StartIndex; i < EndIndex + 1; i++) {
      if (i > StartIndex) {
        t1 = corr_times[i - 1];
      } else {
        /* First part */
        t1 = fwd_start_time;
      }

      if (i == EndIndex || StartIndex == EndIndex) {
        /* Last part */
        t2 = fix_time;
      } else {
        t2 = corr_times[i];
      }

      I1 += corr_dom_fx_ts[i] * (t2 - t1);
    }

    *the_corr = I1 / (fix_time - fwd_start_time);
    return NULL;
  }

  /*	Calibrate IR models */
  /*	Domestic */
  err = old_cpd_calib_diagonal_wrapper(
      dom_yc_name, dom_vc_name, dom_ref_name, get_cash_vol, 0.0, 0, 0, NULL,
      pay_date, NULL, NULL, 0, 0.0, 0.0, dom_swap_freq, dom_swap_basis, 1, 0, 0,
      0, 0, 0, 0, NULL, &dom_lam, 1, 0.0, 0.0, 0.0, &ndomsig, &domsigtime,
      &domsig, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Foreign */
  err = old_cpd_calib_diagonal_wrapper(
      for_yc_name, for_vc_name, for_ref_name, get_cash_vol, 0.0, 0, 0, NULL,
      pay_date, NULL, NULL, 0, 0.0, 0.0, for_swap_freq, for_swap_basis, 1, 0, 0,
      0, 0, 0, 0, NULL, &for_lam, 1, 0.0, 0.0, 0.0, &nforsig, &forsigtime,
      &forsig, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Merge rates ts */

  err = merge_rates_ts(domsigtime, domsig, ndomsig, forsigtime, forsig, nforsig,
                       &mergesigtime, &mergedomsig, &mergeforsig, &nmergesig);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Update correlations */
  mergecorr = (double *)calloc(nmergesig, sizeof(double));
  mergecorrdomfor = (double *)calloc(nmergesig, sizeof(double));
  mergecorrdomfx = (double *)calloc(nmergesig, sizeof(double));
  mergecorrforfx = (double *)calloc(nmergesig, sizeof(double));

  if (!mergecorr || !mergecorrdomfor || !mergecorrdomfx || !mergecorrforfx) {
    err = "Memory allocation failure in quanto_correl";
    goto FREE_RETURN;
  }
  for (i = 0; i < nmergesig; i++) {
    mergecorr[i] =
        interp(corr_mats, corr_dom_yld_idx, num_corr,
               (fix_time > mergesigtime[i] ? fix_time - mergesigtime[i] : 0.0),
               0, mergecorr + i);

    index = Get_Index(mergesigtime[i], corr_times, corr_n_times);
    mergecorrdomfor[i] = corr_dom_for_ts[index];
    mergecorrdomfx[i] = corr_dom_fx_ts[index];
    mergecorrforfx[i] = corr_for_fx_ts[index];
  }

  /*	Calibrate Fx */
  err = Fx3DtsCalibration_corr(
      &fx_fix_time, &pay_time, &fx_ivol_pay, 1, mergesigtime, nmergesig,
      mergedomsig, dom_lam, mergeforsig, for_lam, corr_times, corr_dom_for_ts,
      corr_dom_fx_ts, corr_for_fx_ts, corr_n_times, &temp_fx_vol_curve);
  if (err) {
    goto FREE_RETURN;
  }
  sigma_s = temp_fx_vol_curve[0];
  free(temp_fx_vol_curve);
  temp_fx_vol_curve = NULL;

  /*	Calculate the integrals	*/
  I1 = I2 = I3 = I4 = I5 = 0.0;
  StartIndex = Get_Index(fwd_start_time, mergesigtime, nmergesig);
  EndIndex = Get_Index(fix_time, mergesigtime, nmergesig);

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      t1 = mergesigtime[i - 1];
    } else {
      /* First part */
      t1 = fwd_start_time;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      t2 = fix_time;
    } else {
      t2 = mergesigtime[i];
    }

    I1 += mergecorrdomfx[i] * mergedomsig[i] *
          (exp(dom_lam * t2) - exp(dom_lam * t1));
    I2 += mergecorr[i] * mergedomsig[i] * mergedomsig[i] *
          (exp(dom_lam * t2) - exp(dom_lam * t1) -
           exp(-dom_lam * pay_time) *
               (exp(2.0 * dom_lam * t2) - exp(2.0 * dom_lam * t1)) / 2.0);
    I3 += mergecorrdomfor[i] * mergedomsig[i] * mergeforsig[i] *
          ((exp(dom_lam * t2) - exp(dom_lam * t1)) / dom_lam -
           exp(-for_lam * pay_time) *
               (exp((dom_lam + for_lam) * t2) - exp((dom_lam + for_lam) * t1)) /
               (dom_lam + for_lam));
    I4 += mergedomsig[i] * mergedomsig[i] *
          (exp(2.0 * dom_lam * t2) - exp(2.0 * dom_lam * t1));
  }

  I1 *= sigma_s * (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) /
        dom_lam / dom_lam / (end_time - start_time);
  I2 *= (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
        dom_lam / dom_lam / (end_time - start_time);
  I3 *= (exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
        for_lam / (end_time - start_time);
  I4 *= ((exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
         (end_time - start_time)) *
        ((exp(-dom_lam * start_time) - exp(-dom_lam * end_time)) / dom_lam /
         (end_time - start_time)) /
        2.0 / dom_lam;

  /*	Implied vol of Fx */
  err = Fx3DtsImpliedVol_corr(pay_time, fwd_start_time, fix_time, mergesigtime,
                              nmergesig, mergedomsig, dom_lam, mergeforsig,
                              for_lam, &fx_fix_time, &sigma_s, 1, corr_times,
                              corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts,
                              corr_n_times, &denomivol);
  if (err) {
    goto FREE_RETURN;
  }
  I5 = denomivol * denomivol * fix_time;

  /*	Return result */

  *the_corr = (I1 + I2 - I3) / sqrt(I4 * I5);
  *quanto_adj = (I1 + I2 - I3);

FREE_RETURN:

  if (domsigtime)
    free(domsigtime);
  if (domsig)
    free(domsig);
  if (forsigtime)
    free(forsigtime);
  if (forsig)
    free(forsig);
  if (mergesigtime)
    free(mergesigtime);
  if (mergedomsig)
    free(mergedomsig);
  if (mergeforsig)
    free(mergeforsig);
  if (mergecorr)
    free(mergecorr);
  if (mergecorrdomfor)
    free(mergecorrdomfor);
  if (mergecorrdomfx)
    free(mergecorrdomfx);
  if (mergecorrforfx)
    free(mergecorrforfx);
  if (temp_fx_vol_curve)
    free(temp_fx_vol_curve);

  return err;
}

/*	**********************************************************	*/
/*	**********************************************************	*/
/*	All this part is a second version of the functions above	*/
/*	with a term-structure of coorelation
 */
/*	**********************************************************	*/
/*	**********************************************************	*/

Err Get_FX_StochRate_TermStructures_corr(
    char *underlying, double **sigma_date_dom, double **sigma_dom,
    long *sigma_n_dom, double **tau_date_dom, double **tau_dom, long *tau_n_dom,
    double **sigma_date_for, double **sigma_for, long *sigma_n_for,
    double **tau_date_for, double **tau_for, long *tau_n_for,
    double **sigma_date_fx, double **sigma_fx, long *sigma_n_fx,
    double **correl_date, double **correl_dom_for, double **correl_dom_fx,
    double **correl_for_fx, long *correl_n) {
  SrtUndPtr und;
  char *domname, *forname;
  double *corr = NULL;
  long today;
  int i;

  Err err = NULL;

  *sigma_date_dom = NULL;
  *sigma_dom = NULL;
  *tau_date_dom = NULL;
  *tau_dom = NULL;
  *sigma_date_for = NULL;
  *sigma_for = NULL;
  *tau_date_for = NULL;
  *tau_for = NULL;
  *sigma_date_fx = NULL;
  *sigma_fx = NULL;
  *correl_date = NULL;
  *correl_dom_for = NULL;
  *correl_dom_fx = NULL;
  *correl_for_fx = NULL;

  und = lookup_und(underlying);
  if (!und) {
    err = serror("Couldn't find underlying named %s", underlying);
    goto FREE_RETURN;
  }

  if (get_underlying_type(und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", underlying);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", underlying);
    goto FREE_RETURN;
  }

  domname = get_domname_from_fxund(und);
  err = Get_LGM_TermStructure(domname, sigma_date_dom, sigma_dom, sigma_n_dom,
                              tau_date_dom, tau_dom, tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  forname = get_forname_from_fxund(und);
  err = Get_LGM_TermStructure(forname, sigma_date_for, sigma_for, sigma_n_for,
                              tau_date_for, tau_for, tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = srt_f_display_FX_TermStruct(underlying, sigma_n_fx, sigma_date_fx,
                                    sigma_fx, correl_n, correl_date, &corr);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(und);
  for (i = 0; i < *sigma_n_fx; i++) {
    (*sigma_date_fx)[i] = ((*sigma_date_fx)[i] - today) / 365.0;
  }

  *correl_dom_for = calloc(*correl_n, sizeof(double));
  *correl_dom_fx = calloc(*correl_n, sizeof(double));
  *correl_for_fx = calloc(*correl_n, sizeof(double));

  if (!*correl_dom_for || !*correl_dom_fx || !*correl_for_fx) {
    err = "Memory allocation failure in Get_FX_StochRate_TermStructures_corr";
    goto FREE_RETURN;
  }

  for (i = 0; i < *correl_n; i++) {
    (*correl_date)[i] = ((*correl_date)[i] - today) / 365.0;
    (*correl_dom_for)[i] = corr[i * 3];
    (*correl_dom_fx)[i] = corr[i * 3 + 1];
    (*correl_for_fx)[i] = corr[i * 3 + 2];
  }

FREE_RETURN:

  if (err) {
    if (*sigma_date_dom)
      free(*sigma_date_dom);
    *sigma_date_dom = NULL;
    if (*sigma_dom)
      free(*sigma_dom);
    *sigma_dom = NULL;

    if (*sigma_date_for)
      free(*sigma_date_for);
    *sigma_date_for = NULL;
    if (*sigma_for)
      free(*sigma_for);
    *sigma_for = NULL;

    if (*tau_date_dom)
      free(*tau_date_dom);
    *tau_date_dom = NULL;
    if (*tau_dom)
      free(*tau_dom);
    *tau_dom = NULL;

    if (*tau_date_for)
      free(*tau_date_for);
    *tau_date_for = NULL;
    if (*tau_for)
      free(*tau_for);
    *tau_for = NULL;

    if (*sigma_date_fx)
      free(*sigma_date_fx);
    *sigma_date_fx = NULL;
    if (*sigma_fx)
      free(*sigma_fx);
    *sigma_fx = NULL;

    if (*correl_date)
      free(*correl_date);
    *correl_date = NULL;
    if (*correl_dom_for)
      free(*correl_dom_for);
    *correl_dom_for = NULL;
    if (*correl_dom_fx)
      free(*correl_dom_fx);
    *correl_dom_fx = NULL;
    if (*correl_for_fx)
      free(*correl_for_fx);
    *correl_for_fx = NULL;
  }

  if (corr)
    free(corr);

  return err;
}

/*	Fx Implied volatility */
/*  The rates maturity dates have to be merged ! */
Err Fx3DtsImpliedVol_corr(double opt_maturity, double start_date,
                          double end_date, double *maturity_rates, long nbMat,
                          double *sig_curve_dom, double lda_dom,
                          double *sig_curve_for, double lda_for,
                          double *maturity_fx, double *sig_curve_fx,
                          long nbrMat_fx, double *maturity_corr,
                          double *correl_dom_for, double *correl_dom_fx,
                          double *correl_for_fx, long nbrMat_corr,
                          double *fx_vol) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2, ta, tb;
  double var;
  double var_partial;
  int i, j, k;
  long StartIndex, EndIndex, StartIndex2, EndIndex2, StartIndex3, EndIndex3;
  Err err = NULL;

  if (start_date > end_date) {
    err = "end_date before start_date in Fx3DtsImpliedVol";
    return err;
  }

  if (end_date == 0) {
    (*fx_vol) = 0;
    return err;
  }

  StartIndex = Get_Index(start_date, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(end_date, maturity_fx, nbrMat_fx);

  var = 0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = start_date;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = end_date;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }
      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      StartIndex3 = Get_Index(t1, maturity_corr, nbrMat_corr);
      EndIndex3 = Get_Index(t2, maturity_corr, nbrMat_corr);

      for (k = StartIndex3; k < EndIndex3 + 1; k++) {
        if (k > StartIndex3) {
          ta = maturity_corr[k - 1];
        } else {
          /* First part */
          ta = t1;
        }

        if (k == EndIndex3 || StartIndex3 == EndIndex3) {
          /* Last part */
          tb = t2;
        } else {
          tb = maturity_corr[k];
        }

        err = Partial_Var(opt_maturity, ta, tb, sig_dom, lda_dom, sig_for,
                          lda_for, sig_fx, correl_dom_for[k], correl_dom_fx[k],
                          correl_for_fx[k], &var_partial);

        if (err) {
          return err;
        }

        var += var_partial;
      }
    }
  }

  if (fabs(end_date - start_date) > 1.0e-08) {
    *fx_vol = sqrt(var / (end_date - start_date));
  } else {
    *fx_vol = 0.0;
  }

  return err;
}

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DtsCalibration_corr(double *exercise_opt, double *maturity_opt,
                           double *vol_opt, long nbrOpt, double *maturity_rates,
                           long nbrMat, double *sig_curve_dom, double lda_dom,
                           double *sig_curve_for, double lda_for,
                           double *maturity_corr, double *correl_dom_for,
                           double *correl_dom_fx, double *correl_for_fx,
                           long nbrCorr, double **fx_vol_curve) {
  double a, b, c, a_part, b_part, c_part, c2, delta;
  double sig_dom, sig_for;
  double cumvar, vol;
  double T1, T2, t1, t2, ta, tb;
  double VolMin;
  double val_opt;
  int idxopt, i, j;
  long StartIndex, EndIndex, StartIndex2, EndIndex2;

  Err err = NULL;

  /* loop on the number of options */
  (*fx_vol_curve) = NULL;
  (*fx_vol_curve) = (double *)calloc(nbrOpt, sizeof(double));
  if (!(*fx_vol_curve)) {
    err = "Memory allocation error in Fx3DtsCalibration";
    goto FREE_RETURN;
  }

  for (idxopt = 0; idxopt <= nbrOpt - 1; idxopt++) {
    T2 = exercise_opt[idxopt];
    val_opt = maturity_opt[idxopt];

    /* Get the cumulated variance till the previous maturity T1*/
    if (idxopt > 0) {
      T1 = exercise_opt[idxopt - 1];

      err = Fx3DtsImpliedVol_corr(val_opt, 0, T1, maturity_rates, nbrMat,
                                  sig_curve_dom, lda_dom, sig_curve_for,
                                  lda_for, exercise_opt, *fx_vol_curve, nbrOpt,
                                  maturity_corr, correl_dom_for, correl_dom_fx,
                                  correl_for_fx, nbrCorr, &vol);

      if (err) {
        goto FREE_RETURN;
      }

      cumvar = vol * vol * T1;
    } else {
      T1 = 0;
      cumvar = 0;
    }

    /* Calculate the coefficient of the last cumulative of the variance between
     * T1 and T2 */

    a = b = c2 = 0;
    StartIndex = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex = Get_Index(T2, maturity_rates, nbrMat);

    for (i = StartIndex; i < EndIndex + 1; i++) {
      if (i > StartIndex) {
        t1 = maturity_rates[i - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (i == EndIndex || StartIndex == EndIndex) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[i];
      }

      sig_dom = sig_curve_dom[i];
      sig_for = sig_curve_for[i];

      StartIndex2 = Get_Index(t1, maturity_corr, nbrCorr);
      EndIndex2 = Get_Index(t2, maturity_corr, nbrCorr);

      for (j = StartIndex2; j < EndIndex2 + 1; j++) {
        if (j > StartIndex2) {
          ta = maturity_corr[j - 1];
        } else {
          /* First part */
          ta = t1;
        }

        if (j == EndIndex2 || StartIndex2 == EndIndex2) {
          /* Last part */
          tb = t2;
        } else {
          tb = maturity_corr[j];
        }

        err = Coefs_Partial_Var(val_opt, ta, tb, sig_dom, lda_dom, sig_for,
                                lda_for, correl_dom_for[j], correl_dom_fx[j],
                                correl_for_fx[j], &a_part, &b_part, &c_part);

        if (err) {
          goto FREE_RETURN;
        }

        a += a_part;
        b += b_part;
        c2 += c_part;
      }
    }

    c = c2 + cumvar;
    /* substract value to match */
    c -= vol_opt[idxopt] * vol_opt[idxopt] * T2;

    /* just solve the second order equation */
    delta = b * b - 4 * a * c;
    /* delta < 0 or solutions are negatives */
    if ((delta < 0) || ((c > 0) && (b > 0))) {
      VolMin = sqrt((c2 + cumvar - b * b / (4 * a)) / val_opt) * 100;

      err = serror("Cannot find solution to match the option %d (its vol "
                   "should be > %.2f %%)",
                   idxopt + 1, VolMin);

      goto FREE_RETURN;
    }

    /* if there is two positive solutions we are taking the smallest one */
    /* it would be easier to match the next one */
    if ((c > 0) && (b < 0)) {
      (*fx_vol_curve)[idxopt] = (-b - sqrt(delta)) / (2. * a);
    } else {
      (*fx_vol_curve)[idxopt] = (-b + sqrt(delta)) / (2. * a);
    }
  }

FREE_RETURN:

  if (err) {
    if (*fx_vol_curve) {
      free(*fx_vol_curve);
      *fx_vol_curve = NULL;
    }
  }

  return err;
}

/*	Implied vol direct from underlying */
Err Fx3DImpliedVol_corr(char *fx_underlying, double val_time, double start_time,
                        double end_time, double *vol) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx,
      nb_merge_dates, nb_corr;
  double *sigma_date_dom = NULL, *sigma_dom = NULL;
  double *tau_date_dom = NULL, *tau_dom = NULL;
  double *sigma_date_for = NULL, *sigma_for = NULL;
  double *tau_date_for = NULL, *tau_for = NULL;
  double *sigma_date_fx = NULL, *sigma_fx = NULL;
  double *correl_date = NULL, *correl_dom_for = NULL, *correl_dom_fx = NULL,
         *correl_for_fx = NULL;
  double lda_dom, lda_for;
  double *sig_dom = NULL, *sig_for = NULL, *merge_dates = NULL;

  Err err = NULL;

  err = Get_FX_StochRate_TermStructures_corr(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_date, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx, &nb_corr);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DtsImpliedVol_corr(
      val_time, start_time, end_time, merge_dates, nb_merge_dates, sig_dom,
      lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx,
      correl_date, correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr, vol);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (correl_date) {
    free(correl_date);
  }

  if (correl_dom_for) {
    free(correl_dom_for);
  }

  if (correl_dom_fx) {
    free(correl_dom_fx);
  }

  if (correl_for_fx) {
    free(correl_for_fx);
  }

  return err;
}

/*	Fx calibration direct from underlying */
Err Fx3DCalibration_corr(char *dom_underlying, char *for_underlying,
                         double *maturity_correl, double *correl_dom_for,
                         double *correl_dom_fx, double *correl_for_fx,
                         long nb_correl, double *exercise_opt,
                         double *maturity_opt, double *vol_opt, long nbropt,
                         double **fx_vol_curve) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;
  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  err = Fx3DtsCalibration_corr(
      exercise_opt, maturity_opt, vol_opt, nbropt, merge_dates, nb_merge_dates,
      sig_dom, lda_dom, sig_for, lda_for, maturity_correl, correl_dom_for,
      correl_dom_fx, correl_for_fx, nb_correl, fx_vol_curve);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

/*	Calculate adjustment between Q-Tpay and Q-beta        , such that
        Q-beta expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay expect log ( S (
   Tfix ) / S ( 0 ) ) + adj */
Err Fx3DtsFwdBetaAdjustment_corr(
    double T0,   /*	Forward start date */
    double Tval, /*	Value date of the forward */
    double Tpay, /*	Pay date of the forward */
    double Tfix, /*	Fix date of the forward */
    /*	Model data */
    double *maturity_rates, long nbrMat, double *sig_curve_dom, double lda_dom,
    double *sig_curve_for, double lda_for, double *maturity_fx,
    double *sig_curve_fx, long nbrMat_fx, double *maturity_corr,
    double *correl_dom_for, double *correl_dom_fx, double *correl_for_fx,
    long nbrMat_corr,
    /*	Result */
    double *adjust) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2, ta, tb;
  double adj;
  double adjust_partial;
  int i, j, k;
  long StartIndex, EndIndex, StartIndex2, EndIndex2, StartIndex3, EndIndex3;
  Err err = NULL;

  if (T0 > Tfix) {
    err = "end_date before start_date in Fx3DtsFwdBetaAdjustment";
    return err;
  }

  if (Tfix == 0) {
    (*adjust) = 0;
    return err;
  }

  StartIndex = Get_Index(T0, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(Tfix, maturity_fx, nbrMat_fx);

  adj = 0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = T0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = Tfix;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbrMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }

      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      StartIndex3 = Get_Index(t1, maturity_corr, nbrMat_corr);
      EndIndex3 = Get_Index(t2, maturity_corr, nbrMat_corr);

      for (k = StartIndex3; k < EndIndex3 + 1; k++) {
        if (k > StartIndex3) {
          ta = maturity_corr[k - 1];
        } else {
          /* First part */
          ta = t1;
        }

        if (k == EndIndex3 || StartIndex3 == EndIndex3) {
          /* Last part */
          tb = t2;
        } else {
          tb = maturity_corr[k];
        }

        adjust_partial =
            correl_dom_fx[k] * sig_dom * Etha_Func(lda_dom, Tpay, ta, tb) *
                sig_fx -
            correl_dom_for[k] * sig_dom * sig_for *
                Psi2_Func(lda_for, lda_dom, Tval, Tpay, ta, tb) +
            sig_dom * sig_dom * Psi2_Func(lda_dom, lda_dom, Tval, Tpay, ta, tb);

        adj += adjust_partial;
      }
    }
  }

  *adjust = adj;

  return err;
}

/*	Q beta */
Err Fx3DFwdBeta_corr(char *fx_underlying,
                     /*	Forward fixing time */
                     double Tfix,
                     /*	Forward value time */
                     double Tval,
                     /*	0: expect [log Fx]        , 1: expect [Fx] */
                     int log_flag,
                     /*	Answer */
                     double *expect) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx, nb_corr;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *correl_date = NULL,
         *correl_dom_for = NULL, *correl_dom_fx = NULL, *correl_for_fx = NULL,
         *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL,
         *sigma_date_fx = NULL, *sigma_fx = NULL;

  double lda_dom, lda_for;

  double adjustment, fx_vol, fx_fwd;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  long today, fix_date, val_date;
  double spot_fx;

  Err err = NULL;

  /*	Read model data */

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures_corr(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_date, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx, &nb_corr);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  today = get_today_from_underlying(dom_und);
  fix_date = (long)(today + Tfix * DAYS_IN_YEAR + 1.0e-08);
  val_date = (long)(today + Tval * DAYS_IN_YEAR + 1.0e-08);

  /*	Calculate fwd        , spot and implied vol */

  /*	Forward */
  err = Fx3DFwdFx(fx_underlying, val_date, &fx_fwd);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Spot */
  err = Fx3DFwdFx(fx_underlying, today, &spot_fx);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Implied vol */
  err = Fx3DImpliedVol_corr(fx_underlying, Tval, 0.0, Tfix, &fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculate expectation of the log under Q-Tfix */

  *expect = log(fx_fwd / spot_fx) - 0.5 * fx_vol * fx_vol * Tfix;

  /*	Adjust for Q-beta */

  err = Fx3DtsFwdBetaAdjustment_corr(
      0.0, Tval, Tval, Tfix, merge_dates, nb_merge_dates, sig_dom, lda_dom,
      sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx, correl_date,
      correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr, &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  if (log_flag) {
    *expect = spot_fx * exp(*expect + 0.5 * fx_vol * fx_vol * Tfix);
  }

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (correl_date) {
    free(correl_date);
  }

  if (correl_dom_for) {
    free(correl_dom_for);
  }

  if (correl_dom_fx) {
    free(correl_dom_fx);
  }

  if (correl_for_fx) {
    free(correl_for_fx);
  }

  return err;
}

/*	Calculate adjustment between Q-Tpay1 and Q-Tpay2        , such that
        Q-Tpay2 expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay1 expect log ( S (
   Tfix ) / S ( 0 ) ) + adj */
Err Fx3DtsFwdPayAdjustment_corr(double T0,    /*	Forward start date */
                                double Tval,  /*	Value date of the forward */
                                double Tpay1, /*	Original pay date */
                                double Tpay2, /*	Adjusted pay date */
                                double Tfix,  /*	Fix date of the forward */
                                /*	Model data */
                                double *maturity_rates, long nbrMat,
                                double *sig_curve_dom, double lda_dom,
                                double *sig_curve_for, double lda_for,
                                double *maturity_fx, double *sig_curve_fx,
                                long nbrMat_fx, double *maturity_corr,
                                double *correl_dom_for, double *correl_dom_fx,
                                double *correl_for_fx, long nbrMat_corr,
                                /*	Result */
                                double *adjust) {
  double adj1, adj2;
  Err err = NULL;

  err = Fx3DtsFwdBetaAdjustment_corr(
      T0, Tval, Tpay1, Tfix, maturity_rates, nbrMat, sig_curve_dom, lda_dom,
      sig_curve_for, lda_for, maturity_fx, sig_curve_fx, nbrMat_fx,
      maturity_corr, correl_dom_for, correl_dom_fx, correl_for_fx, nbrMat_corr,
      &adj1);

  if (err) {
    return err;
  }

  err = Fx3DtsFwdBetaAdjustment_corr(
      T0, Tval, Tpay2, Tfix, maturity_rates, nbrMat, sig_curve_dom, lda_dom,
      sig_curve_for, lda_for, maturity_fx, sig_curve_fx, nbrMat_fx,
      maturity_corr, correl_dom_for, correl_dom_fx, correl_for_fx, nbrMat_corr,
      &adj2);

  *adjust = adj1 - adj2;

  return err;
}

Err OptRainbowts3F(double T, double df, double K, double nX, double X0,
                   double nY, double Y0, double *maturity, long nbMat,
                   double *sig_curve_dom1, double lda_dom1,
                   double *sig_curve_for1, double lda_for1,
                   double *sig_curve_fx1, double *sig_curve_dom2,
                   double lda_dom2, double *sig_curve_for2, double lda_for2,
                   double *sig_curve_fx2, double ***correlations, double *res,
                   double *correl) {
  double vol1, vol2, var1, var2, fwd_cor, covar;
  Err err = NULL;

  err = Fx3DtsImpliedVol_corr(
      T, 0, T, maturity, nbMat, sig_curve_dom1, lda_dom1, sig_curve_for1,
      lda_for1, maturity, sig_curve_fx1, nbMat, maturity, correlations[0][1],
      correlations[0][2], correlations[1][2], nbMat, &vol1);

  if (err)
    return err;

  err = Fx3DtsImpliedVol_corr(
      T, 0, T, maturity, nbMat, sig_curve_dom2, lda_dom2, sig_curve_for2,
      lda_for2, maturity, sig_curve_fx2, nbMat, maturity, correlations[3][4],
      correlations[3][5], correlations[4][5], nbMat, &vol2);

  if (err)
    return err;

  err = Fx3DtsFxFwdCov_corr(T, T, 0, T, maturity, nbMat, sig_curve_dom1,
                            lda_dom1, sig_curve_for1, lda_for1, sig_curve_fx1,
                            sig_curve_dom2, lda_dom2, sig_curve_for2, lda_for2,
                            sig_curve_fx2, correlations, &var1, &var2, &covar);

  if (err)
    return err;

  fwd_cor = covar / (vol1 * vol2 * T);

  *res = srt_f_optrainbow(T, df, K, nX, X0, vol1, nY, Y0, vol2, fwd_cor);

  *correl = fwd_cor;

  return err;
}

/*	Implied correlations */
Err Fx3DtsCorrel(double start_mat, double end_mat, double val_mat,
                 double *sig_dates, long nb_sig_dates, double *sig_curve_dom,
                 double lda_dom, double *sig_curve_for, double lda_for,
                 double *sig_curve_fx, double *correl_dom_for_ts,
                 double *correl_dom_fx_ts, double *correl_for_fx_ts,
                 double *dom_vol, double *for_vol, double *fx_vol,
                 double *dom_for_cor, double *dom_fx_cor, double *for_fx_cor) {
  double sig_dom, sig_dom2;
  double sig_for, sig_for2;
  double sig_fx, sig_domfor;
  double T1, T2;
  double var_dom, var_for;
  double x_dom, y_dom, x_domfor;
  double x_for, y_for;
  double mat, sqmat;
  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double correl12, correl13_1, correl13_2, correl13_3, correl23_1, correl23_2,
      correl23_3;
  int i;
  long StartIndex, EndIndex;
  Err err = NULL;

  StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);
  EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
  mat = (end_mat - start_mat);
  sqmat = sqrt(mat);

  var_dom = 0;
  var_for = 0;

  /*	Implied Fx Vol */
  err = Fx3DtsImpliedVol_corr(
      val_mat, start_mat, end_mat, sig_dates, nb_sig_dates, sig_curve_dom,
      lda_dom, sig_curve_for, lda_for, sig_dates, sig_curve_fx, nb_sig_dates,
      sig_dates, correl_dom_for_ts, correl_dom_fx_ts, correl_for_fx_ts,
      nb_sig_dates, fx_vol);

  correl12 = 0;
  correl13_1 = correl13_2 = correl13_3 = 0;
  correl23_1 = correl23_2 = correl23_3 = 0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = sig_dates[i - 1];
    } else {
      /* First part */
      T1 = start_mat;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = end_mat;
    } else {
      T2 = sig_dates[i];
    }

    sig_dom = sig_dom2 = sig_domfor = sig_curve_dom[i];
    sig_dom2 *= sig_dom2;
    sig_for = sig_for2 = sig_curve_for[i];
    sig_for2 *= sig_for2;
    sig_domfor *= sig_for;
    sig_fx = sig_curve_fx[i];
    correl_dom_for = correl_dom_for_ts[i];
    correl_dom_fx = correl_dom_fx_ts[i];
    correl_for_fx = correl_for_fx_ts[i];

    x_dom = Z_Func(lda_dom, T1, T2);
    y_dom = Z_Func(2 * lda_dom, T1, T2);

    x_for = Z_Func(lda_for, T1, T2);
    y_for = Z_Func(2 * lda_for, T1, T2);

    x_domfor = Z_Func(lda_dom + lda_for, T1, T2);

    var_dom += sig_dom2 * y_dom;
    var_for += sig_for2 * y_for;

    /* Correlations */

    /* domestic / foreign */
    correl12 += correl_dom_for * sig_dom * sig_for * x_domfor;

    /* domestic / fx */
    correl13_1 += sig_dom2 * (x_dom - exp(-lda_dom * val_mat) * y_dom);
    correl13_2 += correl_dom_for * sig_domfor *
                  (x_dom - exp(-lda_for * val_mat) * x_domfor);
    correl13_3 += correl_dom_fx * sig_dom * sig_fx * x_dom;

    /* foreign fx */
    correl23_1 += sig_for2 * (x_for - exp(-lda_for * val_mat) * y_for);
    correl23_2 += correl_dom_for * sig_domfor *
                  (x_for - exp(-lda_dom * val_mat) * x_domfor);
    correl23_3 += correl_for_fx * sig_for * sig_fx * x_for;
  }

  *dom_vol = exp(-lda_dom * val_mat) * sqrt(var_dom / mat);
  *for_vol = exp(-lda_for * val_mat) * sqrt(var_for / mat);

  *dom_for_cor = correl12 * exp(-(lda_dom + lda_for) * val_mat) /
                 ((*dom_vol) * (*for_vol) * mat);

  *dom_fx_cor = exp(-lda_dom * val_mat) *
                (correl13_1 / lda_dom - correl13_2 / lda_for + correl13_3) /
                ((*dom_vol) * (*fx_vol) * mat);

  *for_fx_cor = exp(-lda_for * val_mat) *
                (-correl23_1 / lda_for + correl23_2 / lda_dom + correl23_3) /
                ((*for_vol) * (*fx_vol) * mat);

  return err;
}

/*	Implied correlations */
Err Fx3DCorrel(char *fx_underlying, double val_time, double start_time,
               double end_time, double *dom_vol, double *for_vol,
               double *fx_vol, double *dom_for_cor, double *dom_fx_cor,
               double *for_fx_cor) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx,
      nb_merge_dates, nb_corr;
  double *sigma_date_dom = NULL, *sigma_dom = NULL;
  double *tau_date_dom = NULL, *tau_dom = NULL;
  double *sigma_date_for = NULL, *sigma_for = NULL;
  double *tau_date_for = NULL, *tau_for = NULL;
  double *sigma_date_fx = NULL, *sigma_fx = NULL;
  double *correl_date = NULL, *correl_dom_for = NULL, *correl_dom_fx = NULL,
         *correl_for_fx = NULL;
  double lda_dom, lda_for;
  double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL, *corr_dom_for = NULL,
         *corr_dom_fx = NULL, *corr_for_fx = NULL, *merge_dates = NULL;

  Err err = NULL;

  err = Get_FX_StochRate_TermStructures_corr(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_date, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx, &nb_corr);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_fx_corr_ts(
      sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for, sigma_for,
      sigma_n_for, sigma_date_fx, sigma_fx, sigma_n_fx, correl_date,
      correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr, &merge_dates,
      &sig_dom, &sig_for, &sig_fx, &corr_dom_for, &corr_dom_fx, &corr_for_fx,
      &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DtsCorrel(start_time, end_time, val_time, merge_dates,
                     nb_merge_dates, sig_dom, lda_dom, sig_for, lda_for, sig_fx,
                     corr_dom_for, corr_dom_fx, corr_for_fx, dom_vol, for_vol,
                     fx_vol, dom_for_cor, dom_fx_cor, for_fx_cor);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (sig_fx) {
    free(sig_fx);
  }

  if (corr_dom_for) {
    free(corr_dom_for);
  }

  if (corr_dom_fx) {
    free(corr_dom_fx);
  }

  if (corr_for_fx) {
    free(corr_for_fx);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (correl_date) {
    free(correl_date);
  }

  if (correl_dom_for) {
    free(correl_dom_for);
  }

  if (correl_dom_fx) {
    free(correl_dom_fx);
  }

  if (correl_for_fx) {
    free(correl_for_fx);
  }

  return err;
}

Err srt_f_optrainbow3F(double T, double K, double nX, double nY, char *fx_und1,
                       char *fx_und2, double *res, double *correl) {

  double *sig_curve_dom1 = NULL, *sig_curve_for1 = NULL, *sig_curve_fx1 = NULL,
         *sig_curve_dom2 = NULL, *sig_curve_for2 = NULL, *sig_curve_fx2 = NULL;

  double ***correlations = NULL;

  char *dom_name1, *dom_name2, *for_name1, *for_name2;

  double *sigma_date_dom1 = NULL, *sigma_dom1 = NULL, *tau_date_dom1 = NULL,
         *tau_dom1 = NULL, *sigma_date_dom2 = NULL, *sigma_dom2 = NULL,
         *tau_date_dom2 = NULL, *tau_dom2 = NULL, *sigma_date_for1 = NULL,
         *sigma_for1 = NULL, *tau_date_for1 = NULL, *tau_for1 = NULL,
         *sigma_fx1 = NULL, *sigma_date_for2 = NULL, *sigma_for2 = NULL,
         *tau_date_for2 = NULL, *tau_for2 = NULL, *sigma_fx2 = NULL,
         *sigma_date_fx1 = NULL, *sigma_date_fx2 = NULL, *corr1 = NULL,
         *corr_date1 = NULL, *corr2 = NULL, *corr_date2 = NULL,
         *all_dates = NULL;

  double lda_dom1, lda_dom2, lda_for1, lda_for2;
  long sigma_n_dom1, tau_n_dom1, sigma_n_dom2, tau_n_dom2, sigma_n_for1,
      tau_n_for1, sigma_n_for2, tau_n_for2, sigma_n_fx1, sigma_n_fx2,
      nb_corr_date1, nb_corr_date2, nb_dates;

  long i;
  double today;
  long maturity_date;
  char *yc;
  double df;

  double X0, Y0;

  SrtUndPtr und_fx1, und_fx2, und_dom1, und_dom2;
  SrtCorrLstPtr corr_list = NULL;

  Err err = NULL;

  /* get the correlation list */
  corr_list = srt_f_GetTheCorrelationList();

  if (!corr_list) {
    err = "Cannot get the correlation structure";
    goto FREE_RETURN;
  }

  und_fx1 = lookup_und(fx_und1);
  und_fx2 = lookup_und(fx_und2);

  today = get_today_from_underlying(und_fx1);
  maturity_date = (long)(today + 365 * T + 1E-12);

  dom_name1 = get_domname_from_fxund(und_fx1);
  for_name1 = get_forname_from_fxund(und_fx1);
  dom_name2 = get_domname_from_fxund(und_fx2);
  for_name2 = get_forname_from_fxund(und_fx2);

  if (!(!strcmp(dom_name1, dom_name2))) {
    smessage("Be carreful        , Fx underlyings must have the same domestic "
             "underlying");
  }

  und_dom1 = lookup_und(dom_name1);
  und_dom2 = lookup_und(dom_name2);

  yc = get_ycname_from_irund(und_dom1);

  df = swp_f_df(today, maturity_date, yc);

  err = Fx3DFwdFx(fx_und1, maturity_date, &X0);
  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DFwdFx(fx_und2, maturity_date, &Y0);
  if (err) {
    goto FREE_RETURN;
  }

  /* Get all the term structures */
  err = Get_LGM_TermStructure(dom_name1, &sigma_date_dom1, &sigma_dom1,
                              &sigma_n_dom1, &tau_date_dom1, &tau_dom1,
                              &tau_n_dom1);
  if (err) {
    goto FREE_RETURN;
  }

  err = Get_LGM_TermStructure(for_name1, &sigma_date_for1, &sigma_for1,
                              &sigma_n_for1, &tau_date_for1, &tau_for1,
                              &tau_n_for1);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom1, tau_n_dom1, &lda_dom1);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for1, tau_n_for1, &lda_for1);
  if (err) {
    goto FREE_RETURN;
  }

  err = srt_f_display_FX_TermStruct(fx_und1, &sigma_n_fx1, &sigma_date_fx1,
                                    &sigma_fx1, &nb_corr_date1, &corr_date1,
                                    &corr1);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_corr_date1; i++)
    corr_date1[i] = (corr_date1[i] - today) / 365.0;

  for (i = 0; i < sigma_n_fx1; i++) {
    sigma_date_fx1[i] = (sigma_date_fx1[i] - today) / 365.0;
  }

  err = Get_LGM_TermStructure(dom_name2, &sigma_date_dom2, &sigma_dom2,
                              &sigma_n_dom2, &tau_date_dom2, &tau_dom2,
                              &tau_n_dom2);
  if (err) {
    goto FREE_RETURN;
  }

  err = Get_LGM_TermStructure(for_name2, &sigma_date_for2, &sigma_for2,
                              &sigma_n_for2, &tau_date_for2, &tau_for2,
                              &tau_n_for2);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom2, tau_n_dom2, &lda_dom2);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for2, tau_n_for2, &lda_for2);
  if (err) {
    goto FREE_RETURN;
  }

  err = srt_f_display_FX_TermStruct(fx_und2, &sigma_n_fx2, &sigma_date_fx2,
                                    &sigma_fx2, &nb_corr_date2, &corr_date2,
                                    &corr2);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_corr_date2; i++)
    corr_date2[i] = (corr_date2[i] - today) / 365.0;

  for (i = 0; i < sigma_n_fx2; i++) {
    sigma_date_fx2[i] = (sigma_date_fx2[i] - today) / 365.0;
  }

  /* Merge everything */

  all_dates = calloc(sigma_n_dom1, sizeof(double));
  if (!all_dates) {
    err = "Memory allocation failure in srt_f_optrainbow3F (1)";
    goto FREE_RETURN;
  }
  memcpy(all_dates, sigma_date_dom1, sigma_n_dom1 * sizeof(double));
  nb_dates = sigma_n_dom1;

  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_for1, sigma_date_for1);
  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_fx1, sigma_date_fx1);
  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_dom2, sigma_date_dom2);
  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_for2, sigma_date_for2);
  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_fx2, sigma_date_fx2);
  num_f_concat_vector(&nb_dates, &all_dates, nb_corr_date1, corr_date1);
  num_f_concat_vector(&nb_dates, &all_dates, nb_corr_date2, corr_date2);
  num_f_sort_vector(nb_dates, all_dates);
  num_f_unique_vector(&nb_dates, all_dates);

  sig_curve_dom1 = calloc(nb_dates, sizeof(double));
  sig_curve_for1 = calloc(nb_dates, sizeof(double));
  sig_curve_fx1 = calloc(nb_dates, sizeof(double));
  sig_curve_dom2 = calloc(nb_dates, sizeof(double));
  sig_curve_for2 = calloc(nb_dates, sizeof(double));
  sig_curve_fx2 = calloc(nb_dates, sizeof(double));

  correlations = f3tensor(0, 9, 0, 9, 0, nb_dates - 1);

  if (!sig_curve_dom1 || !sig_curve_for1 || !sig_curve_fx1 || !sig_curve_dom2 ||
      !sig_curve_for2 || !sig_curve_fx2 || !correlations) {
    err = "Memory allocation failure in srt_f_optrainbow3F (2)";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_dates; i++) {
    sig_curve_dom1[i] =
        sigma_dom1[Get_Index(all_dates[i], sigma_date_dom1, sigma_n_dom1)];
    sig_curve_for1[i] =
        sigma_for1[Get_Index(all_dates[i], sigma_date_for1, sigma_n_for1)];
    sig_curve_fx1[i] =
        sigma_fx1[Get_Index(all_dates[i], sigma_date_fx1, sigma_n_fx1)];
    sig_curve_dom2[i] =
        sigma_dom2[Get_Index(all_dates[i], sigma_date_dom2, sigma_n_dom2)];
    sig_curve_for2[i] =
        sigma_for2[Get_Index(all_dates[i], sigma_date_for2, sigma_n_for2)];
    sig_curve_fx2[i] =
        sigma_fx2[Get_Index(all_dates[i], sigma_date_fx2, sigma_n_fx2)];

    /* get the correlations */

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, for_name1,
                                       all_dates[i], &(correlations[0][1][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, fx_und1,
                                       all_dates[i], &(correlations[0][2][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, dom_name2,
                                       all_dates[i], &(correlations[0][3][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, for_name2,
                                       all_dates[i], &(correlations[0][4][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, fx_und2,
                                       all_dates[i], &(correlations[0][5][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, for_name1, fx_und1,
                                       all_dates[i], &(correlations[1][2][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, for_name1, dom_name2,
                                       all_dates[i], &(correlations[1][3][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, for_name1, for_name2,
                                       all_dates[i], &(correlations[1][4][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, for_name1, fx_und2,
                                       all_dates[i], &(correlations[1][5][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, fx_und1, dom_name2,
                                       all_dates[i], &(correlations[2][3][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, fx_und1, for_name2,
                                       all_dates[i], &(correlations[2][4][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name2, for_name2,
                                       all_dates[i], &(correlations[3][4][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name2, fx_und2,
                                       all_dates[i], &(correlations[3][5][i]));
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_get_corr_from_CorrList(corr_list, fx_und1, fx_und2,
                                       all_dates[i], &(correlations[2][5][i]));
    if (err) {
      goto FREE_RETURN;
    }

    correlations[3][5][i] = correlations[0][5][i];

    err = srt_f_get_corr_from_CorrList(corr_list, for_name2, fx_und2,
                                       all_dates[i], &(correlations[4][5][i]));
    if (err) {
      goto FREE_RETURN;
    }
  }

  err = OptRainbowts3F(T, df, K, nX, X0, nY, Y0, all_dates, nb_dates,
                       sig_curve_dom1, lda_dom1, sig_curve_for1, lda_for1,
                       sig_curve_fx1, sig_curve_dom2, lda_dom2, sig_curve_for2,
                       lda_for2, sig_curve_fx2, correlations, res, correl);

  if (!(!strcmp(dom_name1, dom_name2))) {
    *res = 0.0;
  }

FREE_RETURN:

  if (sigma_date_dom1) {
    free(sigma_date_dom1);
  }

  if (sigma_dom1) {
    free(sigma_dom1);
  }

  if (tau_date_dom1) {
    free(tau_date_dom1);
  }

  if (tau_dom1) {
    free(tau_dom1);
  }

  if (sigma_date_for1) {
    free(sigma_date_for1);
  }

  if (sigma_for1) {
    free(sigma_for1);
  }

  if (tau_date_for1) {
    free(tau_date_for1);
  }

  if (tau_for1) {
    free(tau_for1);
  }

  if (sigma_date_fx1) {
    free(sigma_date_fx1);
  }

  if (sigma_fx1) {
    free(sigma_fx1);
  }

  if (corr_date1) {
    free(corr_date1);
  }

  if (corr1) {
    free(corr1);
  }

  if (sigma_date_dom2) {
    free(sigma_date_dom2);
  }

  if (sigma_dom2) {
    free(sigma_dom2);
  }

  if (tau_date_dom2) {
    free(tau_date_dom2);
  }

  if (tau_dom2) {
    free(tau_dom2);
  }

  if (sigma_date_for2) {
    free(sigma_date_for2);
  }

  if (sigma_for2) {
    free(sigma_for2);
  }

  if (tau_date_for2) {
    free(tau_date_for2);
  }

  if (tau_for2) {
    free(tau_for2);
  }

  if (sigma_date_fx2) {
    free(sigma_date_fx2);
  }

  if (sigma_fx2) {
    free(sigma_fx2);
  }

  if (corr_date2) {
    free(corr_date2);
  }

  if (corr2) {
    free(corr2);
  }

  if (all_dates) {
    free(all_dates);
  }

  if (sig_curve_dom1) {
    free(sig_curve_dom1);
  }

  if (sig_curve_for1) {
    free(sig_curve_for1);
  }

  if (sig_curve_fx1) {
    free(sig_curve_fx1);
  }

  if (sig_curve_dom2) {
    free(sig_curve_dom2);
  }

  if (sig_curve_for2) {
    free(sig_curve_for2);
  }

  if (sig_curve_fx2) {
    free(sig_curve_fx2);
  }

  if (correlations) {
    free_f3tensor(correlations, 0, 9, 0, 9, 0, nb_dates - 1);
  }

  return err;
}

/*	Q Tpay */
Err Fx3DFwdTpay_corr(char *fx_underlying,
                     /*	Forward fixing time */
                     double Tfix,
                     /*	Forward value time */
                     double Tval,
                     /*	Payment time */
                     double Tpay,
                     /*	0: expect [log Fx]        , 1: expect [Fx] */
                     int log_flag,
                     /*	Answer */
                     double *expect) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx, nb_corr;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL, *correl_date = NULL, *correl_dom_for = NULL,
         *correl_dom_fx = NULL, *correl_for_fx = NULL;

  double lda_dom, lda_for;

  double adjustment, fx_vol, fx_fwd;
  SrtUndPtr dom_und, fx_und;
  char *domname;
  long today, fix_date, val_date, pay_date;
  double spot_fx;

  Err err = NULL;

  /*	Read model data */

  fx_und = lookup_und(fx_underlying);
  if (!fx_und) {
    return serror("Couldn't fin underlying named %s", fx_underlying);
  }

  if (get_underlying_type(fx_und) != FOREX_UND) {
    return serror("Underlying %s is not of type IR", fx_underlying);
  }

  err = Get_FX_StochRate_TermStructures_corr(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_date, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx, &nb_corr);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);

  today = get_today_from_underlying(dom_und);
  fix_date = (long)(today + Tfix * DAYS_IN_YEAR + 1.0e-08);
  val_date = (long)(today + Tval * DAYS_IN_YEAR + 1.0e-08);
  pay_date = (long)(today + Tpay * DAYS_IN_YEAR + 1.0e-08);

  /*	Calculate fwd        , spot and implied vol */

  /*	Forward */
  err = Fx3DFwdFx(fx_underlying, val_date, &fx_fwd);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Spot */
  err = Fx3DFwdFx(fx_underlying, today, &spot_fx);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Implied vol */
  err = Fx3DImpliedVol_corr(fx_underlying, Tval, 0.0, Tfix, &fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculate expectation of the log under Q-Tfix */

  *expect = log(fx_fwd / spot_fx) - 0.5 * fx_vol * fx_vol * Tfix;

  /*	Adjust for Q-Tpay */

  err = Fx3DtsFwdPayAdjustment_corr(
      0.0, Tval, Tval, Tpay, Tfix, merge_dates, nb_merge_dates, sig_dom,
      lda_dom, sig_for, lda_for, sigma_date_fx, sigma_fx, sigma_n_fx,
      correl_date, correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr,
      &adjustment);

  if (err) {
    goto FREE_RETURN;
  }

  *expect += adjustment;

  if (log_flag) {
    *expect = spot_fx * exp(*expect + 0.5 * fx_vol * fx_vol * Tfix);
  }

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (correl_date) {
    free(correl_date);
  }

  if (correl_dom_for) {
    free(correl_dom_for);
  }

  if (correl_dom_fx) {
    free(correl_dom_fx);
  }

  if (correl_for_fx) {
    free(correl_for_fx);
  }

  return err;
}

Err Fx3DGenCorrelation(char *und1_name, char *und2_name, double val_time,
                       double start_time, double end_time, double *vol1,
                       double *vol2, double *correl) {

  double *sig_curve_dom1 = NULL, *sig_curve_for1 = NULL, *sig_curve_fx1 = NULL,
         *sig_curve_dom2 = NULL, *sig_curve_for2 = NULL, *sig_curve_fx2 = NULL;

  double ***correlations = NULL;

  char *dom_name1, *dom_name2, *for_name1, *for_name2;

  double *sigma_date_dom1 = NULL, *sigma_dom1 = NULL, *tau_date_dom1 = NULL,
         *tau_dom1 = NULL, *sigma_date_dom2 = NULL, *sigma_dom2 = NULL,
         *tau_date_dom2 = NULL, *tau_dom2 = NULL, *sigma_date_for1 = NULL,
         *sigma_for1 = NULL, *tau_date_for1 = NULL, *tau_for1 = NULL,
         *sigma_fx1 = NULL, *sigma_date_for2 = NULL, *sigma_for2 = NULL,
         *tau_date_for2 = NULL, *tau_for2 = NULL, *sigma_fx2 = NULL,
         *sigma_date_fx1 = NULL, *sigma_date_fx2 = NULL, *corr1 = NULL,
         *corr_date1 = NULL, *corr2 = NULL, *corr_date2 = NULL,
         *all_dates = NULL;

  double lda_dom1, lda_dom2, lda_for1, lda_for2;
  long sigma_n_dom1, sigma_n_dom2, sigma_n_for1, sigma_n_for2, sigma_n_fx1,
      sigma_n_fx2, nb_corr_date1, nb_corr_date2, nb_dates;

  int is_fx1, is_fx2;

  long i;
  double today;

  SrtUndPtr und1, und2;
  SrtCorrLstPtr corr_list = NULL;

  Err err = NULL;

  /* get the correlation list */
  corr_list = srt_f_GetTheCorrelationList();

  if (!corr_list) {
    err = "Cannot get the correlation structure";
    goto FREE_RETURN;
  }

  /* check type of und1 */
  und1 = lookup_und(und1_name);

  if (!und1) {
    err = serror("Couldn't find underlying named %s", und1_name);
    goto FREE_RETURN;
  }

  switch (get_underlying_type(und1)) {
  case INTEREST_RATE_UND: {
    if (get_mdltype_from_irund(und1) != LGM &&
        get_mdldim_from_irund(und1) != ONE_FAC) {
      err = serror("Underlying named %s must be of type LGM1F", und1_name);
      goto FREE_RETURN;
    }

    is_fx1 = 0;
    dom_name1 = und1_name;

    break;
  }

  case FOREX_UND: {
    if (get_mdltype_from_fxund(und1) != FX_STOCH_RATES) {
      err = serror("Underlying %s is not of type FX Stoch Rates", und1_name);
      goto FREE_RETURN;
    }

    is_fx1 = 1;
    dom_name1 = get_domname_from_fxund(und1);
    for_name1 = get_forname_from_fxund(und1);

    break;
  }

  default: {
    err = ("Underlying named %s must be of type LGM1F of FX_STOCH_RATES",
           und1_name);
    goto FREE_RETURN;
    break;
  }
  }

  /* check type of und2 */
  und2 = lookup_und(und2_name);

  if (!und2) {
    err = serror("Couldn't find underlying named %s", und2_name);
    goto FREE_RETURN;
  }

  switch (get_underlying_type(und2)) {
  case INTEREST_RATE_UND: {
    if (get_mdltype_from_irund(und2) != LGM &&
        get_mdldim_from_irund(und2) != ONE_FAC) {
      err = serror("Underlying named %s must be of type LGM1F", und2_name);
      goto FREE_RETURN;
    }

    is_fx2 = 0;
    dom_name2 = und2_name;

    break;
  }

  case FOREX_UND: {
    if (get_mdltype_from_fxund(und2) != FX_STOCH_RATES) {
      err = serror("Underlying %s is not of type FX Stoch Rates", und2_name);
      goto FREE_RETURN;
    }

    is_fx2 = 1;
    dom_name2 = get_domname_from_fxund(und2);
    for_name2 = get_forname_from_fxund(und2);

    break;
  }

  default: {
    err = ("Underlying named %s must be of type LGM1F of FX_STOCH_RATES",
           und2_name);
    goto FREE_RETURN;
    break;
  }
  }

  /* Get und1 Term Struct */

  today = get_today_from_underlying(und1);

  /* Get all the term structures */
  err = Get_LGM_TermStructure2(dom_name1, &sigma_date_dom1, &sigma_dom1,
                               &sigma_n_dom1, &lda_dom1);

  if (err)
    goto FREE_RETURN;

  lda_dom1 = 1.0 / lda_dom1;

  if (is_fx1) {
    err = Get_LGM_TermStructure2(for_name1, &sigma_date_for1, &sigma_for1,
                                 &sigma_n_for1, &lda_for1);

    if (err)
      goto FREE_RETURN;

    lda_for1 = 1.0 / lda_for1;

    err = srt_f_display_FX_TermStruct(und1_name, &sigma_n_fx1, &sigma_date_fx1,
                                      &sigma_fx1, &nb_corr_date1, &corr_date1,
                                      &corr1);

    if (err)
      goto FREE_RETURN;

    for (i = 0; i < nb_corr_date1; i++)
      corr_date1[i] = (corr_date1[i] - today) * YEARS_IN_DAY;
    for (i = 0; i < sigma_n_fx1; i++)
      sigma_date_fx1[i] = (sigma_date_fx1[i] - today) * YEARS_IN_DAY;
  }

  err = Get_LGM_TermStructure2(dom_name2, &sigma_date_dom2, &sigma_dom2,
                               &sigma_n_dom2, &lda_dom2);

  if (err)
    goto FREE_RETURN;

  lda_dom2 = 1.0 / lda_dom2;

  if (is_fx2) {
    err = Get_LGM_TermStructure2(for_name2, &sigma_date_for2, &sigma_for2,
                                 &sigma_n_for2, &lda_for2);

    if (err)
      goto FREE_RETURN;

    lda_for2 = 1.0 / lda_for2;

    err = srt_f_display_FX_TermStruct(und2_name, &sigma_n_fx2, &sigma_date_fx2,
                                      &sigma_fx2, &nb_corr_date2, &corr_date2,
                                      &corr2);

    if (err)
      goto FREE_RETURN;

    for (i = 0; i < nb_corr_date2; i++)
      corr_date2[i] = (corr_date2[i] - today) * YEARS_IN_DAY;
    for (i = 0; i < sigma_n_fx2; i++)
      sigma_date_fx2[i] = (sigma_date_fx2[i] - today) * YEARS_IN_DAY;
  }

  /* Merge everything */

  all_dates = calloc(sigma_n_dom1, sizeof(double));
  if (!all_dates) {
    err = "Memory allocation failure in srt_f_optrainbow3F (1)";
    goto FREE_RETURN;
  }

  memcpy(all_dates, sigma_date_dom1, sigma_n_dom1 * sizeof(double));
  nb_dates = sigma_n_dom1;

  num_f_concat_vector(&nb_dates, &all_dates, sigma_n_dom2, sigma_date_dom2);

  if (is_fx1) {
    num_f_concat_vector(&nb_dates, &all_dates, sigma_n_for1, sigma_date_for1);
    num_f_concat_vector(&nb_dates, &all_dates, sigma_n_fx1, sigma_date_fx1);
    num_f_concat_vector(&nb_dates, &all_dates, nb_corr_date1, corr_date1);
  }

  if (is_fx2) {
    num_f_concat_vector(&nb_dates, &all_dates, sigma_n_for2, sigma_date_for2);
    num_f_concat_vector(&nb_dates, &all_dates, sigma_n_fx2, sigma_date_fx2);
    num_f_concat_vector(&nb_dates, &all_dates, nb_corr_date2, corr_date2);
  }

  num_f_sort_vector(nb_dates, all_dates);
  num_f_unique_vector(&nb_dates, all_dates);

  sig_curve_dom1 = calloc(nb_dates, sizeof(double));
  sig_curve_dom2 = calloc(nb_dates, sizeof(double));

  if (is_fx1) {
    sig_curve_for1 = calloc(nb_dates, sizeof(double));
    sig_curve_fx1 = calloc(nb_dates, sizeof(double));
  }

  if (is_fx2) {
    sig_curve_for2 = calloc(nb_dates, sizeof(double));
    sig_curve_fx2 = calloc(nb_dates, sizeof(double));
  }

  correlations = f3tensor(0, 9, 0, 9, 0, nb_dates - 1);

  if (!sig_curve_dom1 || !sig_curve_dom2 ||
      (is_fx1 && (!sig_curve_for1 || !sig_curve_fx1)) ||
      (is_fx2 && (!sig_curve_for2 || !sig_curve_fx2)) || !correlations) {
    err = "Memory allocation failure in srt_f_optrainbow3F (2)";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_dates; i++) {
    sig_curve_dom1[i] =
        sigma_dom1[Get_Index(all_dates[i], sigma_date_dom1, sigma_n_dom1)];
    sig_curve_dom2[i] =
        sigma_dom2[Get_Index(all_dates[i], sigma_date_dom2, sigma_n_dom2)];

    err = srt_f_get_corr_from_CorrList(corr_list, dom_name1, dom_name2,
                                       all_dates[i], &(correlations[0][3][i]));
    if (err)
      goto FREE_RETURN;

    if (is_fx1 && is_fx2) {
      err =
          srt_f_get_corr_from_CorrList(corr_list, for_name1, for_name2,
                                       all_dates[i], &(correlations[1][4][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, for_name1, und2_name,
                                       all_dates[i], &(correlations[1][5][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, und1_name, for_name2,
                                       all_dates[i], &(correlations[2][4][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, und1_name, und2_name,
                                       all_dates[i], &(correlations[2][5][i]));
      if (err)
        goto FREE_RETURN;
    }

    if (is_fx1) {
      sig_curve_for1[i] =
          sigma_for1[Get_Index(all_dates[i], sigma_date_for1, sigma_n_for1)];
      sig_curve_fx1[i] =
          sigma_fx1[Get_Index(all_dates[i], sigma_date_fx1, sigma_n_fx1)];

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name1, for_name1,
                                       all_dates[i], &(correlations[0][1][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name1, und1_name,
                                       all_dates[i], &(correlations[0][2][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, for_name1, und1_name,
                                       all_dates[i], &(correlations[1][2][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, for_name1, dom_name2,
                                       all_dates[i], &(correlations[1][3][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, und1_name, dom_name2,
                                       all_dates[i], &(correlations[2][3][i]));
      if (err)
        goto FREE_RETURN;
    }

    if (is_fx2) {
      sig_curve_for2[i] =
          sigma_for2[Get_Index(all_dates[i], sigma_date_for2, sigma_n_for2)];
      sig_curve_fx2[i] =
          sigma_fx2[Get_Index(all_dates[i], sigma_date_fx2, sigma_n_fx2)];

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name2, for_name2,
                                       all_dates[i], &(correlations[3][4][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name2, und2_name,
                                       all_dates[i], &(correlations[3][5][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, for_name2, und2_name,
                                       all_dates[i], &(correlations[4][5][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name1, for_name2,
                                       all_dates[i], &(correlations[0][4][i]));
      if (err)
        goto FREE_RETURN;

      err =
          srt_f_get_corr_from_CorrList(corr_list, dom_name1, und2_name,
                                       all_dates[i], &(correlations[0][5][i]));
      if (err)
        goto FREE_RETURN;
    }
  }

  err = Fx3DtsFxFwdCov_corr(val_time, val_time, start_time, end_time, all_dates,
                            nb_dates, sig_curve_dom1, lda_dom1, sig_curve_for1,
                            lda_for1, sig_curve_fx1, sig_curve_dom2, lda_dom2,
                            sig_curve_for2, lda_for2, sig_curve_fx2,
                            correlations, vol1, vol2, correl);

  if (err)
    goto FREE_RETURN;

  *vol1 = sqrt((*vol1) / (end_time - start_time));
  *vol2 = sqrt((*vol2) / (end_time - start_time));
  *correl = *correl / ((*vol1) * (*vol2) * (end_time - start_time));

FREE_RETURN:

  if (sigma_date_dom1) {
    free(sigma_date_dom1);
  }

  if (sigma_dom1) {
    free(sigma_dom1);
  }

  if (tau_date_dom1) {
    free(tau_date_dom1);
  }

  if (tau_dom1) {
    free(tau_dom1);
  }

  if (sigma_date_for1) {
    free(sigma_date_for1);
  }

  if (sigma_for1) {
    free(sigma_for1);
  }

  if (tau_date_for1) {
    free(tau_date_for1);
  }

  if (tau_for1) {
    free(tau_for1);
  }

  if (sigma_date_fx1) {
    free(sigma_date_fx1);
  }

  if (sigma_fx1) {
    free(sigma_fx1);
  }

  if (corr_date1) {
    free(corr_date1);
  }

  if (corr1) {
    free(corr1);
  }

  if (sigma_date_dom2) {
    free(sigma_date_dom2);
  }

  if (sigma_dom2) {
    free(sigma_dom2);
  }

  if (tau_date_dom2) {
    free(tau_date_dom2);
  }

  if (tau_dom2) {
    free(tau_dom2);
  }

  if (sigma_date_for2) {
    free(sigma_date_for2);
  }

  if (sigma_for2) {
    free(sigma_for2);
  }

  if (tau_date_for2) {
    free(tau_date_for2);
  }

  if (tau_for2) {
    free(tau_for2);
  }

  if (sigma_date_fx2) {
    free(sigma_date_fx2);
  }

  if (sigma_fx2) {
    free(sigma_fx2);
  }

  if (corr_date2) {
    free(corr_date2);
  }

  if (corr2) {
    free(corr2);
  }

  if (all_dates) {
    free(all_dates);
  }

  if (sig_curve_dom1) {
    free(sig_curve_dom1);
  }

  if (sig_curve_for1) {
    free(sig_curve_for1);
  }

  if (sig_curve_fx1) {
    free(sig_curve_fx1);
  }

  if (sig_curve_dom2) {
    free(sig_curve_dom2);
  }

  if (sig_curve_for2) {
    free(sig_curve_for2);
  }

  if (sig_curve_fx2) {
    free(sig_curve_fx2);
  }

  if (correlations) {
    free_f3tensor(correlations, 0, 9, 0, 9, 0, nb_dates - 1);
  }

  return err;
}