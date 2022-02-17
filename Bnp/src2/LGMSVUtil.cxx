#include "LGMSVUtil.h"
#include "Fx3FUtils.h"
#include "LGMSVPDE.h"

#define LGMSV_MINVOL 0.0001
#define LGMSV_MAXVOL 0.1

/* -----------------------------------------------------------------------------------------------
        LGMSVIntMult
                Example p=2
                  Calculate
                        Integral between t1 and t of
                                Integral between t1 and s of (u-t1)^n *
   exp(-lambda[p]*(s-u))*du times exp(-lambda[p-1]*(t-s))*ds
   ------------------------------------------------------------------------------------------------
 */
static double LGMSVIntMult(
    /* Inputs */
    int p,          /* Number of integration */
    double *lambda, /* Array of lambda */

    int n,     /* Power in the last integral */
    double dDt /* = t-t1 */
) {
  double x1, x2;
  int iLevel;

  if (p == 0) {
    /* Zero level of integration */
    return pow(dDt, n);
  } else {
    /* p>0 level of integration */
    if (fabs(lambda[p]) < 1.0e-06) {
      return 1.0 / (n + 1) * LGMSVIntMult(p - 1, lambda, n + 1, dDt);
    } else {
      if (n == 0) {
        x1 = LGMSVIntMult(p - 1, lambda, 0, dDt);
        for (iLevel = 0; iLevel < p; iLevel++) {
          lambda[iLevel] = lambda[iLevel] - lambda[p];
        }
        x2 = LGMSVIntMult(p - 1, lambda, 0, dDt);
        for (iLevel = 0; iLevel < p; iLevel++) {
          lambda[iLevel] = lambda[iLevel] + lambda[p];
        }
        return 1.0 / lambda[p] * (x1 - exp(-lambda[p] * dDt) * x2);
      } else {
        x1 = LGMSVIntMult(p - 1, lambda, n, dDt);
        x2 = LGMSVIntMult(p, lambda, n - 1, dDt);
        return 1.0 / lambda[p] * (x1 - n * x2);
      }
    }
  }
}

/* -----------------------------------------------------------------------------------------------
        LGMSVFuncValue
                Evaluation of
                        Xt  = Xt1*exp(-lambda*(t-t1))+integral between t1 and t
   of fs*exp(-lambda*(s-t1))

                where fs is of the same type as Xt
   ------------------------------------------------------------------------------------------------
 */
double LGMSVFuncValue(
    /* Inputs */
    LGMSVSolFunc Xt, /*  Function to be valued */

    double dDt,     /* = t-t1 */
    double *Lambda, /* Array of lambda */
    int iLevelIntegration /* 0 for valuation */) {
  double dXt1, dLambda, accu;
  LGMSVSolFunc *pft;
  LGMSVSolFunc *pgt;
  LGMSVSolFunc *pht;
  int iNumLevel;

  /* Extract value from structure of Xt */
  dXt1 = Xt.dXt1;       /* Value at t1 */
  dLambda = Xt.dLambda; /* Lambda */

  /* Fill the Lambda array for the iLevelIntegration */
  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++)
    Lambda[iNumLevel] = Lambda[iNumLevel] - dLambda;

  accu = dXt1 * exp(-dLambda * dDt) *
         LGMSVIntMult(iLevelIntegration, Lambda, 0, dDt);

  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++)
    Lambda[iNumLevel] = Lambda[iNumLevel] + dLambda;

  /* Fill the Lambda array for the iLevelIntegration+1 */
  Lambda[iLevelIntegration + 1] = dLambda;

  /* First sub function */
  if (Xt.a != 0) {
    if (Xt.bIsft1 == 1)
      accu += Xt.a * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pft = (LGMSVSolFunc *)(Xt.pft);

      accu += Xt.a * LGMSVFuncValue((*pft), dDt, Lambda, iLevelIntegration + 1);
    }
  }

  /* Second sub function */
  if (Xt.b != 0) {
    if (Xt.bIsgt1 == 1)
      accu += Xt.b * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pgt = (LGMSVSolFunc *)(Xt.pgt);

      accu += Xt.b * LGMSVFuncValue((*pgt), dDt, Lambda, iLevelIntegration + 1);
    }
  }

  /* Second sub function */
  if (Xt.cxx != 0) {
    if (Xt.bIsht1 == 1)
      accu += Xt.cxx * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pht = (LGMSVSolFunc *)(Xt.pht);

      accu +=
          Xt.cxx * LGMSVFuncValue((*pht), dDt, Lambda, iLevelIntegration + 1);
    }
  }

  return accu;
}

/* -----------------------------------------------------------------------------------------------
        LGMSVFuncValue
                Evaluation of
                        Xt  = Xt1*exp(-lambda*(t-t1))+integral between t1 and t
   of fs*exp(-lambda*(s-t1))

                where fs is of the same type as Xt
   ------------------------------------------------------------------------------------------------
 */
void LGMSVFuncValue2(
    /* Inputs */
    LGMSVSolFunc *Xt, /*  Function to be valued */

    double dDt,     /* = t-t1 */
    double *Lambda, /* Array of lambda */
    int iLevelIntegration /* 0 for valuation */, double *res) {
  double dXt1, dLambda, accu;
  LGMSVSolFunc *pft;
  LGMSVSolFunc *pgt;
  LGMSVSolFunc *pht;
  int iNumLevel;
  double temp;

  /* Extract value from structure of Xt */
  dXt1 = Xt->dXt1;       /* Value at t1 */
  dLambda = Xt->dLambda; /* Lambda */

  /* Fill the Lambda array for the iLevelIntegration */
  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++)
    Lambda[iNumLevel] = Lambda[iNumLevel] - dLambda;

  accu = dXt1 * exp(-dLambda * dDt) *
         LGMSVIntMult(iLevelIntegration, Lambda, 0, dDt);

  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++)
    Lambda[iNumLevel] = Lambda[iNumLevel] + dLambda;

  /* Fill the Lambda array for the iLevelIntegration+1 */
  Lambda[iLevelIntegration + 1] = dLambda;

  /* First sub function */
  if (Xt->a != 0) {
    if (Xt->bIsft1 == 1)
      accu += Xt->a * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pft = (LGMSVSolFunc *)(Xt->pft);

      LGMSVFuncValue2(pft, dDt, Lambda, iLevelIntegration + 1, &temp);
      accu += Xt->a * temp;
    }
  }

  /* Second sub function */
  if (Xt->b != 0) {
    if (Xt->bIsgt1 == 1)
      accu += Xt->b * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pgt = (LGMSVSolFunc *)(Xt->pgt);
      LGMSVFuncValue2(pgt, dDt, Lambda, iLevelIntegration + 1, &temp);
      accu += Xt->b * temp;
    }
  }

  /* Second sub function */
  if (Xt->c != 0) {
    if (Xt->bIsht1 == 1)
      accu += Xt->c * LGMSVIntMult(iLevelIntegration + 1, Lambda, 0, dDt);
    else {
      /* Casting */
      pht = (LGMSVSolFunc *)(Xt->pht);
      LGMSVFuncValue2(pht, dDt, Lambda, iLevelIntegration + 1, &temp);
      accu += Xt->c * temp;
    }
  }

  *res = accu;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/

/* -----------------------------------------------------------------------------------------------
        LGMSVIntMult
                Example p=2
                  Calculate
                        Integral between t1 and t of
                                Integral between t1 and s of (u-t1)^n *
   exp(-lambda[p]*(s-u))*du times exp(-lambda[p-1]*(t-s))*ds
   ------------------------------------------------------------------------------------------------
 */
static double LGMSVIntMult1(
    /* Inputs */
    int p,          /* Number of integration */
    double *lambda, /* Array of lambda */

    int n,     /* Power in the last integral */
    double dDt /* = t-t1 */
) {
  double x1, x2;
  int iLevel;

  if (p == 0) {
    /* Zero level of integration */
    return pow(dDt, n);
  } else {
    /* p>0 level of integration */
    if (fabs(lambda[p]) < 1.0e-06) {
      return 1.0 / (n + 1) * LGMSVIntMult1(p - 1, lambda, n + 1, dDt);
    } else {
      if (n == 0) {
        x1 = LGMSVIntMult1(p - 1, lambda, 0, dDt);
        for (iLevel = 0; iLevel < p; iLevel++) {
          lambda[iLevel] = lambda[iLevel] - lambda[p];
        }
        x2 = LGMSVIntMult1(p - 1, lambda, 0, dDt);
        for (iLevel = 0; iLevel < p; iLevel++) {
          lambda[iLevel] = lambda[iLevel] + lambda[p];
        }
        return 1.0 / lambda[p] * (x1 - exp(-lambda[p] * dDt) * x2);
      } else {
        x1 = LGMSVIntMult1(p - 1, lambda, n, dDt);
        x2 = LGMSVIntMult1(p, lambda, n - 1, dDt);
        return 1.0 / lambda[p] * (x1 - n * x2);
      }
    }
  }
}

/* -----------------------------------------------------------------------------------------------
        LGMSVFuncValue
                Evaluation of
                        Xt  = Xt1*exp(-lambda*(t-t1))+integral between t1 and t
   of fs*exp(-lambda*(s-t1))

                where fs is of the same type as Xt
   ------------------------------------------------------------------------------------------------
 */
double LGMSVFuncValue1(
    /* Inputs */
    LGMSVSolFunc1 *Xt, /*  Function to be valued */
    double dDt,        /* = t - t1 */
    double *Lambda,    /* Array of lambda */
    int iLevelIntegration /* 0 for valuation */) {
  double dLambda, accu;
  LGMSVSolFunc1 *pft;
  int iNumLevel, iNumSubFunction;

  /* Extract value from structure of Xt */
  dLambda = Xt->dLambda; /* Lambda */

  /* Fill the Lambda array for the iLevelIntegration */
  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++) {
    Lambda[iNumLevel] = Lambda[iNumLevel] - dLambda;
  }

  accu = Xt->dXt1 * exp(-dLambda * dDt) *
         LGMSVIntMult1(iLevelIntegration, Lambda, 0, dDt);

  for (iNumLevel = 0; iNumLevel <= iLevelIntegration; iNumLevel++) {
    Lambda[iNumLevel] = Lambda[iNumLevel] + dLambda;
  }

  /* Fill the Lambda array for the iLevelIntegration+1 */
  Lambda[iLevelIntegration + 1] = dLambda;

  /* For all sub function */
  for (iNumSubFunction = 0; iNumSubFunction < Xt->iNbFunction;
       iNumSubFunction++) {
    if (fabs(Xt->a[iNumSubFunction]) > 1.0E-12) {
      if (Xt->bIsft1[iNumSubFunction] == 1) {
        accu += Xt->a[iNumSubFunction] *
                LGMSVIntMult1(iLevelIntegration + 1, Lambda, 0, dDt);
      } else {
        /* Casting */
        pft = (LGMSVSolFunc1 *)(Xt->pft[iNumSubFunction]);

        accu += Xt->a[iNumSubFunction] *
                LGMSVFuncValue1(pft, dDt, Lambda, iLevelIntegration + 1);
      }
    }
  }

  return accu;
}

void ConvertTS_LGMSV_to_LGM(long sigma_n, double *sigma_time, double *sigma,
                            double lambda, double Tstar) {
  double t1, t2, val1, val2, lambda2;
  int i;

  lambda2 = 2.0 * lambda;
  t1 = 0.0;
  val1 = exp(-lambda2 * Tstar);

  for (i = 0; i < sigma_n; i++) {
    t2 = sigma_time[i];
    val2 = exp(-lambda2 * (Tstar - t2));
    sigma[i] /= sqrt((val2 - val1) / lambda2 / (t2 - t1));
    t1 = t2;
    val1 = val2;
  }
}

void ConvertAlphaRho_LGM_to_LGMSV(long sigma_n, double *sigma_time,
                                  double tstar, double lambda, double alpha,
                                  double gamma, double rho, double *new_alpha,
                                  double *new_rho) {
  double t1, t2, val1_old, val1_new, val2_old, val2_new, val12_old, val12_new;
  double lambda2, lambda12, lambda22, lambdasum, lambda_ratio, coef_lam;
  int i;

  lambda2 = lambda + gamma;
  lambda12 = 2.0 * lambda;
  lambda22 = 2.0 * lambda2;
  lambdasum = lambda + lambda2;
  lambda_ratio = lambda / lambda2;

  coef_lam = 4.0 * lambda * lambda2 / (lambda + lambda2) / (lambda + lambda2);

  t1 = 0.0;
  val1_old = exp(lambda12 * (t1 - tstar));
  val2_old = exp(lambda22 * (t1 - tstar));
  val12_old = exp(lambdasum * (t1 - tstar));

  for (i = 0; i < sigma_n; i++) {
    t2 = sigma_time[i];
    val1_new = exp(lambda12 * (t2 - tstar));
    val2_new = exp(lambda22 * (t2 - tstar));
    val12_new = exp(lambdasum * (t2 - tstar));

    new_alpha[i] = alpha * sqrt(lambda_ratio * (val2_new - val2_old) /
                                (val1_new - val1_old));
    new_rho[i] = rho * sqrt(coef_lam * (val12_new - val12_old) *
                            (val12_new - val12_old) / (val1_new - val1_old) /
                            (val2_new - val2_old));

    t1 = t2;
    val1_old = val1_new;
    val2_old = val2_new;
    val12_old = val12_new;
  }
}

void ConvertTS_LGM_to_LGMSV(long sigma_n, double *sigma_time, double *sigma,
                            double lambda, double Tstar,
                            /* extra for 2 Factor */
                            int one2F, double alpha, double gamma, double rho,
                            double *new_alpha, double *new_rho) {
  double t1, t2, val1_old, val1_new, lambda12;
  int i;

  lambda12 = 2.0 * lambda;

  t1 = 0.0;
  val1_old = exp(-lambda12 * Tstar);

  for (i = 0; i < sigma_n; i++) {
    t2 = sigma_time[i];
    val1_new = exp(-lambda12 * (Tstar - t2));
    sigma[i] *= sqrt((val1_new - val1_old) / lambda12 / (t2 - t1));
    t1 = t2;
    val1_old = val1_new;
  }

  if (one2F == 2) {
    ConvertAlphaRho_LGM_to_LGMSV(sigma_n, sigma_time, Tstar, lambda, alpha,
                                 gamma, rho, new_alpha, new_rho);
  }
}

Err Get_LGMSV_TermStructure(char *underlying, double **sigma_time,
                            double **sigma, double **alpha, double **rho,
                            double **lambdaeps, double *Tstar, long *sigma_n,
                            double *fixed_tau,
                            /* extra 2F */
                            int *one2F, double **alpha2F_sv, double *gamma2F_sv,
                            double **rho2F_sv)

{
  SrtUndPtr und;
  TermStruct *ts;
  long today;
  int i;

  double *sigma2 = NULL, *beta = NULL, *beta2 = NULL, *rhoTS = NULL,
         *TauDates = NULL, *Tau = NULL, *Tau2 = NULL;

  long nbTau;
  double tau;

  double alpha2F, gamma2F, rho2F;

  SrtMdlDim mdl_dim;

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

  if (get_mdltype_from_irund(und) != LGM_STOCH_VOL) {
    err = serror("Underlying %s is not of type LGMSV", underlying);
    goto FREE_RETURN;
  }

  ts = get_ts_from_irund(und);
  today = get_today_from_underlying(und);

  err = get_underlying_mdldim(und, &mdl_dim);
  if (err)
    goto FREE_RETURN;

  if (mdl_dim == ONE_FAC) {
    *one2F = 1;
  } else if (mdl_dim == TWO_FAC) {
    *one2F = 2;
  } else {
    err = "LGM dimension cannot be more than 2 for LGMSV";
    goto FREE_RETURN;
  }

  if (*one2F == 1) {
    err = srt_f_display_IRM_OneFac_TermStruct(ts, sigma_time, sigma, &beta,
                                              alpha, rho, lambdaeps, sigma_n,
                                              &TauDates, &Tau, &nbTau);
  } else {
    /* to be implemented */
    err = srt_f_display_IRM_TwoFac_TermStruct(ts, sigma_time, sigma, &beta,
                                              &sigma2, &beta2, &rhoTS, sigma_n,
                                              &TauDates, &Tau, &Tau2, &nbTau);

    alpha2F = sigma2[0] / *sigma[0];
    gamma2F = 1.0 / Tau2[0] - 1.0 / Tau[0];
    rho2F = rhoTS[0];

    for (i = 1; i < *sigma_n; i++) {
      if (fabs(sigma2[i] / (*sigma)[i] - alpha2F) < 1.0E-08 ||
          fabs(1.0 / Tau2[i] - 1.0 / Tau[i] - gamma2F) < 1.0E-08 ||
          fabs(rhoTS[i] - rho2F) < 1.0E-08) {
        err = "Alpha2F        , Gamma2F and Rho2F must be constant in LGMSV2F";
        goto FREE_RETURN;
      }
    }
  }

  if (err)
    goto FREE_RETURN;

  /* for the moment TStar is not saved */
  *Tstar = LGMSV_Tstar;

  /* Now transform dates into maturities and check tau	*/
  for (i = 0; i < *sigma_n; i++) {
    (*sigma_time)[i] = ((*sigma_time)[i] - today) / 365.0;
    (*alpha)[i] *= 2.0;
    (*lambdaeps)[i] *= 2.0;
  }

  tau = Tau[0];

  for (i = 1; i < nbTau; i++) {
    if (fabs(Tau[i] - tau) > 1.0E-08) {
      err = "Tau must be constant in the LGMSV term struct";
      goto FREE_RETURN;
    }
  }

  (*fixed_tau) = tau;

  if (*one2F == 2) {
    *alpha2F_sv = calloc(*sigma_n, sizeof(double));
    *rho2F_sv = calloc(*sigma_n, sizeof(double));

    if (!*alpha2F_sv || !*rho2F_sv) {
      err = "Memory allocation faillure in Get_LGMSV_TermStructure";
      goto FREE_RETURN;
    }
  } else {
    *alpha2F_sv = NULL;
    *rho2F_sv = NULL;
  }

  /* Conversion of the volatility */
  ConvertTS_LGM_to_LGMSV(*sigma_n, *sigma_time, *sigma, 1.0 / tau, *Tstar,
                         *one2F, alpha2F, gamma2F, rho2F, *alpha2F_sv,
                         *rho2F_sv);

FREE_RETURN:

  if (beta)
    free(beta);
  if (TauDates)
    free(TauDates);
  if (Tau)
    free(Tau);

  return err;
}

void Convert_Tstar_LGMSV(double lambda, long sigma_n, double *sigma_time,
                         double *sigma, double *alphaeps, double *lameps,
                         double *rhoeps, double init_tstar, double new_tstar) {
  double adjust_mean, adjust_vol;
  int i;

  adjust_vol = exp(-lambda * (new_tstar - init_tstar));
  adjust_mean = (1.0 - adjust_vol) / lambda;

  for (i = 0; i < sigma_n; i++) {
    lameps[i] += adjust_mean * alphaeps[i] * rhoeps[i] * sigma[i];
    sigma[i] *= adjust_vol;
  }
}

Err SrtInitIRMSVUnd(char *undName, /* und name */
                    char *ycname,  /* mkt name */

                    /* dimension */
                    int one2F,

                    /* volatility and tau */
                    double **sigma_datas, int num_sigma, int sigma_col,
                    double **tau_datas, int num_tau, int tau_col,

                    /* LGM 2F parameters */
                    double alpha, double gamma, double rho,

                    /* Stoch vol parameters */
                    double **smile_datas, int num_smile, int smile_col,

                    /* T* of the model */
                    double tstar) {
  int bCleanUpUndFlag = 1;

  SrtUndPtr pUndPtr;
  SrtIrDesc *pIrPtr;
  SrtUndListPtr und_list;
  SrtCurvePtr pYieldCurve;

  irm_sv *irmsv = NULL;

  char *ccy;
  long today;
  int i, j;

  Err err = NULL;

  /* first check the inputs */
  if (sigma_col != 2) {
    if (sigma_col != 1 || num_sigma > 1) {
      err = "sigma input must have two columns or one single value";
      goto FREE_RETURN;
    }
  }

  if (tau_col != 2) {
    if (tau_col != 1 || num_tau > 1) {
      err = "tau input must have two columns or one single value";
      goto FREE_RETURN;
    }
  }

  if ((one2F == 1 && smile_col != 4) || (one2F == 2 && smile_col != 5)) {
    err = "smile input must have 4 or 5 columns";
    goto FREE_RETURN;
  }

  /* Allocate the irmsv structure */
  irmsv = calloc(1, sizeof(irm_sv));

  /* Get the yield curve and the spot date from the yield curve name */
  pYieldCurve = lookup_curve(ycname);

  if (!pYieldCurve) {
    err = "Cannot find Yield Curve";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(pYieldCurve);
  ccy = get_curve_ccy(pYieldCurve);

  /* Create the new General underlying */

  pUndPtr = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));

  pUndPtr->underl_type = INTEREST_RATE_UND;
  strcpy(pUndPtr->underl_name, undName);
  strupper(pUndPtr->underl_name);
  strip_white_space(pUndPtr->underl_name);
  strcpy(pUndPtr->underl_lbl, "LGMSV_UND");
  pUndPtr->underl_ccy = ccy;

  /* Create the IR underlying */

  pIrPtr = calloc(1, sizeof(SrtIrDesc));

  strcpy(pIrPtr->mdl_lbl, "LGMSV_UND");
  pIrPtr->mdl_type = LGM_STOCH_VOL;

  if (one2F == 1) {
    pIrPtr->mdl_dim = ONE_FAC;
  } else if (one2F == 2) {
    pIrPtr->mdl_dim = TWO_FAC;
  } else {
    err = "LGMSV cannot accept more than two factors";
    goto FREE_RETURN;
  }

  strcpy(pIrPtr->yc_name, pYieldCurve->curve_name);
  pIrPtr->spec = irmsv;

  /* Init Structure */

  irmsv->one2F = one2F;
  irmsv->today = today;

  if (one2F == 2) {
    irmsv->alpha = alpha;
    irmsv->gamma = gamma;
    irmsv->rho = rho;
  } else {
    irmsv->alpha = 0.0;
    irmsv->gamma = 0.0;
    irmsv->rho = 0.0;
  }

  if (fabs(tstar) < 1.0E-08) {
    tstar = LGMSV_Tstar;
  }

  irmsv->tstar = tstar;

  /* volatility */

  if (sigma_col == 2) {
    i = 0;
    while (i < num_sigma && sigma_datas[0][i] < today) {
      i++;
    }

    if (i == num_sigma) {
      err = "No current sigma values in TS";
      goto FREE_RETURN;
    }

    for (j = i; j < num_sigma; j++) {
      //			if (!sigma_datas[0][j] || !sigma_datas[1][j])
      if (!sigma_datas[0][j]) {
        break;
      }
    }

    num_sigma = j - i;

    if (num_sigma == 0) {
      err = "No current sigma values in TS";
      goto FREE_RETURN;
    }

    irmsv->num_sigma = num_sigma;
    irmsv->sigma_dates = calloc(irmsv->num_sigma, sizeof(long));
    irmsv->sigma_times = calloc(irmsv->num_sigma, sizeof(double));
    irmsv->sigma = calloc(irmsv->num_sigma, sizeof(double));

    if (!irmsv->sigma_dates || !irmsv->sigma_times || !irmsv->sigma) {
      err = "Memory allocation faillure in SrtInitIRMSVUnd";
      goto FREE_RETURN;
    }

    for (j = 0; j < num_sigma; j++) {
      irmsv->sigma_dates[j] = DTOL(sigma_datas[0][i + j]);
      irmsv->sigma_times[j] = (sigma_datas[0][i + j] - today) * YEARS_IN_DAY;
      irmsv->sigma[j] = sigma_datas[1][i + j];
    }
  } else {
    irmsv->num_sigma = num_sigma;
    irmsv->sigma_dates = calloc(irmsv->num_sigma, sizeof(long));
    irmsv->sigma_times = calloc(irmsv->num_sigma, sizeof(double));
    irmsv->sigma = calloc(irmsv->num_sigma, sizeof(double));

    if (!irmsv->sigma_dates || !irmsv->sigma_times || !irmsv->sigma) {
      err = "Memory allocation faillure in SrtInitIRMSVUnd";
      goto FREE_RETURN;
    }

    irmsv->sigma_dates[0] = today + 365;
    irmsv->sigma_times[0] = 1.0;
    irmsv->sigma[0] = sigma_datas[0][0];
  }

  /* tau */
  if (tau_col == 2) {
    i = 0;
    while (i < num_tau && tau_datas[0][i] < today) {
      i++;
    }

    if (i == num_tau) {
      err = "No current tau values in TS";
      goto FREE_RETURN;
    }

    for (j = i; j < num_tau; j++) {
      if (!tau_datas[0][j] || !tau_datas[1][j]) {
        break;
      }
    }

    num_tau = j - i;

    if (num_tau == 0) {
      err = "No current tau values in TS";
      goto FREE_RETURN;
    }

    irmsv->num_tau = num_tau;
    irmsv->tau_dates = calloc(irmsv->num_tau, sizeof(long));
    irmsv->tau_times = calloc(irmsv->num_tau, sizeof(double));
    irmsv->tau = calloc(irmsv->num_tau, sizeof(double));

    if (!irmsv->tau_dates || !irmsv->tau_times || !irmsv->tau) {
      err = "Memory allocation faillure in SrtInitIRMSVUnd";
      goto FREE_RETURN;
    }

    for (j = 0; j < num_tau; j++) {
      irmsv->tau_dates[j] = DTOL(tau_datas[0][i + j]);
      irmsv->tau_times[j] = (tau_datas[0][i + j] - today) * YEARS_IN_DAY;
      irmsv->tau[j] = tau_datas[1][i + j];
    }
  } else {
    irmsv->num_tau = num_tau;
    irmsv->tau_dates = calloc(irmsv->num_tau, sizeof(long));
    irmsv->tau_times = calloc(irmsv->num_tau, sizeof(double));
    irmsv->tau = calloc(irmsv->num_tau, sizeof(double));

    if (!irmsv->tau_dates || !irmsv->tau_times || !irmsv->tau) {
      err = "Memory allocation faillure in SrtInitIRMSVUnd";
      goto FREE_RETURN;
    }

    irmsv->tau_dates[0] = irmsv->sigma_dates[0];
    irmsv->tau_times[0] = irmsv->sigma_times[0];
    irmsv->tau[0] = tau_datas[0][0];
  }

  /* SmileSV */
  i = 0;
  while (i < num_smile && smile_datas[0][i] < today) {
    i++;
  }

  if (i == num_smile) {
    err = "No current smile values in TS";
    goto FREE_RETURN;
  }

  for (j = i; j < num_smile; j++) {
    if (!smile_datas[0][j]) // || !smile_datas[1][j])
    {
      break;
    }
  }

  num_smile = j - i;

  if (num_smile == 0) {
    err = "No current smile values in TS";
    goto FREE_RETURN;
  }

  irmsv->num_smile = num_smile;
  irmsv->smile_dates = calloc(irmsv->num_smile, sizeof(long));
  irmsv->smile_times = calloc(irmsv->num_smile, sizeof(double));
  irmsv->alphaSV = calloc(irmsv->num_smile, sizeof(double));
  irmsv->lamSV = calloc(irmsv->num_smile, sizeof(double));
  irmsv->rhoSV = calloc(irmsv->num_smile, sizeof(double));

  if (!irmsv->smile_dates || !irmsv->smile_times || !irmsv->alphaSV ||
      !irmsv->lamSV || !irmsv->rhoSV) {
    err = "Memory allocation faillure in SrtInitIRMSVUnd";
    goto FREE_RETURN;
  }

  if (one2F == 2) {
    irmsv->rho2SV = calloc(irmsv->num_smile, sizeof(double));
    if (!irmsv->rho2SV) {
      err = "Memory allocation faillure in SrtInitIRMSVUnd";
      goto FREE_RETURN;
    }
  } else {
    irmsv->rho2SV = NULL;
  }

  if (num_smile == 1) {
    irmsv->smile_dates[0] = irmsv->sigma_dates[0];
    irmsv->smile_times[0] = irmsv->sigma_times[0];

    irmsv->alphaSV[0] = smile_datas[1][i];
    irmsv->lamSV[0] = smile_datas[2][i];
    irmsv->rhoSV[0] = smile_datas[3][i];

    if (one2F == 2) {
      irmsv->rho2SV[0] = smile_datas[4][i];
    }
  } else {
    for (j = 0; j < num_smile; j++) {
      irmsv->smile_dates[j] = DTOL(smile_datas[0][i + j]);
      irmsv->smile_times[j] = (smile_datas[0][i + j] - today) * YEARS_IN_DAY;

      irmsv->alphaSV[j] = smile_datas[1][i + j];
      irmsv->lamSV[j] = smile_datas[2][i + j];
      irmsv->rhoSV[j] = smile_datas[3][i + j];

      if (one2F == 2) {
        irmsv->rho2SV[j] = smile_datas[4][i + j];
      }
    }
  }

  pUndPtr->spec_desc = pIrPtr;

  /* Put the underlying into the depot */

  und_list = get_underlying_list();

  err = srt_f_lstins(und_list, pUndPtr->underl_name, 0.0, OBJ_PTR_UND,
                     (void *)pUndPtr, &irm_sv_free_und_struct,
                     &(pUndPtr->underl_ticker));

FREE_RETURN:

  if (err) {
    irm_sv_free_struct(irmsv);
  }

  return err;
}

void irm_sv_free_struct(IRM_SV irm_sv) {
  if (irm_sv) {
    if (irm_sv->sigma_dates)
      free(irm_sv->sigma_dates);
    irm_sv->sigma_dates = NULL;

    if (irm_sv->sigma_times)
      free(irm_sv->sigma_times);
    irm_sv->sigma_times = NULL;

    if (irm_sv->sigma)
      free(irm_sv->sigma);
    irm_sv->sigma = NULL;

    if (irm_sv->tau_dates)
      free(irm_sv->tau_dates);
    irm_sv->tau_dates = NULL;

    if (irm_sv->tau_times)
      free(irm_sv->tau_times);
    irm_sv->tau_times = NULL;

    if (irm_sv->tau)
      free(irm_sv->tau);
    irm_sv->tau = NULL;

    if (irm_sv->smile_dates)
      free(irm_sv->smile_dates);
    irm_sv->smile_dates = NULL;

    if (irm_sv->smile_times)
      free(irm_sv->smile_times);
    irm_sv->smile_times = NULL;

    if (irm_sv->alphaSV)
      free(irm_sv->alphaSV);
    irm_sv->alphaSV = NULL;

    if (irm_sv->lamSV)
      free(irm_sv->lamSV);
    irm_sv->lamSV = NULL;

    if (irm_sv->rhoSV)
      free(irm_sv->rhoSV);
    irm_sv->rhoSV = NULL;

    if (irm_sv->rho2SV)
      free(irm_sv->rho2SV);
    irm_sv->rho2SV = NULL;
  }
}

Err irm_sv_free_und_struct(SrtUndPtr pUndDesc) {
  SrtIrDesc *pSrtIrPtr;

  pSrtIrPtr = (SrtIrDesc *)(pUndDesc->spec_desc);
  irm_sv_free_struct(pSrtIrPtr->spec);
  free(pSrtIrPtr->spec);
  free(pUndDesc->spec_desc);
  free(pUndDesc);
  pUndDesc = NULL;
  return NULL;
}

Err irm_sv_get_struct_from_und(char *und, irm_sv **irmsv) {
  Err err = NULL;
  SrtUndPtr pUndPtr;
  SrtIrDesc *pIrPtr;

  // Get the underlying through its name and check it exists
  // Check on the underlying type
  pUndPtr = lookup_und(und);
  if (!pUndPtr) {
    err = "Undefined underlying";
    return err;
  }

  if (get_underlying_type(pUndPtr) != INTEREST_RATE_UND) {
    return serror("Underlying %s is not of type IR", und);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(pUndPtr) != LGM_STOCH_VOL) {
    err = serror("Underlying %s is not of type LGMSV", und);
    goto FREE_RETURN;
  }

  // Extract the information from the underlying
  pIrPtr = (SrtIrDesc *)(pUndPtr->spec_desc);
  *irmsv = (irm_sv *)(pIrPtr->spec);

FREE_RETURN:

  return err;
}

Err irm_sv_get_term_struct(char *und, int *one2F, int *num_sigma,
                           double **sigma_times, double **sigma, int *num_tau,
                           double **tau_times, double **tau, double *alpha,
                           double *gamma, double *rho, int *num_smile,
                           double **smile_times, double **alphaSV,
                           double **lamSV, double **rhoSV, double **rho2SV,
                           double *tstar) {
  Err err = NULL;
  irm_sv *irmsv;

  *sigma_times = NULL;
  *sigma = NULL;
  *tau_times = NULL;
  *tau = NULL;
  *smile_times = NULL;
  *alphaSV = NULL;
  *lamSV = NULL;
  *rhoSV = NULL;
  *rho2SV = NULL;

  err = irm_sv_get_struct_from_und(und, &irmsv);

  if (err)
    goto FREE_RETURN;

  *one2F = irmsv->one2F;
  *num_sigma = irmsv->num_sigma;
  *num_tau = irmsv->num_tau;
  *alpha = irmsv->alpha;
  *gamma = irmsv->gamma;
  *rho = irmsv->rho;
  *num_smile = irmsv->num_smile;
  *tstar = irmsv->tstar;

  /* Memory allocation */
  *sigma_times = calloc(*num_sigma, sizeof(double));
  *sigma = calloc(*num_sigma, sizeof(double));
  *tau_times = calloc(*num_tau, sizeof(double));
  *tau = calloc(*num_tau, sizeof(double));
  *smile_times = calloc(*num_smile, sizeof(double));
  *alphaSV = calloc(*num_smile, sizeof(double));
  *lamSV = calloc(*num_smile, sizeof(double));
  *rhoSV = calloc(*num_smile, sizeof(double));

  if (!*sigma_times || !*sigma_times || !*sigma_times || !*sigma_times ||
      !*sigma_times || !*sigma_times || !*sigma_times || !*sigma_times) {
    err = "Memomry allocation faillure in irm_sv_get_term_struct";
    goto FREE_RETURN;
  }

  memcpy(*sigma_times, irmsv->sigma_times, *num_sigma * sizeof(double));
  memcpy(*sigma, irmsv->sigma, *num_sigma * sizeof(double));
  memcpy(*tau_times, irmsv->tau_times, *num_tau * sizeof(double));
  memcpy(*tau, irmsv->tau, *num_tau * sizeof(double));

  memcpy(*smile_times, irmsv->smile_times, *num_smile * sizeof(double));
  memcpy(*alphaSV, irmsv->alphaSV, *num_smile * sizeof(double));
  memcpy(*lamSV, irmsv->lamSV, *num_smile * sizeof(double));
  memcpy(*rhoSV, irmsv->rhoSV, *num_smile * sizeof(double));

  if (*one2F == 2) {
    *rho2SV = calloc(*num_smile, sizeof(double));

    if (!*rho2SV) {
      err = "Memomry allocation faillure in irm_sv_get_term_struct";
      goto FREE_RETURN;
    }

    memcpy(*rho2SV, irmsv->rho2SV, *num_smile * sizeof(double));
  }

FREE_RETURN:

  if (err) {
    if (*sigma_times)
      free(*sigma_times);
    if (*sigma)
      free(*sigma);
    if (*tau_times)
      free(*tau_times);
    if (*tau)
      free(*tau);
    if (*smile_times)
      free(*smile_times);
    if (*alphaSV)
      free(*alphaSV);
    if (*lamSV)
      free(*lamSV);
    if (*rhoSV)
      free(*rhoSV);
    if (*rho2SV)
      free(*rho2SV);
  }

  return err;
}

Err irm_sv_get_term_struct_date(
    char *und, int *one2F, int *num_sigma, long **sigma_dates, double **sigma,
    int *num_tau, long **tau_dates, double **tau, double *alpha, double *gamma,
    double *rho, int *num_smile, long **smile_dates, double **alphaSV,
    double **lamSV, double **rhoSV, double **rho2SV, double *tstar) {
  Err err = NULL;
  irm_sv *irmsv;

  *sigma_dates = NULL;
  *sigma = NULL;
  *tau_dates = NULL;
  *tau = NULL;
  *smile_dates = NULL;
  *alphaSV = NULL;
  *lamSV = NULL;
  *rhoSV = NULL;
  *rho2SV = NULL;

  err = irm_sv_get_struct_from_und(und, &irmsv);

  if (err)
    goto FREE_RETURN;

  *one2F = irmsv->one2F;
  *num_sigma = irmsv->num_sigma;
  *num_tau = irmsv->num_tau;
  *alpha = irmsv->alpha;
  *gamma = irmsv->gamma;
  *rho = irmsv->rho;
  *num_smile = irmsv->num_smile;
  *tstar = irmsv->tstar;

  /* Memory allocation */
  *sigma_dates = calloc(*num_sigma, sizeof(double));
  *sigma = calloc(*num_sigma, sizeof(double));
  *tau_dates = calloc(*num_tau, sizeof(double));
  *tau = calloc(*num_tau, sizeof(double));
  *smile_dates = calloc(*num_smile, sizeof(double));
  *alphaSV = calloc(*num_smile, sizeof(double));
  *lamSV = calloc(*num_smile, sizeof(double));
  *rhoSV = calloc(*num_smile, sizeof(double));

  if (!*sigma_dates || !*sigma_dates || !*sigma_dates || !*sigma_dates ||
      !*sigma_dates || !*sigma_dates || !*sigma_dates || !*sigma_dates) {
    err = "Memomry allocation faillure in irm_sv_get_term_struct";
    goto FREE_RETURN;
  }

  memcpy(*sigma_dates, irmsv->sigma_dates, *num_sigma * sizeof(long));
  memcpy(*sigma, irmsv->sigma, *num_sigma * sizeof(double));
  memcpy(*tau_dates, irmsv->tau_dates, *num_tau * sizeof(long));
  memcpy(*tau, irmsv->tau, *num_tau * sizeof(double));

  memcpy(*smile_dates, irmsv->smile_dates, *num_smile * sizeof(long));
  memcpy(*alphaSV, irmsv->alphaSV, *num_smile * sizeof(double));
  memcpy(*lamSV, irmsv->lamSV, *num_smile * sizeof(double));
  memcpy(*rhoSV, irmsv->rhoSV, *num_smile * sizeof(double));

  if (*one2F == 2) {
    *rho2SV = calloc(*num_smile, sizeof(double));

    if (!*rho2SV) {
      err = "Memomry allocation faillure in irm_sv_get_term_struct";
      goto FREE_RETURN;
    }

    memcpy(*rho2SV, irmsv->rho2SV, *num_smile * sizeof(double));
  }

FREE_RETURN:

  if (err) {
    if (*sigma_dates)
      free(*sigma_dates);
    if (*sigma)
      free(*sigma);
    if (*tau_dates)
      free(*tau_dates);
    if (*tau)
      free(*tau);
    if (*smile_dates)
      free(*smile_dates);
    if (*alphaSV)
      free(*alphaSV);
    if (*lamSV)
      free(*lamSV);
    if (*rhoSV)
      free(*rhoSV);
    if (*rho2SV)
      free(*rho2SV);
  }

  return err;
}

void init_NULL_LGMSV_model(LGMSV_MODEL model) {
  if (model) {
    model->dPWTime = NULL;
    model->dSigma = NULL;
    model->dAlpha = NULL;
    model->dLambdaEps = NULL;
    model->dLvlEps = NULL;
    model->dRho = NULL;
    model->dLGMAlpha = NULL;
    model->dLGMRho = NULL;
    model->dRho2 = NULL;
  }
}

Err init_LGMSV_model(LGMSV_MODEL model, long lToday, int iOne2F, int iNbPWTime,
                     double dLambdaX, double *dPWTime, double *dSigma,
                     double *dAlpha, double *dLambdaEps, double *dLvlEps,
                     double *dRho, double dTStar, double dLGMAlpha,
                     double dLGMGamma, double dLGMRho, double *dRho2) {
  Err err = NULL;

  init_NULL_LGMSV_model(model);

  model->iNbPWTime = iNbPWTime;
  model->dLambdaX = dLambdaX;

  if (fabs(dLambdaX) > 1.0E-08) {
    model->dTau = 1.0 / dLambdaX;
  } else {
    err = "Invalid Mean reversion in init_LGMSV_model";
    goto FREE_RETURN;
  }

  model->iOne2F = iOne2F;
  model->lToday = lToday;
  model->dTStar = dTStar;
  model->dInitTStar = dTStar;
  model->lTStarDate =
      (long)(model->lToday + model->dTStar * DAYS_IN_YEAR + 1.0e-08);

  if (model->iOne2F == 2) {
    model->dInitLGMAlpha = dLGMAlpha;
    model->dLGMGamma = dLGMGamma;
    model->dInitLGMRho = dLGMRho;

    model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;

    if (fabs(model->dLambdaX2) > 1.0E-08) {
      model->dTau2 = 1.0 / model->dLambdaX2;
    } else {
      err = "Invalid Gamma in init_LGMSV_model";
      goto FREE_RETURN;
    }

    model->dLGMAlpha = calloc(iNbPWTime, sizeof(double));
    model->dLGMRho = calloc(iNbPWTime, sizeof(double));
    model->dRho2 = calloc(iNbPWTime, sizeof(double));

    if (!model->dLGMAlpha || !model->dLGMRho) {
      err = "Memory allocation faillure in init_LGMSV_model";
      goto FREE_RETURN;
    }

    ConvertAlphaRho_LGM_to_LGMSV(iNbPWTime, dPWTime, model->dTStar,
                                 model->dLambdaX, dLGMAlpha, dLGMGamma, dLGMRho,
                                 model->dLGMAlpha, model->dLGMRho);

    memcpy(model->dRho2, dRho2, iNbPWTime * sizeof(double));
  } else {
    model->dLGMAlpha = NULL;
    model->dLGMGamma = 0.0;
    model->dLGMRho = NULL;
    model->dRho2 = NULL;
  }

  model->dPWTime = calloc(iNbPWTime, sizeof(double));
  model->dSigma = calloc(iNbPWTime, sizeof(double));
  model->dAlpha = calloc(iNbPWTime, sizeof(double));
  model->dLambdaEps = calloc(iNbPWTime, sizeof(double));
  model->dLvlEps = calloc(iNbPWTime, sizeof(double));
  model->dRho = calloc(iNbPWTime, sizeof(double));

  if (!dPWTime || !dSigma || !dAlpha || !dLambdaEps || !dLvlEps || !dRho) {
    err = "Memory allocation faillure in init_LGMSV_model";
    goto FREE_RETURN;
  }

  memcpy(model->dPWTime, dPWTime, iNbPWTime * sizeof(double));
  memcpy(model->dSigma, dSigma, iNbPWTime * sizeof(double));
  memcpy(model->dAlpha, dAlpha, iNbPWTime * sizeof(double));
  memcpy(model->dLambdaEps, dLambdaEps, iNbPWTime * sizeof(double));
  memcpy(model->dLvlEps, dLvlEps, iNbPWTime * sizeof(double));
  memcpy(model->dRho, dRho, iNbPWTime * sizeof(double));

  /* Check the correlation matrix */
  err = Check_Model_Corr_Matrix(model);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (err)
    free_LGMSV_model(model);

  return err;
}

void free_LGMSV_model(LGMSV_MODEL model) {
  if (model) {
    if (model->dPWTime)
      free(model->dPWTime);
    model->dPWTime = NULL;
    if (model->dSigma)
      free(model->dSigma);
    model->dSigma = NULL;
    if (model->dAlpha)
      free(model->dAlpha);
    model->dAlpha = NULL;
    if (model->dLambdaEps)
      free(model->dLambdaEps);
    model->dLambdaEps = NULL;
    if (model->dLvlEps)
      free(model->dLvlEps);
    model->dLvlEps = NULL;
    if (model->dRho)
      free(model->dRho);
    model->dRho = NULL;

    if (model->dRho2)
      free(model->dRho2);
    model->dRho2 = NULL;
    if (model->dLGMAlpha)
      free(model->dLGMAlpha);
    model->dLGMAlpha = NULL;
    if (model->dLGMRho)
      free(model->dLGMRho);
    model->dLGMRho = NULL;

    model->iNbPWTime = 0;
  }
}

Err fill_LGMSV_model_from_irm_sv(irm_sv *irmsv, LGMSV_MODEL model) {
  Err err = NULL;
  int i, index;

  /* initialisation */
  model->dPWTime = NULL;
  model->dSigma = NULL;
  model->dAlpha = NULL;
  model->dLambdaEps = NULL;
  model->dLvlEps = NULL;
  model->dRho = NULL;
  model->dLGMAlpha = NULL;
  model->dLGMRho = NULL;
  model->dRho2 = NULL;

  /* Get the Tau */
  model->dTau = irmsv->tau[0];

  for (i = 1; i < irmsv->num_tau; i++) {
    if (fabs(irmsv->tau[i] - model->dTau) > 1.0E-08) {
      err = "LGMSV cannot accept Term Structure of Tau for the moment";
      goto FREE_RETURN;
    }
  }

  if (fabs(model->dTau) < 1.0E-08) {
    err = "Tau cannot be equal to 0.0 in LGMSV";
    goto FREE_RETURN;
  }

  model->dLambdaX = 1.0 / model->dTau;

  /* Get the volatility and smile parameters */
  model->iNbPWTime = irmsv->num_sigma;
  model->dPWTime = calloc(irmsv->num_sigma, sizeof(double));

  if (!model->dPWTime) {
    err = "Memory allocation faillure in fill_LGMSV_model_from_irm_sv";
    goto FREE_RETURN;
  }

  memcpy(model->dPWTime, irmsv->sigma_times, irmsv->num_sigma * sizeof(double));

  for (i = 0; i < irmsv->num_smile; i++) {
    num_f_add_number(&(model->iNbPWTime), &(model->dPWTime),
                     irmsv->smile_times[i]);
  }

  num_f_sort_vector(model->iNbPWTime, model->dPWTime);
  num_f_unique_vector(&(model->iNbPWTime), model->dPWTime);

  model->dSigma = calloc(model->iNbPWTime, sizeof(double));
  model->dAlpha = calloc(model->iNbPWTime, sizeof(double));
  model->dLambdaEps = calloc(model->iNbPWTime, sizeof(double));
  model->dLvlEps = calloc(model->iNbPWTime, sizeof(double));
  model->dRho = calloc(model->iNbPWTime, sizeof(double));

  if (!model->dSigma || !model->dAlpha || !model->dLambdaEps ||
      !model->dLvlEps || !model->dRho) {
    err = "Memory allocation faillure in fill_LGMSV_model_from_irm_sv";
    goto FREE_RETURN;
  }

  for (i = 0; i < model->iNbPWTime; i++) {
    model->dSigma[i] = irmsv->sigma[Get_Index(
        model->dPWTime[i], irmsv->sigma_times, irmsv->num_sigma)];
    index = Get_Index(model->dPWTime[i], irmsv->smile_times, irmsv->num_smile);

    model->dAlpha[i] = irmsv->alphaSV[index] * 2.0;
    model->dLambdaEps[i] = irmsv->lamSV[index] * 2.0;
    model->dLvlEps[i] = model->dLambdaEps[i];
    model->dRho[i] = irmsv->rhoSV[index];
  }

  model->lToday = irmsv->today;
  model->dTStar = irmsv->tstar;
  model->dInitTStar = irmsv->tstar;

  model->lTStarDate =
      (long)(model->lToday + model->dTStar * DAYS_IN_YEAR + 0.5);

  model->iOne2F = irmsv->one2F;

  if (model->iOne2F == 2) {
    model->dRho2 = calloc(model->iNbPWTime, sizeof(double));
    model->dLGMAlpha = calloc(model->iNbPWTime, sizeof(double));
    model->dLGMRho = calloc(model->iNbPWTime, sizeof(double));

    if (!model->dRho2 || !model->dLGMAlpha || !model->dLGMRho) {
      err = "Memory allocation faillure in fill_LGMSV_model_from_irm_sv";
      goto FREE_RETURN;
    }

    for (i = 0; i < model->iNbPWTime; i++) {
      model->dRho2[i] = irmsv->rho2SV[Get_Index(
          model->dPWTime[i], irmsv->smile_times, irmsv->num_smile)];
    }

    model->dInitLGMAlpha = irmsv->alpha;
    model->dLGMGamma = irmsv->gamma;
    model->dInitLGMRho = irmsv->rho;

    model->dLambdaX2 = model->dLambdaX + model->dLGMGamma;

    if (fabs(model->dLambdaX2) < 1.0E-08) {
      err = "Second Lambda cannot be 0.0 in LGMSV 2 Factor";
      goto FREE_RETURN;
    }

    model->dTau2 = 1.0 / model->dLambdaX2;
  }

  /* conversion of the vol        , alpha and rho */
  /* Conversion of the volatility */
  ConvertTS_LGM_to_LGMSV(model->iNbPWTime, model->dPWTime, model->dSigma,
                         model->dLambdaX, model->dTStar, model->iOne2F,
                         model->dInitLGMAlpha, model->dLGMGamma,
                         model->dInitLGMRho, model->dLGMAlpha, model->dLGMRho);

FREE_RETURN:

  if (err) {
    free_LGMSV_model(model);
  }

  return err;
}

Err Get_LGMSV_model(char *undname, LGMSV_MODEL model) {
  irm_sv *irmsv;
  Err err = NULL;

  err = irm_sv_get_struct_from_und(undname, &irmsv);

  if (err)
    return err;

  err = fill_LGMSV_model_from_irm_sv(irmsv, model);

  return err;
}

Err Check_Model_Input_Corr_Matrix(int iOne2F, int iNbPWTime, double *dLGMRho,
                                  double *dRho, double *dRho2) {
  Err err = NULL;

  double temp;
  int i;

  if (iOne2F == 2) {
    for (i = 0; i < iNbPWTime; i++) {
      temp = dLGMRho[i] * dLGMRho[i] + dRho[i] * dRho[i] + dRho2[i] * dRho2[i] -
             2.0 * dLGMRho[i] * dRho[i] * dRho2[i];

      if (temp >= 1.0) {
        err =
            "LGMSV Correlation matrix is not positive definite        , change "
            "Rho2 for example";
        return err;
      }
    }
  }

  return err;
}

Err Check_Model_Corr_Matrix(LGMSV_MODEL model) {
  Err err = NULL;

  err =
      Check_Model_Input_Corr_Matrix(model->iOne2F, model->iNbPWTime,
                                    model->dLGMRho, model->dRho, model->dRho2);

  return err;
}

void Convert_Tstar_model(LGMSV_MODEL model, double new_tstar) {
  double adjust_mean, adjust_vol, adjust_mean2, adjust_vol2, adjust_alpha;
  int i;

  adjust_vol = exp(-model->dLambdaX * (new_tstar - model->dTStar));
  adjust_mean = (1.0 - adjust_vol) / model->dLambdaX;

  if (model->iOne2F == 1) {
    for (i = 0; i < model->iNbPWTime; i++) {
      model->dLambdaEps[i] +=
          adjust_mean * model->dAlpha[i] * model->dRho[i] * model->dSigma[i];
      model->dSigma[i] *= adjust_vol;
    }
  } else if (model->iOne2F == 2) {
    adjust_vol2 = exp(-model->dLambdaX2 * (new_tstar - model->dTStar));
    adjust_mean2 = (1.0 - adjust_vol2) / model->dLambdaX2;
    adjust_alpha = adjust_vol2 / adjust_vol;

    for (i = 0; i < model->iNbPWTime; i++) {
      model->dLambdaEps[i] +=
          model->dAlpha[i] * model->dSigma[i] *
          (model->dRho[i] * adjust_mean +
           model->dRho2[i] * adjust_mean2 * model->dLGMAlpha[i]);

      model->dSigma[i] *= adjust_vol;
      model->dLGMAlpha[i] *= adjust_alpha;
    }
  }

  /* update TStar */
  model->dTStar = new_tstar;
  model->lTStarDate =
      (long)(model->lToday + DAYS_IN_YEAR * model->dTStar + 1.0E-08);
}
