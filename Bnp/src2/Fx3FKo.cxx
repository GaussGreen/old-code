/* ==========================================================================
   FILE_NAME:	Fx3FKo.cxx

   PURPOSE:		Closed Form pricing of KO Power Duals

   DATE:		08/28/00
   ========================================================================== */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err Fx3FKoDig(char *undName, int nKo, long *koDates, double *koLvl,
              int *ab, /*	1: Above        , -1: Below */
              long payDate, long npth, double *ans) {
  double *koMat = NULL, payMat;

  char *domYcName;
  SrtUndPtr und;
  long today;

  double spotFx, df, *fwdFx = NULL, *logBar = NULL, **covar = NULL;

  int i, j;

  Err err = NULL;

  /*	Allocation */

  koMat = (double *)calloc(nKo, sizeof(double));
  fwdFx = (double *)calloc(nKo, sizeof(double));
  logBar = (double *)calloc(nKo, sizeof(double));
  covar = dmatrix(0, nKo - 1, 0, nKo - 1);

  if (!koMat || !fwdFx || !logBar || !covar) {
    err = "Memory allocation error in Fx3FKoDig";
    goto FREE_RETURN;
  }

  /*	Read data outside loop */
  und = lookup_und(undName);
  if (!und) {
    err = "Underlying doesn't exist in Fx3FKoDig";
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(und);

  err = Fx3DFwdFx(undName, today, &spotFx);
  if (err) {
    goto FREE_RETURN;
  }

  payMat = (payDate - today) * YEARS_IN_DAY;

  domYcName = get_discname_from_fxund(und);
  df = swp_f_df(today, payDate, domYcName);

  /*	Read data */
  for (i = 0; i < nKo; i++) {
    koMat[i] = (koDates[i] - today) * YEARS_IN_DAY;
    if (i == 0) {
      if (koDates[i] <= today) {
        err = "Barrier dates must be > today in Fx3FKoDig";
        goto FREE_RETURN;
      }
    } else {
      if (koDates[i] <= koDates[i - 1]) {
        err = "Barrier dates must be increasing in Fx3FKoDig";
        goto FREE_RETURN;
      }
    }

    logBar[i] = log(koLvl[i] / spotFx);

    err = Fx3DFwdTpay(undName, koMat[i], koMat[i], payMat, 0, &(fwdFx[i]));
    if (err) {
      goto FREE_RETURN;
    }
  }

  for (i = 0; i < nKo; i++) {
    for (j = i; j < nKo; j++) {
      err = Fx3DCovar(undName, koMat[i], koMat[i], koMat[j], koMat[j],
                      &(covar[i][j]));
      if (err) {
        goto FREE_RETURN;
      }
      covar[j][i] = covar[i][j];
    }
  }

  /*	Main calculation */
  *ans = df * num_f_multinorm(nKo, fwdFx, covar, logBar, ab, npth);

  /*	Free and go */

FREE_RETURN:

  if (koMat)
    free(koMat);
  if (fwdFx)
    free(fwdFx);
  if (logBar)
    free(logBar);
  if (covar)
    free_dmatrix(covar, 0, nKo - 1, 0, nKo - 1);

  return err;
}

Err Fx3FKoFwd(char *undName, int nKo, long *koDates, double *koLvl,
              int *ab, /*	1: Above        , -1: Below */
              long fixDate, long payDate, long npth, double *ans) {
  double *koMat = NULL, fixMat, payMat;

  char *domYcName;
  SrtUndPtr und;
  long today;

  double spotFx, paidFwdFx, ivol, df, *fwdFx = NULL, *logBar = NULL,
                                      **covar = NULL;

  int i, j;

  Err err = NULL;

  /*	Allocation */

  koMat = (double *)calloc(nKo, sizeof(double));
  fwdFx = (double *)calloc(nKo, sizeof(double));
  logBar = (double *)calloc(nKo, sizeof(double));
  covar = dmatrix(0, nKo - 1, 0, nKo - 1);

  if (!koMat || !fwdFx || !logBar || !covar) {
    err = "Memory allocation error in Fx3FKoFwd";
    goto FREE_RETURN;
  }

  /*	Read data outside loop */
  und = lookup_und(undName);
  if (!und) {
    err = "Underlying doesn't exist in Fx3FKoFwd";
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(und);

  err = Fx3DFwdFx(undName, today, &spotFx);
  if (err) {
    goto FREE_RETURN;
  }

  payMat = (payDate - today) * YEARS_IN_DAY;
  fixMat = (fixDate - today) * YEARS_IN_DAY;

  domYcName = get_discname_from_fxund(und);
  df = swp_f_df(today, payDate, domYcName);

  err = Fx3DFwdTpay(undName, fixMat, fixMat, payMat, 0, &paidFwdFx);
  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DImpliedVol(undName, fixMat, 0.0, fixMat, &ivol);
  if (err) {
    goto FREE_RETURN;
  }
  paidFwdFx = spotFx * exp(paidFwdFx + 0.5 * ivol * ivol * fixMat);

  /*	Read data */
  for (i = 0; i < nKo; i++) {
    koMat[i] = (koDates[i] - today) * YEARS_IN_DAY;
    if (i == 0) {
      if (koDates[i] <= today) {
        err = "Barrier dates must be > today in Fx3FKoFwd";
        goto FREE_RETURN;
      }
    } else {
      if (koDates[i] <= koDates[i - 1]) {
        err = "Barrier dates must be increasing in Fx3FKoFwd";
        goto FREE_RETURN;
      }
    }

    logBar[i] = log(koLvl[i] / spotFx);

    err = Fx3DFwdQf(undName, koMat[i], koMat[i], fixMat, fixMat, payMat, 0,
                    &(fwdFx[i]));
    if (err) {
      goto FREE_RETURN;
    }
  }

  if (koDates[nKo - 1] > fixDate) {
    err = "Fixing date must be after all barrier dates in Fx3FKoFwd";
    goto FREE_RETURN;
  }

  for (i = 0; i < nKo; i++) {
    for (j = i; j < nKo; j++) {
      err = Fx3DCovar(undName, koMat[i], koMat[i], koMat[j], koMat[j],
                      &(covar[i][j]));

      if (err) {
        goto FREE_RETURN;
      }
      covar[j][i] = covar[i][j];
    }
  }

  /*	Main calculation */
  *ans = df * paidFwdFx * num_f_multinorm(nKo, fwdFx, covar, logBar, ab, npth);

  /*	Free and go */

FREE_RETURN:

  if (koMat)
    free(koMat);
  if (fwdFx)
    free(fwdFx);
  if (logBar)
    free(logBar);
  if (covar)
    free_dmatrix(covar, 0, nKo - 1, 0, nKo - 1);

  return err;
}

Err Fx3FKoOption(char *undName, int nKo, long *koDates, double *koLvl,
                 int *ab,                    /*	1: Above        , -1: Below */
                 double strike, int callPut, /*	1: Call        , -1: Put */
                 long fixDate, long payDate, long npth, double *ans) {
  long *newKoDates = NULL;
  double *newKoLvl = NULL;
  int *newAb = NULL;
  double koFwd, koDig;
  Err err = NULL;

  newKoDates = (long *)calloc(nKo + 1, sizeof(long));
  newKoLvl = (double *)calloc(nKo + 1, sizeof(double));
  newAb = (int *)calloc(nKo + 1, sizeof(int));

  if (!newKoDates || !newKoLvl || !newAb) {
    err = "Memory allocation error in Fx3FKoOption";
    goto FREE_RETURN;
  }

  memcpy(newKoDates, koDates, nKo * sizeof(long));
  memcpy(newKoLvl, koLvl, nKo * sizeof(double));
  memcpy(newAb, ab, nKo * sizeof(int));

  if (fixDate == koDates[nKo - 1]) {
    if (ab[nKo - 1] == 1 && callPut == 1) {
      newKoLvl[nKo - 1] = (strike > koLvl[nKo - 1] ? strike : koLvl[nKo - 1]);
    } else if (ab[nKo - 1] == -1 && callPut == -1) {
      newKoLvl[nKo - 1] = (strike < koLvl[nKo - 1] ? strike : koLvl[nKo - 1]);
    } else if (ab[nKo - 1] == -1 && callPut == 1) {
      if (strike > koLvl[nKo - 1]) {
        *ans = 0.0;
        goto FREE_RETURN;
      } else {
        nKo++;
        fixDate++;
        newKoDates[nKo - 1] = fixDate;
        newKoLvl[nKo - 1] = strike;
        newAb[nKo - 1] = callPut;
      }
    } else {
      if (strike < koLvl[nKo - 1]) {
        *ans = 0.0;
        goto FREE_RETURN;
      } else {
        nKo++;
        fixDate++;
        newKoDates[nKo - 1] = fixDate;
        newKoLvl[nKo - 1] = strike;
        newAb[nKo - 1] = callPut;
      }
    }
  } else {
    nKo++;
    newKoDates[nKo - 1] = fixDate;
    newKoLvl[nKo - 1] = strike;
    newAb[nKo - 1] = callPut;
  }

  err = Fx3FKoFwd(undName, nKo, newKoDates, newKoLvl, newAb, fixDate, payDate,
                  npth, &koFwd);
  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3FKoDig(undName, nKo, newKoDates, newKoLvl, newAb, payDate, npth,
                  &koDig);
  if (err) {
    goto FREE_RETURN;
  }

  *ans = callPut * (koFwd - strike * koDig);

FREE_RETURN:

  if (newKoDates)
    free(newKoDates);
  if (newKoLvl)
    free(newKoLvl);
  if (newAb)
    free(newAb);

  return err;
}