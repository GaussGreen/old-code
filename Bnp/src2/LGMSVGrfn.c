#include "LGMSVgrfn.h"
#include "LGMSVcalib.h"
#include "math.h"
#include "srt_h_all.h"

/****************************************************************
 *			MODEL UNDER QTstar g(t)	piecewise constant *
 ****************************************************************/

Err payoff_lgmsv_pde(double evt_date, double evt_time, void *func_parm,

                     /* Market data	*/
                     void *yc,

                     /* Model data	*/
                     double lamx, double dTstar, /* In time */

                     /* Grid data	*/
                     int lphi, int uphi, int lx, int ux, int leps, int ueps,
                     double *phi, double *ftTstar,

                     /* Vector of results to be updated */
                     int nprod,
                     /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
                     double ****prod_val) {
  GRFNPARMLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpphi, tmpx;
  static int i, j, k, l, m;

  static Err err = NULL;

  double dDerBeta, dBetaDerBeta, dBTTstar, dBeta, dBeta2, coef;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV)func_parm;
  global = total->global;
  local = total->local;

  /*	Pre calculations		*/
  if (total->num_df > 0) {
    for (k = 0; k < total->num_df; k++) {
      total->dff[k] = swp_f_df(evt_date, total->df_dts[k], (char *)yc);
      total->gam[k] = (1.0 - exp(-lamx * total->df_tms[k])) / lamx;
      total->gam2[k] = 0.5 * total->gam[k] * total->gam[k];
    }
  }

  /* Calculation of beta(t  ,Tstar)*D/DTstar ( beta(t  ,Tstar)) and D/DTstar (
   * beta(t  ,Tstar)) */
  dDerBeta = exp(-lamx * (dTstar - evt_time));
  dBeta = (1 - dDerBeta) / lamx;
  dBetaDerBeta = dDerBeta * dBeta;
  dBeta2 = dBeta * dBeta;

  /* Multiplication of the old payoff by B(EventDate  ,Tstar) to  obtain the PV
   * at the event date */
  dBTTstar =
      swp_f_df(evt_date, evt_date + (dTstar - evt_time) * 365.0, (char *)yc);
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      tmpx = (ftTstar[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      coef = exp(-0.5 * dBeta2 * tmpphi - dBeta * tmpx);
      for (l = leps; l <= ueps; l++)
        for (m = 0; m < nprod; m++)
          prod_val[i][j][l][m] *= dBTTstar * coef;
    }
  }

  /*	Evaluation */
  /*	Adding the Cash Flows at the event date */
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      /*	Re-transform into x */
      tmpx = (ftTstar[j] - dBetaDerBeta * tmpphi) / dDerBeta;

      /* ... */
      for (k = 0; k < total->num_df; k++) {
        local->evt->df[0][k] = total->dff[k] * exp(-total->gam[k] * tmpx -
                                                   total->gam2[k] * tmpphi);
      }

      for (l = leps; l <= ueps; l++) {
        err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL,
                             prod_val[i][j][l], &temp);
        if (err) {
          return err;
        }
      }
    }
  }
  /* Division of the new payoff by B(EventDate  ,Tstar)  */
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      tmpx = (ftTstar[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      coef = exp(-0.5 * dBeta2 * tmpphi - dBeta * tmpx);
      for (l = leps; l <= ueps; l++)
        for (m = 0; m < nprod; m++)
          prod_val[i][j][l][m] /= dBTTstar * coef;
    }
  };

  return err;
}

Err payoff_lgmsv_rec_swaption(
    double evt_date, double evt_time, void *func_parm,

    /* Market data	*/
    void *yc,

    /* Model data	*/
    double lamx, double dTstar, /* In time */

    /* Grid data	*/
    int lphi, int uphi, int lx, int ux, int leps, int ueps, double *phi,
    double *x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
    double ****prod_val) {

  static LGMSV_SWAP_PARAM param;

  static int i, j, k;
  static double temp, tmpphi, tmpx, level, floating, iv, df;
  static double dDerBeta, dBetaDerBeta, dBTTstar, dBeta, dBeta2, coef;
  static int l, m;

  param = (LGMSV_SWAP_PARAM)(func_parm);

  /* Calculation of beta(t  ,Tstar)*D/DTstar ( beta(t  ,Tstar)) and D/DTstar (
   * beta(t  ,Tstar)) */
  dDerBeta = exp(-lamx * (dTstar - evt_time));
  dBeta = (1 - dDerBeta) / lamx;
  dBetaDerBeta = dDerBeta * dBeta;
  dBeta2 = dBeta * dBeta;
  dBTTstar =
      swp_f_df(evt_date, evt_date + (dTstar - evt_time) * 365.0, (char *)yc);

  /*	Evaluation			*/
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      /*	Re-transform into x */
      tmpx = (x[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      /* ... */

      /* First Coupon */

      k = param->start_index;

      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);

      floating = df;
      level = 0.0;

      k++;

      /* Middle Coupon */
      for (k; k < param->ncpn - 1; k++) {
        df = param->cpn_df[k] *
             exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
        level += df * param->cpn_cvg[k];
      }

      /* Last Coupon */
      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
      level += df * param->cpn_cvg[k];

      floating -= df;

      iv = max(param->strike * level - floating, 0.0);

      for (k = leps; k <= ueps; k++) {
        prod_val[i][j][k][0] = iv;
      }
    }
  }

  /* Division of the new payoff by B(EventDate  ,Tstar)  */
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      tmpx = (x[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      coef = exp(-0.5 * dBeta2 * tmpphi - dBeta * tmpx);
      for (l = leps; l <= ueps; l++)
        for (m = 0; m < nprod; m++)
          prod_val[i][j][l][m] /= dBTTstar * coef;
    }
  };

  return NULL;
}

/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/************************************************************
 *						MODEL UNDER QBeta
 **
 ************************************************************/
Err payoff_lgmsv_pde_QBeta(double evt_date, double evt_time, void *func_parm,

                           /* Market data	*/
                           void *yc,

                           /* Model data	*/
                           double lamx,

                           /* Grid data	*/
                           int lphi, int uphi, int lx, int ux, int leps,
                           int ueps, double *phi, double *x,

                           /* Vector of results to be updated */
                           int nprod,
                           /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
                           double ****prod_val) {
  GRFNPARMLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpphi, tmpx;
  static int i, j, k, l, m;

  static Err err = NULL;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV)func_parm;
  global = total->global;
  local = total->local;

  /*	Pre calculations		*/
  if (total->num_df > 0) {
    for (k = 0; k < total->num_df; k++) {
      total->dff[k] = swp_f_df(evt_date, total->df_dts[k], (char *)yc);
      total->gam[k] = (1.0 - exp(-lamx * total->df_tms[k])) / lamx;
      total->gam2[k] = 0.5 * total->gam[k] * total->gam[k];
    }
  }

  /*	Evaluation			*/
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      /*	Re-transform into x */
      tmpx = x[j];
      /* ... */
      for (k = 0; k < total->num_df; k++) {
        local->evt->df[0][k] = total->dff[k] * exp(-total->gam[k] * tmpx -
                                                   total->gam2[k] * tmpphi);
      }

      for (l = leps; l <= ueps; l++) {
        err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL,
                             prod_val[i][j][l], &temp);
        if (err) {
          return err;
        }
      }
    }
  }

  return err;
}

Err payoff_lgmsv_rec_swaption_QBeta(
    double evt_date, double evt_time, void *func_parm,

    /* Market data	*/
    void *yc,

    /* Model data	*/
    double lamx,

    /* Grid data	*/
    int lphi, int uphi, int lx, int ux, int leps, int ueps, double *phi,
    double *x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
    double ****prod_val) {

  static LGMSV_SWAP_PARAM param;

  static int i, j, k;
  static double temp, tmpphi, tmpx, level, floating, iv, df;

  param = (LGMSV_SWAP_PARAM)(func_parm);

  /*	Evaluation			*/
  for (i = lphi; i <= uphi; i++) {
    tmpphi = phi[i];
    for (j = lx; j <= ux; j++) {
      /*	Re-transform into x */
      tmpx = x[j];
      /* ... */

      /* First Coupon */

      k = param->start_index;

      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);

      floating = df;
      level = 0.0;

      k++;

      /* Middle Coupon */
      for (k; k < param->ncpn - 1; k++) {
        df = param->cpn_df[k] *
             exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
        level += df * param->cpn_cvg[k];
      }

      /* Last Coupon */
      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
      level += df * param->cpn_cvg[k];

      floating -= df;

      iv = max(param->strike * level - floating, 0.0);

      for (k = leps; k <= ueps; k++) {
        prod_val[i][j][k][0] = iv;
      }
    }
  }

  return NULL;
}

/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/********************************************************************************************
 *			MODEL UNDER QTstar u(t)=g(t)exp(-lambdaX*(Tstar-t))
 *piecewise constant		    *
 ********************************************************************************************/

Err payoff_lgmsv_pde_UtPieceWise(
    double evt_date, double evt_time, void *func_parm,

    /* Market data	*/
    void *yc,

    /* Model data	*/
    double lamx, double dTstar, /* In time */
    double alpha, double rho, double ut,

    /* Grid data	*/
    int lpsi, int upsi, int lx, int ux, int lz, int uz, double *psi, double *x,
    double *z,

    /* Vector of results to be updated */
    int *nprod,
    /* 4 dimensions : Psit  ,Xt  ,Zt  ,Product	*/
    double ****prod_val) {
  GRFNPARMLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmppsi, tmpf, tmp;
  static int k, j, m, i, l;

  static Err err = NULL;

  static double dDerBeta, dBeta, dBetaDivDerBeta, dBetaDivDerBeta2_div2;
  static double dDerBetalamx, dLnBTTstar, dUtRhoDivAlpha, dlnCoef, coef,
      coef_tmp;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV)func_parm;
  global = total->global;
  local = total->local;

  /* Calculation of beta(t  ,Tstar)*D/DTstar ( beta(t  ,Tstar)) and D/DTstar (
   * beta(t  ,Tstar)) */
  dDerBeta = exp(-lamx * (dTstar - evt_time));
  dBeta = (1.0 - dDerBeta) / lamx;

  dDerBetalamx = lamx * dDerBeta;

  dBetaDivDerBeta = dBeta / dDerBeta;
  dBetaDivDerBeta2_div2 = dBetaDivDerBeta * dBetaDivDerBeta / 2.0;
  dUtRhoDivAlpha = ut * rho / alpha;

  /* Calculation of Log(B(EventDate  ,Tstar)) */
  dLnBTTstar = log(swp_f_df(
      evt_date, evt_date + (dTstar - evt_time) * DAYS_IN_YEAR, (char *)yc));

  /*	Pre calculations for calculation of the B(t  ,T) */
  if (total->num_df > 0) {
    for (m = 0; m < total->num_df; m++) {
      total->dff[m] = log(swp_f_df(evt_date, total->df_dts[m], (char *)yc));
      total->gam[m] = (1.0 - exp(-lamx * total->df_tms[m])) / dDerBetalamx;
      total->gam2[m] = total->gam[m] * (dBetaDivDerBeta - 0.5 * total->gam[m]);
    }
  }

  /*	Evaluation */
  /*	Adding the Cash Flows at the event date */

  for (i = lz; i <= uz; i++) {
    tmp = dUtRhoDivAlpha * (z[i] - 1.0);

    for (j = lx; j <= ux; j++) {
      tmpf = x[j] + tmp;
      coef_tmp = dLnBTTstar - dBetaDivDerBeta * tmpf;

      for (k = lpsi; k <= upsi; k++) {
        /* Multiplication of the old payoff by B(EventDate  ,Tstar) to  obtain
         * the PV at the event date */
        coef = exp(coef_tmp + dBetaDivDerBeta2_div2 * psi[k]);

        for (l = 0; l < *nprod; l++)
          prod_val[i][j][k][l] *= coef;

        /* Calculation of the usefull Df for the evt date */
        for (m = 0; m < total->num_df; m++) {
          local->evt->df[0][m] = exp(total->dff[m] - total->gam[m] * tmpf +
                                     total->gam2[m] * psi[k]);
        }

        err = FIRSTEvalEvent(global, local, *nprod, 2, NULL, NULL,
                             prod_val[i][j][k], &temp);
        if (err) {
          return err;
        }

        /* Division of the new payoff by B(EventDate  ,Tstar)  */
        for (l = 0; l < *nprod; l++)
          prod_val[i][j][k][l] /= coef;
      }
    }
  }

  return err;
}

Err payoff_lgmsv_FFT(double evt_date, double evt_time, void *func_parm,

                     /* Market data	*/
                     void *yc,

                     /* Model data	*/
                     double lamx, double dTstar, /* In time */

                     /* Grid data	*/
                     int iNbPhi, int iNbft,

                     int iIndexPhiMean, int iIndexft0,

                     double dPhiStep, double dftStep, double dPhitMean,

                     /* Vector of results to be updated */
                     int nprod,
                     /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
                     double ***prod_val) {
  GRFNPARMLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  double temp, tmpphi, tmpx, tmpft;
  int iNumPhi, iNumft, k;
  int nbPhiHalf, nbftHalf;
  double JumpPhi, Jumpft;

  Err err = NULL;

  double dDerBeta, dBetaDerBeta, dBTTstar, dBeta, dBeta2, coef;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV)func_parm;
  global = total->global;
  local = total->local;

  /*	Pre calculations		*/
  if (total->num_df > 0) {
    for (k = 0; k < total->num_df; k++) {
      total->dff[k] = swp_f_df(evt_date, total->df_dts[k], (char *)yc);
      total->gam[k] = (1.0 - exp(-lamx * total->df_tms[k])) / lamx;
      total->gam2[k] = 0.5 * total->gam[k] * total->gam[k];
    }
  }

  /* Calculation of beta(t  ,Tstar)*D/DTstar ( beta(t  ,Tstar)) and D/DTstar (
   * beta(t  ,Tstar)) */
  dDerBeta = exp(-lamx * (dTstar - evt_time));
  dBeta = (1 - dDerBeta) / lamx;
  dBetaDerBeta = dDerBeta * dBeta;
  dBeta2 = dBeta * dBeta;

  /*	Evaluation 	*/

  dBTTstar =
      swp_f_df(evt_date, evt_date + (dTstar - evt_time) * 365.0, (char *)yc);

  nbPhiHalf = iNbPhi - iIndexPhiMean + 1;
  nbftHalf = iNbft - iIndexft0 + 1;
  JumpPhi = iNbPhi * dPhiStep;
  Jumpft = iNbft * dftStep;

  tmpphi = dPhitMean;

  for (iNumPhi = 1; iNumPhi <= iNbPhi; iNumPhi++) {
    tmpft = 0.0;

    for (iNumft = 1; iNumft <= iNbft; iNumft++) {
      tmpx = (tmpft - dBetaDerBeta * tmpphi) / dDerBeta;
      coef = dBTTstar * exp(-0.5 * dBeta2 * tmpphi - dBeta * tmpx);

      for (k = 0; k < total->num_df; k++) {
        local->evt->df[0][k] = total->dff[k] * exp(-total->gam[k] * tmpx -
                                                   total->gam2[k] * tmpphi);
      }

      err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL,
                           prod_val[iNumPhi][iNumft], &temp);
      if (err) {
        return err;
      }

      for (k = 0; k < nprod; k++) {
        prod_val[iNumPhi][iNumft][k] /= coef;
      }

      tmpft += dftStep;

      if (iNumft == nbftHalf) {
        tmpft -= Jumpft;
      }
    }

    tmpphi += dPhiStep;

    if (iNumPhi == nbPhiHalf) {
      tmpphi -= JumpPhi;
    }
  }

  return err;
}

Err payoff_lgmsv_rec_swaption_UtPieceWise(
    double evt_date, double evt_time, void *func_parm,

    /* Market data	*/
    void *yc,

    /* Model data	*/
    double lamx, double dTstar, /* In time */

    /* Grid data	*/
    int lpsi, int upsi, int lf, int uf, int lz, int uz, double *psi, double *x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Psit  ,Xt  ,Epst  ,Product	*/
    double ****prod_val) {

  static LGMSV_SWAP_PARAM param;

  static int i, j, k;
  static double temp, tmpphi, tmpx, level, floating, iv, df;
  static double dDerBeta, dBetaDerBeta, dBTTstar, dBeta, dBeta2, coef;
  static int l, m;
  static double dPsi2Phi;

  param = (LGMSV_SWAP_PARAM)(func_parm);

  /* Calculation of beta(t  ,Tstar)*D/DTstar ( beta(t  ,Tstar)) and D/DTstar (
   * beta(t  ,Tstar)) */
  dDerBeta = exp(-lamx * (dTstar - evt_time));
  dBeta = (1 - dDerBeta) / lamx;
  dBetaDerBeta = dDerBeta * dBeta;
  dBeta2 = dBeta * dBeta;

  dPsi2Phi = 1.0 / dDerBeta / dDerBeta;

  dBTTstar =
      swp_f_df(evt_date, evt_date + (dTstar - evt_time) * 365.0, (char *)yc);

  /*	Evaluation			*/
  for (i = lpsi; i <= upsi; i++) {
    tmpphi = psi[i] * dPsi2Phi;
    for (j = lf; j <= uf; j++) {
      /*	Re-transform into x */
      tmpx = (x[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      /* ... */

      /* First Coupon */

      k = param->start_index;

      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);

      floating = df;
      level = 0.0;

      k++;

      /* Middle Coupon */
      for (k; k < param->ncpn - 1; k++) {
        df = param->cpn_df[k] *
             exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
        level += df * param->cpn_cvg[k];
      }

      /* Last Coupon */
      df = param->cpn_df[k] *
           exp(-param->gam[k] * tmpx - param->gam2[k] * tmpphi);
      level += df * param->cpn_cvg[k];

      floating -= df;

      iv = max(param->strike * level - floating, 0.0);

      for (k = lz; k <= uz; k++) {
        prod_val[i][j][k][0] = iv;
      }
    }
  }

  /* Division of the new payoff by B(EventDate  ,Tstar)  */
  for (i = lpsi; i <= upsi; i++) {
    tmpphi = psi[i] * dPsi2Phi;
    for (j = lf; j <= uf; j++) {
      tmpx = (x[j] - dBetaDerBeta * tmpphi) / dDerBeta;
      coef = exp(-0.5 * dBeta2 * tmpphi - dBeta * tmpx);
      for (l = lz; l <= uz; l++)
        for (m = 0; m < nprod; m++)
          prod_val[i][j][l][m] /= dBTTstar * coef;
    }
  };

  return NULL;
}

/****************************************************************
 *						Monte Carlo
 **
 ****************************************************************/
Err payoff_lgmsv_mc(long path_index, double evt_date, double evt_time,
                    void *func_parm,

                    double ft, double psi, double v,

                    /* Vector of results to be updated */
                    int nprod,
                    /* Result	*/
                    double *prod_val, int *stop_path) {
  GRFNPARMLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, df;
  static int k;

  static Err err = NULL;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV)func_parm;
  global = total->global;
  local = total->local;

  memset(prod_val, 0, nprod * sizeof(double));

  /*	Evaluation */
  for (k = 0; k < total->num_df; k++) {
    local->evt->df[0][k] =
        exp(total->dff[k] - total->gam[k] * ft - total->gam2[k] * psi);
  }

  local->smp.und[0].sv[PHI] = psi;
  local->smp.und[0].sv[SIGMA] = v;

  err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val, &temp);

  return err;
}

Err payoff_lgmsv2F_mc(long path_index, double evt_date, double evt_time,
                      void *func_parm,

                      double ft1, double ft2, double psi1, double psi2,
                      double psi12, double v,

                      /* Vector of results to be updated */
                      int nprod,
                      /* Result	*/
                      double *prod_val, int *stop_path) {
  GRFNPARMLGMSV2F total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpx1, tmpx2, df;
  static int k;

  static Err err = NULL;

  /*	Get the event		*/
  total = (GRFNPARMLGMSV2F)func_parm;
  global = total->global;
  local = total->local;

  memset(prod_val, 0, nprod * sizeof(double));

  /*	Evaluation */
  for (k = 0; k < total->num_df; k++) {
    local->evt->df[0][k] =
        exp(total->dff[k] - total->gam1[k] * ft1 - total->gam2[k] * ft2 -
            total->gam1_2[k] * psi1 - total->gam2_2[k] * psi2 -
            total->gam12[k] * psi12);
  }

  local->smp.und[0].sv[PHI1] = psi1;
  local->smp.und[0].sv[PHI2] = psi2;
  local->smp.und[0].sv[CROSSPHI] = psi12;
  local->smp.und[0].sv[SIGMA] = v;

  err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val, &temp);

  return err;
}