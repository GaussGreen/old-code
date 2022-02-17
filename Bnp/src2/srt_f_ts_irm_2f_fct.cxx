/* ============================================================================

   FILENAME:	SRT_F_TS_IRM_2F_FCT.C

   PURPOSE:		Retrieve all useful quantities at any date from
                                        a given Term Struct.

   TWO FACTOR EQUIVALENT OF SRT_F_TS_IRM_FCT.C

   ============================================================================
 */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts_irm.h"

/* ---------------------------------------------------------------------------
 */

/* Tau */
Err find_2f_tau(double time, TermStruct *l, double *tau1, double *tau2) {
  SrtLst *ls;

  /* Goto first element after time */
  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }

  /* If not below the end of the ts */
  if (ls != NULL) {
    *tau1 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[0].tau;
    *tau2 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[1].tau;
  } else
  /* Take last */
  {
    *tau1 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[0].tau;
    *tau2 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[1].tau;
  }
  return NULL;
} /* END Err find_2f_tau(...) */

/* ---------------------------------------------------------------------------
 */

/* Sigma: same principle */
Err find_2f_sig(double time, TermStruct *l, double *sig1, double *sig2) {
  SrtLst *ls;

  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }

  if (ls != NULL) {
    *sig1 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[0].sig;
    *sig2 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[1].sig;
  } else {
    *sig1 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[0].sig;
    *sig2 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[1].sig;
  }

  return NULL;

} /* END of Err find_2f_sig (...) */

/* Rho: idem */
Err find_2f_rho(double time, TermStruct *l, double *rho) {
  Err err = NULL;
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }

  if (ls != NULL) {
    *rho = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->rho;
  } else {
    *rho = ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->rho;
  }

  return err;
} /* END of Err find_2f_rho (...) */

/* ------------------------------------------------------------------------- */

/* Finds var/covar between time1 and time2 from term struct */
Err get_2f_sig2_rho_interp(double time1, double time2, TermStruct *l,
                           double *sig1_s, double *sig2_s, double *rho) {
  SrtLst *ls1, *ls2;
  TwoFacIrmTermStructVal *tsval1, *tsval2;
  double var1, var2, covar, sig1, sig2, corr, prev_time;
  int ind;

  /* Time2 should always be >= time1 */
  if (time2 <= time1) {
    return serror(" error in get_2f_sig2_rho_interp: time2 <= time1");
  }

  /* Top of the term struct */
  ls1 = l->head;
  /* If stationary: simply return the values */
  if (l->head == l->tail) {
    sig1 = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->exp_ts[0].sig;
    *sig1_s = sig1 * sig1;
    sig2 = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->exp_ts[1].sig;
    *sig2_s = sig2 * sig2;
    *rho = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->rho;
    return NULL;
  }

  /* If nonstationary term structure: interpolate */

  /* Positions ls1 right after time1 */
  while ((ls1 != NULL) &&
         (((TwoFacIrmTermStructVal *)ls1->element->val.pval)->time < time1)) {
    ls1 = ls1->next;
  }
  /* If time1 is after last step        , positions ls1 at tail */
  if (ls1 == NULL) {
    ls1 = l->tail;
  }

  /* Positions ls2 right after time2        , count number of steps in ind */
  ls2 = ls1;
  ind = 0;
  while ((ls2 != NULL) &&
         (((TwoFacIrmTermStructVal *)ls2->element->val.pval)->time < time2)) {
    ls2 = ls2->next;
    ind++;
  }
  /* If time2 is after last step        , positions ls2 at tail */
  if (ls2 == NULL) {
    ls2 = l->tail;
    ind--;
  }

  /* No step jump: return as if stationary */
  if (ind == 0) {
    sig1 = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->exp_ts[0].sig;
    *sig1_s = sig1 * sig1;
    sig2 = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->exp_ts[1].sig;
    *sig2_s = sig2 * sig2;
    *rho = ((TwoFacIrmTermStructVal *)ls1->element->val.pval)->rho;
    return NULL;
  } else
  /* Interpolate */
  {
    /* Sig1 and sig2 after time2 */
    sig1 = ((TwoFacIrmTermStructVal *)ls2->element->val.pval)->exp_ts[0].sig;
    sig2 = ((TwoFacIrmTermStructVal *)ls2->element->val.pval)->exp_ts[1].sig;
    corr = ((TwoFacIrmTermStructVal *)ls2->element->val.pval)->rho;

    /* Cumulative var/cov between time before time2 and previous timestep */
    tsval2 = (TwoFacIrmTermStructVal *)ls2->previous->element->val.pval;
    var1 = sig1 * sig1 * (time2 - tsval2->time);
    var2 = sig2 * sig2 * (time2 - tsval2->time);
    covar = sig1 * sig2 * corr * (time2 - tsval2->time);

    prev_time = time1;

    /* Add cumulative var/cov between time1 and timestep previous to time2 */
    while (ind > 0) {
      /* Time previous to time1 */
      tsval1 = (TwoFacIrmTermStructVal *)ls1->element->val.pval;
      sig1 = tsval1->exp_ts[0].sig;
      sig2 = tsval1->exp_ts[1].sig;
      corr = tsval1->rho;
      var1 += sig1 * sig1 * (tsval1->time - prev_time);
      var2 += sig2 * sig2 * (tsval1->time - prev_time);
      covar += sig1 * sig2 * corr * (tsval1->time - prev_time);
      prev_time = tsval1->time;
      ls1 = ls1->next;
      ind--;
    }

    *sig1_s = var1 / (time2 - time1);
    *sig2_s = var2 / (time2 - time1);
    *rho = covar / sqrt(var1 * var2);
  }
  return NULL;
} /* Err get_2f_sig2_rho_interp(...) */

/* ------------------------------------------------------------------------------
 */

/* Get Beta in a CHEY_BETA2F model */
Err find_tf_beta(double time, TermStruct *l, double *beta1, double *beta2) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }

  if (ls != NULL) {
    *beta1 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[0].beta;
    *beta2 = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->exp_ts[1].beta;
  } else {
    *beta1 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[0].beta;
    *beta2 =
        ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->exp_ts[1].beta;
  }

  return NULL;

} /* END of double find_tf_beta(...) */

/* ------------------------------------------------------------------------------
 */

/* Get Omega in a MIXED_BETA2F model */
double find_tf_omega(double time, TermStruct *l) {
  double omega;
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }

  if (ls != NULL) {
    omega = ((TwoFacIrmTermStructVal *)ls->element->val.pval)->omega;
  } else {
    omega = ((TwoFacIrmTermStructVal *)l->tail->element->val.pval)->omega;
  }

  return omega;

} /* END of double find_tf_beta2(...) */

/* ------------------------------------------------------------------------------
 */

/* Interpolate F */
Err get_2f_F_funcs(double time, TermStruct *l, SrtTFTSVec *F) {
  SrtLst *ls;
  TwoFacIrmTermStructVal *tsval, *tsval_p;
  Err err = NULL;
  double lambda, temp0, dt;
  long i;

  /* Poisition on step after time */
  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }
  if (ls == NULL) {
    ls = l->tail;
  }
  tsval = (TwoFacIrmTermStructVal *)ls->element->val.pval;

  /* Compute dt */
  if (ls != l->head) {
    tsval_p = (TwoFacIrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  dt = time;

  /* Implement chain rule */
  for (i = 0; i <= 1; i++) {
    lambda = 1.0 / tsval->exp_ts[i].tau;
    temp0 = exp(-lambda * dt);
    (*F)[i] = tsval->exp_ts[i].F * temp0;
  }

  return err;
} /* END of Err get_2f_F_funcs(...) */

/* ------------------------------------------------------------------------------
 */

/* Interpolate Psi */
Err get_2f_Psi_funcs(double time, TermStruct *l, SrtTFTSVec *Psi) {
  SrtLst *ls;
  TwoFacIrmTermStructVal *tsval, *tsval_p;
  double temp0, temp1, lambda, dt;
  int i;
  Err err = NULL;

  /* Position on step after time */
  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }
  if (ls == NULL) {
    ls = l->tail;
  }
  tsval = (TwoFacIrmTermStructVal *)ls->element->val.pval;

  /* Compute dt */
  if (ls != l->head) {
    tsval_p = (TwoFacIrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  dt = time;

  /* Implement chain rule */
  for (i = 0; i <= 1; i++) {
    lambda = 1.0 / tsval->exp_ts[i].tau;
    temp0 = exp(-lambda * dt);
    if (fabs(lambda) > EPS) {
      temp1 = (1.0 - temp0) / lambda;
    } else {
      temp1 = dt;
    }
    (*Psi)[i] = tsval->exp_ts[i].Psi + tsval->exp_ts[i].F * temp1;
  }
  return err;

} /* END of Err get_2f_Psi_funcs(...) */

/* -------------------------------------------------------------------------- */
/* Interpolate G */
Err get_2f_G_funcs(double time, TermStruct *l, SrtTFTSMat *G) {
  SrtLst *ls;
  double temp0[2], temp2[2], tempc1, lambda[2], var[2], covar, dt;
  int i;
  TwoFacIrmTermStructVal *tsval, *tsval_p;

  /* Position on step after time */
  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }
  if (ls == NULL) {
    ls = l->tail;
  }
  tsval = (TwoFacIrmTermStructVal *)ls->element->val.pval;

  /* Compute dt */
  if (ls != l->head) {
    tsval_p = (TwoFacIrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  dt = time;

  for (i = 0; i <= 1; i++) {
    var[i] = pow(tsval->exp_ts[i].sig, 2.0);
    lambda[i] = 1.0 / tsval->exp_ts[i].tau;
    temp0[i] = exp(-lambda[i] * dt);
    if (fabs(lambda[i]) > EPS) {
      temp2[i] = (exp(2 * lambda[i] * dt) - 1.0) / (2 * lambda[i]);
    } else {
      temp2[i] = dt;
    }
  }

  covar = tsval->rho * tsval->exp_ts[0].sig * tsval->exp_ts[1].sig;
  if (fabs(lambda[0] + lambda[1]) > EPS) {
    tempc1 =
        (exp((lambda[0] + lambda[1]) * dt) - 1.0) / (lambda[0] + lambda[1]);
  } else {
    tempc1 = dt;
  }

  for (i = 0; i <= 1; i++) {
    (*G)[i][i] =
        tsval->reb_ts[i][i].G + temp2[i] * var[i] / pow(tsval->exp_ts[i].F, 2);
  }
  (*G)[0][1] = (*G)[1][0] = tsval->reb_ts[0][1].G + tempc1 * covar /
                                                        tsval->exp_ts[0].F /
                                                        tsval->exp_ts[1].F;

  return NULL;

} /* END of get_2f_G_funcs(...) */

/* ------------------------------------------------------------------------- */
/* Interpolate H */
Err get_2f_H_funcs(double time, TermStruct *l, SrtTFTSMat *H) {
  SrtLst *ls;
  double temp0[2], temp1[2], tempc2[2][2], lambda[2], var[2], covar, dt;
  int i, j;
  TwoFacIrmTermStructVal *tsval, *tsval_p;

  /* Position on step after time */
  ls = l->head;
  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }
  if (ls == NULL) {
    ls = l->tail;
  }
  tsval = (TwoFacIrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (TwoFacIrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  dt = time;

  for (i = 0; i <= 1; i++) {
    var[i] = pow(tsval->exp_ts[i].sig, 2.0);
    lambda[i] = 1.0 / tsval->exp_ts[i].tau;
    temp0[i] = exp(-lambda[i] * dt);
    if (fabs(lambda[i]) > EPS) {
      temp1[i] = (1.0 - temp0[i]) / lambda[i];
    } else {
      temp1[i] = dt;
    }
  }

  covar = tsval->rho * tsval->exp_ts[0].sig * tsval->exp_ts[1].sig;
  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 1; j++) {
      if (fabs(lambda[i] + lambda[j]) > EPS) {
        tempc2[i][j] = (exp(-(lambda[i] + lambda[j]) * dt) - 1.0) /
                       (lambda[i] + lambda[j]);
      } else {
        tempc2[i][j] = -dt;
      }
    }
  }

  /* Implement chain rule */
  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 1; j++) {
      (*H)[i][j] =
          temp0[i] * tsval->reb_ts[i][j].H +
          temp0[i] * tsval->exp_ts[i].F * tsval->exp_ts[j].F *
              tsval->reb_ts[i][j].G * temp1[j] +
          (i == j ? var[i] : covar) * (temp1[i] + tempc2[i][j]) / lambda[j];
    }
  }
  return NULL;

} /* END of get_2f_H_funcs(...) */

/* ------------------------------------------------------------------------- */
/* Interpolate J */
/* NOT CHECKED */
Err get_2f_J_funcs(double time, TermStruct *l, SrtTFTSVec *J) {
  SrtLst *ls;
  double temp, temp2;
  TwoFacIrmTermStructVal *tsval, *tsval_p;
  Err err = NULL;

  ls = l->head;

  while ((ls != NULL) &&
         (((TwoFacIrmTermStructVal *)ls->element->val.pval)->time < time)) {
    ls = ls->next;
  }
  if (ls == NULL) {
    ls = l->tail;
  }

  tsval = (TwoFacIrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (TwoFacIrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  temp = exp(time / tsval->exp_ts[0].tau);
  if (fabs(1.0 / tsval->exp_ts[0].tau) > EPS) {
    temp2 = (temp - 1.0) * tsval->exp_ts[0].tau;
  } else {
    temp2 = time;
  }

  *J[0] =
      tsval->exp_ts[0].J + tsval->exp_ts[0].sig * temp2 / tsval->exp_ts[0].F;

  temp = exp(time / tsval->exp_ts[1].tau);
  if (fabs(1.0 / tsval->exp_ts[1].tau) > EPS) {
    temp2 = (temp - 1.0) * tsval->exp_ts[1].tau;
  } else {
    temp2 = time;
  }

  *J[1] =
      tsval->exp_ts[1].J + tsval->exp_ts[1].sig * temp2 / tsval->exp_ts[1].F;

  return err;

} /* END of get_2f_J_funcs(...) */

/* ----------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------- */

/* Computes the forward LGM cumulative volatility between two dates t1 and t2 ,
   defined by:
                CumVol(t        , T1        , T2) = [Psi(T2)-Psi(T1)]^2 * G(t)
   ------------------------------------------------------------------------- */
Err make_2f_lgm_cum_vols(void *tts, double fixing_time, double start_time,
                         double end_time, double *Cum_Vol3, double *Cum_Vol4,
                         double *Var) {
  SrtTFTSVec t1_Psi, t2_Psi;
  SrtTFTSMat t1_G;
  double Sig1, Sig2, corr;
  Err err;
  TermStruct *ts;

  ts = tts;

  err = get_2f_Psi_funcs(start_time, ts, &t1_Psi);
  err = get_2f_Psi_funcs(end_time, ts, &t2_Psi);
  err = get_2f_G_funcs(fixing_time, ts, &t1_G);

  corr = t1_G[1][0] / sqrt(t1_G[0][0] * t1_G[1][1]);

  if (corr < -1.00) {
    corr = -1.00;
  } else if (corr > 1.00) {
    corr = 1.00;
  }

  Sig1 = (t2_Psi[0] - t1_Psi[0]) * sqrt(t1_G[0][0]);
  Sig2 = (t2_Psi[1] - t1_Psi[1]) * sqrt(t1_G[1][1]);

  *Cum_Vol3 = fabs(sqrt((1.00 - corr) / 2.00) * (Sig2 - Sig1));
  *Cum_Vol4 = fabs(sqrt((1.00 + corr) / 2.00) * (Sig1 + Sig2));
  *Var = (*Cum_Vol3) * (*Cum_Vol3) + (*Cum_Vol4) * (*Cum_Vol4);

  return err;

} /* END of Err make_2f_lgm_cum_vols(...) */

/* ------------------------------------------------------------------------------
 */
