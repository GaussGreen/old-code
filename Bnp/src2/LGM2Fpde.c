/* ==========================================================================
   FILE_NAME:	LGM2Fpde.c

   PURPOSE:		PDE implementation of the 2 Factor LGM model.
                                Discretisation is LOD.

   DATE:		02/28/01

   AUTHOR:		L.C.
   ========================================================================== */

#define NSTD_LGM 7.0
#define NB_CHECK 10

#include "LGM2Fpde.h"
#include "Fx3FUtils.h"
#include "LGMSVUtil.h"
#include "math.h"

static double H_func(double lam, double t1, double t2) {
  return exp(lam * t2) - exp(lam * t1);
}

static void disc_normal_center(double *x, long nstp, double xmin, double stdn,
                               long *index, double *meshx);
static void disc_sqrt_center(double *x, long nstpx, double *time, long nstpt,
                             int nb_check, double *varmax, long *index,
                             double *coefx);

/*	Function to evaluate the expecations
        of the variables r1_dim1 and r3
        This function suppose that all the vol dates
        are included in the time discretisation			*/

static void LGM2FExpectations(int nstept, double *time, double lam1,
                              double *sig_time, double *sig1, int nb_sig,
                              double alpha, double gamma, double rho,
                              double *fwd1, double *fwd1_mid, double *fwd3,
                              double *fwd3_mid, double *var1, double *var2,
                              double *phi1, double *phi2, double *phi12) {
  double lam2, vol2, lam21, lam22, lam12;
  double t1, t2, ta, tb;
  double coef11, coef12, coef13, coef21, coef22, coef23, fact1, fact2, rhoalpha;

  double H1, H2, H21, H22, H12;
  double st1;

  int i, j, nb_sig_minus1;

  lam2 = lam1 + gamma;
  lam21 = 2.0 * lam1;
  lam22 = 2.0 * lam2;
  lam12 = (lam1 + lam2);
  rhoalpha = rho * alpha;
  nb_sig_minus1 = nb_sig - 1;

  fact1 = 1.0 / alpha / sqrt(1.0 - rho * rho);
  fact2 = rho / sqrt(1.0 - rho * rho);

  /* constant for expectation r1_dim1 */
  coef11 = (lam2 + rho * alpha * lam1) / (lam1 * lam1 * lam2);
  coef12 = -1.0 / (2.0 * lam1 * lam1);
  coef13 = -rhoalpha / (lam2 * lam12);

  /* constant for expectation r3 */
  coef21 = alpha * (rho * lam2 + alpha * lam1) / (lam1 * lam2 * lam2);
  coef22 = -alpha * alpha / (2.0 * lam2 * lam2);
  coef23 = -rhoalpha / (lam1 * lam12);

  /* initialisation */
  H1 = H2 = H21 = H22 = H12 = 0.0;
  t1 = 0.0;
  j = 0;
  vol2 = sig1[0] * sig1[0];

  phi1[0] = 0.0;
  phi2[0] = 0.0;
  phi12[0] = 0.0;

  fwd1[0] = 0;
  fwd3[0] = 0;

  for (i = 1; i < nstept; i++) {
    t2 = 0.5 * (t1 + time[i]);

    /* first from t1 to tmid */

    ta = t1;
    tb = sig_time[j];

    st1 = 0.0;

    while (tb < t2 && j < nb_sig_minus1) {
      H1 += vol2 * H_func(lam1, ta, tb);
      H2 += vol2 * H_func(lam2, ta, tb);
      H21 += vol2 * H_func(lam21, ta, tb);
      H22 += vol2 * H_func(lam22, ta, tb);
      H12 += vol2 * H_func(lam12, ta, tb);
      st1 += vol2 * (tb - ta);

      j++;
      vol2 = sig1[j] * sig1[j];
      ta = tb;
      tb = sig_time[j];
    }

    H1 += vol2 * H_func(lam1, ta, t2);
    H2 += vol2 * H_func(lam2, ta, t2);
    H21 += vol2 * H_func(lam21, ta, t2);
    H22 += vol2 * H_func(lam22, ta, t2);
    H12 += vol2 * H_func(lam12, ta, t2);
    st1 += vol2 * (t2 - ta);

    var1[i - 1] = st1;
    var2[i - 1] = st1;

    phi1[i] = exp(-lam21 * t2) * H21;
    phi2[i] = exp(-lam22 * t2) * H22;
    phi12[i] = exp(-lam12 * t2) * H12;

    fwd1_mid[i - 1] =
        coef11 * exp(-lam1 * t2) * H1 + coef12 * phi1[i] + coef13 * phi12[i];
    fwd3_mid[i - 1] = fact1 * (coef21 * exp(-lam2 * t2) * H2 +
                               coef22 * phi2[i] + coef23 * phi12[i]) -
                      fact2 * fwd1_mid[i - 1];

    phi1[i] /= lam21;
    phi2[i] *= alpha * alpha / lam22;
    phi12[i] *= alpha / lam12;

    t1 = t2;

    t2 = time[i];

    /* first from tmid to t2 */

    ta = t1;
    tb = sig_time[j];

    while (tb < t2 && j < nb_sig_minus1) {
      H1 += vol2 * H_func(lam1, ta, tb);
      H2 += vol2 * H_func(lam2, ta, tb);
      H21 += vol2 * H_func(lam21, ta, tb);
      H22 += vol2 * H_func(lam22, ta, tb);
      H12 += vol2 * H_func(lam12, ta, tb);
      st1 += vol2 * (tb - ta);

      j++;
      vol2 = sig1[j] * sig1[j];
      ta = tb;
      tb = sig_time[j];
    }

    H1 += vol2 * H_func(lam1, ta, t2);
    H2 += vol2 * H_func(lam2, ta, t2);
    H21 += vol2 * H_func(lam21, ta, t2);
    H22 += vol2 * H_func(lam22, ta, t2);
    H12 += vol2 * H_func(lam12, ta, t2);
    st1 += vol2 * (t2 - ta);

    var1[i - 1] = st1;
    var2[i - 1] = st1;

    phi1[i] = exp(-lam21 * t2) * H21;
    phi2[i] = exp(-lam22 * t2) * H22;
    phi12[i] = exp(-lam12 * t2) * H12;

    fwd1[i] =
        coef11 * exp(-lam1 * t2) * H1 + coef12 * phi1[i] + coef13 * phi12[i];
    fwd3[i] = fact1 * (coef21 * exp(-lam2 * t2) * H2 + coef22 * phi2[i] +
                       coef23 * phi12[i]) -
              fact2 * fwd1[i];

    phi1[i] /= lam21;
    phi2[i] *= alpha * alpha / lam22;
    phi12[i] *= alpha / lam12;

    t1 = t2;
  }
}

static Err LGM2FExpectations2(int nstept, double *time, double *sigma,
                              double *sigma_time, int nb_sigma, double *lambda,
                              double *lambda_time, int nb_lambda, double alpha,
                              double gamma, double rho, double *fwd1,
                              double *fwd3, double *var, double *quad_var,
                              double *var1_max, double *var3_max,
                              double *avg_lam, double *phi1, double *phi2,
                              double *phi12) {
  double t1, t2, ta, tb, dt;
  double alpha2, gamma2, fact1, fact2;
  double st1, lamtemp, sig1, lam1;
  double var1, var2, var3, var1max, var3max;
  double rho2;

  LGMSVSolFunc ExpectR1, ExpectR2, Phi1, Phi2, Phi12;

  double vector[20];
  int i, j, nb_sig_minus1;

  double *ts_time = NULL, *sig = NULL, *lam = NULL;

  int nb_ts;

  Err err = NULL;

  memset(vector, 0, 20);

  /* First merge the two term structure */

  ts_time = calloc(nb_sigma, sizeof(double));
  rho2 = 1.0 - rho * rho;
  fact1 = 1.0 / alpha / sqrt(1.0 - rho * rho);
  fact2 = rho / sqrt(1.0 - rho * rho);

  if (!ts_time) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  memcpy(ts_time, sigma_time, nb_sigma * sizeof(double));
  nb_ts = nb_sigma;

  for (i = 0; i < nb_lambda; i++) {
    num_f_add_number(&nb_ts, &ts_time, lambda_time[i]);
  }

  num_f_sort_vector(nb_ts, ts_time);
  num_f_unique_vector(&nb_ts, ts_time);

  sig = calloc(nb_ts, sizeof(double));
  lam = calloc(nb_ts, sizeof(double));

  if (!sig || !lam) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_ts; i++) {
    sig[i] = sigma[Get_Index(ts_time[i], sigma_time, nb_sigma)];
    lam[i] = lambda[Get_Index(ts_time[i], lambda_time, nb_lambda)];
  }

  nb_sig_minus1 = nb_ts - 1;
  alpha2 = alpha * alpha;
  gamma2 = 2.0 * gamma;
  fact1 = 1.0 / alpha / sqrt(1.0 - rho * rho);
  fact2 = -rho / sqrt(1.0 - rho * rho);

  /* Initialisation */
  phi1[0] = 0.0;
  phi2[0] = 0.0;
  phi12[0] = 0.0;

  fwd1[0] = 0;
  fwd3[0] = 0;

  quad_var[0] = 0.0;

  ExpectR1.a = 0.0;
  ExpectR1.bIsft1 = 1;
  ExpectR1.b = 1.0;
  ExpectR1.bIsgt1 = 0;
  ExpectR1.pgt = &Phi1;
  ExpectR1.c = rho;
  ExpectR1.bIsht1 = 0;
  ExpectR1.pht = &Phi12;

  ExpectR2.a = 0.0;
  ExpectR2.bIsft1 = 1;
  ExpectR2.b = 1.0;
  ExpectR2.bIsgt1 = 0;
  ExpectR2.pgt = &Phi2;
  ExpectR2.c = rho;
  ExpectR2.bIsht1 = 0;
  ExpectR2.pht = &Phi12;

  Phi1.bIsft1 = 1;
  Phi1.b = 0.0;
  Phi1.bIsgt1 = 1;
  Phi1.c = 0.0;
  Phi1.bIsht1 = 1;

  Phi2.bIsft1 = 1;
  Phi2.b = 0.0;
  Phi2.bIsgt1 = 1;
  Phi2.c = 0.0;
  Phi2.bIsht1 = 1;

  Phi12.bIsft1 = 1;
  Phi12.b = 0.0;
  Phi12.bIsgt1 = 1;
  Phi12.c = 0.0;
  Phi12.bIsht1 = 1;

  lam1 = lam[0];
  sig1 = sig[0];

  ExpectR1.dLambda = lam1;
  ExpectR2.dLambda = ExpectR1.dLambda + gamma;
  Phi1.a = sig1 * sig1;
  Phi1.dLambda = 2.0 * lam1;
  Phi2.a = alpha2 * Phi1.a;
  Phi2.dLambda = Phi1.dLambda + gamma2;
  Phi12.a = alpha * Phi1.a;
  Phi12.dLambda = Phi1.dLambda + gamma;

  t1 = 0.0;
  j = 0;

  var1 = 0.0;
  var2 = 0.0;

  var1max = 0.0;
  var3max = 0.0;

  for (i = 1; i < nstept; i++) {
    t2 = time[i];
    dt = t2 - t1;

    ta = t1;
    tb = ts_time[j];
    st1 = 0.0;
    lamtemp = 0.0;

    ExpectR1.dXt1 = fwd1[i - 1];
    ExpectR2.dXt1 = fwd3[i - 1];
    Phi1.dXt1 = phi1[i - 1];
    Phi2.dXt1 = phi2[i - 1];
    Phi12.dXt1 = phi12[i - 1];

    while (tb < t2 && j < nb_sig_minus1) {
      dt = (tb - ta);

      /* Calculation at intermedary time tb */
      ExpectR1.dXt1 = LGMSVFuncValue(ExpectR1, dt, vector, 0);
      ExpectR2.dXt1 = LGMSVFuncValue(ExpectR2, dt, vector, 0);
      Phi1.dXt1 = LGMSVFuncValue(Phi1, dt, vector, 0);
      Phi2.dXt1 = LGMSVFuncValue(Phi2, dt, vector, 0);
      Phi12.dXt1 = LGMSVFuncValue(Phi12, dt, vector, 0);

      st1 += Phi1.a * dt;
      lamtemp += lam1 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = lam[j];
      sig1 = sig[j];

      /* Update the model parameters */
      ExpectR1.dLambda = lam1;
      ExpectR2.dLambda = ExpectR1.dLambda + gamma;
      Phi1.a = sig1 * sig1;
      Phi1.dLambda = 2.0 * lam1;
      Phi2.a = alpha2 * Phi1.a;
      Phi2.dLambda = Phi1.dLambda + gamma2;
      Phi12.a = alpha * Phi1.a;
      Phi12.dLambda = Phi1.dLambda + gamma;
    }

    dt = (t2 - ta);

    fwd1[i] = LGMSVFuncValue(ExpectR1, dt, vector, 0);
    fwd3[i] = LGMSVFuncValue(ExpectR2, dt, vector, 0);
    phi1[i] = LGMSVFuncValue(Phi1, dt, vector, 0);
    phi2[i] = LGMSVFuncValue(Phi2, dt, vector, 0);
    phi12[i] = LGMSVFuncValue(Phi12, dt, vector, 0);

    st1 += Phi1.a * dt;
    lamtemp += lam1 * dt;

    var[i - 1] = st1;
    avg_lam[i - 1] = lamtemp / (t2 - t1);
    quad_var[i] = quad_var[i - 1] + st1;

    if (phi1[i] > var1max) {
      var1max = phi1[i];
    }

    var1_max[i] = var1max;

    var3 = 1.0 / rho2 *
           (phi2[i] / alpha / alpha + rho * rho * phi1[i] -
            2.0 * rho * rho / alpha * phi12[i]);

    if (var3 > var3max) {
      var3max = var3;
    }

    var3_max[i] = var3max;

    t1 = t2;
  }

  for (i = 1; i < nstept; i++) {
    fwd3[i] = fact2 * fwd1[i] + fact1 * fwd3[i];
  }

FREE_RETURN:

  if (sig)
    free(sig);
  if (lam)
    free(lam);
  if (ts_time)
    free(ts_time);

  return err;
}

static Err LGM2FExpectations3(int nstept, double *time, double *sigma,
                              double *sigma_time, int nb_sigma, double *lambda,
                              double *lambda_time, int nb_lambda, double alpha,
                              double gamma, double rho, double *fwd1,
                              double *fwd2, double *var1, double *avg_lam,
                              double *phi1, double *phi2, double *phi12,
                              double *std1_mid, double *std3_mid,
                              double *std1der_mid, double *std3der_mid) {
  double t1, t2, ta, tb, dt;
  double alpha2, gamma2;
  double st1, lamtemp, sig1, lam1;
  double rho2, std3_c1, std3_c2, std3_c3;

  LGMSVSolFunc ExpectR1_, *ExpectR1 = &ExpectR1_, ExpectR2_,
                          *ExpectR2 = &ExpectR2_, Phi1_, *Phi1 = &Phi1_, Phi2_,
                          *Phi2 = &Phi2_, Phi12_, *Phi12 = &Phi12_;

  double vector[20];
  int i, j, nb_sig_minus1;

  double *ts_time = NULL, *sig = NULL, *lam = NULL;

  int nb_ts;

  Err err = NULL;

  memset(vector, 0, 20);

  /* First merge the two term structure */

  ts_time = calloc(nb_sigma, sizeof(double));

  if (!ts_time) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  memcpy(ts_time, sigma_time, nb_sigma * sizeof(double));
  nb_ts = nb_sigma;

  for (i = 0; i < nb_lambda; i++) {
    num_f_add_number(&nb_ts, &ts_time, lambda_time[i]);
  }

  num_f_sort_vector(nb_ts, ts_time);
  num_f_unique_vector(&nb_ts, ts_time);

  sig = calloc(nb_ts, sizeof(double));
  lam = calloc(nb_ts, sizeof(double));

  if (!sig || !lam) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_ts; i++) {
    sig[i] = sigma[Get_Index(ts_time[i], sigma_time, nb_sigma)];
    lam[i] = lambda[Get_Index(ts_time[i], lambda_time, nb_lambda)];
  }

  nb_sig_minus1 = nb_ts - 1;
  alpha2 = alpha * alpha;
  gamma2 = 2.0 * gamma;
  rho2 = 1.0 - rho * rho;
  std3_c1 = 1.0 / (rho2 * alpha * alpha);
  std3_c2 = rho * rho / rho2;
  std3_c3 = -2.0 * rho * rho / alpha / rho2;

  /* Initialisation */
  phi1[0] = 0.0;
  phi2[0] = 0.0;
  phi12[0] = 0.0;

  fwd1[0] = 0;
  fwd2[0] = 0;

  ExpectR1->a = 0.0;
  ExpectR1->bIsft1 = 1;
  ExpectR1->b = 1.0;
  ExpectR1->bIsgt1 = 0;
  ExpectR1->pgt = Phi1;
  ExpectR1->c = rho;
  ExpectR1->bIsht1 = 0;
  ExpectR1->pht = Phi12;

  ExpectR2->a = 0.0;
  ExpectR2->bIsft1 = 1;
  ExpectR2->b = 1.0;
  ExpectR2->bIsgt1 = 0;
  ExpectR2->pgt = Phi2;
  ExpectR2->c = rho;
  ExpectR2->bIsht1 = 0;
  ExpectR2->pht = Phi12;

  Phi1->bIsft1 = 1;
  Phi1->b = 0.0;
  Phi1->bIsgt1 = 1;
  Phi1->c = 0.0;
  Phi1->bIsht1 = 1;

  Phi2->bIsft1 = 1;
  Phi2->b = 0.0;
  Phi2->bIsgt1 = 1;
  Phi2->c = 0.0;
  Phi2->bIsht1 = 1;

  Phi12->bIsft1 = 1;
  Phi12->b = 0.0;
  Phi12->bIsgt1 = 1;
  Phi12->c = 0.0;
  Phi12->bIsht1 = 1;

  lam1 = lam[0];
  sig1 = sig[0];

  ExpectR1->dLambda = lam1;
  ExpectR2->dLambda = ExpectR1->dLambda + gamma;
  Phi1->a = sig1 * sig1;
  Phi1->dLambda = 2.0 * lam1;
  Phi2->a = alpha2 * Phi1->a;
  Phi2->dLambda = Phi1->dLambda + gamma2;
  Phi12->a = alpha * Phi1->a;
  Phi12->dLambda = Phi1->dLambda + gamma;

  t1 = 0.0;
  j = 0;

  for (i = 1; i < nstept; i++) {
    /* First from t1 to tmid */

    t2 = (time[i] + t1) / 2.0;
    dt = t2 - t1;

    ta = t1;
    tb = ts_time[j];
    st1 = 0.0;
    lamtemp = 0.0;

    ExpectR1->dXt1 = fwd1[i - 1];
    ExpectR2->dXt1 = fwd2[i - 1];
    Phi1->dXt1 = phi1[i - 1];
    Phi2->dXt1 = phi2[i - 1];
    Phi12->dXt1 = phi12[i - 1];

    while (tb < t2 && j < nb_sig_minus1) {
      dt = (tb - ta);

      /* Calculation at intermedary time tb */
      LGMSVFuncValue2(ExpectR1, dt, vector, 0, &(ExpectR1->dXt1));
      LGMSVFuncValue2(ExpectR2, dt, vector, 0, &(ExpectR2->dXt1));
      LGMSVFuncValue2(Phi1, dt, vector, 0, &(Phi1->dXt1));
      LGMSVFuncValue2(Phi2, dt, vector, 0, &(Phi2->dXt1));
      LGMSVFuncValue2(Phi12, dt, vector, 0, &(Phi12->dXt1));

      st1 += Phi1->a * dt;
      lamtemp += lam1 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = lam[j];
      sig1 = sig[j];

      /* Update the model parameters */
      ExpectR1->dLambda = lam1;
      ExpectR2->dLambda = ExpectR1->dLambda + gamma;
      Phi1->a = sig1 * sig1;
      Phi1->dLambda = 2.0 * lam1;
      Phi2->a = alpha2 * Phi1->a;
      Phi2->dLambda = Phi1->dLambda + gamma2;
      Phi12->a = alpha * Phi1->a;
      Phi12->dLambda = Phi1->dLambda + gamma;
    }

    dt = (t2 - ta);

    LGMSVFuncValue2(ExpectR1, dt, vector, 0, &(ExpectR1->dXt1));
    LGMSVFuncValue2(ExpectR2, dt, vector, 0, &(ExpectR2->dXt1));
    LGMSVFuncValue2(Phi1, dt, vector, 0, &(Phi1->dXt1));
    LGMSVFuncValue2(Phi2, dt, vector, 0, &(Phi2->dXt1));
    LGMSVFuncValue2(Phi12, dt, vector, 0, &(Phi12->dXt1));

    st1 += Phi1->a * dt;
    lamtemp += lam1 * dt;

    std1_mid[i - 1] = sqrt(Phi1->dXt1);
    std3_mid[i - 1] = sqrt(std3_c1 * Phi2->dXt1 + std3_c2 * Phi1->dXt1 +
                           std3_c3 * Phi12->dXt1);

    std1der_mid[i - 1] =
        (Phi1->a - Phi1->dLambda * std1_mid[i - 1] * std1_mid[i - 1]) /
        (2.0 * std1_mid[i - 1]);
    std3der_mid[i - 1] = (Phi1->a - std3_c1 * Phi2->dLambda * Phi2->dXt1 -
                          std3_c2 * Phi1->dLambda * Phi1->dXt1 -
                          std3_c3 * Phi12->dLambda * Phi12->dXt1) /
                         (2.0 * std3_mid[i - 1]);

    /* Then From tmid to t2 */

    t1 = t2;
    t2 = time[i];
    dt = t2 - t1;

    ta = t1;
    tb = ts_time[j];

    while (tb < t2 && j < nb_sig_minus1) {
      dt = (tb - ta);

      /* Calculation at intermedary time tb */
      LGMSVFuncValue2(ExpectR1, dt, vector, 0, &(ExpectR1->dXt1));
      LGMSVFuncValue2(ExpectR2, dt, vector, 0, &(ExpectR2->dXt1));
      LGMSVFuncValue2(Phi1, dt, vector, 0, &(Phi1->dXt1));
      LGMSVFuncValue2(Phi2, dt, vector, 0, &(Phi2->dXt1));
      LGMSVFuncValue2(Phi12, dt, vector, 0, &(Phi12->dXt1));

      st1 += Phi1->a * dt;
      lamtemp += lam1 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = lam[j];
      sig1 = sig[j];

      /* Update the model parameters */
      ExpectR1->dLambda = lam1;
      ExpectR2->dLambda = ExpectR1->dLambda + gamma;
      Phi1->a = sig1 * sig1;
      Phi1->dLambda = 2.0 * lam1;
      Phi2->a = alpha2 * Phi1->a;
      Phi2->dLambda = Phi1->dLambda + gamma2;
      Phi12->a = alpha * Phi1->a;
      Phi12->dLambda = Phi1->dLambda + gamma;
    }

    dt = (t2 - ta);

    LGMSVFuncValue2(ExpectR1, dt, vector, 0, &(fwd1[i]));
    LGMSVFuncValue2(ExpectR2, dt, vector, 0, &(fwd2[i]));
    LGMSVFuncValue2(Phi1, dt, vector, 0, &(phi1[i]));
    LGMSVFuncValue2(Phi2, dt, vector, 0, &(phi2[i]));
    LGMSVFuncValue2(Phi12, dt, vector, 0, &(phi12[i]));

    /* */

    st1 += Phi1->a * dt;
    lamtemp += lam1 * dt;

    var1[i - 1] = st1;
    avg_lam[i - 1] = lamtemp / 2.0 / (t2 - t1);

    t1 = t2;
  }

FREE_RETURN:

  if (sig)
    free(sig);
  if (lam)
    free(lam);
  if (ts_time)
    free(ts_time);

  return err;
}

/*	Main function without Tau TS */
/*	---------------------------- */
Err lgm2f_adi(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nsteps,

    /*	Model data		*/
    double lam1, double *sig_time, double *sig1, int nb_sig, double alpha,
    double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *r1_dim2,
                       double **r2,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res) {

  Err err = NULL;

  int i, j, step, index_x, index_z;
  int nstepx, nstepz;
  double meshx, meshz, dt;
  double mu_r1i, mu_r3i;
  double lam2;
  double fact1, rho2, rhoq, rhoalpha, rhoalpha_plus1, rho2alpha;
  double r_temp;

  double std1, std3, r2_temp;

  double *r1_bar = NULL, *r3_bar = NULL, *r1_dim1 = NULL, **r1_dim2 = NULL,
         **r2 = NULL, ***values = NULL, ***values_p1 = NULL,
         ***values_temp = NULL, **mu_r1 = NULL, **mu_r3 = NULL, **var_r1 = NULL,
         **var_r3 = NULL, **r = NULL, **r_init = NULL, *fwd1 = NULL,
         *fwd1_mid = NULL, *fwd3 = NULL, *fwd3_mid = NULL, *var1 = NULL,
         *var3 = NULL, *phi1 = NULL, *phi2 = NULL, *phi12 = NULL;

  int lx, ux, lz, uz, tem;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = &pdestr;

  t1 = clock();

  /* Constant	*/
  lam2 = lam1 + gamma;

  rho2 = sqrt(1.0 - rho * rho);
  rhoalpha = rho * alpha;
  rhoalpha_plus1 = 1.0 + rhoalpha;
  rho2alpha = alpha * rho2;
  fact1 = -rho * gamma / rho2;
  rhoq = rho / rho2;

  /*	Allocations	of time vectors */
  phi1 = dvector(0, nstept - 1);
  phi2 = dvector(0, nstept - 1);
  phi12 = dvector(0, nstept - 1);
  var1 = dvector(0, nstept - 1);
  var3 = dvector(0, nstept - 1);
  fwd1 = dvector(0, nstept - 1);
  fwd3 = dvector(0, nstept - 1);
  fwd1_mid = dvector(0, nstept - 1);
  fwd3_mid = dvector(0, nstept - 1);

  if (!pde || !var1 || !var3 || !phi1 || !phi2 || !phi12 || !fwd1 || !fwd3 ||
      !fwd1_mid || !fwd3_mid) {
    err = "Memory allocation error (1) in lgm2fpde";
    goto FREE_RETURN;
  }

  /*	Calculate the expecations of r1 and r3 */
  LGM2FExpectations(nstept, time, lam1, sig_time, sig1, nb_sig, alpha, gamma,
                    rho, fwd1, fwd1_mid, fwd3, fwd3_mid, var1, var3, phi1, phi2,
                    phi12);

  /*	Calculation of the number of steps in each direction: since local
   * volatility			*/
  /*	is the same  , the mesh has to be the same  , but the number of steps
   * has to be adjusted	*/

  std1 = sqrt(phi1[nstept - 1]);
  std3 = 1.0 / rho2 *
         sqrt(phi2[nstept - 1] / alpha / alpha + rho * rho * phi1[nstept - 1] -
              2.0 * rho * rho / alpha * phi12[nstept - 1]);

  nstepx = (int)(nsteps * sqrt(std1 / std3) + 0.5);
  nstepz = (int)(nsteps * sqrt(std3 / std1) + 0.5);

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepx / 2)) * 2 + 1;
  nstepz = ((int)(nstepz / 2)) * 2 + 1;

  /*	we want at least three points in each directions
   */
  if (nstepx < 3) {
    nstepx = 3;
    nstepz = (int)(nsteps * nsteps / 3.0 + 0.5);
    nstepz = ((int)(nstepz / 2)) * 2 + 1;

    if (nstepz < 3) {
      nstepz = 3;
    }
  }

  if (nstepz < 3) {
    nstepz = 3;
    nstepx = (int)(nsteps * nsteps / 3.0 + 0.5);
    nstepx = ((int)(nstepx / 2)) * 2 + 1;

    if (nstepx < 3) {
      nstepx = 3;
    }
  }

  /*	corresponding index to the 0 value of x and z	*/
  index_x = (nstepx - 1) / 2;
  index_z = (nstepz - 1) / 2;

  meshx = 2.0 * NSTD_LGM * sqrt(phi1[nstept - 1]) / (nstepx - 1);
  meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);

  /*	Allocations of space vectors					*/
  r1_bar = dvector(0, nstepx - 1);
  r3_bar = dvector(0, nstepz - 1);
  r1_dim1 = dvector(0, nstepx - 1);
  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  mu_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  mu_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r_init = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r1_dim2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 ||
      !var_r3 || !r || !r_init || !r1_bar || !r1_dim1 || !r3_bar) {
    err = "Memory allocation error (2) in lgm2fpde";
    goto FREE_RETURN;
  }

  t2 = clock();

  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  smessage("Phase 2 -convolution  , stept: %d stepx: %d stepz: %d", nstept,
           nstepx, nstepz);

  t1 = clock();

  /*	Then discretise space in the orthogonal system r1 / r3	*/

  r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
  r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

  for (i = 1; i < nstepx; i++) {
    r1_bar[i] = r1_bar[i - 1] + meshx;
  }

  for (j = 1; j < nstepz; j++) {
    r3_bar[j] = r3_bar[j - 1] + meshz;
  }

  /*	Corresponding r1 and r2					*/
  for (i = 0; i < nstepx; i++) {
    r1_dim1[i] = r1_bar[i] + fwd1[nstept - 1];
    r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[nstept - 1];
    r_temp = rhoalpha_plus1 * r1_bar[i];

    for (j = 0; j < nstepz; j++) {
      r1_dim2[i][j] = r1_dim1[i];
      r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
      r_init[i][j] = r_temp + rho2alpha * r3_bar[j];
    }
  }

  /*	Final payoff valuation					*/
  if (!eval_evt[nstept - 1]) {
    err = "No event at last step in lgm2f_pde";
    goto FREE_RETURN;
  }

  /*	Eval payoff */
  err = payoff_func(date[nstept - 1], time[nstept - 1],
                    func_parm_tab[nstept - 1], yc, &lam1, date, 1, gamma, rho,
                    phi1[nstept - 1], phi2[nstept - 1], phi12[nstept - 1], 0,
                    nstepx - 1, 0, nstepz - 1, r1_dim1, r2, nprod, values_p1);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Initialize the CNPDE_TEMP_2D		*/

  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

  if (!pde) {
    err = "Memory allocation error (3) in lgm2fpde";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  /*	now do the backward pde					*/

  for (step = nstept - 2; step >= 0; step--) {
    dt = time[step + 1] - time[step];

    r_temp = ifr[step] + rhoalpha_plus1 * fwd1_mid[step] +
             rho2alpha * fwd3_mid[step];

    for (i = lx; i <= ux; i++) {
      mu_r1i = -lam1 * r1_bar[i] * dt;
      mu_r3i = fact1 * r1_bar[i];

      for (j = lz; j <= uz; j++) {
        mu_r1[i][j] = mu_r1i;
        mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;

        var_r1[i][j] = var1[step];
        var_r3[i][j] = var3[step];

        r[i][j] = (r_init[i][j] + r_temp) * dt;
      }
    }

    /*	convolve							*/

    num_f_pde_one_step_backward_2f_adi(
        pde, nstepx, r1_bar, nstepz, r3_bar, 0, nprod - 1, values_p1, mu_r1,
        mu_r3, var_r1, var_r3, r, values, lx, ux, lz, uz);

    /*	Eval payoff */
    if (eval_evt[step]) {

      /*	Corresponding r1 and r2					*/
      for (i = lx; i <= ux; i++) {
        r1_dim1[i] = r1_bar[i] + fwd1[step];
        r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[step];

        for (j = lz; j <= uz; j++) {
          r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
        }
      }

      err =
          payoff_func(date[step], time[step], func_parm_tab[step], yc, &lam1,
                      date, 1, gamma, rho, phi1[step], phi2[step], phi12[step],
                      lx, ux, lz, uz, r1_dim1, r2, nprod, values);
      if (err) {
        goto FREE_RETURN;
      }
    }

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;

    /* new indexes: we cut the PDE after NSTD_LGM number of standard deviations
     */
    /* but we need at least three points to do the pde
     */

    tem = (int)(NSTD_LGM * sqrt(phi1[step]) / meshx + 0.5);
    ux = min(nstepx - 1, index_x + tem);
    lx = max(0, index_x - tem);

    if (ux - lx < 2) {
      tem += 1;
      ux = min(nstepx - 1, index_x + tem);
      lx = max(0, index_x - tem);
    }

    tem = (int)(NSTD_LGM / rho2 *
                    sqrt(phi2[step] / alpha / alpha + rho * rho * phi1[step] -
                         2.0 * rho * rho / alpha * phi12[step]) /
                    meshx +
                0.5);
    uz = min(nstepz - 1, index_z + tem);
    lz = max(0, index_z - tem);

    if (uz - lz < 2) {
      tem += 1;
      uz = min(nstepz - 1, index_z + tem);
      lz = max(0, index_z - tem);
    }
  }

  /* copy the result					*/
  for (i = 0; i < nprod; i++) {
    res[i] = values_p1[index_x][index_z][i];
  }

  t2 = clock();

  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /*	allocation (1)		*/

  if (phi1)
    free_dvector(phi1, 0, nstept - 1);
  if (phi2)
    free_dvector(phi2, 0, nstept - 1);
  if (phi12)
    free_dvector(phi12, 0, nstept - 1);
  if (var1)
    free_dvector(var1, 0, nstept - 1);
  if (var3)
    free_dvector(var3, 0, nstept - 1);
  if (fwd1)
    free_dvector(fwd1, 0, nstept - 1);
  if (fwd3)
    free_dvector(fwd3, 0, nstept - 1);
  if (fwd1_mid)
    free_dvector(fwd1_mid, 0, nstept - 1);
  if (fwd3_mid)
    free_dvector(fwd3_mid, 0, nstept - 1);

  /*	allocation (2)		*/
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  if (r1_bar)
    free_dvector(r1_bar, 0, nstepx - 1);
  if (r3_bar)
    free_dvector(r3_bar, 0, nstepz - 1);
  if (r1_dim1)
    free_dvector(r1_dim1, 0, nstepx - 1);
  if (r1_dim2)
    free_dmatrix(r1_dim2, 0, nstepx - 1, 0, nstepz - 1);
  if (r2)
    free_dmatrix(r2, 0, nstepx - 1, 0, nstepz - 1);
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (mu_r1)
    free_dmatrix(mu_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (mu_r3)
    free_dmatrix(mu_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r1)
    free_dmatrix(var_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r3)
    free_dmatrix(var_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);
  if (r_init)
    free_dmatrix(r_init, 0, nstepx - 1, 0, nstepz - 1);

  return err;
}

#if 0

Err	 lgm2fTau_adi_new(	
					/*	Time data		*/
					int			nstept  ,					
					double		*time  ,
					double		*date  ,

					/*	Discretisation	*/
					int			nsteps  ,					
										
					/*	Model data		*/

					double		*pdLambda_time  ,
					double		*pdLambda  ,
					int				nb_lambda  ,
					
					double		*sig_time  ,
					double		*sig1  ,
					int			nb_sig  ,
					double		alpha  ,
					double		gamma  ,
					double		rho  ,					

					/*	Product data */
					void		**func_parm_tab  , 
					int			*eval_evt  ,
					
					/*	Market data */
					double		*ifr  ,					
					char		*yc  ,
					
					/*	Payoff function */
					Err (*payoff_func)(/* Event */
										double	evt_date  ,
										double	evt_time  ,
										void	*func_parm  ,
										
										/* Market data	*/										
										void	*yc  ,
										
										/* Model data	*/
										double	*lam  ,
										double	*ts_time  ,
										int		nb_ts  ,
										double	gamma  ,										
										double	rho  ,
										double	phi1  ,
										double	phi2  ,
										double	phi12  ,
										
										/* Gride data	*/
										int		l1  ,
										int		u1  ,
										int		l2  ,
										int		u2  ,
										double	*r1_dim2  ,
										double	**r2  ,
																
										/* Vector of results to be updated */
										int		nprod  ,
										double	***prod_val
										)  ,
					/*	Result */
					int			nprod  , 
					double		*res)
{

Err				err = NULL;
				
int				i  , j  , step  , index_x  , index_z;
int				nstepx  , nstepz;
double			meshx  , meshz  , dt;
double			mu_r1i  , mu_r3i;
//double			lam2;
double			fact1  , rho2  , rhoq  , rhoalpha  , rhoalpha_plus1  , rho2alpha;
double			r_temp;

double			std1  , std3  , r2_temp;

double			*r1_bar			= NULL  ,
				*r3_bar			= NULL  ,
				*r1_dim1		= NULL  ,
				**r1_dim2		= NULL  ,
				**r2			= NULL  ,
				***values		= NULL  ,
				***values_p1	= NULL  ,
				***values_temp	= NULL  ,
				**mu_r1			= NULL  ,
				**mu_r3			= NULL  ,
				**var_r1		= NULL  ,
				**var_r3		= NULL  ,
				**r				= NULL  ,
				**r_init		= NULL  ,
				*fwd1			= NULL  ,
				*fwd1_mid		= NULL  ,
				*fwd3			= NULL  ,
				*fwd3_mid		= NULL  ,
				*var1			= NULL  ,
				*var3			= NULL  ,
				*phi1			= NULL  ,
				*phi2			= NULL  ,
				*phi12			= NULL;

int				lx  , ux  , lz  , uz  , tem;

clock_t			t1  , t2;

CNPDE_TEMP_2D_ADI	pdestr  , *pde = &pdestr;

	t1 = clock();

	/* Constant	*/
	//lam2 = lam1 + gamma;
	
	rho2 = sqrt(1.0 - rho * rho);	
	rhoalpha = rho * alpha;	
	rhoalpha_plus1 = 1.0 + rhoalpha;
	rho2alpha = alpha * rho2;
	fact1 = - rho * gamma / rho2;
	rhoq = rho / rho2;
								
	/*	Allocations	of time vectors						*/
	phi1 = dvector(0  , nstept-1);
	phi2 = dvector(0  , nstept-1);
	phi12 = dvector(0  , nstept-1);
	var1 = dvector(0  , nstept-1);
	var3 = dvector(0  , nstept-1);
	fwd1 = dvector(0  , nstept-1);
	fwd3 = dvector(0  , nstept-1);
	fwd1_mid = dvector(0  , nstept-1);
	fwd3_mid = dvector(0  , nstept-1);

	if (!pde || !var1 || !var3 || !phi1 || !phi2 || !phi12 || !fwd1 || !fwd3 || !fwd1_mid || !fwd3_mid)
	{
		err = "Memory allocation error (1) in lgm2fpde";
		goto FREE_RETURN;
	}
		
	/*	Calculate the expecations of r1 and r3 */
	LGM2FExpectations2(	nstept  ,
						time  ,
						
						sig_time  ,
						sig1  ,
						nb_sig  ,
						
						pdLambda  ,
						pdLambda_time  ,
						nb_lambda  ,

						alpha  ,
						gamma  ,
						rho  ,
						fwd1  ,
						//fwd1_mid  ,
						fwd3  ,
						//fwd3_mid  ,
						var1  ,
						var3  ,
						phi1  ,
						phi2  ,
						phi12
						);

	/*	Calculation of the number of steps in each direction: since local volatility			*/
	/*	is the same  , the mesh has to be the same  , but the number of steps has to be adjusted	*/

	std1 = sqrt(phi1[nstept-1]);
	std3 = 1.0 / rho2 * sqrt( phi2[nstept-1] / alpha / alpha 
							 + rho * rho * phi1[nstept-1] 
							 - 2.0 * rho * rho / alpha * phi12[nstept-1]);

	nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
	nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
		
	/*	nstep has to be a odd nuber			*/
	nstepx = ((int) (nstepx / 2)) * 2 + 1;
	nstepz = ((int) (nstepz / 2)) * 2 + 1;

	/*	we want at least three points in each directions										*/
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}

	/*	corresponding index to the 0 value of x and z	*/
	index_x = (nstepx - 1) / 2;
	index_z = (nstepz - 1) / 2;
	
	meshx = 2.0 * NSTD_LGM * sqrt(phi1[nstept-1]) / (nstepx - 1);
	meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);
		
	/*	Allocations of space vectors					*/
	r1_bar = dvector(0  , nstepx-1);
	r3_bar = dvector(0  , nstepz-1);
	r1_dim1 = dvector(0  , nstepx-1);
	values = f3tensor(0  , nstepx-1  , 0  , nstepz-1  , 0  , nprod-1);
	values_p1 = f3tensor(0  , nstepx-1  , 0  , nstepz-1  , 0  , nprod-1);
	mu_r1 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	mu_r3 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	var_r1 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	var_r3 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	r = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	r_init = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	r1_dim2 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);
	r2 = dmatrix(0  , nstepx-1  , 0  , nstepz-1);

	if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r || !r_init ||
		!r1_bar || !r1_dim1 || !r3_bar)
	{
		err = "Memory allocation error (2) in lgm2fpde";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing  , time in sec: %.2f"  , (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution  , stept: %d stepx: %d stepz: %d"  , nstept  , nstepx  , nstepz);

	t1 = clock();


	/*	Then discretise space in the orthogonal system r1 / r3	*/
						
	r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
	r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

	for (i=1; i<nstepx; i++)
	{
		r1_bar[i] = r1_bar[i-1] + meshx;
	}

	for (j=1; j<nstepz; j++)
	{
		r3_bar[j] = r3_bar[j-1] + meshz;
	}

	/*	Corresponding r1 and r2					*/
	for (i=0; i<nstepx; i++)
	{
		r1_dim1[i] = r1_bar[i] + fwd1[nstept-1];		
		r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[nstept-1];
		r_temp = rhoalpha_plus1 * r1_bar[i];

		for (j=0; j<nstepz; j++)
		{
			r1_dim2[i][j] = r1_dim1[i];
			r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
			r_init[i][j] = r_temp + rho2alpha * r3_bar[j];
		}
	}


	/*	Final payoff valuation					*/
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in lgm2f_pde";
		goto FREE_RETURN;
	}

	/*	Eval payoff */
	err = payoff_func (	date[nstept-1]  ,
						time[nstept-1]  ,
						func_parm_tab[nstept-1]  ,
						yc  ,
						&lam1  ,
						date  ,
						1  ,
						gamma  ,
						rho  ,
						phi1[nstept-1]  ,
						phi2[nstept-1]  ,
						phi12[nstept-1]  ,
						0  ,
						nstepx-1  ,
						0  ,
						nstepz-1  ,
						r1_dim1  ,						
						r2  ,
						nprod  ,
						values_p1
						);

	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Initialize the CNPDE_TEMP_2D		*/

	num_f_pde_init_2d_adi(	pde  ,
							nstepx  ,
							nstepz  ,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (3) in lgm2fpde";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;

	/*	now do the backward pde					*/

	for (step=nstept-2; step>=0; step--)
	{
		dt = time[step+1] - time[step];
						
		r_temp = ifr[step] + rhoalpha_plus1 * fwd1_mid[step] + rho2alpha * fwd3_mid[step];

		for (i=lx; i<=ux; i++)
		{							
			mu_r1i = -lam1 * r1_bar[i] * dt;
			mu_r3i = fact1 * r1_bar[i];

			for (j=lz; j<=uz; j++)
			{
				mu_r1[i][j] = mu_r1i;
				mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;
				
				var_r1[i][j] = var1[step];
				var_r3[i][j] = var3[step];
												
				r[i][j] = (r_init[i][j] + r_temp) * dt;
			}
		}

		/*	convolve							*/

		num_f_pde_one_step_backward_2f_adi(	pde  ,
											nstepx  ,
											r1_bar  ,
											nstepz  ,
											r3_bar  ,
											0  ,
											nprod-1  ,
											values_p1  ,
											mu_r1  ,
											mu_r3  ,
											var_r1  ,
											var_r3  ,
											r  ,
											values  ,
											lx  ,
											ux  ,
											lz  ,
											uz);
		
		/*	Eval payoff */
		if (eval_evt[step])
		{

			/*	Corresponding r1 and r2					*/
			for (i=lx; i<=ux; i++)
			{
				r1_dim1[i] = r1_bar[i] + fwd1[step];		
				r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[step];			

				for (j=lz; j<=uz; j++)
				{					
					r2[i][j] = r2_temp + rho2alpha * r3_bar[j];				
				}
			}

			err = payoff_func (	date[step]  ,
								time[step]  ,
								func_parm_tab[step]  ,						
								yc  ,
								&lam1  ,
								date  ,
								1  ,
								gamma  ,
								rho  ,
								phi1[step]  ,
								phi2[step]  ,
								phi12[step]  ,
								lx  ,
								ux  ,
								lz  ,
								uz  ,
								r1_dim1  ,
								r2  ,
								nprod  ,
								values
								);
			if (err)
			{
				goto FREE_RETURN;
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;

		/* new indexes: we cut the PDE after NSTD_LGM number of standard deviations */
		/* but we need at least three points to do the pde						*/

		tem = (int) (NSTD_LGM * sqrt(phi1[step]) / meshx + 0.5);
		ux = min(nstepx-1  , index_x + tem);
		lx = max(0  , index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1  , index_x + tem);
			lx = max(0  , index_x - tem);
		}

		tem = (int) (NSTD_LGM / rho2 * sqrt(phi2[step] / alpha / alpha
									  + rho * rho * phi1[step]
									  - 2.0 * rho * rho / alpha * phi12[step]) / meshx + 0.5);			
		uz = min(nstepz-1  , index_z + tem);
		lz = max(0  , index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1  , index_z + tem);
			lz = max(0  , index_z - tem);
		}		
	}

	/* copy the result					*/
	for (i=0; i<nprod; i++)
	{
		res[i] = values_p1[index_x][index_z][i];
	}

	t2 = clock();

	smessage ("Phase 2 -convolution  , time in sec: %.2f"  , (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

	/*	allocation (1)		*/

	if (phi1) free_dvector(phi1  , 0  , nstept-1);
	if (phi2) free_dvector(phi2  , 0  , nstept-1);
	if (phi12) free_dvector(phi12  , 0  , nstept-1);
	if (var1) free_dvector(var1  , 0  , nstept-1);
	if (var3) free_dvector(var3  , 0  , nstept-1);
	if (fwd1) free_dvector(fwd1  , 0  , nstept-1);
	if (fwd3) free_dvector(fwd3  , 0  , nstept-1);
	if (fwd1_mid) free_dvector(fwd1_mid  , 0  , nstept-1);
	if (fwd3_mid) free_dvector(fwd3_mid  , 0  , nstept-1);


	/*	allocation (2)		*/
	if (pde) num_f_pde_free_2d_adi(pde  , nstepx  , nstepz  , nprod);

	if (r1_bar) free_dvector(r1_bar  , 0  , nstepx-1);
	if (r3_bar) free_dvector(r3_bar  , 0  , nstepz-1);
	if (r1_dim1) free_dvector(r1_dim1  , 0  , nstepx-1);
	if (r1_dim2) free_dmatrix(r1_dim2  , 0  , nstepx-1  , 0  , nstepz-1);
	if (r2) free_dmatrix(r2  , 0  , nstepx-1  , 0  , nstepz-1);
	if (values) free_f3tensor(values  , 0  , nstepx-1  , 0  , nstepz-1  , 0  , nprod-1);
	if (values_p1) free_f3tensor(values_p1  , 0  , nstepx-1  , 0  , nstepz-1  , 0  , nprod-1);
	if (mu_r1) free_dmatrix(mu_r1  , 0  , nstepx-1  , 0  , nstepz-1);
	if (mu_r3) free_dmatrix(mu_r3  , 0  , nstepx-1  , 0  , nstepz-1);
	if (var_r1) free_dmatrix(var_r1  , 0  , nstepx-1  , 0  , nstepz-1);
	if (var_r3) free_dmatrix(var_r3  , 0  , nstepx-1  , 0  , nstepz-1);
	if (r) free_dmatrix(r  , 0  , nstepx-1  , 0  , nstepz-1);
	if (r_init) free_dmatrix(r_init  , 0  , nstepx-1  , 0  , nstepz-1);
			
	return err;
}

#endif // if 0

/*	Main function with Tau TS */
/*	------------------------- */
Err lgm2fTau_adi(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nsteps,
    int disc_method, /* 0: linear distretisation  , 1: normal discretisation */

    /*	Model data		*/
    double *sig, double *sig_time, int nb_sig, double *lam, double *lam_time,
    int nb_lam, double alpha, double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *r1_dim2,
                       double **r2,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res) {

  Err err = NULL;

  int i, j, step, index_x, index_z;
  int nstepx, nstepz;
  double meshx, meshz, dt;
  double mu_r1i, mu_r3i;
  double lam1, lam2;
  double fact1, rho2, rhoq, rhoalpha, rhoalpha_plus1, rho2alpha;
  double r_temp;

  double std1, std3, r2_temp;

  double *r1_bar = NULL, *r3_bar = NULL, *r1_dim1 = NULL, **r1_dim2 = NULL,
         **r2 = NULL, ***values = NULL, ***values_p1 = NULL,
         ***values_temp = NULL, **mu_r1 = NULL, **mu_r3 = NULL, **var_r1 = NULL,
         **var_r3 = NULL, **r = NULL, **r_init = NULL, *fwd1 = NULL,
         *fwd3 = NULL, *var1 = NULL, *quad_var = NULL, *avg_lam = NULL,
         *phi1 = NULL, *phi2 = NULL, *phi12 = NULL, *var1_max = NULL,
         *var3_max = NULL;

  int lx, ux, lz, uz, tem;
  double coefx, coefz;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = &pdestr;

  t1 = clock();

  /* Constant	*/
  rho2 = sqrt(1.0 - rho * rho);
  rhoalpha = rho * alpha;
  rhoalpha_plus1 = 1.0 + rhoalpha;
  rho2alpha = alpha * rho2;
  fact1 = -rho * gamma / rho2;
  rhoq = rho / rho2;

  /*	Allocations	of time vectors */
  phi1 = dvector(0, nstept - 1);
  phi2 = dvector(0, nstept - 1);
  phi12 = dvector(0, nstept - 1);
  var1 = dvector(0, nstept - 1);
  quad_var = dvector(0, nstept - 1);
  avg_lam = dvector(0, nstept - 1);
  fwd1 = dvector(0, nstept - 1);
  fwd3 = dvector(0, nstept - 1);
  var1_max = dvector(0, nstept - 1);
  var3_max = dvector(0, nstept - 1);

  if (!pde || !var1 || !avg_lam || !phi1 || !phi2 || !phi12 || !fwd1 || !fwd3 ||
      !quad_var || !var1_max || !var3_max) {
    err = "Memory allocation error (1) in lgm2fpde";
    goto FREE_RETURN;
  }

  /*	Calculate the expecations of r1 and r3 */
  err =
      LGM2FExpectations2(nstept, time, sig, sig_time, nb_sig, lam, lam_time,
                         nb_lam, alpha, gamma, rho, fwd1, fwd3, var1, quad_var,
                         var1_max, var3_max, avg_lam, phi1, phi2, phi12);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calculation of the number of steps in each direction: since local
   * volatility			*/
  /*	is the same  , the mesh has to be the same  , but the number of steps
   * has to be adjusted	*/

  /*
  std1 = sqrt(phi1[nstept-1]);
  std3 = 1.0 / rho2 * sqrt( phi2[nstept-1] / alpha / alpha
                                                   + rho * rho * phi1[nstept-1]
                                                   - 2.0 * rho * rho / alpha *
  phi12[nstept-1]);
  */

  /*
  std1 = sqrt(quad_var[nstept-1]);
  std3 = std1;
  */

  std1 = sqrt(var1_max[nstept - 1]);
  std3 = sqrt(var3_max[nstept - 1]);

  if (nsteps == (((int)(nsteps / 2)) * 2 + 1)) {
    disc_method = 2;
  } else {
    disc_method = 0;
  }

  nstepx = (int)(nsteps * sqrt(std1 / std3) + 0.5);
  nstepz = (int)(nsteps * sqrt(std3 / std1) + 0.5);

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepx / 2)) * 2 + 1;
  nstepz = ((int)(nstepz / 2)) * 2 + 1;

  /*	we want at least three points in each directions
   */
  if (nstepx < 3) {
    nstepx = 3;
    nstepz = (int)(nsteps * nsteps / 3.0 + 0.5);
    nstepz = ((int)(nstepz / 2)) * 2 + 1;

    if (nstepz < 3) {
      nstepz = 3;
    }
  }

  if (nstepz < 3) {
    nstepz = 3;
    nstepx = (int)(nsteps * nsteps / 3.0 + 0.5);
    nstepx = ((int)(nstepx / 2)) * 2 + 1;

    if (nstepx < 3) {
      nstepx = 3;
    }
  }

  /*	corresponding index to the 0 value of x and z	*/
  index_x = (nstepx - 1) / 2;
  index_z = (nstepz - 1) / 2;

  meshx = 2.0 * NSTD_LGM * std1 / (nstepx - 1);
  meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);

  /*	Allocations of space vectors					*/
  r1_bar = dvector(0, nstepx - 1);
  r3_bar = dvector(0, nstepz - 1);
  r1_dim1 = dvector(0, nstepx - 1);
  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  mu_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  mu_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r_init = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r1_dim2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 ||
      !var_r3 || !r || !r_init || !r1_bar || !r1_dim1 || !r3_bar) {
    err = "Memory allocation error (2) in lgm2fpde";
    goto FREE_RETURN;
  }

  t2 = clock();

  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  smessage("Phase 2 -convolution  , stept: %d stepx: %d stepz: %d", nstept,
           nstepx, nstepz);

  t1 = clock();

  /*	Then discretise space in the orthogonal system r1 / r3	*/

  r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
  r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

  if (disc_method == 1) {
    disc_normal_center(r1_bar, nstepx, r1_bar[0], std1, &index_x, &meshx);
    disc_normal_center(r3_bar, nstepz, r3_bar[0], std3, &index_z, &meshz);
  } else if (disc_method == 0) {
    for (i = 1; i < nstepx; i++) {
      r1_bar[i] = r1_bar[i - 1] + meshx;
    }

    for (j = 1; j < nstepz; j++) {
      r3_bar[j] = r3_bar[j - 1] + meshz;
    }
  } else {
    disc_sqrt_center(r1_bar, nstepx, time, nstept, NB_CHECK, var1_max, &index_x,
                     &coefx);
    disc_sqrt_center(r3_bar, nstepz, time, nstept, NB_CHECK, var3_max, &index_z,
                     &coefz);
  }

  /*	Corresponding r1 and r2					*/
  for (i = 0; i < nstepx; i++) {
    r1_dim1[i] = r1_bar[i] + fwd1[nstept - 1];
    r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[nstept - 1];
    r_temp = rhoalpha_plus1 * r1_bar[i];

    for (j = 0; j < nstepz; j++) {
      r1_dim2[i][j] = r1_dim1[i];
      r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
      r_init[i][j] = r_temp + rho2alpha * r3_bar[j];
    }
  }

  /*	Final payoff valuation					*/
  if (!eval_evt[nstept - 1]) {
    err = "No event at last step in lgm2f_pde";
    goto FREE_RETURN;
  }

  /*	Eval payoff */
  err =
      payoff_func(date[nstept - 1], time[nstept - 1], func_parm_tab[nstept - 1],
                  yc, lam, lam_time, nb_lam, gamma, rho, phi1[nstept - 1],
                  phi2[nstept - 1], phi12[nstept - 1], 0, nstepx - 1, 0,
                  nstepz - 1, r1_dim1, r2, nprod, values_p1);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Initialize the CNPDE_TEMP_2D		*/

  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

  if (!pde) {
    err = "Memory allocation error (3) in lgm2fpde";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  /*	now do the backward pde					*/

  for (step = nstept - 2; step >= 0; step--) {
    dt = time[step + 1] - time[step];

    r_temp = ifr[step] + rhoalpha_plus1 * fwd1[step] + rho2alpha * fwd3[step];

    lam1 = avg_lam[step];
    lam2 = lam1 + gamma;

    for (i = lx; i <= ux; i++) {
      mu_r1i = -lam1 * r1_bar[i] * dt;
      mu_r3i = fact1 * r1_bar[i];

      for (j = lz; j <= uz; j++) {
        mu_r1[i][j] = mu_r1i;
        mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;

        var_r1[i][j] = var1[step];
        var_r3[i][j] = var1[step];

        r[i][j] = (r_init[i][j] + r_temp) * dt;
      }
    }

    /*	convolve							*/

    num_f_pde_one_step_backward_2f_adi(
        pde, nstepx, r1_bar, nstepz, r3_bar, 0, nprod - 1, values_p1, mu_r1,
        mu_r3, var_r1, var_r3, r, values, lx, ux, lz, uz);

    /*	Eval payoff */
    if (eval_evt[step]) {

      /*	Corresponding r1 and r2					*/
      for (i = lx; i <= ux; i++) {
        r1_dim1[i] = r1_bar[i] + fwd1[step];
        r2_temp = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[step];

        for (j = lz; j <= uz; j++) {
          r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
        }
      }

      err =
          payoff_func(date[step], time[step], func_parm_tab[step], yc, lam,
                      lam_time, nb_lam, gamma, rho, phi1[step], phi2[step],
                      phi12[step], lx, ux, lz, uz, r1_dim1, r2, nprod, values);

      if (err) {
        goto FREE_RETURN;
      }
    }

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;

    /* new indexes: we cut the PDE after NSTD_LGM number of standard deviations
     */
    /* but we need at least three points to do the pde
     */

    if (disc_method == 1) {
      tem = (int)((norm(NSTD_LGM * sqrt(var1_max[step]) / std1) - 0.5) / meshx +
                  0.5);
    } else if (disc_method == 0) {
      tem = (int)(NSTD_LGM * sqrt(var1_max[step]) / meshx + 0.5);
    } else {
      tem = (int)(coefx * sqrt(time[step]) / 2.0 + 0.5);
    }

    ux = min(nstepx - 1, index_x + tem);
    lx = max(0, index_x - tem);

    while (ux - lx < 2) {
      tem += 1;
      ux = min(nstepx - 1, index_x + tem);
      lx = max(0, index_x - tem);
    }

    if (disc_method == 1) {
      tem = (int)((norm(NSTD_LGM * sqrt(var3_max[step]) / std3) - 0.5) / meshz +
                  0.5);
    } else if (disc_method == 0) {
      tem = (int)(NSTD_LGM * sqrt(var3_max[step]) / meshx + 0.5);
    } else {
      tem = (int)(coefz * sqrt(time[step]) / 2.0 + 0.5);
    }

    uz = min(nstepz - 1, index_z + tem);
    lz = max(0, index_z - tem);

    while (uz - lz < 2) {
      tem += 1;
      uz = min(nstepz - 1, index_z + tem);
      lz = max(0, index_z - tem);
    }
  }

  /* copy the result					*/
  for (i = 0; i < nprod; i++) {
    res[i] = values_p1[index_x][index_z][i];
  }

  //	fclose(stream);

  t2 = clock();

  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /*	allocation (1)		*/

  if (phi1)
    free_dvector(phi1, 0, nstept - 1);
  if (phi2)
    free_dvector(phi2, 0, nstept - 1);
  if (phi12)
    free_dvector(phi12, 0, nstept - 1);
  if (var1)
    free_dvector(var1, 0, nstept - 1);
  if (avg_lam)
    free_dvector(avg_lam, 0, nstept - 1);
  if (fwd1)
    free_dvector(fwd1, 0, nstept - 1);
  if (fwd3)
    free_dvector(fwd3, 0, nstept - 1);
  if (quad_var)
    free_dvector(quad_var, 0, nstept - 1);
  if (var1_max)
    free_dvector(var1_max, 0, nstept - 1);
  if (var3_max)
    free_dvector(var3_max, 0, nstept - 1);

  /*	allocation (2)		*/
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  if (r1_bar)
    free_dvector(r1_bar, 0, nstepx - 1);
  if (r3_bar)
    free_dvector(r3_bar, 0, nstepz - 1);
  if (r1_dim1)
    free_dvector(r1_dim1, 0, nstepx - 1);
  if (r1_dim2)
    free_dmatrix(r1_dim2, 0, nstepx - 1, 0, nstepz - 1);
  if (r2)
    free_dmatrix(r2, 0, nstepx - 1, 0, nstepz - 1);
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (mu_r1)
    free_dmatrix(mu_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (mu_r3)
    free_dmatrix(mu_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r1)
    free_dmatrix(var_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r3)
    free_dmatrix(var_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);
  if (r_init)
    free_dmatrix(r_init, 0, nstepx - 1, 0, nstepz - 1);

  return err;
}

/*	Other function with Tau TS but specially adapted for neg Tau */
/*	------------------------------------------------------------ */
Err lgm2fTau2_adi(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nsteps,
    int disc_method, /* 0: linear distretisation  , 1: normal discretisation */

    /*	Model data		*/
    double *sig, double *sig_time, int nb_sig, double *lam, double *lam_time,
    int nb_lam, double alpha, double gamma, double rho,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,

                       /* Model data	*/
                       double *lam, double *ts_time, int nb_ts, double gamma,
                       double rho, double phi1, double phi2, double phi12,

                       /* Gride data	*/
                       int l1, int u1, int l2, int u2, double *r1_dim2,
                       double **r2,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res) {

  Err err = NULL;

  int i, j, step, index_x, index_z;
  int nstepx, nstepz;
  double meshx, meshz, dt;
  double mu_r1i, mu_r3i;
  double lam1, lam2;
  double fact1, rho2, rhoq, rhoalpha, rhoalpha_plus1, rho2alpha, std3_c1,
      std3_c2, std3_c3;
  double r_temp;

  double std1, std3, r2_temp, std1_next, std3_next, std1_mean, std3_mean,
      var1_temp, var3_temp, fact1_temp, rhoalpha_plus1_temp;

  double *r1_bar = NULL, *r3_bar = NULL, *r1_dim1 = NULL, **r2 = NULL,
         ***values = NULL, ***values_p1 = NULL, ***values_temp = NULL,
         **mu_r1 = NULL, **mu_r3 = NULL, **var_r1 = NULL, **var_r3 = NULL,
         **r = NULL, *fwd1 = NULL, *fwd2 = NULL, *std1_mid = NULL,
         *std3_mid = NULL, *std1der_mid = NULL, *std3der_mid = NULL,
         *var1 = NULL, *avg_lam = NULL, *phi1 = NULL, *phi2 = NULL,
         *phi12 = NULL;

  int lx, ux, lz, uz;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = &pdestr;

  t1 = clock();

  /* Constant	*/
  rho2 = sqrt(1.0 - rho * rho);
  rhoalpha = rho * alpha;
  rhoalpha_plus1 = 1.0 + rhoalpha;
  rho2alpha = alpha * rho2;
  fact1 = -rho * gamma / rho2;
  rhoq = rho / rho2;
  std3_c1 = 1.0 / (alpha * alpha * rho2 * rho2);
  std3_c2 = rho * rho / rho2 / rho2;
  std3_c3 = -2.0 * rho * rho / rho2 / rho2 / alpha;

  /*	Allocations	of time vectors */
  phi1 = dvector(0, nstept - 1);
  phi2 = dvector(0, nstept - 1);
  phi12 = dvector(0, nstept - 1);
  var1 = dvector(0, nstept - 1);
  avg_lam = dvector(0, nstept - 1);
  fwd1 = dvector(0, nstept - 1);
  fwd2 = dvector(0, nstept - 1);
  std1_mid = dvector(0, nstept - 1);
  std3_mid = dvector(0, nstept - 1);
  std1der_mid = dvector(0, nstept - 1);
  std3der_mid = dvector(0, nstept - 1);

  if (!pde || !var1 || !avg_lam || !phi1 || !phi2 || !phi12 || !fwd1 || !fwd2 ||
      !std1_mid || !std3_mid || !std1der_mid || !std3der_mid) {
    err = "Memory allocation error (1) in lgm2fpde";
    goto FREE_RETURN;
  }

  /*	Calculate the expecations of r1 and r3 */

  err = LGM2FExpectations3(nstept, time, sig, sig_time, nb_sig, lam, lam_time,
                           nb_lam, alpha, gamma, rho, fwd1, fwd2, var1, avg_lam,
                           phi1, phi2, phi12, std1_mid, std3_mid, std1der_mid,
                           std3der_mid);

  if (err) {
    goto FREE_RETURN;
  }

  nstepx = nsteps;
  nstepz = nsteps;

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepx / 2)) * 2 + 1;
  nstepz = ((int)(nstepz / 2)) * 2 + 1;

  /*	corresponding index to the 0 value of x and z	*/
  index_x = (nstepx - 1) / 2;
  index_z = (nstepz - 1) / 2;

  meshx = 2.0 * NSTD_LGM / (nstepx - 1);
  meshz = 2.0 * NSTD_LGM / (nstepz - 1);

  /*	Allocations of space vectors					*/
  r1_bar = dvector(0, nstepx - 1);
  r3_bar = dvector(0, nstepz - 1);
  r1_dim1 = dvector(0, nstepx - 1);
  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  mu_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  mu_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r1 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  var_r3 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  if (!r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 ||
      !r || !r1_bar || !r1_dim1 || !r3_bar) {
    err = "Memory allocation error (2) in lgm2fpde";
    goto FREE_RETURN;
  }

  t2 = clock();

  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  smessage("Phase 2 -convolution  , stept: %d stepx: %d stepz: %d", nstept,
           nstepx, nstepz);

  t1 = clock();

  /*	Then discretise space in the orthogonal system r1 / r3	*/
  r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
  r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

  if (disc_method) {
    disc_normal_center(r1_bar, nstepx, r1_bar[0], 1.0, &index_x, &meshx);
    disc_normal_center(r3_bar, nstepz, r3_bar[0], 1.0, &index_z, &meshz);
  } else {
    for (i = 1; i < nstepx; i++) {
      r1_bar[i] = r1_bar[i - 1] + meshx;
    }

    for (j = 1; j < nstepz; j++) {
      r3_bar[j] = r3_bar[j - 1] + meshz;
    }
  }

  /*	Corresponding r1 and r2					*/

  std1 = sqrt(phi1[nstept - 1]);
  std3 = 1.0 / rho2 *
         sqrt(phi2[nstept - 1] / alpha / alpha + rho * rho * phi1[nstept - 1] -
              2.0 * rho * rho / alpha * phi12[nstept - 1]);

  for (i = 0; i < nstepx; i++) {
    r1_dim1[i] = r1_bar[i] * std1 + fwd1[nstept - 1];
    r2_temp = rhoalpha * r1_bar[i] * std1 + fwd2[nstept - 1];

    for (j = 0; j < nstepz; j++) {
      r2[i][j] = r2_temp + rho2alpha * r3_bar[j] * std3;
    }
  }

  /*	Final payoff valuation					*/
  if (!eval_evt[nstept - 1]) {
    err = "No event at last step in lgm2f_pde";
    goto FREE_RETURN;
  }

  /*	Eval payoff */
  err =
      payoff_func(date[nstept - 1], time[nstept - 1], func_parm_tab[nstept - 1],
                  yc, lam, lam_time, nb_lam, gamma, rho, phi1[nstept - 1],
                  phi2[nstept - 1], phi12[nstept - 1], 0, nstepx - 1, 0,
                  nstepz - 1, r1_dim1, r2, nprod, values_p1);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Initialize the CNPDE_TEMP_2D		*/
  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

  if (!pde) {
    err = "Memory allocation error (3) in lgm2fpde";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  std1_next = std1;
  std3_next = std3;

  /*	now do the backward pde					*/

  for (step = nstept - 2; step >= 0; step--) {
    dt = time[step + 1] - time[step];

    r_temp = ifr[step] + fwd1[step] + fwd2[step];
    std1 = sqrt(phi1[step]);
    std3 = sqrt(std3_c1 * phi2[step] + std3_c2 * phi1[step] +
                std3_c3 * phi12[step]);

    std1_mean = std1_mid[step];
    std3_mean = std3_mid[step];

    lam1 = avg_lam[step] + std1der_mid[step] / std1_mean;
    lam2 = avg_lam[step] + gamma + std3der_mid[step] / std3_mean;

    var1_temp = var1[step] / std1_mean / std1_mean;
    var3_temp = var1[step] / std3_mean / std3_mean;

    fact1_temp = fact1 * std1_mean / std3_mean;
    rhoalpha_plus1_temp = rhoalpha_plus1 * std1_mean;

    for (i = lx; i <= ux; i++) {
      mu_r1i = -lam1 * r1_bar[i] * dt;
      mu_r3i = fact1_temp * r1_bar[i];

      r2_temp = r_temp + rhoalpha_plus1_temp * r1_bar[i];

      for (j = lz; j <= uz; j++) {
        mu_r1[i][j] = mu_r1i;
        mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;

        var_r1[i][j] = var1_temp;
        var_r3[i][j] = var3_temp;

        r[i][j] = (r2_temp + rho2alpha * r3_bar[j] * std3_mean) * dt;
      }
    }

    /*	convolve							*/

    num_f_pde_one_step_backward_2f_adi(
        pde, nstepx, r1_bar, nstepz, r3_bar, 0, nprod - 1, values_p1, mu_r1,
        mu_r3, var_r1, var_r3, r, values, lx, ux, lz, uz);

    /*	Eval payoff */
    if (eval_evt[step]) {

      /*	Corresponding r1 and r2					*/
      for (i = lx; i <= ux; i++) {
        r1_dim1[i] = r1_bar[i] * std1 + fwd1[step];
        r2_temp = rhoalpha * r1_bar[i] * std1 + fwd2[step];

        for (j = lz; j <= uz; j++) {
          r2[i][j] = r2_temp + rho2alpha * r3_bar[j] * std3;
        }
      }

      err =
          payoff_func(date[step], time[step], func_parm_tab[step], yc, lam,
                      lam_time, nb_lam, gamma, rho, phi1[step], phi2[step],
                      phi12[step], lx, ux, lz, uz, r1_dim1, r2, nprod, values);
      if (err) {
        goto FREE_RETURN;
      }
    }

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;

    std1_next = std1;
    std3_next = std3;
  }

  /* copy the result					*/
  for (i = 0; i < nprod; i++) {
    res[i] = values_p1[index_x][index_z][i];
  }

  t2 = clock();

  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /*	allocation (1)		*/

  if (phi1)
    free_dvector(phi1, 0, nstept - 1);
  if (phi2)
    free_dvector(phi2, 0, nstept - 1);
  if (phi12)
    free_dvector(phi12, 0, nstept - 1);
  if (var1)
    free_dvector(var1, 0, nstept - 1);
  if (avg_lam)
    free_dvector(avg_lam, 0, nstept - 1);
  if (fwd1)
    free_dvector(fwd1, 0, nstept - 1);
  if (fwd2)
    free_dvector(fwd2, 0, nstept - 1);
  if (std1_mid)
    free_dvector(std1_mid, 0, nstept - 1);
  if (std3_mid)
    free_dvector(std3_mid, 0, nstept - 1);
  if (std1der_mid)
    free_dvector(std1der_mid, 0, nstept - 1);
  if (std3der_mid)
    free_dvector(std3der_mid, 0, nstept - 1);

  /*	allocation (2)		*/
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  if (r1_bar)
    free_dvector(r1_bar, 0, nstepx - 1);
  if (r3_bar)
    free_dvector(r3_bar, 0, nstepz - 1);
  if (r1_dim1)
    free_dvector(r1_dim1, 0, nstepx - 1);
  if (r2)
    free_dmatrix(r2, 0, nstepx - 1, 0, nstepz - 1);
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (mu_r1)
    free_dmatrix(mu_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (mu_r3)
    free_dmatrix(mu_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r1)
    free_dmatrix(var_r1, 0, nstepx - 1, 0, nstepz - 1);
  if (var_r3)
    free_dmatrix(var_r3, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);

  return err;
}

void disc_normal_center(double *x, long nstp, double xmin, double stdn,
                        long *index, double *meshx) {
  long i;
  double d;

  (*index) = (nstp - 1) / 2;
  d = norm(xmin / stdn);
  *meshx = (0.5 - d) / (*index);
  d += (*meshx);

  x[0] = xmin;
  for (i = 1; i < nstp - 1; i++) {
    x[i] = stdn * inv_cumnorm_fast(d);
    d += (*meshx);
  }
  x[nstp - 1] = -xmin;
}

void disc_sqrt_center(double *x, long nstpx, double *time, long nstpt,
                      int nb_check, double *varmax, long *index,
                      double *coefx) {
  long i, j, k, ind, nbp, nbp2, nbtot;
  double coef, dx, dt, t;

  FILE *stream = fopen("C:\\toto.txt", "w+");

  ind = (nstpx - 1) / 2;
  (*index) = ind;

  coef = nstpx / sqrt(time[nstpt - 1]);

  dt = time[nstpt - 1] / (nb_check - 1);
  j = 1;
  nbtot = 0;

  k = ind;
  x[k] = 0.0;

  t = dt - 1.0E-08;

  while (j < nstpt) {
    while (j < nstpt && time[j] < t) {
      j++;
    }

    nbp = (long)(coef * sqrt(time[j]) + 1.0E-10);
    nbp = ((int)(nbp / 2)) * 2 + 1;
    nbp2 = (nbp - 1) / 2 - nbtot;

    dx = (NSTD_LGM * sqrt(varmax[j]) - x[k]) / nbp2;

    fprintf(stream, "%f	%d	%f	%f \n", time[j], nbp2, dx,
            sqrt(varmax[j]));

    if (dx > 0) {
      for (i = 1; i <= nbp2; i++) {
        x[k + i] = x[k + i - 1] + dx;
      }

      nbtot = (nbp - 1) / 2;
      k = ind + nbtot;
    }

    j++;
    t += dt;
  }

  for (i = 1; i <= ind; i++) {
    x[ind - i] = -x[ind + i];
  }

  fclose(stream);

  *coefx = coef;
}
