// FxSVOptMC.c : MC european option pricing in the FXSV model

#include "FxSVOptMC.h"
#include "Fx3FCalib.h"
#include "RandomGen.h"
#include "math.h"
#include "srt_h_all.h"

Err OptFXSVMC(char *fxundname, double gamma, double alpha, double **rho,
              long mat, double *K, int nK, long npaths, int nsteps, double *res,
              double *std) {
  Err err = NULL;
  SrtUndPtr und, dom_und, for_und;
  char *dom_yc, *for_yc, fxund[80];
  long j, today, spot_date, ex_date, seed = -12345678;
  int i, k, m;
  char *domname, *forname;
  double *sig_d = NULL, *sigtms_d = NULL;
  double *sig_f = NULL, *sigtms_f = NULL;
  double *sig_x = NULL, *sigtms_x = NULL;
  double *times = NULL;
  int nsig_d, nsig_f, nsig_x, ntimes;
  double lam_d, lam_f, ExTime, MatTime, ffx0;
  int idx_d = 0, idx_f = 0, idx_x = 0;
  SRandomGen rg;
  double **sv = NULL, **chol = NULL, br[4], cbr[4], *pvs = NULL;
  double dt, sqrtdt, Gamma_d, Gamma_f, V12, V123, V123_ln, pv1, pv2;

  memset(&rg, 0, sizeof(SRandomGen));

  strcpy(fxund, fxundname);
  strupper(fxund);
  rem_tick_string(fxund, fxund);

  und = lookup_und(fxund);
  if (!und) {
    err = serror("Couldn't find underlying named %s", fxund);
    goto FREE_RETURN;
  }
  if (get_underlying_type(und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", fxund);
    goto FREE_RETURN;
  }
  if (get_mdltype_from_fxund(und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", fxund);
    goto FREE_RETURN;
  }

  domname = get_domname_from_fxund(und);
  err = Get_LGM_TermStructure2(domname, &sigtms_d, &sig_d, &nsig_d, &lam_d);
  if (err)
    goto FREE_RETURN;
  lam_d = 1.0 / lam_d;

  forname = get_forname_from_fxund(und);
  err = Get_LGM_TermStructure2(forname, &sigtms_f, &sig_f, &nsig_f, &lam_f);
  if (err)
    goto FREE_RETURN;
  lam_f = 1.0 / lam_f;

  err = srt_f_display_FXBS_TermStruct(fxund, &nsig_x, &sigtms_x, &sig_x);
  if (err)
    goto FREE_RETURN;

  today = get_today_from_underlying(und);
  for (i = 0; i < nsig_x; i++)
    sigtms_x[i] = (sigtms_x[i] - today) * YEARS_IN_DAY;

  ex_date = add_unit(mat, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
  ExTime = (ex_date - today) * YEARS_IN_DAY;
  if (ExTime < 1e-5) {
    err = serror("Exercise date <= today");
    goto FREE_RETURN;
  }
  MatTime = (mat - today) * YEARS_IN_DAY;

  ntimes = nsig_d + nsig_f + nsig_x + 2;
  times = (double *)calloc(ntimes, sizeof(double));
  if (!times) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }
  memcpy(times, sigtms_d, nsig_d * sizeof(double));
  memcpy(times + nsig_d, sigtms_f, nsig_f * sizeof(double));
  memcpy(times + nsig_d + nsig_f, sigtms_x, nsig_x * sizeof(double));
  times[ntimes - 2] = 0.0;
  times[ntimes - 1] = ExTime;

  num_f_sort_vector(ntimes, times);
  num_f_unique_vector(&ntimes, times);
  while (times[ntimes - 1] > ExTime + 1e-5)
    ntimes--;
  num_f_fill_vector_newalgo(&ntimes, &times, (long)(nsteps * ExTime + 1e-5));

  dom_und = lookup_und(domname);
  for_und = lookup_und(forname);
  dom_yc = get_ycname_from_irund(dom_und);
  for_yc = get_ycname_from_irund(for_und);

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  ffx0 = get_spot_from_fxund(und) * swp_f_df(spot_date, mat, for_yc) /
         swp_f_df(spot_date, mat, dom_yc);

  err = ABS_Init(&rg, seed, npaths, 4, 1);
  if (err)
    goto FREE_RETURN;

  sv = dmatrix(0, npaths - 1, 0, 2);
  chol = dmatrix(0, 3, 0, 3);
  pvs = (double *)calloc(npaths, sizeof(double));
  if (!sv || !chol || !pvs) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  // Initialize state variables at time 0 for all paths
  for (j = 0; j < npaths; j++) {
    sv[j][0] = sv[j][2] = log(ffx0);
    sv[j][1] = 1.0;
  }

  // Cholesky decomposition
  err = choldc(4, rho, chol);
  if (err)
    goto FREE_RETURN;

  // Move forward through time
  for (i = 0; i < ntimes - 1; i++) {
    dt = times[i + 1] - times[i];
    sqrtdt = sqrt(dt);
    while (idx_d < nsig_d - 1 && sigtms_d[idx_d] < times[i] + 1e-5)
      idx_d++;
    while (idx_f < nsig_f - 1 && sigtms_f[idx_f] < times[i] + 1e-5)
      idx_f++;
    while (idx_x < nsig_x - 1 && sigtms_x[idx_x] < times[i] + 1e-5)
      idx_x++;

    Gamma_d =
        -sig_d[idx_d] * (1.0 - exp(-lam_d * (MatTime - times[i]))) / lam_d;
    Gamma_f =
        -sig_f[idx_f] * (1.0 - exp(-lam_f * (MatTime - times[i]))) / lam_f;

    V12 = Gamma_f * Gamma_f / 2.0 + Gamma_d * Gamma_d / 2.0 -
          rho[0][1] * Gamma_f * Gamma_d;
    V123_ln = sig_x[idx_x] * sig_x[idx_x] / 2.0 +
              rho[1][2] * sig_x[idx_x] * Gamma_f -
              rho[0][2] * sig_x[idx_x] * Gamma_d + V12;

    for (j = 0; j < npaths; j++) {
      // Generate independent brownian increments
      for (k = 0; k < 4; k++) {
        err = rg.Gauss(&rg, &br[k]);
        if (err)
          goto FREE_RETURN;
        br[k] *= sqrtdt;
      }
      // Correlate the increments
      memset(cbr, 0, 4 * sizeof(double));
      for (k = 0; k < 4; k++)
        for (m = 0; m <= k; m++)
          cbr[k] += chol[k][m] * br[m];

      // Variance:
      V123 = sig_x[idx_x] * sig_x[idx_x] * sv[j][1] * sv[j][1] / 2.0 +
             rho[1][2] * sig_x[idx_x] * sv[j][1] * Gamma_f -
             rho[0][2] * sig_x[idx_x] * sv[j][1] * Gamma_d + V12;

      // Calculate state variables at the next time step
      sv[j][0] += -V123 * dt + sig_x[idx_x] * sv[j][1] * cbr[2] +
                  Gamma_f * cbr[1] - Gamma_d * cbr[0];
      sv[j][1] += gamma * (1.0 - sv[j][1]) * dt + alpha * cbr[3];
      sv[j][2] += -V123_ln * dt + sig_x[idx_x] * cbr[2] + Gamma_f * cbr[1] -
                  Gamma_d * cbr[0];
    }
  }

  // Calculate PVs for all paths and the expectation for all strikes
  for (k = 0; k < nK; k++) {
    res[k] = 0.0;
    for (j = 0; j < npaths; j++) {
      pv1 = exp(sv[j][0]) - K[k];
      if (K[k] > ffx0)
        pv1 = -pv1;
      if (pv1 < 0.0)
        pv1 = 0.0;

      pv2 = exp(sv[j][2]) - K[k];
      if (K[k] > ffx0)
        pv2 = -pv2;
      if (pv2 < 0.0)
        pv2 = 0.0;

      pvs[j] = pv1 - pv2;
      res[k] += pvs[j] / npaths;
    }
    // Calculate std
    std[k] = 0.0;
    for (j = 0; j < npaths; j++)
      std[k] += (pvs[j] - res[k]) * (pvs[j] - res[k]) / (npaths - 1);
    std[k] = sqrt(std[k] / npaths);

    res[k] *= swp_f_df(today, mat, dom_yc);
    std[k] *= swp_f_df(today, mat, dom_yc);

    err = Fx3DFxOption(fxund, K[k], ex_date, mat, mat, "CALL", &pv1);
    if (err)
      goto FREE_RETURN;

    res[k] += pv1;
  }

FREE_RETURN:
  free(sig_d);
  free(sigtms_d);
  free(sig_f);
  free(sigtms_f);
  free(sig_x);
  free(sigtms_x);
  free(times);

  ABS_Free(&rg);
  if (sv)
    free_dmatrix(sv, 0, npaths - 1, 0, 2);
  if (chol)
    free_dmatrix(chol, 0, 3, 0, 3);
  free(pvs);

  return err;
}
