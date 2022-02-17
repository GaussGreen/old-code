/* ==========================================================================
   FILE_NAME:	Fx3FMc.cxx

   PURPOSE:		Monte Carlo Three Factors

   DATE:		08/28/00
   ========================================================================== */

#include "MCEBOptimisation.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "stdio.h"

/*	Main function */
/*	------------- */
Err mc_main_3dBetaDLMfx(
    /*	Time data */
    long nb_paths, int nb_col, double *time, double *date, long nb_dates,

    /*	Model Data */
    double *dom_fwd, double *dom_std, double *for_fwd_const,
    double *for_fwd_lin, double *for_std, double *fx_fwd, double *fx_std,

    double *bond_pay_const, double *bond_pay_lin,

    double *dom_for_cov, double *dom_fx_cov, double *for_fx_cov,

    /*	Product data */
    void **func_parm_tab,

    /* do PECS adjustment */
    int do_pecs,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path
                    or NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,
                       /* Market data */
                       double Xdom, double Yfor, double Zfx,
                       /* Results */
                       int nb_col, double *res, int *stop_path),
    /*	Results */
    double **res) {
  int stop_path;
  long i, j, k, m;
  double df;
  double X, Y, Z;
  clock_t time1, time2;
  double *temp = NULL, *res_evt = NULL, *sum_price = NULL, *sum_2price = NULL,
         **matrix = NULL, ***save_values = NULL, **infos = NULL;

  int optim_today;

  Err err = NULL;

  time1 = clock();

  /* nb_paths has to be odd */
  nb_paths = 2 * ((long)(nb_paths / 2)) + 1;

  temp = dvector(0, nb_col - 1);
  res_evt = dvector(0, nb_col - 1);
  sum_price = dvector(0, nb_col - 1);
  sum_2price = dvector(0, nb_col - 1);
  matrix = dmatrix(0, nb_paths - 1, 0, 3 * (nb_dates - 1) - 1);

  if (!temp || !res_evt || !sum_price || !sum_2price || !matrix) {
    err = "Memory allocation failure in mc_main_3dBetaDLMfx";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(nb_paths, nb_dates, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  for (k = 0; k < nb_col; k++) {
    sum_price[k] = sum_2price[k] = 0.0;
  }

  /* fill the Brownian matrix */
  err = balsam_generation(nb_paths, 3 * (nb_dates - 1), matrix);

  if (err)
    goto FREE_RETURN;

  time2 = clock();
  smessage("Phase 2 -BalSam generation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  if (do_pecs) {
    adjust_emp_covar(matrix, nb_paths, 3 * (nb_dates - 1));
  }

  time2 = clock();
  smessage("Phase 3 -PECS adjustment        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  err = correlate_variable_covar(
      &(dom_std[1]), &(for_std[1]), &(fx_std[1]), &(dom_for_cov[1]),
      &(dom_fx_cov[1]), &(for_fx_cov[1]), nb_paths, nb_dates - 1, matrix);

  if (err)
    goto FREE_RETURN;

  time2 = clock();
  smessage("Phase 4 -correlation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  for (i = 0; i < nb_paths; i++) {
    X = Y = Z = 0;

    memset(temp, 0, nb_col * sizeof(double));
    stop_path = 0;

    if (init_func) {
      init_func();
    }

    /*	If today is an event date */
    if (func_parm_tab[0]) {
      df = exp(bond_pay_const[0]);

      err = payoff_func(date[0], time[0], func_parm_tab[0], X, Y, Z, nb_col,
                        res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < nb_col; k++) {
        temp[k] += res_evt[k] / df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, df, params);
      }
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < nb_dates; j++) {
      X += dom_fwd[j] + matrix[i][m];
      Y += for_fwd_const[j] + for_fwd_lin[j] * Z + matrix[i][m + 1];
      Z += fx_fwd[j] + matrix[i][m + 2];

      df = exp(bond_pay_const[j] + bond_pay_lin[j] * X);

      err = payoff_func(date[j], time[j], func_parm_tab[j], X, Y, Z, nb_col,
                        res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < nb_col; k++) {
        temp[k] += res_evt[k] / df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[j], res_evt, i, df, params);
      }

      m += 3;
    }

    for (k = 0; k < nb_col; k++) {
      sum_price[k] += temp[k] / nb_paths;
      sum_2price[k] += temp[k] * temp[k] / nb_paths;
    }

    if (do_optimisation && params->iKnockInCol) {
      /* we recopy in the col pay the pv of the column */
      for (j = 0; j < nb_dates; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              temp[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < nb_col; k++) {
    static double temp;

    res[k][0] = sum_price[k];
    temp = (sum_2price[k] - sum_price[k] * sum_price[k]) / nb_paths;
    if (temp > 1.0e-08) {
      res[k][1] = sqrt(temp);
    } else {
      res[k][1] = 0.0;
    }
  }

  time2 = clock();
  smessage("Phase 5 -evaluation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    time1 = clock();

    /* free the useless memory */
    if (matrix) {
      free_dmatrix(matrix, 0, nb_paths - 1, 0, 3 * (nb_dates - 1) - 1);
      matrix = NULL;
    }

    if (time[0] < 1.0E-08 && optimise[0]) {
      optimise[0] = 0;
      optim_today = 1;
    } else {
      optim_today = 0;
    }

    err = find_and_optimise_boundary(save_values, nb_dates, nb_paths, optimise,
                                     params, &(res[nb_col][0]),
                                     &(res[nb_col][1]));

    if (err)
      goto FREE_RETURN;

    if (optim_today) {
      if (params->iIsKO) {
        if (res[nb_col][0] < 0.0) {
          res[nb_col][0] = 0.0;
          res[nb_col][1] = 0.0;
        }
      } else {
        if (save_values[0][params->iNbIndex][0] > res[nb_col][0]) {
          res[nb_col][0] = save_values[0][params->iNbIndex][0];
          res[nb_col][1] = 0.0;
        }
      }

      optimise[0] = 1;
    }

    time2 = clock();
    smessage("Phase 6 -optimisation        , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (temp) {
    free_dvector(temp, 0, nb_col - 1);
  }

  if (res_evt) {
    free_dvector(res_evt, 0, nb_col - 1);
  }

  if (sum_price) {
    free_dvector(sum_price, 0, nb_col - 1);
  }

  if (sum_2price) {
    free_dvector(sum_2price, 0, nb_col - 1);
  }

  if (matrix) {
    free_dmatrix(matrix, 0, nb_paths - 1, 0, 3 * (nb_dates - 1) - 1);
  }

  mceb_free_savevalues_for_GRFN(save_values, nb_paths, nb_dates, params);

  if (infos) {
    free_dmatrix(infos, 0, nb_dates - 1, 0, 2 + 2 * nb_dates);
  }

  return err;
}