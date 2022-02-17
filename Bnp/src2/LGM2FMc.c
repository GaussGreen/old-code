/* ==========================================================================
   FILE_NAME:	Fx3FMc.c

   PURPOSE:		Monte Carlo Three Factors

   DATE:		08/28/00
   ========================================================================== */

#include "Fx3FUtils.h"
#include "MCEBOptimisation.h"
#include "math.h"
#include "srt_h_all.h"

/*	Main function under QTfinal */
/*	------------- */
Err mc_main_lgm2f(
    /*	Time data */
    long nb_paths, int nb_col, double *time, double *date, long nb_dates,

    /* Model datas */
    int jumping_num, double *dom_fwd1, double *dom_fwd2, double *dom_exp1,
    double *dom_exp2, double *dom_phi1, double *dom_phi2, double *dom_phi12,
    double *dom_gam1_fwd, double *dom_gam2_fwd, double *dom_bond_pay,
    double *dom_gam1_pay, double *dom_gam2_pay, double ***covar,

    /*	Product data */
    void **func_parm_tab,

    /* do PECS adjustment */
    int do_pecs,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or
       NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date, double evt_time, void *func_parm,
        /* Market data */
        double R1D, double R2D,
        /* Results */
        int nb_col, double *res, int *stop_path),
    /*	Results */
    double **res) {
  int stop_path;
  long i, j, k, l, m;
  double df;
  double R1D, R2D;

  double *temp = NULL, *res_evt = NULL, *sum_price = NULL, *sum_2price = NULL,
         **matrix = NULL, **chol = NULL, ***save_values = NULL;

  double *matrixi;

  double tempc[2];
  double sum;
  int optim_today;
  clock_t time1, time2;

  Err err = NULL;

  time1 = clock();

  /*	nb_paths has to be odd */
  nb_paths = 2 * ((long)(nb_paths / 2)) + 1;

  temp = dvector(0, nb_col - 1);
  res_evt = dvector(0, nb_col - 1);
  sum_price = dvector(0, nb_col - 1);
  sum_2price = dvector(0, nb_col - 1);
  matrix = dmatrix(0, nb_paths - 1, 0, 2 * (nb_dates - 1) - 1);
  chol = dmatrix(0, 1, 0, 1);

  if (!temp || !res_evt || !sum_price || !sum_2price || !matrix || !covar) {
    err = "Memory allocation (1) failure in mc_main_lgm2f";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(nb_paths, nb_dates, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  /* Initialisation */
  memset(sum_price, 0, nb_col * sizeof(double));
  memset(sum_2price, 0, nb_col * sizeof(double));

  /* fill the Brownian matrix */
  err = balsam_generation(nb_paths, 2 * (nb_dates - 1), matrix);
  if (err) {
    goto FREE_RETURN;
  }
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  time1 = clock();
  if (do_pecs) {
    adjust_emp_covar(matrix, nb_paths, 2 * (nb_dates - 1));
  }
  time2 = clock();
  smessage("Phase 2 -PECS adjustment  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  time1 = clock();

  /* Adjust the  matrix to have the right covariance matrix */
  m = 0;
  for (j = 1; j < nb_dates; j++) {
    /* Find the Cholesky of the covar */
    nr_choldc(2, covar[j], chol);

    for (i = 0; i < nb_paths; i++) {
      for (l = 0; l < 2; l++) {
        sum = 0.0;
        for (k = 0; k < 2; k++) {
          sum += chol[l][k] * matrix[i][m + k];
        }

        tempc[l] = sum;
      }

      for (l = 0; l < 2; l++) {
        matrix[i][m + l] = tempc[l];
      }
    }

    m += 2;
  }

  time2 = clock();
  smessage("Phase 3 -correlation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  time1 = clock();
  /* Do the simulation */
  for (i = 0; i < nb_paths; i++) {
    R1D = R2D = 0.0;
    df = 1.0;

    matrixi = matrix[i];

    memset(temp, 0, nb_col * sizeof(double));
    stop_path = 0;

    if (init_func) {
      init_func();
    }

    /*	If today is an event date */
    if (func_parm_tab[0]) {
      err = payoff_func(date[0], time[0], func_parm_tab[0], R1D, R2D, nb_col,
                        res_evt, &stop_path);

      if (err) {
        goto FREE_RETURN;
      }

      if (jumping_num) {
        for (k = 0; k < nb_col; k++) {
          temp[k] += res_evt[k] * df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, df,
                                         params);
        }

        df *= dom_bond_pay[0];
      } else {
        df = dom_bond_pay[0];

        for (k = 0; k < nb_col; k++) {
          temp[k] += res_evt[k] / df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, df,
                                         params);
        }
      }
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < nb_dates; j++) {
      R1D = R1D * dom_exp1[j] + dom_fwd1[j] + matrixi[m];
      R2D = R2D * dom_exp2[j] + dom_fwd2[j] + matrixi[m + 1];

      err = payoff_func(date[j], time[j], func_parm_tab[j], R1D, R2D, nb_col,
                        res_evt, &stop_path);

      if (err) {
        goto FREE_RETURN;
      }

      /* Discount the payoff */
      if (jumping_num) {
        for (k = 0; k < nb_col; k++) {
          temp[k] += res_evt[k] * df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[j], res_evt, i, df,
                                         params);
        }

        df *= dom_bond_pay[j] *
              exp(-dom_gam1_pay[j] * R1D - dom_gam2_pay[j] * R2D);
      } else {
        df = dom_bond_pay[j] *
             exp(-dom_gam1_pay[j] * R1D - dom_gam2_pay[j] * R2D);

        for (k = 0; k < nb_col; k++) {
          temp[k] += res_evt[k] / df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[j], res_evt, i, df,
                                         params);
        }
      }

      m += 2;
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

  /* Save results */
  for (k = 0; k < nb_col; k++) {
    res[k][0] = sum_price[k];
    res[k][1] = (sum_2price[k] - sum_price[k] * sum_price[k]) / nb_paths;

    if (fabs(res[k][1]) < 1.0E-10) {
      res[k][1] = 0.0;
    } else {
      res[k][1] = sqrt(res[k][1]);
    }
  }

  time2 = clock();
  smessage("Phase 4 -evaluation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    time1 = clock();

    /* free the useless memory */
    if (matrix) {
      free_dmatrix(matrix, 0, nb_paths - 1, 0, 2 * (nb_dates - 1) - 1);
      matrix = NULL;
    }

    if (chol) {
      free_dmatrix(chol, 0, 1, 0, 1);
      chol = NULL;
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
    smessage("Phase 5 -optimisation  , time in sec: %.2f",
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
    free_dmatrix(matrix, 0, nb_paths - 1, 0, 2 * (nb_dates - 1) - 1);
  }

  if (chol) {
    free_dmatrix(chol, 0, 2, 0, 2);
  }

  mceb_free_savevalues_for_GRFN(save_values, nb_paths, nb_dates, params);

  return err;
}
