/* ==========================================================================
   FILE_NAME:	Fx3FMc.cxx

   PURPOSE:		Monte Carlo Three Factors

   DATE:		08/28/00
   ========================================================================== */

#include "MCEBOptimisation.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

/*	Cholesky decomposition which disregards errors */
static void nr_choldc_tolerance(int n, /* 0.. n-1*/
                                double **cov, double **chol) {
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      for (sum = cov[i][j], k = i - 1; k >= 0; k--) {
        sum -= cov[i][k] * cov[j][k];
      }
      if (i == j) {
        if (sum > 0.0) {
          chol[i][i] = sqrt(sum);
        } else {
          chol[i][i] = 0.0;
        }
      } else {
        cov[j][i] = chol[j][i] = sum / chol[i][i];
      }
    }
  }
}

static void correlate_underlyings(double ***cov, long nb_paths, int n,
                                  int nb_dates, double **matrix) {
  long i, j, l, k;
  double **chol;
  double **temp;
  double sum;

  chol = dmatrix(0, n - 1, 0, n - 1);
  temp = dmatrix(0, nb_paths - 1, 0, n - 1);

  for (l = 0; l < nb_dates; l++) {
    nr_choldc_tolerance(n, cov[l + 1], chol);

    for (i = 0; i < nb_paths; i++) {
      for (j = 0; j < n; j++) {
        sum = 0.0;
        for (k = 0; k <= j; k++) {
          sum += chol[j][k] * matrix[i][l * n + k];
        }
        temp[i][j] = sum;
      }
    }
    for (i = 0; i < nb_paths; i++) {
      for (j = 0; j < n; j++) {
        matrix[i][l * n + j] = temp[i][j];
      }
    }
  }

  if (chol)
    free_dmatrix(chol, 0, n - 1, 0, n - 1);
  if (temp)
    free_dmatrix(temp, 0, nb_paths - 1, 0, n - 1);

  chol = NULL;
  temp = NULL;
}

/*	Main function */
/*	------------- */
Err mc_main_multi_3dfx(
    /*	Time data */
    long nb_paths, long nb_col, double *time, long *date, long nb_dates,
    int do_pecs, int do_optimisation, int *optimise, MCEBPARAMS params,
    LINK_UND link_und, void **func_parm_tab,
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date, double evt_time, void *func_parm,
        /* Market data */
        LINK_UND link, double *sv,
        /* Results */
        int nb_col, double *res),
    double **res) {
  long i, j, k, l;
  double t1, t2;
  double df;

  double *temp = NULL, *res_evt = NULL, *sum_price = NULL, *sum_2price = NULL,
         **matrix = NULL, *sv = NULL, ***save_values = NULL;

  clock_t time1, time2;
  int optim_today;
  Err err = NULL;

  time1 = clock();

  /* nb_paths has to be odd */
  nb_paths = 2 * ((long)(nb_paths / 2)) + 1;

  temp = dvector(0, nb_col - 1);
  res_evt = dvector(0, nb_col - 1);
  sum_price = dvector(0, nb_col - 1);
  sum_2price = dvector(0, nb_col - 1);
  matrix =
      dmatrix(0, nb_paths - 1, 0, (link_und->num_und) * (nb_dates - 1) - 1);
  sv = dvector(0, link_und->num_und - 1);

  if (!temp || !res_evt || !sum_price || !sum_2price || !matrix || !sv) {
    err = "Memory allocation failure in mc_main_multi_3dfx";
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
  err =
      balsam_generation(nb_paths, (link_und->num_und) * (nb_dates - 1), matrix);
  if (err) {
    goto FREE_RETURN;
  }
  time2 = clock();
  smessage("Phase 2 -BalSam generation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  if (do_pecs) {
    adjust_emp_covar(matrix, nb_paths, (link_und->num_und) * (nb_dates - 1));
  }
  time2 = clock();
  smessage("Phase 3 -PECS adjustment        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  correlate_underlyings(link_und->covariance, nb_paths, link_und->num_und,
                        nb_dates - 1, matrix);
  if (err) {
    goto FREE_RETURN;
  }
  time2 = clock();
  smessage("Phase 4 -correlation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  for (i = 0; i < nb_paths; i++) {
    t1 = 0;

    memset(sv, 0, link_und->num_und * sizeof(double));
    memset(temp, 0, nb_col * sizeof(double));

    /*	If today is an event date */
    if (func_parm_tab[0]) {
      df = link_und->dom_bond_pay[0];

      err = payoff_func(date[0], time[0], func_parm_tab[0], link_und, sv,
                        nb_col, res_evt);

      if (err) {
        goto FREE_RETURN;
      }

      for (k = 0; k < nb_col; k++) {
        temp[k] += res_evt[k] / df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, df, params);
      }
    }

    for (j = 1; j < nb_dates; j++) {
      t2 = time[j];

      /* get next fx underlyings */
      for (l = 0; l < link_und->num_und; l++) {
        if (link_und->type[l] == 2) {
          sv[l] += link_und->fwd[l][j] +
                   link_und->beta[link_und->for_forex[l]][j - 1] *
                       sv[link_und->for_forex[l]] -
                   link_und->beta[link_und->dom_forex[l]][j - 1] *
                       sv[link_und->dom_forex[l]] +
                   matrix[i][(link_und->num_und) * (j - 1) + l];
        }
      }

      /* get next LGM underlyings */
      for (l = 0; l < link_und->num_und; l++) {
        if (link_und->type[l] < 2) {
          sv[l] = sv[l] * exp(-link_und->lambda[l] * (t2 - t1)) +
                  link_und->fwd[l][j] +
                  matrix[i][(link_und->num_und) * (j - 1) + l];
        }
      }

      df = link_und->dom_bond_pay[j] *
           exp(link_und->dom_beta_pay[j] *
               (-0.5 * link_und->dom_beta_pay[j] *
                    link_und->phi[link_und->dom_index][j] +
                sv[link_und->dom_index]));

      err = payoff_func(date[j], t2, func_parm_tab[j], link_und, sv, nb_col,
                        res_evt);

      if (err) {
        goto FREE_RETURN;
      }

      for (k = 0; k < nb_col; k++) {
        temp[k] += res_evt[k] / df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[j], res_evt, i, df, params);
      }

      t1 = t2;
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
    res[k][0] = sum_price[k];
    res[k][1] = sqrt((sum_2price[k] - sum_price[k] * sum_price[k]) / nb_paths);
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
    temp = NULL;
  }

  if (res_evt) {
    free_dvector(res_evt, 0, nb_col - 1);
    res_evt = NULL;
  }

  if (sum_price) {
    free_dvector(sum_price, 0, nb_col - 1);
    sum_price = NULL;
  }

  if (sum_2price) {
    free_dvector(sum_2price, 0, nb_col - 1);
    sum_2price = NULL;
  }

  /* unknown error here only in Release version... */
  if (matrix) {
    free_dmatrix(matrix, 0, nb_paths - 1, 0,
                 (link_und->num_und) * (nb_dates - 1) - 1);
    matrix = NULL;
  }

  if (sv) {
    free_dvector(sv, 0, link_und->num_und - 1);
    sv = NULL;
  }

  mceb_free_savevalues_for_GRFN(save_values, nb_paths, nb_dates, params);

  return err;
}
