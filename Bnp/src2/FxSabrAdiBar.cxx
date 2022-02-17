/* ==========================================================================
   FILE_NAME:	LGM2FSabr.cxx

   PURPOSE:		ADI implementation of the Sabr Fx model.
                                Discretisation is ADI.

   DATE:		03/16/01

   AUTHOR:		L.C.
   ========================================================================== */

#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "FxSabrGrfn.h"
#include "math.h"
#include "opfnctns.h"

#define NSTD_FX 5.0
#define NSTD_VOL 5.0

#define MINNODE 25
#define MAXNODE2 50

Err FxSabr_adi_bar(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int bar_col,
    int *is_bar, int *is_up,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       long today, double spot_fx, /*	The cash one */
                       void *dom_yc, void *for_yc,

                       /* Grid data	*/
                       int l1, int u1, int l2, int u2, double *x,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res) {
  Err err = NULL;

  long today;
  int i, j, k, step;
  int index_x, index_z, nstepx, nstepz;
  double fux, flx, dt, dft;
  double alpha2, alpharho, alpha2rho2, alpharhobeta5;
  double temp;
  double std1, std3, spotBeta;
  double z0;

  double *exp_h = NULL, *expect_z = NULL, *expect_z2 = NULL, *drift_z = NULL,
         *x = NULL, *x_bar = NULL, *z = NULL, ***values = NULL,
         ***values_p1 = NULL, ***values_temp = NULL, **mux = NULL, **muz = NULL,
         **muzinit = NULL, **sigBeta = NULL, **varx = NULL, **varz = NULL,
         **varxinit = NULL, **varzinit = NULL, **r = NULL;

  double *varxi, *varzi, *muzi, **valuesi, *valuesij;

  double exp_z, const_sig, sigBetaij, const1, const2;
  double const_sigi, const_varx, const_muz1, const_muz2, const_varz;
  double adj_lambda, adj_lambda1, adj_lambda2;

  int lx, ux, lz, uz;

  double sig, varxf;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr1, *pde1 = NULL, pdestr2, *pde2 = NULL, pdestr3,
                             *pde3 = NULL, pdestr4, *pde4 = NULL;

  int ms1, me1, ms2, me2, ms3, me3, ms4, me4, n1, n2, n3, n4, adj_index,
      change_lim;

  double logx, logmax, logmin, save1, save2, save3;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  if (fabs(beta - 1.0) < 1.0E-08) {
    beta = 1.0 - 1.0E-08;
  }

  t1 = clock();

  /* Constant	calculations */
  today = (long)(date[0] + 1.0e-06);

  spotBeta = pow(spot_fx, beta - 1.0);

  alpha2 = alpha * alpha;
  alpharho = alpha * rho;
  alpharhobeta5 = 0.5 * alpharho * beta;
  alpha2rho2 = alpha * alpha * (1.0 - rho * rho);
  const_sig = alpharho / (1.0 - beta);

  z0 = rho * alpha / (1.0 - beta) * (1.0 / spotBeta - 1) - sig0;

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepfx / 2)) * 2 + 1;
  nstepz = ((int)(nstepvol / 2)) * 2 + 1;

  /*	we want at least three points in each directions
   */
  if (nstepx < 3) {
    nstepx = 3;
  }
  if (nstepz < 3) {
    nstepz = 3;
  }

  /*	Memory allocations */

  x = calloc(nstepx, sizeof(double));
  z = dvector(0, nstepz - 1);

  exp_h = dvector(0, nstp - 1);
  expect_z = dvector(0, nstp - 1);
  drift_z = dvector(0, nstp - 1);

  if (!exp_h || !expect_z || !drift_z || !x || !z) {
    err = "Memory allocation error (1) in FxSabr_adi";
    goto FREE_RETURN;
  }

  /*	Calculate the standard deviation of X and the expectaion of Z	*/

  const1 = -0.5 * alpharho * beta * sig0 * sig0 * spotBeta;
  const2 = -alpharho * (1.0 / spotBeta - 1.0) / (1.0 - beta);

  temp = 0.0;
  exp_z = z0;
  varxf = 0.0;

  for (i = 0; i < nstp - 1; i++) {
    dt = (time[i + 1] - time[i]) / 2.0;
    sig = exp(temp);
    varxf += sig * sig * 2.0 * dt;

    exp_z += const1 * sig * exp(alpha2 * time[i]) / (drift[i] + alpha2) *
                 (exp((drift[i] + alpha2) * dt) - 1.0) -
             const2 / sig * (exp(-drift[i] * dt) - 1.0);

    expect_z[i] = exp_z;

    temp += drift[i] * dt;
    sig = exp(temp);

    exp_h[i] = sig;

    drift_z[i] = (const1 * sig * exp(alpha2 * (time[i] + dt)) +
                  const2 / sig * drift[i]) *
                 2.0 * dt;

    exp_z += const1 * sig * exp(alpha2 * (time[i] + dt)) / (drift[i] + alpha2) *
                 (exp((drift[i] + alpha2) * dt) - 1.0) -
             const2 / sig * (exp(-drift[i] * dt) - 1.0);

    temp += drift[i] * dt;
  }

  std1 = sig0 * sqrt(varxf) * spotBeta;
  std3 = sig0 * sqrt(alpha2rho2 * time[nstp - 1]);

  /*	Then discretise space in the orthogonal system x / z */

  /*	Lognormal bounds of X */

  fux = find_beta_lim(spot_fx, sig0 * sqrt(varxf), beta, 0.001);

  flx = find_beta_lim(spot_fx, sig0 * sqrt(varxf), beta, 0.999);

  /*
  flx = max(exp(log(spot_fx) -std1 * std1 / 2.0 - NSTD_FX * std1)        ,
  -1.0); fux = exp(log(flx)  + 2.0 * NSTD_FX * std1);
  */

  /*	Bound for the calculation of max and min drift and var */
  logmax = log(spot_fx) - std1 * std1 / 2.0 + floorstd * std1;
  logmin = log(spot_fx) - std1 * std1 / 2.0 - floorstd * std1;

  logmax = log(fux);
  logmin = log(flx);

  /*	Discretisation of x and z */

  disc_linleft_logright(x, nstepx, spot_fx, flx, fux, std1, NSTD_FX, &index_x);
  disc_normal_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3,
                     NSTD_VOL, &index_z);

  /*	Precalculate variances and expectations				*/

  /* first allocate memory */
  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  mux = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muzinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  sigBeta = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varx = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varxinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varzinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  if (!values || !values_p1 || !mux || !muz || !muzinit || !sigBeta || !varx ||
      !varxinit || !varz || !varzinit || !r) {
    err = "Memory allocation error (2) in FxSabr_adi";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstepx; i++) {
    logx = min(max(log(x[i]), logmin), logmax);
    temp = exp((beta - 1.0) * logx);

    sigBeta[i][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
    muzinit[i][0] = alpharhobeta5 * temp;
    varxinit[i][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

    for (j = 0; j < nstepz; j++) {
      mux[i][j] = 0.0;
      r[i][j] = 0.0;
    }
  }

  /* save precalculations at fwd point */
  save1 = sigBeta[index_x][0];
  save2 = muzinit[index_x][0];
  save3 = varxinit[index_x][0];

  /*	Final payoff valuation					*/
  if (!eval_evt[nstp - 1]) {
    err = "No event at last step in FxSabr_Adi";
    goto FREE_RETURN;
  }

  if (err) {
    goto FREE_RETURN;
  }

  /*	Initialize the 4 CNPDE_TEMP_2D		*/

  ms1 = 0;
  me1 = bar_col - 1;
  n1 = me1 - ms1 + 1;

  if (n1 > 0) {
    pde1 = &pdestr1;

    num_f_pde_init_2d_adi(pde1, nstepx, nstepz, n1);
  }

  pde2 = &pdestr2;
  ms2 = bar_col;
  me2 = bar_col;
  n2 = 1;

  num_f_pde_init_2d_adi(pde2, nstepx, nstepz, n2);

  ms3 = bar_col + 1;
  me3 = nprod - 1;
  n3 = me3 - ms3 + 1;

  if (n3 > 0) {
    pde3 = &pdestr3;

    num_f_pde_init_2d_adi(pde3, nstepx, nstepz, n3);
  }

  /* fourth one is a regular PDE on everything */
  pde4 = &pdestr4;
  ms4 = 0;
  me4 = nprod - 1;
  n4 = nprod;

  num_f_pde_init_2d_adi(pde4, nstepx, nstepz, n4);

  if ((!pde1 && n1) || !pde2 || (!pde3 && n3) || !pde4) {
    err = "Memory allocation error (2) in FxSabr_Adi";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  /* Check that the barrier are inside the gride */
  for (i = 0; i < nstp; i++) {
    if (is_bar[i] && (bar_lvl[i] < x[1] || bar_lvl[i] > x[nstepx - 2])) {
      /* remove the barrier */
      is_bar[i] = 0;
    }
  }

  /* Initialize the barrier */

  if (is_bar[nstp - 1]) {
    if (is_up[nstp - 1]) {
      while (x[ux] > bar_lvl[nstp - 1] - 1.0E-08) {
        ux--;
      }

      ux++;

      /* update precalculations */
      x[ux] = bar_lvl[nstp - 1];

      logx = min(max(log(x[ux]), logmin), logmax);
      temp = exp((beta - 1.0) * logx);
      sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
      muzinit[ux][0] = alpharhobeta5 * temp;
      varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

      if (err) {
        goto FREE_RETURN;
      }

    } else {
      while (x[lx] < bar_lvl[nstp - 1] + 1.0E-08) {
        lx++;
      }

      lx--;

      /* update precalculations */
      x[lx] = bar_lvl[nstp - 1];

      logx = min(max(log(x[lx]), logmin), logmax);
      temp = exp((beta - 1.0) * logx);
      sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
      muzinit[lx][0] = alpharhobeta5 * temp;
      varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

      if (err) {
        goto FREE_RETURN;
      }
    }
  }

  /*	Eval payoff */
  err = payoff_func(date[nstp - 1], time[nstp - 1], func_parm_tab[nstp - 1],
                    today, spot_fx, dom_yc, for_yc, 0, nstepx - 1, 0,
                    nstepz - 1, x, nprod, values_p1);

  /*	now do the backward pde					*/

  for (step = nstp - 2; step >= 0; step--) {

    dt = time[step + 1] - time[step];

    adj_lambda1 = lambda / (sig0 * exp_h[step]) * dt;
    adj_lambda2 = lambda * dt;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];

      const_sigi = sigBeta[i][0] - expect_z[step] * exp_h[step];
      const_varx = varxinit[i][0];
      const_muz1 = -muzinit[i][0] / exp_h[step];
      const_muz2 =
          -sigBeta[i][0] / exp_h[step] * drift[step] * dt - drift_z[step];
      const_varz = alpha2rho2 / exp_h[step] / exp_h[step];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j] * exp_h[step], 1.0E-08);
        adj_lambda = adj_lambda1 * sigBetaij - adj_lambda2;

        sigBetaij *= sigBetaij * dt;

        varxi[j] = const_varx * sigBetaij;
        muzi[j] = const_muz1 * sigBetaij + const_muz2 + adj_lambda;
        varzi[j] = const_varz * sigBetaij;
      }
    }

    /*	convolve the 3 PDE					*/

    if (is_bar[step + 1]) {
      /* do the three PDE seperatly */

      /*	first one is a regular PDE */

      if (n1) {
        num_f_pde_one_step_backward_2f_adi(
            pde1, nstepx, x, nstepz, z, ms1, me1, values_p1, mux, muz, varx,
            varz, r, values, 0, nstepx - 1, 0, nstepz - 1);
      }

      /*	second one with barrier */
      num_f_pde_one_step_backward_2f_adi_bar(
          pde2, nstepx, x, nstepz, z, ms2, me2, values_p1, mux, muz, varx, varz,
          r, values, lx, ux, lz, uz, is_up[step + 1], !is_up[step + 1]);

      /*	third one is a regular PDE */
      if (n3) {
        num_f_pde_one_step_backward_2f_adi(
            pde3, nstepx, x, nstepz, z, ms3, me3, values_p1, mux, muz, varx,
            varz, r, values, 0, nstepx - 1, 0, nstepz - 1);
      }
    } else {
      /* everything is normal */
      num_f_pde_one_step_backward_2f_adi(pde4, nstepx, x, nstepz, z, ms4, me4,
                                         values_p1, mux, muz, varx, varz, r,
                                         values, 0, nstepx - 1, 0, nstepz - 1);
    }

    /*  Apply discounting	*/

    dft = swp_f_df(date[step], date[step + 1], dom_yc);

    for (i = 0; i < nstepx; i++) {
      valuesi = values[i];

      for (j = lz; j <= uz; j++) {
        valuesij = valuesi[j];

        for (k = 0; k < nprod; k++) {
          valuesij[k] *= dft;
        }
      }
    }

    /*	Now we have to define the new gride for next step */

    if (is_bar[step]) {
      if (is_up[step]) {
        if (!is_up[step + 1]) {
          change_lim = 1;

          ux = nstepx - 1;

          while (x[ux] > bar_lvl[step] - 1.0E-08) {
            ux--;
          }

          ux++;
        } else {
          change_lim = 0;
        }

        if (bar_lvl[step] < bar_lvl[step + 1]) {
          /* decreasing barrier */

          if (ux == 0 || bar_lvl[step] >= x[ux - 1]) {
            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, ux, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we do not change the grid */
            x[ux] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (ux == index_x || lx == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              if (ux == index_x) {
                adj_index = ux;
              } else {
                adj_index = lx;
              }

              Interpolate_Payoff(values, x, adj_index, nstepx, lz, uz, nprod,
                                 bar_col, spot_fx, 1);

              x[index_x] = spot_fx;
              sigBeta[index_x][0] = save1;
              muzinit[index_x][0] = save2;
              varxinit[index_x][0] = save3;
            }

            /* we change of point */
            while (ux >= 0 && x[ux] > bar_lvl[step]) {
              ux--;
            }
            ux++;

            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, ux, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we adjust the grid */
            x[ux] = bar_lvl[step];
          }

          /* adjust precalculations */
          logx = min(max(log(x[ux]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[ux][0] = alpharhobeta5 * temp;
          varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        } else {
          /* increasing barrier */

          if (ux == nstepx - 1 || bar_lvl[step] <= x[ux + 1]) {
            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, ux, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we do not change the grid */
            x[ux] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (ux == index_x || lx == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              if (ux == index_x) {
                adj_index = ux;
              } else {
                adj_index = lx;
              }

              Interpolate_Payoff(values, x, adj_index, nstepx, lz, uz, nprod,
                                 bar_col, spot_fx, 1);

              x[index_x] = spot_fx;
              sigBeta[index_x][0] = save1;
              muzinit[index_x][0] = save2;
              varxinit[index_x][0] = save3;
            }

            /* we change of point */
            while (ux < nstepx && x[ux] < bar_lvl[step]) {
              ux++;
            }
            ux--;

            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, ux, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we adjust the grid */
            x[ux] = bar_lvl[step];
          }

          /* adjust precalculations */
          logx = min(max(log(x[ux]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[ux][0] = alpharhobeta5 * temp;
          varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        }

        if (change_lim) {
          lx = 0;
        }
      } else {
        if (is_up[step + 1]) {
          change_lim = 1;

          lx = 0;

          while (x[lx] < bar_lvl[step] + 1.0E-08) {
            lx++;
          }

          lx--;
        } else {
          change_lim = 0;
        }

        if (bar_lvl[step] < bar_lvl[step + 1]) {
          /* decreasing barrier */

          if (lx == 0 || bar_lvl[step] >= x[lx - 1]) {
            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, lx, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we do not change the grid */
            x[lx] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (lx == index_x || ux == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              if (ux == index_x) {
                adj_index = ux;
              } else {
                adj_index = lx;
              }

              Interpolate_Payoff(values, x, adj_index, nstepx, lz, uz, nprod,
                                 bar_col, spot_fx, 1);

              x[index_x] = spot_fx;
              sigBeta[index_x][0] = save1;
              muzinit[index_x][0] = save2;
              varxinit[index_x][0] = save3;
            }

            /* we change of point */
            while (lx >= 0 && x[lx] > bar_lvl[step]) {
              lx--;
            }
            lx++;

            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, lx, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we adjust the grid */
            x[lx] = bar_lvl[step];
          }

          /* adjust precalculations */
          logx = min(max(log(x[lx]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[lx][0] = alpharhobeta5 * temp;
          varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        } else {
          /* increasing barrier */

          if (lx == nstepx - 1 || bar_lvl[step] <= x[lx + 1]) {
            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, lx, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we do not change the grid */
            x[lx] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (lx == index_x || ux == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              if (ux == index_x) {
                adj_index = ux;
              } else {
                adj_index = lx;
              }

              Interpolate_Payoff(values, x, adj_index, nstepx, lz, uz, nprod,
                                 bar_col, spot_fx, 1);

              x[index_x] = spot_fx;
              sigBeta[index_x][0] = save1;
              muzinit[index_x][0] = save2;
              varxinit[index_x][0] = save3;
            }

            /* we change of point */
            while (lx < nstepx && x[lx] < bar_lvl[step]) {
              lx++;
            }
            lx--;

            /* Interpolate the payoff for the new point */
            Interpolate_Payoff(values, x, lx, nstepx, lz, uz, nprod, bar_col,
                               bar_lvl[step], 0);

            /* we adjust the grid */
            x[lx] = bar_lvl[step];
          }

          /* adjust precalculations */
          logx = min(max(log(x[lx]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[lx][0] = alpharhobeta5 * temp;
          varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        }

        if (change_lim) {
          ux = nstepx - 1;
        }
      }
    }

    /*	Eval payoff */
    if (eval_evt[step]) {
      err = payoff_func(date[step], time[step], func_parm_tab[step], today,
                        spot_fx, dom_yc, for_yc, 0, nstepx - 1, lz, uz, x,
                        nprod, values);
      if (err) {
        goto FREE_RETURN;
      }
    }

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;
  }

  if (fabs(x[index_x] - spot_fx) > 1.0E-08) {
    /* we need to interpolate */
    Interpolate_Payoff(values, x, index_x, nstepx, index_z, index_z, nprod,
                       bar_col, spot_fx, 1);
  }

  /* copy the result					*/
  for (k = 0; k < nprod; k++) {
    res[k] = values_p1[index_x][index_z][k];
  }

  t2 = clock();

  smessage("Convolution        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  if (exp_h)
    free_dvector(exp_h, 0, nstp - 1);
  if (expect_z)
    free_dvector(expect_z, 0, nstp - 1);
  if (drift_z)
    free_dvector(drift_z, 0, nstp - 1);

  if (pde1)
    num_f_pde_free_2d_adi(pde1, nstepx, nstepz, n1);
  if (pde2)
    num_f_pde_free_2d_adi(pde2, nstepx, nstepz, n2);
  if (pde3)
    num_f_pde_free_2d_adi(pde3, nstepx, nstepz, n3);
  if (pde4)
    num_f_pde_free_2d_adi(pde4, nstepx, nstepz, n4);

  if (x)
    free(x);
  if (z)
    free_dvector(z, 0, nstepz - 1);
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (mux)
    free_dmatrix(mux, 0, nstepx - 1, 0, nstepz - 1);
  if (muz)
    free_dmatrix(muz, 0, nstepx - 1, 0, nstepz - 1);
  if (muzinit)
    free_dmatrix(muzinit, 0, nstepx - 1, 0, nstepz - 1);
  if (sigBeta)
    free_dmatrix(sigBeta, 0, nstepx - 1, 0, nstepz - 1);
  if (varx)
    free_dmatrix(varx, 0, nstepx - 1, 0, nstepz - 1);
  if (varxinit)
    free_dmatrix(varxinit, 0, nstepx - 1, 0, nstepz - 1);
  if (varz)
    free_dmatrix(varz, 0, nstepx - 1, 0, nstepz - 1);
  if (varzinit)
    free_dmatrix(varzinit, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);

  return err;
}

static void Interpolate_Payoff(double ***values, double *x, int ux, int nstepx,
                               int lz, int uz, int nprod, int bar_col,
                               double barrier, int interp_bar) {
  static double d1, d2, coef1, coef2, coef3;
  static int j, k;
  static double **valuesi1, **valuesi2, **valuesi3;

  if (ux + 1 < nstepx && ux > 0) {
    valuesi1 = values[ux + 1];
    valuesi2 = values[ux];
    valuesi3 = values[ux - 1];

    /* adjust payoff */

    if (nprod != 1) {
      d1 = x[ux + 1] - x[ux];
      d2 = x[ux] - x[ux - 1];

      coef1 = (barrier * (barrier - x[ux] - x[ux - 1]) + x[ux] * x[ux - 1]) /
              (d1 * (d1 + d2));
      coef2 = (barrier * (barrier - x[ux + 1] - x[ux - 1]) +
               x[ux + 1] * x[ux - 1]) /
              (-d1 * d2);

      coef3 = 1.0 - coef1 - coef2;
    }

    for (j = lz; j <= uz; j++) {
      /* knock out payoff */
      if (interp_bar) {
        for (k = 0; k < nprod; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] +
                           coef3 * valuesi3[j][k];
        }
      } else {
        valuesi2[j][bar_col] = 0.0;

        /* other columns        , quadratic interpolation */
        for (k = 0; k < bar_col; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] +
                           coef3 * valuesi3[j][k];
        }
        for (k = bar_col + 1; k < nprod; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] +
                           coef3 * valuesi3[j][k];
        }
      }
    }
  } else if (ux == nstepx - 1) {
    coef2 = (barrier - x[ux - 1]) / (x[ux] - x[ux - 1]);
    coef3 = 1.0 - coef2;

    valuesi2 = values[ux];
    valuesi3 = values[ux - 1];

    for (j = lz; j <= uz; j++) {
      if (interp_bar) {
        for (k = 0; k < nprod; k++) {
          valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
        }
      } else {
        /* knock out payoff */
        valuesi2[j][bar_col] = 0.0;

        /* other columns        , quadratic interpolation */
        for (k = 0; k < bar_col; k++) {
          valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
        }
        for (k = bar_col + 1; k < nprod; k++) {
          valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
        }
      }
    }
  } else {
    coef1 = (barrier - x[ux]) / (x[ux + 1] - x[ux]);
    coef2 = 1.0 - coef1;

    valuesi1 = values[ux + 1];
    valuesi2 = values[ux];

    for (j = lz; j <= uz; j++) {
      if (interp_bar) {
        for (k = 0; k < nprod; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
        }
      } else {
        /* knock out payoff */
        valuesi2[j][bar_col] = 0.0;

        /* other columns        , quadratic interpolation */
        for (k = 0; k < bar_col; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
        }
        for (k = bar_col + 1; k < nprod; k++) {
          valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
        }
      }
    }
  }
}

Err FxSabr_KOOption(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    double strike, int is_call, /* 1 Call        , 0: Put */
    double *bar_lvl, int is_up, /* 1 Up        , 0: Down */
    int is_cvx,                 /* 1 use 1/Fx        , 0: use Fx */

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Result */
    double *res) {
  Err err = NULL;

  long today;
  int i, j, k, step;
  int index_x, index_z, nstepx, nstepz;
  double fux, flx, dt, dft;
  double alpha2, alpharho, alpha2rho2, alpharhobeta5;
  double temp;
  double std1, std3, spotBeta;
  double z0, bar_ex;

  double *exp_h = NULL, *expect_z = NULL, *expect_z2 = NULL, *drift_z = NULL,
         *x = NULL, *z = NULL, ***values = NULL, ***values_p1 = NULL,
         ***values_temp = NULL, **mux = NULL, **muz = NULL, **muzinit = NULL,
         **sigBeta = NULL, **varx = NULL, **varz = NULL, **varxinit = NULL,
         **varzinit = NULL, **r = NULL;

  double *varxi, *varzi, *muzi, **valuesi, *valuesij;

  double exp_z, const_sig, sigBetaij, const1, const2;
  double const_sigi, const_varx, const_muz1, const_muz2, const_varz;
  double adj_lambda, adj_lambda1, adj_lambda2;

  int lx, ux, lz, uz;

  double sig, varxf;

  long nprod;
  double df_ratio, strike2, pay;
  long mat_date;
  int do_bar;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = NULL;

  double logx, logmax, logmin, save1, save2, save3, coef;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  if (fabs(beta - 1.0) < 1.0E-08) {
    beta = 1.0 - 1.0E-08;
  }

  t1 = clock();

  /* Constant	calculations */
  today = (long)(date[0] + 1.0e-06);

  spotBeta = pow(spot_fx, beta - 1.0);

  alpha2 = alpha * alpha;
  alpharho = alpha * rho;
  alpharhobeta5 = 0.5 * alpharho * beta;
  alpha2rho2 = alpha * alpha * (1.0 - rho * rho);
  const_sig = alpharho / (1.0 - beta);

  nprod = 1;

  z0 = rho * alpha / (1.0 - beta) * (1.0 / spotBeta - 1) - sig0;

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepfx / 2)) * 2 + 1;
  nstepz = ((int)(nstepvol / 2)) * 2 + 1;

  /*	we want at least three points in each directions
   */
  if (nstepx < 3) {
    nstepx = 3;
  }
  if (nstepz < 3) {
    nstepz = 3;
  }

  /*	Memory allocations */

  x = calloc(nstepx + 3, sizeof(double));
  z = dvector(0, nstepz - 1);

  exp_h = dvector(0, nstp - 1);
  expect_z = dvector(0, nstp - 1);
  drift_z = dvector(0, nstp - 1);

  if (!exp_h || !expect_z || !drift_z || !x || !z) {
    err = "Memory allocation error (1) in FxSabr_adi";
    goto FREE_RETURN;
  }

  /*	Calculate the standard deviation of X and the expectaion of Z	*/

  const1 = -0.5 * alpharho * beta * sig0 * sig0 * spotBeta;
  const2 = -alpharho * (1.0 / spotBeta - 1.0) / (1.0 - beta);

  temp = 0.0;
  exp_z = z0;
  varxf = 0.0;

  for (i = 0; i < nstp - 1; i++) {
    dt = (time[i + 1] - time[i]) / 2.0;
    sig = exp(temp);
    varxf += sig * sig * 2.0 * dt;

    exp_z += const1 * sig * exp(alpha2 * time[i]) / (drift[i] + alpha2) *
                 (exp((drift[i] + alpha2) * dt) - 1.0) -
             const2 / sig * (exp(-drift[i] * dt) - 1.0);

    expect_z[i] = exp_z;

    temp += drift[i] * dt;
    sig = exp(temp);

    exp_h[i] = sig;

    drift_z[i] = (const1 * sig * exp(alpha2 * (time[i] + dt)) +
                  const2 / sig * drift[i]) *
                 2.0 * dt;

    exp_z += const1 * sig * exp(alpha2 * (time[i] + dt)) / (drift[i] + alpha2) *
                 (exp((drift[i] + alpha2) * dt) - 1.0) -
             const2 / sig * (exp(-drift[i] * dt) - 1.0);

    temp += drift[i] * dt;
  }

  std1 = sig0 * sqrt(varxf) * spotBeta;
  std3 = sig0 * sqrt(alpha2rho2 * time[nstp - 1]);

  /*	Then discretise space in the orthogonal system x / z */

  /*	Lognormal bounds of X */
  flx = max(exp(log(spot_fx) - std1 * std1 / 2.0 - NSTD_FX * std1), -1.0);
  fux = exp(log(flx) + 2.0 * NSTD_FX * std1);

  /*	Bound for the calculation of max and min drift and var */
  logmax = log(spot_fx) - std1 * std1 / 2.0 + floorstd * std1;
  logmin = log(spot_fx) - std1 * std1 / 2.0 - floorstd * std1;

  /*	Discretisation of x and z */

  do_bar = 1;

  mat_date = add_unit((long)(date[nstp - 1] + 1.0E-08), 2, SRT_BDAY,
                      MODIFIED_SUCCEEDING);

  df_ratio =
      swp_f_df(today, mat_date, for_yc) / swp_f_df(today, mat_date, dom_yc);
  if (is_cvx) {
    df_ratio = 1.0 / df_ratio;
  }

  strike2 = strike / df_ratio;

  if (is_up) {
    bar_ex = bar_lvl[0];
    for (i = 1; i < nstp; i++) {
      if (bar_lvl[i] > bar_ex) {
        bar_ex = bar_lvl[i];
      }
    }

    if (bar_ex > fux) {
      bar_ex = fux;
      do_bar = 0;
    }

    if (strike2 > 1.0E-08) {
      /* add the strike */
      if (is_cvx) {
        disc_linleft_logright_strike(x, nstepx, spot_fx, flx, bar_ex, std1,
                                     NSTD_FX, 1.0 / strike2, &index_x);
      } else {
        disc_linleft_logright_strike(x, nstepx, spot_fx, flx, bar_ex, std1,
                                     NSTD_FX, strike2, &index_x);
      }
    } else {
      disc_linleft_logright(x, nstepx, spot_fx, flx, bar_ex, std1, NSTD_FX,
                            &index_x);
    }

  } else {
    bar_ex = bar_lvl[0];
    for (i = 1; i < nstp; i++) {
      if (bar_lvl[i] < bar_ex) {
        bar_ex = bar_lvl[i];
      }
    }

    if (bar_ex < flx) {
      bar_ex = flx;
      do_bar = 0;
    }

    if (strike2 > 1.0E-08) {
      /* add the strike */
      if (is_cvx) {
        disc_linleft_logright_strike(x, nstepx, spot_fx, bar_ex, fux, std1,
                                     NSTD_FX, 1.0 / strike2, &index_x);
      } else {
        disc_linleft_logright_strike(x, nstepx, spot_fx, bar_ex, fux, std1,
                                     NSTD_FX, strike2, &index_x);
      }
    } else {
      disc_linleft_logright(x, nstepx, spot_fx, bar_ex, fux, std1, NSTD_FX,
                            &index_x);
    }
  }

  disc_normal_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3,
                     NSTD_VOL, &index_z);

  /*	Precalculate variances and expectations				*/

  /* first allocate memory */
  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  mux = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muzinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  sigBeta = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varx = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varxinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varzinit = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  if (!values || !values_p1 || !mux || !muz || !muzinit || !sigBeta || !varx ||
      !varxinit || !varz || !varzinit || !r) {
    err = "Memory allocation error (2) in FxSabr_adi";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstepx; i++) {
    logx = min(max(log(x[i]), logmin), logmax);
    temp = exp((beta - 1.0) * logx);

    sigBeta[i][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
    muzinit[i][0] = alpharhobeta5 * temp;
    varxinit[i][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

    for (j = 0; j < nstepz; j++) {
      mux[i][j] = 0.0;
      r[i][j] = 0.0;
    }
  }

  save1 = sigBeta[index_x][0];
  save2 = muzinit[index_x][0];
  save3 = varxinit[index_x][0];

  /* Up and Out Call */
  if (is_call && is_up) {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] > strike2 && x[i] < bar_lvl[nstp - 1] - 1.0E-08) {
          pay = (1.0 / x[i] - strike2) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] > strike2 && x[i] < bar_lvl[nstp - 1] - 1.0E-08) {
          pay = (x[i] - strike2) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    }
  }
  /* Up and Out Put */
  else if (!is_call && is_up) {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] < strike2 && x[i] < bar_lvl[nstp - 1] - 1.0E-08) {
          pay = (strike2 - 1.0 / x[i]) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] < strike2 && x[i] < bar_lvl[nstp - 1] - 1.0E-08) {
          pay = (strike2 - x[i]) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    }
  }
  /* Down and Out Call */
  else if (is_call && !is_up) {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] > strike2 && x[i] > bar_lvl[nstp - 1] + 1.0E-08) {
          pay = (1.0 / x[i] - strike2) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] > strike2 && x[i] > bar_lvl[nstp - 1] + 1.0E-08) {
          pay = (x[i] - strike2) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    }
  } else
  /* Down and Out Put */
  {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] < strike2 && x[i] > bar_lvl[nstp - 1] + 1.0E-08) {
          pay = (strike2 - 1.0 / x[i]) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] < strike2 && x[i] > bar_lvl[nstp - 1] + 1.0E-08) {
          pay = (strike2 - x[i]) * df_ratio;

          for (j = 0; j < nstepz; j++) {
            values_p1[i][j][0] = pay;
          }
        }
      }
    }
  }

  /*	Initialize the CNPDE_TEMP_2D		*/

  pde = &pdestr;

  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

  if (!pde) {
    err = "Memory allocation error (2) in FxSabr_Adi";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  /* Initialize the barrier */

  if (do_bar) {
    if (is_up) {
      while (x[ux] > bar_lvl[nstp - 1] - 1.0E-08) {
        ux--;
      }

      ux++;

      /* update precalculations */
      x[ux] = bar_lvl[nstp - 1];

      logx = min(max(log(x[ux]), logmin), logmax);
      temp = exp((beta - 1.0) * logx);
      sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
      muzinit[ux][0] = alpharhobeta5 * temp;
      varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

      /* update payoff */

      for (j = 0; j < nstepz - 1; j++) {
        values[ux][j][0] = 0.0;
      }
    } else {
      while (x[lx] < bar_lvl[nstp - 1] + 1.0E-08) {
        lx++;
      }

      lx--;

      /* update precalculations */
      x[lx] = bar_lvl[nstp - 1];

      logx = min(max(log(x[lx]), logmin), logmax);
      temp = exp((beta - 1.0) * logx);
      sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
      muzinit[lx][0] = alpharhobeta5 * temp;
      varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */

      /* update payoff */

      for (j = 0; j < nstepz - 1; j++) {
        values[lx][j][0] = 0.0;
      }
    }
  }

  /*	now do the backward pde					*/

  for (step = nstp - 2; step >= 0; step--) {

    dt = time[step + 1] - time[step];

    adj_lambda1 = lambda / (sig0 * exp_h[step]) * dt;
    adj_lambda2 = lambda * dt;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];

      const_sigi = sigBeta[i][0] - expect_z[step] * exp_h[step];
      const_varx = varxinit[i][0];
      const_muz1 = -muzinit[i][0] / exp_h[step];
      const_muz2 =
          -sigBeta[i][0] / exp_h[step] * drift[step] * dt - drift_z[step];
      const_varz = alpha2rho2 / exp_h[step] / exp_h[step];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j] * exp_h[step], 1.0E-08);
        adj_lambda = adj_lambda1 * sigBetaij - adj_lambda2;

        sigBetaij *= sigBetaij * dt;

        varxi[j] = const_varx * sigBetaij;
        muzi[j] = const_muz1 * sigBetaij + const_muz2 + adj_lambda;
        varzi[j] = const_varz * sigBetaij;
      }
    }

    /*	convolve							*/

    if (do_bar) {
      num_f_pde_one_step_backward_2f_adi_bar(
          pde, nstepx, x, nstepz, z, 0, nprod - 1, values_p1, mux, muz, varx,
          varz, r, values, lx, ux, lz, uz, is_up, !is_up);
    } else {
      num_f_pde_one_step_backward_2f_adi(pde, nstepx, x, nstepz, z, 0,
                                         nprod - 1, values_p1, mux, muz, varx,
                                         varz, r, values, lx, ux, lz, uz);
    }

    /*  Apply discounting	*/

    dft = swp_f_df(date[step], date[step + 1], dom_yc);

    for (i = lx; i <= ux; i++) {
      valuesi = values[i];

      for (j = lz; j <= uz; j++) {
        valuesij = valuesi[j];

        for (k = 0; k < nprod; k++) {
          valuesij[k] *= dft;
        }
      }
    }

    /*	Now we have to define the new gride for next step */

    if (do_bar) {
      if (is_up) {
        if (bar_lvl[step] < bar_lvl[step + 1]) {
          /* decreasing barrier */

          if (ux == 0 || bar_lvl[step] >= x[ux - 1]) {
            /* we do not change the grid */
            x[ux] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (ux == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              coef = 1.0 / (x[ux - 1] - x[ux]) / (x[ux - 1] - x[ux + 1]);
              coef *= (spot_fx * (spot_fx - (x[ux] + x[ux + 1])) +
                       x[ux] * x[ux + 1]);

              for (j = lz; j <= uz; j++) {
                values[ux][j][0] = coef * values[ux - 1][j][0];
              }

              x[ux] = spot_fx;
              sigBeta[ux][0] = save1;
              muzinit[ux][0] = save2;
              varxinit[ux][0] = save3;
            }

            /* we change of point */
            while (ux >= 0 && x[ux] > bar_lvl[step]) {
              ux--;
            }
            ux++;

            /* we adjust the grid */
            x[ux] = bar_lvl[step];

            /* adjust payoff */
            for (j = lz; j <= uz; j++) {
              values[ux][j][0] = 0.0;
            }
          }

          /* adjust precalculations */
          logx = min(max(log(x[ux]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[ux][0] = alpharhobeta5 * temp;
          varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        } else {
          /* increasing barrier */

          if (ux == nstepx - 1 || bar_lvl[step] <= x[ux + 1]) {
            /* we do not change the grid */
            x[ux] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (ux == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              coef = 1.0 / (x[ux - 1] - x[ux]) / (x[ux - 1] - x[ux + 1]);
              coef *= (spot_fx * (spot_fx - (x[ux] + x[ux + 1])) +
                       x[ux] * x[ux + 1]);

              for (j = lz; j <= uz; j++) {
                values[ux][j][0] = coef * values[ux - 1][j][0];
              }

              x[ux] = spot_fx;
              sigBeta[ux][0] = save1;
              muzinit[ux][0] = save2;
              varxinit[ux][0] = save3;
            }

            /* we change of point */
            while (ux < nstepx && x[ux] < bar_lvl[step]) {
              ux++;
            }
            ux--;

            /* we adjust the grid */
            x[ux] = bar_lvl[step];

            /* adjust payoff */
            for (j = lz; j <= uz; j++) {
              values[ux][j][0] = 0.0;
            }
          }

          /* adjust precalculations */
          logx = min(max(log(x[ux]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[ux][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[ux][0] = alpharhobeta5 * temp;
          varxinit[ux][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        }
      } else {
        if (bar_lvl[step] < bar_lvl[step + 1]) {
          /* decreasing barrier */

          if (lx == 0 || bar_lvl[step] >= x[lx - 1]) {
            /* we do not change the grid */
            x[lx] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (lx == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              coef = 1.0 / (x[lx + 1] - x[lx]) / (x[lx + 1] - x[lx - 1]);
              coef *= (spot_fx * (spot_fx - (x[lx] + x[lx - 1])) +
                       x[lx] * x[lx - 1]);

              for (j = lz; j <= uz; j++) {
                values[lx][j][0] = coef * values[lx + 1][j][0];
              }

              x[lx] = spot_fx;
              sigBeta[lx][0] = save1;
              muzinit[lx][0] = save2;
              varxinit[lx][0] = save3;
            }

            /* we change of point */
            while (lx >= 0 && x[lx] > bar_lvl[step]) {
              lx--;
            }
            lx++;

            /* we adjust the grid */
            x[lx] = bar_lvl[step];

            /* adjust payoff */
            for (j = lz; j <= uz; j++) {
              values[lx][j][0] = 0.0;
            }
          }

          /* adjust precalculations */
          logx = min(max(log(x[lx]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[lx][0] = alpharhobeta5 * temp;
          varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        } else {
          /* increasing barrier */

          if (lx == nstepx - 1 || bar_lvl[step] <= x[lx + 1]) {
            /* we do not change the grid */
            x[lx] = bar_lvl[step];
          } else {
            /* special case where the barrier point was the one of the fwd ! */
            if (lx == index_x) {
              /* we adjust the payoff by quadratic interpolation */

              coef = 1.0 / (x[lx + 1] - x[lx]) / (x[lx + 1] - x[lx - 1]);
              coef *= (spot_fx * (spot_fx - (x[lx] + x[lx - 1])) +
                       x[lx] * x[lx - 1]);

              for (j = lz; j <= uz; j++) {
                values[lx][j][0] = coef * values[lx + 1][j][0];
              }

              x[lx] = spot_fx;
              sigBeta[lx][0] = save1;
              muzinit[lx][0] = save2;
              varxinit[lx][0] = save3;
            }

            /* we change of point */
            while (lx < nstepx && x[lx] < bar_lvl[step]) {
              lx++;
            }
            lx--;

            /* we adjust the grid */
            x[lx] = bar_lvl[step];

            /* adjust payoff */
            for (j = lz; j <= uz; j++) {
              values[lx][j][0] = 0.0;
            }
          }

          /* adjust precalculations */
          logx = min(max(log(x[lx]), logmin), logmax);
          temp = exp((beta - 1.0) * logx);

          sigBeta[lx][0] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
          muzinit[lx][0] = alpharhobeta5 * temp;
          varxinit[lx][0] = exp(2.0 * beta * logx); /* z^(2.0 * beta) */
        }
      }
    }

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;
  }

  /* copy the result					*/
  for (k = 0; k < nprod; k++) {
    res[k] = values_p1[index_x][index_z][k];
  }

  t2 = clock();

  smessage("Convolution        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  if (exp_h)
    free_dvector(exp_h, 0, nstp - 1);
  if (expect_z)
    free_dvector(expect_z, 0, nstp - 1);
  if (drift_z)
    free_dvector(drift_z, 0, nstp - 1);

  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  if (x)
    free(x);
  if (z)
    free_dvector(z, 0, nstepz - 1);
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (mux)
    free_dmatrix(mux, 0, nstepx - 1, 0, nstepz - 1);
  if (muz)
    free_dmatrix(muz, 0, nstepx - 1, 0, nstepz - 1);
  if (muzinit)
    free_dmatrix(muzinit, 0, nstepx - 1, 0, nstepz - 1);
  if (sigBeta)
    free_dmatrix(sigBeta, 0, nstepx - 1, 0, nstepz - 1);
  if (varx)
    free_dmatrix(varx, 0, nstepx - 1, 0, nstepz - 1);
  if (varxinit)
    free_dmatrix(varxinit, 0, nstepx - 1, 0, nstepz - 1);
  if (varz)
    free_dmatrix(varz, 0, nstepx - 1, 0, nstepz - 1);
  if (varzinit)
    free_dmatrix(varzinit, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);

  return err;
}