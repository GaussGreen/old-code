/* ==========================================================================
   FILE_NAME:	LGM2FSabr.c

   PURPOSE:		ADI implementation of the Sabr Fx model.
                                Discretisation is ADI.

   DATE:		03/16/01

   AUTHOR:		L.C.
   ========================================================================== */

#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "FxSabrGrfn.h"
#include "FxSabrSLAdi.h"
#include "math.h"
#include "opfnctns.h"

#define NSTD_FX 5.0
#define NSTD_VOL 5.0

#define MINNODE 25
#define MAXNODE2 50

Err FxSabrSL_KOOption(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double beta, double rho,
    double lambda,

    double floorstd,

    /*	Product data */
    double strike, int is_call,  /* 1 Call  , 0: Put */
    int is_american, int is_cvx, /* 1 use 1/Fx  , 0: use Fx */
    int is_digital, /* 1: digital payoff  , 0  , regular option payoff */
    double *bar_lvl_up, double *bar_lvl_down, double rebate_up,
    double rebate_down,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,

    /*	Result */
    double *res,

    /* Additional informations */
    int calc_greeks,
    double *greeks) /* array 6 * nprod containing delta  , gamma  , theta  ,
                       vega  , volga and vanna */
{
  Err err = NULL;

  long today;
  int i, j, k, step;
  int index_x, index_z, nstepx, nstepz;
  double fux, flx, dt, dft;

  double std1, std3;

  double *expect_z = NULL, *drift_z = NULL, *x = NULL, *z = NULL,
         ***values = NULL, ***values_p1 = NULL, ***values_temp = NULL,
         **mux = NULL, **muz = NULL, *sigBeta = NULL, **varx = NULL,
         **varz = NULL, *varxinit = NULL, **r = NULL;

  double *varxi, *varzi, *muzi, **valuesi, *valuesij;

  double const_sigi, const_varz, const_muz1, const_muz2, const_varxi, sigBetaij;
  double alpharho_a, alpharhoa5;

  int lx, ux, lz, uz;

  double bar_ex_up, bar_ex_down;

  long nprod;
  double df_ratio, strike2, pay;
  long mat_date;
  int do_bar_up, do_bar_down;

  double a, b;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = NULL;

  double dxu, dxd, dzu, dzd;
  double logx, logmax, logmin, save1, save3, coef;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  if (fabs(beta - 1.0) < 1.0E-08) {
    beta = 1.0 - 1.0E-08;
  }

  t1 = clock();

  /* Constant	calculations */
  a = beta * exp((beta - 1.0) * log(spot_fx));
  b = (1.0 - beta) * exp(beta * log(spot_fx));

  alpharho_a = alpha * rho / a;
  alpharhoa5 = 0.5 * alpha * rho * a;
  const_varz = alpha * alpha * (1.0 - rho * rho);
  today = (long)(date[0] + 1.0e-06);

  nprod = 1;

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
  expect_z = dvector(0, nstp - 1);
  drift_z = dvector(0, nstp - 1);

  if (!expect_z || !drift_z) {
    err = "Memory allocation error (1) in FxSabrSL_KOOption";
    goto FREE_RETURN;
  }

  /*	Calculate the standard deviation of X and the expectaion of Z	*/
  err = FxSabrSLPrecalculations(nstp, time, sig0, drift, alpha, a, b, rho,
                                lambda, spot_fx, expect_z, drift_z, &std1);

  /*	Boundaries */
  if (fabs(lambda) < 1.0E-10) {
    std1 *= a * exp(alpha * sqrt(time[nstp - 1]));
  } else {
    std1 *= a * exp(alpha * sqrt((1.0 - exp(-2.0 * lambda * time[nstp - 1])) /
                                 (2.0 * lambda)));
  }

  flx = ((a * spot_fx + b) * exp(-0.5 * std1 * std1 - NSTD_FX * std1) - b) / a;
  fux = ((a * spot_fx + b) * exp(-0.5 * std1 * std1 + NSTD_FX * std1) - b) / a;

  /*	Bound for the calculation of max and min drift and var */

  logmin =
      ((a * spot_fx + b) * exp(-0.5 * std1 * std1 - floorstd * std1) - b) / a;
  logmax =
      ((a * spot_fx + b) * exp(-0.5 * std1 * std1 + floorstd * std1) - b) / a;

  std1 /= a;

  std3 = sig0 * sqrt(const_varz * time[nstp - 1]);

  /*	Discretisation of x and z */

  mat_date = add_unit((long)(date[nstp - 1] + 1.0E-08), 2, SRT_BDAY,
                      MODIFIED_SUCCEEDING);

  df_ratio =
      swp_f_df(today, mat_date, for_yc) / swp_f_df(today, mat_date, dom_yc);
  if (is_cvx) {
    df_ratio = 1.0 / df_ratio;
  }

  strike2 = strike / df_ratio;

  if (fabs(strike2) < 1.0E-10) {
    strike2 = 1.0E-10;
  }

  if (is_digital) {
    df_ratio = swp_f_df(date[nstp - 1], mat_date, dom_yc);
  } else {
    df_ratio *= swp_f_df(date[nstp - 1], mat_date, dom_yc);
  }

  bar_ex_up = bar_lvl_up[0];
  bar_ex_down = bar_lvl_down[0];

  for (i = 1; i < nstp; i++) {
    if (bar_lvl_up[i] > bar_ex_up) {
      bar_ex_up = bar_lvl_up[i];
    }

    if (bar_lvl_down[i] < bar_ex_down) {
      bar_ex_down = bar_lvl_down[i];
    }
  }

  if (bar_ex_up > fux) {
    bar_ex_up = fux;
    do_bar_up = 0;
  } else {
    do_bar_up = 1;
  }

  if (bar_ex_down < flx) {
    bar_ex_down = flx;
    do_bar_down = 0;
  } else {
    do_bar_down = 1;
  }

  /*	Discretisation of space in the orthogonal system x / z */

  /* first allocate memory */

  x = calloc(nstepx, sizeof(double));
  z = dvector(0, nstepz - 1);

  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);

  mux = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varx = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  sigBeta = dvector(0, nstepx - 1);
  varxinit = dvector(0, nstepx - 1);

  if (!x || !z || !values || !values_p1 || !mux || !muz || !sigBeta || !varx ||
      !varxinit || !varz || !r) {
    err = "Memory allocation error (2) in FxSabrSL_KOOption";
    goto FREE_RETURN;
  }

  /* add the strike */
  if (is_cvx) {
    disc_SL_strike(x, nstepx, spot_fx, bar_ex_down, bar_ex_up, a, b, std1,
                   NSTD_FX, 1.0 / strike2, 0, &index_x);
  } else {
    disc_SL_strike(x, nstepx, spot_fx, bar_ex_down, bar_ex_up, a, b, std1,
                   NSTD_FX, strike2, 0, &index_x);
  }

  disc_linleft_linright_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3,
                               std3, NSTD_VOL, &index_z);

  /*	Precalculate variances and expectations				*/
  for (i = 0; i < nstepx; i++) {
    logx = min(max(x[i], logmin), logmax);
    logx = a * logx + b;

    varxinit[i] = logx * logx;

    logx = log(logx);
    sigBeta[i] = alpharho_a * logx;
  }

  save1 = sigBeta[index_x];
  save3 = varxinit[index_x];

  /* Update Payoff */

  /* Down / Up and Out Call */
  if (is_call) {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] >= strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          pay = (1.0 / x[i] - strike2) * df_ratio;
        } else if (x[i] > bar_lvl_up[nstp - 1] - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_up[nstp - 1] + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        for (j = 0; j < nstepz; j++) {
          values_p1[i][j][0] = pay;
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] >= strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          if (is_digital) {
            if (fabs(x[i] - strike2) > 1.0E-08) {
              pay = df_ratio;
            } else {
              pay = 0.5 * df_ratio;
            }
          } else {
            pay = (x[i] - strike2) * df_ratio;
          }
        } else if (x[i] > bar_lvl_up[nstp - 1] - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down[nstp - 1] + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        for (j = 0; j < nstepz; j++) {
          values_p1[i][j][0] = pay;
        }
      }
    }
  }
  /* Down / Up and Out Put */
  else {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] <= strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          pay = (strike2 - 1.0 / x[i]) * df_ratio;
        } else if (x[i] > bar_lvl_up[nstp - 1] - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down[nstp - 1] + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        for (j = 0; j < nstepz; j++) {
          values_p1[i][j][0] = pay;
        }
      }
    } else {
      for (i = 0; i < nstepx; i++) {
        if (x[i] <= strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          if (is_digital) {
            if (fabs(x[i] - strike2) > 1.0E-08) {
              pay = df_ratio;
            } else {
              pay = 0.5 * df_ratio;
            }
          } else {
            pay = (strike2 - x[i]) * df_ratio;
          }
        } else if (x[i] > bar_lvl_up[nstp - 1] - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down[nstp - 1] + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        for (j = 0; j < nstepz; j++) {
          values_p1[i][j][0] = pay;
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

  /* Initialize the barrier */

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  if (do_bar_up) {
    ux = nstepx - 1;
    while (x[ux] > bar_lvl_up[nstp - 1] - 1.0E-08) {
      ux--;
    }

    ux++;

    /* update precalculations */
    x[ux] = bar_lvl_up[nstp - 1];
    logx = a * min(max(x[ux], logmin), logmax) + b;
    varxinit[ux] = logx * logx;
    logx = log(logx);
    sigBeta[ux] = alpharho_a * logx;

    /* adjust payoff */
    for (j = lz; j <= uz; j++) {
      values[ux][j][0] = rebate_up;
    }
  }

  if (do_bar_down) {
    lz = 0;
    while (x[lx] < bar_lvl_down[nstp - 1] + 1.0E-08) {
      lx++;
    }

    lx--;

    /* update precalculations */
    x[lx] = bar_lvl_down[nstp - 1];
    logx = a * min(max(x[lx], logmin), logmax) + b;
    varxinit[lx] = logx * logx;
    logx = log(logx);
    sigBeta[lx] = alpharho_a * logx;

    /* adjust payoff */
    for (j = lz; j <= uz; j++) {
      values[lx][j][0] = rebate_down;
    }
  }

  /*	now do the backward pde					*/

  for (step = nstp - 2; step >= 0; step--) {

    dt = time[step + 1] - time[step];

    const_muz1 = (drift[step] - lambda) * dt;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];

      const_sigi = sigBeta[i] - expect_z[step];
      const_varxi = varxinit[i];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j], 1.0E-08);

        const_muz2 = const_muz1 * sigBetaij;

        sigBetaij *= sigBetaij * dt;

        varxi[j] = const_varxi * sigBetaij;
        muzi[j] = -alpharhoa5 * sigBetaij - const_muz2 - drift_z[step];
        varzi[j] = const_varz * sigBetaij;
      }
    }

    /*	convolve							*/

    dft = swp_f_df(date[step], date[step + 1], dom_yc);

    num_f_pde_one_step_backward_2f_adi_bar(
        pde, nstepx, x, nstepz, z, 0, nprod - 1, values_p1, mux, muz, varx,
        varz, r, values, lx, ux, lz, uz, do_bar_up, do_bar_down);

    /*  Apply discounting	*/

    if (is_american) {
      df_ratio = swp_f_df(today, date[step], for_yc) /
                 swp_f_df(today, date[step], dom_yc);
      if (is_cvx) {
        df_ratio = 1.0 / df_ratio;
      }

      strike2 = strike / df_ratio;

      if (fabs(strike2) < 1.0E-10) {
        strike2 = 1.0E-10;
      }

      /* American option */
      for (i = lx; i <= ux; i++) {
        valuesi = values[i];

        if (is_call) {
          if (is_cvx) {
            if (1.0 / x[i] > strike2 && (x[i] < bar_lvl_up[step] - 1.0E-08) &&
                (x[i] > bar_lvl_down[step] + 1.0E-08)) {
              pay = (1.0 / x[i] - strike2) * df_ratio;
            } else if (x[i] > bar_lvl_up[step] - 1.0E-08) {
              pay = rebate_up;
            } else if (x[i] < bar_lvl_down[step] + 1.0E-08) {
              pay = rebate_down;
            } else {
              pay = 0.0;
            }
          } else {
            if (x[i] > strike2 && (x[i] < bar_lvl_up[step] - 1.0E-08) &&
                (x[i] > bar_lvl_down[step] + 1.0E-08)) {
              if (is_digital) {
                if (fabs(x[i] - strike2) > 1.0E-08) {
                  pay = 1.0;
                } else {
                  pay = 0.5;
                }
              } else {
                pay = (x[i] - strike2) * df_ratio;
              }
            } else if (x[i] > bar_lvl_up[step] - 1.0E-08) {
              pay = rebate_up;
            } else if (x[i] < bar_lvl_down[step] + 1.0E-08) {
              pay = rebate_down;
            } else {
              pay = 0.0;
            }
          }
        } else {
          if (is_cvx) {
            if (1.0 / x[i] < strike2 && (x[i] < bar_lvl_up[step] - 1.0E-08) &&
                (x[i] > bar_lvl_down[step] + 1.0E-08)) {
              pay = (strike2 - 1.0 / x[i]) * df_ratio;
            } else if (x[i] > bar_lvl_up[step] - 1.0E-08) {
              pay = rebate_up;
            } else if (x[i] < bar_lvl_down[step] + 1.0E-08) {
              pay = rebate_down;
            } else {
              pay = 0.0;
            }
          } else {
            if (x[i] < strike2 && (x[i] < bar_lvl_up[step] - 1.0E-08) &&
                (x[i] > bar_lvl_down[step] + 1.0E-08)) {
              if (is_digital) {
                if (fabs(x[i] - strike2) > 1.0E-08) {
                  pay = 1.0;
                } else {
                  pay = 0.5;
                }
              } else {
                pay = (strike2 - x[i]) * df_ratio;
              }
            } else if (x[i] > bar_lvl_up[step] - 1.0E-08) {
              pay = rebate_up;
            } else if (x[i] < bar_lvl_down[step] + 1.0E-08) {
              pay = rebate_down;
            } else {
              pay = 0.0;
            }
          }
        }

        for (j = lz; j <= uz; j++) {
          valuesij = valuesi[j];

          for (k = 0; k < nprod; k++) {
            valuesij[k] = max(valuesij[k] * dft, pay);
          }
        }
      }
    } else {
      for (i = lx + do_bar_down; i <= ux - do_bar_up; i++) {
        valuesi = values[i];

        for (j = lz; j <= uz; j++) {
          valuesij = valuesi[j];

          for (k = 0; k < nprod; k++) {
            valuesij[k] *= dft;
          }
        }
      }
    }

    /*	Now we have to define the new gride for next step */

    if (do_bar_up) {
      if (bar_lvl_up[step] < bar_lvl_up[step + 1]) {
        /* decreasing barrier */

        if (ux == 0 || bar_lvl_up[step] >= x[ux - 1]) {
          /* we do not change the grid */
          x[ux] = bar_lvl_up[step];
        } else {
          /* special case where the barrier point was the one of the fwd ! */
          if (ux == index_x) {
            /* we adjust the payoff by quadratic interpolation */

            coef = 1.0 / (x[ux - 1] - x[ux]) / (x[ux - 1] - x[ux + 1]);
            coef *=
                (spot_fx * (spot_fx - (x[ux] + x[ux + 1])) + x[ux] * x[ux + 1]);

            for (j = lz; j <= uz; j++) {
              values[ux][j][0] = coef * values[ux - 1][j][0];
            }

            x[ux] = spot_fx;
            sigBeta[ux] = save1;
            varxinit[ux] = save3;
          }

          /* we change of point */
          while (ux >= 0 && x[ux] > bar_lvl_up[step]) {
            ux--;
          }
          ux++;

          /* we adjust the grid */
          x[ux] = bar_lvl_up[step];

          /* adjust payoff */
          for (j = lz; j <= uz; j++) {
            values[ux][j][0] = rebate_up;
          }
        }

        /* adjust precalculations */
        logx = a * min(max(x[ux], logmin), logmax) + b;
        varxinit[ux] = logx * logx;
        logx = log(logx);
        sigBeta[ux] = alpharho_a * logx;
      } else {
        /* increasing barrier */

        if (ux == nstepx - 1 || bar_lvl_up[step] <= x[ux + 1]) {
          /* we do not change the grid */
          x[ux] = bar_lvl_up[step];
        } else {
          /* special case where the barrier point was the one of the fwd ! */
          if (ux == index_x) {
            /* we adjust the payoff by quadratic interpolation */

            coef = 1.0 / (x[ux - 1] - x[ux]) / (x[ux - 1] - x[ux + 1]);
            coef *=
                (spot_fx * (spot_fx - (x[ux] + x[ux + 1])) + x[ux] * x[ux + 1]);

            for (j = lz; j <= uz; j++) {
              values[ux][j][0] = coef * values[ux - 1][j][0];
            }

            x[ux] = spot_fx;
            sigBeta[ux] = save1;
            varxinit[ux] = save3;
          }

          /* we change of point */
          while (ux < nstepx && x[ux] < bar_lvl_up[step]) {
            ux++;
          }
          ux--;

          /* we adjust the grid */
          x[ux] = bar_lvl_up[step];

          /* adjust payoff */
          for (j = lz; j <= uz; j++) {
            values[ux][j][0] = rebate_up;
          }
        }

        /* adjust precalculations */
        logx = a * min(max(x[ux], logmin), logmax) + b;
        varxinit[ux] = logx * logx;
        logx = log(logx);
        sigBeta[ux] = alpharho_a * logx;
      }
    }

    if (do_bar_down) {
      if (bar_lvl_down[step] < bar_lvl_down[step + 1]) {
        /* decreasing barrier */

        if (lx == 0 || bar_lvl_down[step] >= x[lx - 1]) {
          /* we do not change the grid */
          x[lx] = bar_lvl_down[step];
        } else {
          /* special case where the barrier point was the one of the fwd ! */
          if (lx == index_x) {
            /* we adjust the payoff by quadratic interpolation */

            coef = 1.0 / (x[lx + 1] - x[lx]) / (x[lx + 1] - x[lx - 1]);
            coef *=
                (spot_fx * (spot_fx - (x[lx] + x[lx - 1])) + x[lx] * x[lx - 1]);

            for (j = lz; j <= uz; j++) {
              values[lx][j][0] = coef * values[lx + 1][j][0];
            }

            x[lx] = spot_fx;
            sigBeta[lx] = save1;
            varxinit[lx] = save3;
          }

          /* we change of point */
          while (lx >= 0 && x[lx] > bar_lvl_down[step]) {
            lx--;
          }
          lx++;

          /* we adjust the grid */
          x[lx] = bar_lvl_down[step];

          /* adjust payoff */
          for (j = lz; j <= uz; j++) {
            values[lx][j][0] = rebate_down;
          }
        }

        /* adjust precalculations */
        logx = a * min(max(x[lx], logmin), logmax) + b;
        varxinit[lx] = logx * logx;
        logx = log(logx);
        sigBeta[lx] = alpharho_a * logx;
      } else {
        /* increasing barrier */

        if (lx == nstepx - 1 || bar_lvl_down[step] <= x[lx + 1]) {
          /* we do not change the grid */
          x[lx] = bar_lvl_down[step];
        } else {
          /* special case where the barrier point was the one of the fwd ! */
          if (lx == index_x) {
            /* we adjust the payoff by quadratic interpolation */

            coef = 1.0 / (x[lx + 1] - x[lx]) / (x[lx + 1] - x[lx - 1]);
            coef *=
                (spot_fx * (spot_fx - (x[lx] + x[lx - 1])) + x[lx] * x[lx - 1]);

            for (j = lz; j <= uz; j++) {
              values[lx][j][0] = coef * values[lx + 1][j][0];
            }

            x[lx] = spot_fx;
            sigBeta[lx] = save1;
            varxinit[lx] = save3;
          }

          /* we change of point */
          while (lx < nstepx && x[lx] < bar_lvl_down[step]) {
            lx++;
          }
          lx--;

          /* we adjust the grid */
          x[lx] = bar_lvl_down[step];

          /* adjust payoff */
          for (j = lz; j <= uz; j++) {
            values[lx][j][0] = rebate_down;
          }
        }

        /* adjust precalculations */
        logx = a * min(max(x[lx], logmin), logmax) + b;
        varxinit[lx] = logx * logx;
        logx = log(logx);
        sigBeta[lx] = alpharho_a * logx;
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

  /* compute greeks if needed */
  if (calc_greeks) {
    dxu = x[index_x + 1] - x[index_x];
    dxd = x[index_x] - x[index_x - 1];

    dzu = -(z[index_z + 1] - z[index_z]);
    dzd = -(z[index_z] - z[index_z - 1]);

    dt = time[1] - time[0];

    greeks[0] = (values_p1[index_x + 1][index_z][0] -
                 values_p1[index_x - 1][index_z][0]) /
                (dxu + dxd);
    greeks[1] =
        2.0 *
        (values_p1[index_x + 1][index_z][0] * dxd +
         values_p1[index_x - 1][index_z][0] * dxu - (dxd + dxu) * res[0]) /
        (dxu * dxd * (dxu + dxd));
    greeks[2] =
        -(values[index_x][index_z][0] - values_p1[index_x][index_z][0]) / dt /
        DAYS_IN_YEAR;

    greeks[3] = (values_p1[index_x][index_z + 1][0] -
                 values_p1[index_x][index_z - 1][0]) /
                (dzu + dzd);
    greeks[4] =
        2.0 *
        (values_p1[index_x][index_z + 1][0] * dzd +
         values_p1[index_x][index_z - 1][0] * dzu - (dzu + dzd) * res[0]) /
        (dzu * dzd * (dzu + dzd));

    greeks[5] = ((values_p1[index_x + 1][index_z + 1][0] -
                  values_p1[index_x - 1][index_z + 1][0]) /
                     (dxu + dxd) -
                 (values_p1[index_x + 1][index_z - 1][0] -
                  values_p1[index_x - 1][index_z - 1][0]) /
                     (dxu + dxd)) /
                (dzu + dzd);
  }

  t2 = clock();

  smessage("Convolution  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /* Allocation 1 */
  if (expect_z)
    free_dvector(expect_z, 0, nstp - 1);
  if (drift_z)
    free_dvector(drift_z, 0, nstp - 1);

  if (x)
    free(x);
  if (z)
    free_dvector(z, 0, nstepz - 1);

  /* Allocation 2 */
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);

  if (mux)
    free_dmatrix(mux, 0, nstepx - 1, 0, nstepz - 1);
  if (muz)
    free_dmatrix(muz, 0, nstepx - 1, 0, nstepz - 1);
  if (varx)
    free_dmatrix(varx, 0, nstepx - 1, 0, nstepz - 1);
  if (varz)
    free_dmatrix(varz, 0, nstepx - 1, 0, nstepz - 1);
  if (r)
    free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);

  if (varxinit)
    free_dvector(varxinit, 0, nstepx - 1);
  if (sigBeta)
    free_dvector(sigBeta, 0, nstepx - 1);

  /* Allocation 3 */
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  return err;
}