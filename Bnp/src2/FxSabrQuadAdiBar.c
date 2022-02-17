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
#include "FxSabrQuadAdi.h"
#include "FxSabrSLAdi.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrgeneric.h"

#define NSTD_FX 5.0
#define NSTD_VOL 5.0

#define MINNODE 25
#define MAXNODE2 50

Err FxSabrQuad_KOOption(
    /*	Time data		*/
    int nstp, double *time, double *date,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double a, double b, double c,
    double rho, double lambda,

    double floorstd,

    /*	Product data */
    long settlmt_date, double strike, int is_call, /* 1 Call  , 0: Put */
    int is_american, int is_cvx,                   /* 1 use 1/Fx  , 0: use Fx */
    int is_digital, /* 1: digital payoff  , 0  , regular option payoff */
    double *bar_lvl_up, double *bar_lvl_down, double rebate_up,
    double rebate_down,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,
    int eod_fix_flag, /*	EOD Fixing Flag 0: I  , 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I  , 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I  , 1: E */

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

  double std1, std3, stdln;

  double *expect_z = NULL, *drift_z = NULL, *x = NULL, *z = NULL,
         ***values = NULL, ***values_p1 = NULL, ***values_temp = NULL,
         **mux = NULL, **muz = NULL, *sigBeta = NULL, **varx = NULL,
         **varz = NULL, *varxinit = NULL, *muzinit = NULL, **r = NULL;

  double *varxi, *varzi, *muzi, **valuesi, *valuesij;

  double const_sigi, const_varz, const_muz, const_muz1, const_muz2, const_varxi,
      sigBetaij;
  double alpharho_a, alpharhoa5;

  int lx, ux, lz, uz;

  double bar_ex_up, bar_ex_down;

  long nprod;
  double df_ratio, strike2, pay;
  int do_bar_up, do_bar_down;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = NULL;

  double dxu, dxd, dzu, dzd;
  double logx, logmax, logmin, save1, save2, save3, coef;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  t1 = clock();

  /* Constant	calculations */

  alpharho_a = alpha * rho;
  alpharhoa5 = 0.5 * alpha * rho;
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
    err = "Memory allocation error (1) in FxSabrQuad_KOOption";
    goto FREE_RETURN;
  }

  /*	Calculate the standard deviation of X and the expectaion of Z	*/
  err = FxSabrQuadPrecalculations(nstp, time, sig0, drift, alpha, a, b, c, rho,
                                  lambda, spot_fx, expect_z, drift_z, &std1);

  /*	Boundaries */
  std1 /= sqrt(time[nstp - 1]);

  stdln = op_sabrgen(spot_fx, spot_fx, time[nstp - 1], std1, alpha, a, b, c,
                     rho, vol_quadra);

  stdln *= sqrt(time[nstp - 1]);

  flx = spot_fx * exp(-0.5 * stdln * stdln - NSTD_FX * stdln);
  fux = spot_fx * exp(-0.5 * stdln * stdln + NSTD_FX * stdln);

  /*	Bound for the calculation of max and min drift and var */

  logmin = spot_fx * exp(-0.5 * stdln * stdln - floorstd * stdln);
  logmax = spot_fx * exp(-0.5 * stdln * stdln + floorstd * stdln);

  std1 /= a;

  std3 = sig0 * sqrt(const_varz * time[nstp - 1]);

  /*	Discretisation of x and z */

  df_ratio = swp_f_df(today, settlmt_date, for_yc) /
             swp_f_df(today, settlmt_date, dom_yc);
  if (is_cvx) {
    df_ratio = 1.0 / df_ratio;
  }

  strike2 = strike / df_ratio;

  dft = swp_f_df(date[nstp - 1], settlmt_date, dom_yc);
  df_ratio *= dft;

  if (fabs(strike2) < 1.0E-10) {
    strike2 = 1.0E-10;
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
  muzinit = dvector(0, nstepx - 1);

  if (!x || !z || !values || !values_p1 || !mux || !muz || !sigBeta || !varx ||
      !varxinit || !muzinit || !varz || !r) {
    err = "Memory allocation error (2) in FxSabrQuad_KOOption";
    goto FREE_RETURN;
  }

  /* add the strike */
  if (is_cvx) {
    /*
    disc_SL_strike(x  , nstepx  , spot_fx  , bar_ex_down  , bar_ex_up  , 0.9999
    , 0.00001  , stdln  , NSTD_FX  , 1.0 / strike2  , 0  , &index_x);
    */

    disc_SL_center_linleft_logright(x, nstepx, spot_fx, bar_ex_down, bar_ex_up,
                                    0.9999, 0.00001, stdln, NSTD_FX, &index_x);
  } else {
    disc_SL_strike(x, nstepx, spot_fx, bar_ex_down, bar_ex_up, 0.9999, 0.00001,
                   stdln, NSTD_FX, strike2, 0, &index_x);
    /*
    disc_linleft_logright_center_strike(x  , nstepx  , spot_fx  , bar_ex_down  ,
    bar_ex_up  , stdln  , NSTD_FX  , strike2  , &index_x);
    */
  }

  disc_linleft_linright_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3,
                               std3, NSTD_VOL, &index_z);

  /*	Precalculate variances and expectations				*/
  for (i = 0; i < nstepx; i++) {
    logx = min(max(x[i], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[i] = logx * logx;

    logx = min(max(x[i], logmin), logmax);
    muzinit[i] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[i] = alpharho_a * vol_quadra(logx, a, b, c, 3);
  }

  save1 = sigBeta[index_x];
  save2 = muzinit[index_x];
  save3 = varxinit[index_x];

  /* Update Payoff */

  /* Down / Up and Out Call */
  if (is_call) {
    if (is_cvx) {
      for (i = 0; i < nstepx; i++) {
        if (1.0 / x[i] > strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
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
        if (x[i] > strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          if (is_digital) {
            pay = 1.0;
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
        if (1.0 / x[i] < strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
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
        if (x[i] < strike2 && (x[i] < bar_lvl_up[nstp - 1] - 1.0E-08) &&
            (x[i] > bar_lvl_down[nstp - 1] + 1.0E-08)) {
          if (is_digital) {
            pay = 1.0;
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
    err = "Memory allocation error (2) in FxSabrQuad_Adi";
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
    logx = min(max(x[ux], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[ux] = logx * logx;

    logx = min(max(x[ux], logmin), logmax);
    muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);

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
    logx = min(max(x[lx], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[lx] = logx * logx;

    logx = min(max(x[lx], logmin), logmax);
    muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);

    /* adjust payoff */
    for (j = lz; j <= uz; j++) {
      values[lx][j][0] = rebate_down;
    }
  }

  /*	now do the backward pde					*/

  for (step = nstp - 2; step >= 0; step--) {
    if (step == 0) {
      if (eod_ex_flag) {
        is_american = 0;
      }
      if (eod_fix_flag) {
        do_bar_up = 0;
        do_bar_down = 0;
      }
      if (eod_pay_flag) {
        rebate_up = 0.0;
        rebate_down = 0.0;
      }
    }

    dt = time[step + 1] - time[step];

    const_muz1 = (drift[step] - lambda) * dt;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];
      const_muz = muzinit[i];

      const_sigi = sigBeta[i] - expect_z[step];
      const_varxi = varxinit[i];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j], 1.0E-08);

        const_muz2 = const_muz1 * sigBetaij;

        sigBetaij *= sigBetaij * dt;

        varxi[j] = const_varxi * sigBetaij;
        muzi[j] = -const_muz * sigBetaij - const_muz2 - drift_z[step];
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
                pay = 1.0;
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
                pay = 1.0;
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
            muzinit[ux] = save2;
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
        logx = min(max(x[ux], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[ux] = logx * logx;

        logx = min(max(x[ux], logmin), logmax);
        muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[ux] = save2;
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
        logx = min(max(x[ux], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[ux] = logx * logx;

        logx = min(max(x[ux], logmin), logmax);
        muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[lx] = save2;
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
        logx = min(max(x[lx], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[lx] = logx * logx;

        logx = min(max(x[lx], logmin), logmax);
        muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[lx] = save2;
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
        logx = min(max(x[lx], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[lx] = logx * logx;

        logx = min(max(x[lx], logmin), logmax);
        muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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

    dzu = -(z[index_z + 1] - z[index_z]) / sig0;
    dzd = -(z[index_z] - z[index_z - 1]) / sig0;

    dt = time[1] - time[0];

    greeks[0] = (values_p1[index_x + 1][index_z][0] -
                 values_p1[index_x - 1][index_z][0]) /
                (dxu + dxd);
    greeks[1] =
        2.0 *
        (values_p1[index_x + 1][index_z][0] * dxd +
         values_p1[index_x - 1][index_z][0] * dxu - (dxd + dxu) * res[0]) /
        (dxu * dxd * (dxu + dxd));
    greeks[2] = (values[index_x][index_z][0] - values_p1[index_x][index_z][0]) /
                dt / DAYS_IN_YEAR;

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
  if (muzinit)
    free_dvector(muzinit, 0, nstepx - 1);

  /* Allocation 3 */
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  return err;
}

Err KO_Payoff(
    double event_date,

    long today, char *dom_yc, char *for_yc,

    double notional, long settlmt_date, double strike,
    int is_call,    /* 1 Call  , 0: Put */
    int is_cvx,     /* 1 use 1/Fx  , 0: use Fx */
    int is_digital, /* 1: digital payoff  , 0  , regular option payoff */
    double bar_lvl_up, double bar_lvl_down, double rebate_up,
    double rebate_down,

    double *x, long start_x, long end_x, long start_z, long end_z,
    double ***values_p1) {
  int i, j;
  double pay, df_ratio, dft, strike2;
  Err err = NULL;

  df_ratio = swp_f_df(today, settlmt_date, for_yc) /
             swp_f_df(today, settlmt_date, dom_yc);
  if (is_cvx) {
    df_ratio = 1.0 / df_ratio;
  }

  strike2 = strike / df_ratio;

  dft = swp_f_df(event_date, settlmt_date, dom_yc);
  df_ratio *= dft;

  if (fabs(strike2) < 1.0E-10) {
    strike2 = 1.0E-10;
  }

  /* Down / Up and Out Call */
  if (is_call) {
    if (is_cvx) {
      for (i = start_x; i <= end_x; i++) {
        if (1.0 / x[i] > strike2 && (x[i] < bar_lvl_up - 1.0E-08) &&
            (x[i] > bar_lvl_down + 1.0E-08)) {
          pay = (1.0 / x[i] - strike2) * df_ratio;
        } else if (x[i] > bar_lvl_up - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_up + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        pay *= notional;

        for (j = start_z; j <= end_z; j++) {
          values_p1[i][j][0] += pay;
        }
      }
    } else {
      for (i = start_x; i <= end_x; i++) {
        if (x[i] > strike2 && (x[i] < bar_lvl_up - 1.0E-08) &&
            (x[i] > bar_lvl_down + 1.0E-08)) {
          if (is_digital) {
            pay = 1.0;
          } else {
            pay = (x[i] - strike2) * df_ratio;
          }
        } else if (x[i] > bar_lvl_up - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        pay *= notional;

        for (j = start_z; j <= end_z; j++) {
          values_p1[i][j][0] += pay;
        }
      }
    }
  }
  /* Down / Up and Out Put */
  else {
    if (is_cvx) {
      for (i = start_x; i <= end_x; i++) {
        if (1.0 / x[i] < strike2 && (x[i] < bar_lvl_up - 1.0E-08) &&
            (x[i] > bar_lvl_down + 1.0E-08)) {
          pay = (strike2 - 1.0 / x[i]) * df_ratio;
        } else if (x[i] > bar_lvl_up - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        pay *= notional;

        for (j = start_z; j <= end_z; j++) {
          values_p1[i][j][0] += pay;
        }
      }
    } else {
      for (i = start_x; i <= end_x; i++) {
        if (x[i] < strike2 && (x[i] < bar_lvl_up - 1.0E-08) &&
            (x[i] > bar_lvl_down + 1.0E-08)) {
          if (is_digital) {
            pay = 1.0;
          } else {
            pay = (strike2 - x[i]) * df_ratio;
          }
        } else if (x[i] > bar_lvl_up - 1.0E-08) {
          pay = rebate_up;
        } else if (x[i] < bar_lvl_down + 1.0E-08) {
          pay = rebate_down;
        } else {
          pay = 0.0;
        }

        pay *= notional;

        for (j = start_z; j <= end_z; j++) {
          values_p1[i][j][0] += pay;
        }
      }
    }
  }

  return err;
}

Err FxSabrQuad_KO_MultiOption(
    /*	Time data		*/
    int nstp, double *time, double *date, int *eval_evt,

    /*	Discretisation	*/
    int nstepfx, int nstepvol,

    /*	Model data		*/
    double sig0, double *drift, double alpha, double a, double b, double c,
    double rho, double lambda,

    double floorstd,

    /*	Product data */
    int nb_product, double *notional, long *exercise_date, long *settlmt_date,
    double *strike, int *is_call, /* 1 Call  , 0: Put */
    int *is_cvx,                  /* 1 use 1/Fx  , 0: use Fx */
    int *is_digital, /* 1: digital payoff  , 0  , regular option payoff */
    double *bar_lvl_up, double *bar_lvl_down, double *rebate_up,
    double *rebate_down,

    /*	Market data */
    double spot_fx, /*	The cash one */
    char *dom_yc, char *for_yc,
    int eod_fix_flag, /*	EOD Fixing Flag 0: I  , 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I  , 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I  , 1: E */

    /*	Result */
    double *res,

    /* Additional informations */
    int calc_greeks,
    double *greeks) /* array 6 * nprod containing delta  , gamma  , theta  ,
                       vega  , volga and vanna */
{
  Err err = NULL;

  long today;
  int i, j, k, step, prod, m;
  int index_x, index_z, nstepx, nstepz;
  double fux, flx, dt, dft;

  double std1, std3, stdln;

  double *expect_z = NULL, *drift_z = NULL, *x = NULL, *z = NULL,
         ***values = NULL, ***values_p1 = NULL, ***values_temp = NULL,
         **mux = NULL, **muz = NULL, *sigBeta = NULL, **varx = NULL,
         **varz = NULL, *varxinit = NULL, *muzinit = NULL, **r = NULL;

  double *varxi, *varzi, *muzi, **valuesi, *valuesij;

  double const_sigi, const_varz, const_muz, const_muz1, const_muz2, const_varxi,
      sigBetaij;
  double alpharho_a, alpharhoa5;

  int lx, ux, lz, uz;

  double bar_ex_up, bar_ex_down;

  long nprod;
  int do_bar_up, do_bar_down;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = NULL;

  double dxu, dxd, dzu, dzd;
  double logx, logmax, logmin, save1, save2, save3, coef;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  t1 = clock();

  /* Constant	calculations */

  alpharho_a = alpha * rho;
  alpharhoa5 = 0.5 * alpha * rho;
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
    err = "Memory allocation error (1) in FxSabrQuad_KOOption";
    goto FREE_RETURN;
  }

  /*	Calculate the standard deviation of X and the expectaion of Z	*/
  err = FxSabrQuadPrecalculations(nstp, time, sig0, drift, alpha, a, b, c, rho,
                                  lambda, spot_fx, expect_z, drift_z, &std1);

  /*	Boundaries */
  std1 /= sqrt(time[nstp - 1]);

  stdln = op_sabrgen(spot_fx, spot_fx, time[nstp - 1], std1, alpha, a, b, c,
                     rho, vol_quadra);

  stdln *= sqrt(time[nstp - 1]);

  flx = spot_fx * exp(-0.5 * stdln * stdln - NSTD_FX * stdln);
  fux = spot_fx * exp(-0.5 * stdln * stdln + NSTD_FX * stdln);

  /*	Bound for the calculation of max and min drift and var */

  logmin = spot_fx * exp(-0.5 * stdln * stdln - floorstd * stdln);
  logmax = spot_fx * exp(-0.5 * stdln * stdln + floorstd * stdln);

  std1 /= a;

  std3 = sig0 * sqrt(const_varz * time[nstp - 1]);

  /*	Discretisation of x and z */
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
  muzinit = dvector(0, nstepx - 1);

  if (!x || !z || !values || !values_p1 || !mux || !muz || !sigBeta || !varx ||
      !varxinit || !muzinit || !varz || !r) {
    err = "Memory allocation error (2) in FxSabrQuad_KOOption";
    goto FREE_RETURN;
  }

  if (is_cvx[nb_product - 1]) {
    disc_SL_center_linleft_logright(x, nstepx, spot_fx, bar_ex_down, bar_ex_up,
                                    0.9999, 0.00001, stdln, NSTD_FX, &index_x);
  } else {
    disc_SL_strike(x, nstepx, spot_fx, bar_ex_down, bar_ex_up, 0.9999, 0.00001,
                   stdln, NSTD_FX, bar_ex_up * 1000, 0, &index_x);
  }

  disc_linleft_linright_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3,
                               std3, NSTD_VOL, &index_z);

  /*	Precalculate variances and expectations				*/
  for (i = 0; i < nstepx; i++) {
    logx = min(max(x[i], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[i] = logx * logx;

    logx = min(max(x[i], logmin), logmax);
    muzinit[i] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[i] = alpharho_a * vol_quadra(logx, a, b, c, 3);
  }

  save1 = sigBeta[index_x];
  save2 = muzinit[index_x];
  save3 = varxinit[index_x];

  /* Update payoff */

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  step = nstp - 1;
  prod = nb_product - 1;

  for (i = lx; i <= ux; i++) {
    for (j = lz; j <= uz; j++) {
      values_p1[i][j][0] = 0.0;
    }
  }

  for (m = 0; m < eval_evt[step]; m++) {
    err = KO_Payoff(date[step], today, dom_yc, for_yc, notional[prod],
                    settlmt_date[prod], strike[prod], is_call[prod],
                    is_cvx[prod], is_digital[prod], bar_lvl_up[step],
                    bar_lvl_down[step], rebate_up[prod], rebate_down[prod], x,
                    lx, ux, lz, uz, values_p1);

    if (err)
      goto FREE_RETURN;

    prod--;
  }

  /*	Initialize the CNPDE_TEMP_2D		*/

  pde = &pdestr;

  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

  if (!pde) {
    err = "Memory allocation error (2) in FxSabrQuad_Adi";
    goto FREE_RETURN;
  }

  /* Initialize the barrier */

  if (do_bar_up) {
    ux = nstepx - 1;
    while (x[ux] > bar_lvl_up[nstp - 1] - 1.0E-08) {
      ux--;
    }

    ux++;

    /* update precalculations */
    logx = min(max(x[ux], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[ux] = logx * logx;

    logx = min(max(x[ux], logmin), logmax);
    muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);

    /* adjust payoff */
    for (j = lz; j <= uz; j++) {
      values[ux][j][0] = rebate_up[nb_product - 1];
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
    logx = min(max(x[lx], logmin), logmax);
    logx = logx * (a * logx + b) + c;
    varxinit[lx] = logx * logx;

    logx = min(max(x[lx], logmin), logmax);
    muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
    sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);

    /* adjust payoff */
    for (j = lz; j <= uz; j++) {
      values[lx][j][0] = rebate_down[nb_product - 1];
    }
  }

  /*	now do the backward pde					*/

  for (step = nstp - 2; step >= 0; step--) {
    if (step == 0) {
      if (eod_ex_flag) {
      }
      if (eod_fix_flag) {
        do_bar_up = 0;
        do_bar_down = 0;
      }
      if (eod_pay_flag) {
        rebate_up[0] = 0.0;
        rebate_down[0] = 0.0;
      }
    }

    dt = time[step + 1] - time[step];

    const_muz1 = (drift[step] - lambda) * dt;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];
      const_muz = muzinit[i];

      const_sigi = sigBeta[i] - expect_z[step];
      const_varxi = varxinit[i];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j], 1.0E-08);

        const_muz2 = const_muz1 * sigBetaij;

        sigBetaij *= sigBetaij * dt;

        varxi[j] = const_varxi * sigBetaij;
        muzi[j] = -const_muz * sigBetaij - const_muz2 - drift_z[step];
        varzi[j] = const_varz * sigBetaij;
      }
    }

    /*	convolve							*/

    dft = swp_f_df(date[step], date[step + 1], dom_yc);

    num_f_pde_one_step_backward_2f_adi_bar(
        pde, nstepx, x, nstepz, z, 0, nprod - 1, values_p1, mux, muz, varx,
        varz, r, values, lx, ux, lz, uz, do_bar_up, do_bar_down);

    /*  Apply discounting	*/

    for (i = lx + do_bar_down; i <= ux - do_bar_up; i++) {
      valuesi = values[i];

      for (j = lz; j <= uz; j++) {
        valuesij = valuesi[j];

        for (k = 0; k < nprod; k++) {
          valuesij[k] *= dft;
        }
      }
    }

    if (eval_evt[step]) {
      for (m = 0; m < eval_evt[step]; m++) {
        err = KO_Payoff(date[step], today, dom_yc, for_yc, notional[prod],
                        settlmt_date[prod], strike[prod], is_call[prod],
                        is_cvx[prod], is_digital[prod], bar_lvl_up[step],
                        bar_lvl_down[step], rebate_up[prod], rebate_down[prod],
                        x, lx, ux, lz, uz, values);

        if (err)
          goto FREE_RETURN;

        prod--;
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
            muzinit[ux] = save2;
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
            values[ux][j][0] = rebate_up[prod + 1];
          }
        }

        /* adjust precalculations */
        logx = min(max(x[ux], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[ux] = logx * logx;

        logx = min(max(x[ux], logmin), logmax);
        muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[ux] = save2;
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
            values[ux][j][0] = rebate_up[prod + 1];
          }
        }

        /* adjust precalculations */
        logx = min(max(x[ux], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[ux] = logx * logx;

        logx = min(max(x[ux], logmin), logmax);
        muzinit[ux] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[ux] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[lx] = save2;
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
            values[lx][j][0] = rebate_down[prod + 1];
          }
        }

        /* adjust precalculations */
        logx = min(max(x[lx], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[lx] = logx * logx;

        logx = min(max(x[lx], logmin), logmax);
        muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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
            muzinit[lx] = save2;
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
            values[lx][j][0] = rebate_down[prod + 1];
          }
        }

        /* adjust precalculations */
        logx = min(max(x[lx], logmin), logmax);
        logx = logx * (a * logx + b) + c;
        varxinit[lx] = logx * logx;

        logx = min(max(x[lx], logmin), logmax);
        muzinit[lx] = alpharhoa5 * (2.0 * a * logx + b);
        sigBeta[lx] = alpharho_a * vol_quadra(logx, a, b, c, 3);
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

    dzu = -(z[index_z + 1] - z[index_z]) / sig0;
    dzd = -(z[index_z] - z[index_z - 1]) / sig0;

    dt = time[1] - time[0];

    greeks[0] = (values_p1[index_x + 1][index_z][0] -
                 values_p1[index_x - 1][index_z][0]) /
                (dxu + dxd);
    greeks[1] =
        2.0 *
        (values_p1[index_x + 1][index_z][0] * dxd +
         values_p1[index_x - 1][index_z][0] * dxu - (dxd + dxu) * res[0]) /
        (dxu * dxd * (dxu + dxd));
    greeks[2] = (values[index_x][index_z][0] - values_p1[index_x][index_z][0]) /
                dt / DAYS_IN_YEAR;

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
  if (muzinit)
    free_dvector(muzinit, 0, nstepx - 1);

  /* Allocation 3 */
  if (pde)
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

  return err;
}