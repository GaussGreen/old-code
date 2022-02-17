/* ==========================================================================
   FILE_NAME:	LGM2FSabr.cxx

   PURPOSE:		ADI implementation of the Sabr Fx model.
                                Discretisation is ADI.

   DATE:		03/16/01

   AUTHOR:		L.C.
   ========================================================================== */

#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "FxSabrGrfn.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrcalib.h"

#define NSTD_FX 5.0
#define NSTD_VOL 5.0

#define MINNODE 25
#define MAXNODE2 50

#define EPSALPHA 0.02
#define EPSRHO 0.05

#define EPSALPHAMIN 0.005
#define EPSRHOMIN 0.005

#define NBITERMAX 50

#define DEFAULT_PERC 0.005
#define LIMIT_DOWN_P 0.0005

Err op_sabrSL_adi(
    double forward, double *strike, int nb_strike, double maturity, double disc,
    int call_put_var, /*	0:	(F - K)+
                                        1:	(K - F)+
                                        2:	F^2	*/
    double sigma_beta, double alpha, double beta, double rho, double lambda,
    double floorstd, int nstp, int nstepx, int nstepz, double *res,
    /* Additional informations */
    int calc_greeks,
    double **greeks, /* array 6 * nprod containing delta        , gamma
                        , theta        , vega        , volga and vanna */

    /* For calibration purpose */
    int calc_at_point, int column, double target, double *vol,
    double *res_at_point) {
  Err err = NULL;

  int i, j, k, step;
  int index_x, index_z;
  double fux, flx, dt, t;
  double std1, std3;

  double *x = NULL, *z = NULL, ***values = NULL, ***values_p1 = NULL,
         ***values_temp = NULL, **mux = NULL, **muz = NULL, **varx = NULL,
         **varz = NULL, *varxinit = NULL, *sigBeta = NULL, **r = NULL;

  double *varxi, *varzi, *muzi;

  double expect_z, drift_z, z0;
  double const_sigi, const_varz, const_varxi, sigBetaij;
  double alpharho_a, alpharhoa5, alpha2, lambda2, mean, lambdadt, sigma_beta2;
  double const_expect1, const_expect2, const_muz, pay;

  int lx, ux, lz, uz;

  double a, b;

  double dxu, dxd, dzu, dzd, vega, coef, is_up, newz;

  clock_t t1, t2;

  CNPDE_TEMP_2D_ADI pdestr, *pde = NULL;

  double logx, logmax, logmin;

  /*	For stability reasons */

  if (alpha < 1.0e-05) {
    alpha = 1.0E-05;
  }

  if (fabs(beta - 1.0) < 1.0E-08) {
    beta = 1.0 - 1.0E-08;
  }

  t1 = clock();

  /* Constant	calculations */

  a = beta * exp((beta - 1.0) * log(forward));
  b = (1.0 - beta) * exp(beta * log(forward));

  /*	Precalculate dt */
  dt = maturity / (nstp - 1);

  alpharho_a = alpha * rho / a;
  alpharhoa5 = 0.5 * alpha * rho * a * dt;
  const_varz = alpha * alpha * (1.0 - rho * rho) * dt;
  z0 = alpharho_a * log(a * forward + b) - sigma_beta;
  lambda2 = 2.0 * lambda;
  alpha2 = alpha * alpha;
  sigma_beta2 = sigma_beta * sigma_beta;
  lambdadt = lambda * dt;

  mean = (lambda2 - alpha2);
  if (fabs(mean) < 1.0E-10) {
    mean = 1.0E-10;
  }

  const_expect1 = rho * alpha * a * sigma_beta2 * lambda2 / mean;

  const_expect2 = rho * alpha * a * sigma_beta2 * alpha2 / mean / mean;

  /*	nstep has to be a odd nuber			*/
  nstepx = ((int)(nstepx / 2)) * 2 + 1;
  nstepz = ((int)(nstepz / 2)) * 2 + 1;

  /*	we want at least three points in each directions */
  if (nstepx < 3) {
    nstepx = 3;
  }
  if (nstepz < 3) {
    nstepz = 3;
  }

  /*	Memory allocations */

  x = calloc(nstepx, sizeof(double));
  z = dvector(0, nstepz - 1);

  if (!x || !z) {
    err = "Memory allocation error (1) in op_sabrSL_adi";
    goto FREE_RETURN;
  }

  /*	Boundaries */

  if (fabs(lambda) < 1.0E-10) {
    std1 = a * sigma_beta * exp(alpha * sqrt(maturity)) * sqrt(maturity);
  } else {
    std1 = a * sigma_beta *
           exp(alpha *
               sqrt((1.0 - exp(-2.0 * lambda * maturity)) / (2.0 * lambda))) *
           sqrt(maturity);
  }

  flx = ((a * forward + b) * exp(-0.5 * std1 * std1 - NSTD_FX * std1) - b) / a;
  fux = ((a * forward + b) * exp(-0.5 * std1 * std1 + NSTD_FX * std1) - b) / a;

  /*	Bound for the calculation of max and min drift and var */

  logmin =
      ((a * forward + b) * exp(-0.5 * std1 * std1 - floorstd * std1) - b) / a;
  logmax =
      ((a * forward + b) * exp(-0.5 * std1 * std1 + floorstd * std1) - b) / a;

  std3 = sigma_beta * alpha * sqrt((1.0 - rho * rho) * maturity);

  /*	Discretisation of space in the orthogonal system x / z */
  if (nb_strike > 1) {
    if (floorstd > 5) {
      disc_linleft_logright_center(x, nstepx, a * forward + b, a * flx + b,
                                   a * fux + b, std1, NSTD_FX, &index_x);
    } else {
      disc_SL_center(x, nstepx, forward, flx, fux, a, b, std1 / a, NSTD_FX,
                     &index_x);
    }
  } else {
    if (floorstd > 5) {
      disc_linleft_logright_center_strike(x, nstepx, a * forward + b,
                                          a * flx + b, a * fux + b, std1,
                                          NSTD_FX, a * strike[0] + b, &index_x);
    } else {
      disc_SL_strike(x, nstepx, forward, flx, fux, a, b, std1 / a, NSTD_FX,
                     strike[0], 1, &index_x);
    }
  }

  if (floorstd > 5) {
    /* adjust */
    for (i = 0; i < nstepx; i++) {
      x[i] = (x[i] - b) / a;
    }
  }

  /*
  disc_normal_center(z        , nstepz        , 0        , -NSTD_VOL * std3 ,
  NSTD_VOL * std3        , std3        , NSTD_VOL        , &index_z);
  */

  disc_linleft_linright_center(z, nstepz, 0.0, -NSTD_VOL * std3,
                               NSTD_VOL * std3, std3, NSTD_VOL, &index_z);

  /* first allocate memory */

  values = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nb_strike - 1);
  values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nb_strike - 1);

  mux = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  muz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varx = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  varz = dmatrix(0, nstepx - 1, 0, nstepz - 1);
  r = dmatrix(0, nstepx - 1, 0, nstepz - 1);

  sigBeta = dvector(0, nstepx - 1);
  varxinit = dvector(0, nstepx - 1);

  if (!x || !z || !values || !values_p1 || !mux || !muz || !sigBeta || !varx ||
      !varxinit || !varz || !r) {
    err = "Memory allocation error (2) in op_sabrSL_adi";
    goto FREE_RETURN;
  }

  /*	Precalculate variances and expectations				*/
  for (i = 0; i < nstepx; i++) {
    logx = min(max(x[i], logmin), logmax);
    logx = a * logx + b;
    varxinit[i] = logx * logx * dt;

    logx = log(logx);
    sigBeta[i] = alpharho_a * logx;
  }

  /*	Final payoff */
  for (k = 0; k < nb_strike; k++) {
    for (i = 0; i < nstepx; i++) {
      if (call_put_var == 0) {
        pay = max(0, x[i] - strike[k]);
      } else if (call_put_var == 1) {
        pay = max(0, strike[k] - x[i]);
      } else if (call_put_var == 2) {
        pay = x[i] * x[i];
      }

      for (j = 0; j < nstepz; j++) {
        values_p1[i][j][k] = pay;
      }
    }
  }

  /*	Initialize the CNPDE_TEMP_2D		*/

  pde = &pdestr;

  num_f_pde_init_2d_adi(pde, nstepx, nstepz, nb_strike);

  if (!pde) {
    err = "Memory allocation error (3) in op_sabrSL_adi";
    goto FREE_RETURN;
  }

  lx = 0;
  ux = nstepx - 1;
  lz = 0;
  uz = nstepz - 1;

  t = maturity + dt / 2.0;

  /*	now do the backward pde					*/
  for (step = nstp - 2; step >= 0; step--) {
    t -= dt;

    drift_z = exp(-mean * t);
    expect_z = z0 - const_expect1 * t + const_expect2 * (1.0 - drift_z);
    drift_z = sigma_beta2 * (lambda2 - alpha2 * drift_z) / mean;

    for (i = lx; i <= ux; i++) {
      varxi = varx[i];
      varzi = varz[i];
      muzi = muz[i];

      const_sigi = sigBeta[i] - expect_z;
      const_varxi = varxinit[i];

      for (j = lz; j <= uz; j++) {
        /* reconstruction of the vol */
        sigBetaij = max(const_sigi - z[j], 1.0E-08);

        const_muz = lambdadt * (sigBetaij - sigma_beta);
        sigBetaij *= sigBetaij;

        varxi[j] = const_varxi * sigBetaij;
        muzi[j] = -alpharhoa5 * (sigBetaij - drift_z) + const_muz;
        varzi[j] = const_varz * sigBetaij;
      }
    }

    /*	convolve							*/
    num_f_pde_one_step_backward_2f_adi(pde, nstepx, x, nstepz, z, 0,
                                       nb_strike - 1, values_p1, mux, muz, varx,
                                       varz, r, values, lx, ux, lz, uz);

    values_temp = values_p1;
    values_p1 = values;
    values = values_temp;
  }

  /* copy the result					*/
  for (k = 0; k < nb_strike; k++) {
    res[k] = values_p1[index_x][index_z][k];
  }

  /* compute greeks if needed */
  if (calc_greeks) {
    dxu = x[index_x + 1] - x[index_x];
    dxd = x[index_x] - x[index_x - 1];

    dzu = -(z[index_z + 1] - z[index_z]);
    dzd = -(z[index_z] - z[index_z - 1]);

    for (k = 0; k < nb_strike; k++) {
      greeks[0][k] = (values_p1[index_x + 1][index_z][k] -
                      values_p1[index_x - 1][index_z][k]) /
                     (dxu + dxd);
      greeks[1][k] =
          2.0 *
          (values_p1[index_x + 1][index_z][k] * dxd +
           values_p1[index_x - 1][index_z][k] * dxu - (dxd + dxu) * res[0]) /
          (dxu * dxd * (dxu + dxd));
      greeks[2][k] =
          (values[index_x][index_z][k] - values_p1[index_x][index_z][k]) / dt /
          DAYS_IN_YEAR;

      greeks[3][k] = (values_p1[index_x][index_z + 1][k] -
                      values_p1[index_x][index_z - 1][k]) /
                     (dzu + dzd);
      greeks[4][k] =
          2.0 *
          (values_p1[index_x][index_z + 1][k] * dzd +
           values_p1[index_x][index_z - 1][k] * dzu - (dzu + dzd) * res[0]) /
          (dzu * dzd * (dzu + dzd));

      greeks[5][k] = ((values_p1[index_x + 1][index_z + 1][k] -
                       values_p1[index_x - 1][index_z + 1][k]) /
                          (dxu + dxd) -
                      (values_p1[index_x + 1][index_z - 1][k] -
                       values_p1[index_x - 1][index_z - 1][k]) /
                          (dxu + dxd)) /
                     (dzu + dzd);
    }
  }

  if (calc_at_point) {
    /* look on z for the closest value to target */

    if ((values_p1[index_x][index_z + 1][column] - res[column]) /
            (z[index_z + 1] - z[index_z]) >
        0.0) {
      /* increasing function of index_z */
      is_up = 1;
    } else {
      is_up = -1;
    }

    if ((target - res[column]) * is_up > 0.0) {
      i = index_z;

      while (target < values_p1[index_x][i][column] && i < nstepz) {
        i++;
      }
      if (i == nstepz) {
        i = nstepz - 1;
      }

      vega =
          (values_p1[index_x][i][column] - values_p1[index_x][i - 1][column]) /
          (z[i] - z[i - 1]);

      newz = (target - values_p1[index_x][i - 1][column]) / vega + z[i - 1];

      *vol = sigma_beta - newz;

      coef = (newz - z[i - 1]) / (z[i] - z[i - 1]);

      for (k = 0; k < nb_strike; k++) {
        res_at_point[k] = coef * values_p1[index_x][i][k] +
                          (1.0 - coef) * values_p1[index_x][i - 1][k];
      }
    } else {
      i = index_z;

      while (target > values_p1[index_x][i][column] && i >= 0) {
        i--;
      }
      if (i == -1) {
        i = 0;
      }

      vega =
          (values_p1[index_x][i + 1][column] - values_p1[index_x][i][column]) /
          (z[i + 1] - z[i]);

      newz = (target - values_p1[index_x][i][column]) / vega + z[i];

      *vol = sigma_beta - newz;

      coef = (newz - z[i]) / (z[i + 1] - z[i]);

      for (k = 0; k < nb_strike; k++) {
        res_at_point[k] = coef * values_p1[index_x][i + 1][k] +
                          (1.0 - coef) * values_p1[index_x][i][k];
      }
    }
  }

  t2 = clock();

  smessage("Convolution        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /* Allocation 1 */
  if (x)
    free(x);
  if (z)
    free_dvector(z, 0, nstepz - 1);

  /* Allocation 2 */
  if (values)
    free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nb_strike - 1);
  if (values_p1)
    free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nb_strike - 1);

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
    num_f_pde_free_2d_adi(pde, nstepx, nstepz, nb_strike);

  return err;
}

Err op_sabrSL_calib_adi(double forward, double strike, double maturity,
                        double tgt_vol, SrtDiffusionType input_vol_type,
                        double alpha, double beta, double rho, double lambda,
                        int nt, int nx, int nz, int nbIter, double precision,
                        double floor_std, double *guess, double *res) {

  Err err = NULL;

  int j, k, l, nb_iter;
  double bsvol, vega;
  double premium, premium2, premium_tgt;
  double error, errorv;
  double sig0, newSig0;
  double shift_fwd;

  double coefa, coefb;
  double **res_iter = NULL;
  double shift_strike;

  time_t t1, t2;

  t1 = clock();

  /*	For stability reasons */
  if (alpha < 1.0E-05) {
    alpha = 1.0E-05;
  }

  coefa = beta * exp((beta - 1.0) * log(forward));
  coefb = (1.0 - beta) * exp(beta * log(forward));

  shift_fwd = coefa * forward + coefb;

  res_iter = dmatrix(0, nbIter, 0, 1);

  if (!res_iter) {
    err = "Memory allocation error (1) in op_sabrSL_calib_adi";
    goto FREE_RETURN;
  }

  strike = forward;
  shift_strike = shift_fwd;

  smessage("Starting Calibration of SabrSL");

  for (j = 0; j <= nbIter; j++) {
    res_iter[j][0] = 0.0;
    res_iter[j][1] = 0.0;
  }

  nb_iter = 0;

  if (input_vol_type == SRT_LOGNORMAL) {
    premium_tgt = srt_f_optblksch(forward, strike, tgt_vol, maturity, 1.0,
                                  SRT_CALL, PREMIUM);
  } else {
    err = "Input type shoulb be lognormal";
    goto FREE_RETURN;
  }

  err = srt_f_optimpvol(premium_tgt * coefa, shift_fwd, shift_strike, maturity,
                        1.0, SRT_CALL, SRT_LOGNORMAL, &sig0);

  sig0 /= coefa;

  /* use the SABR guess */
  sig0 = op_sabrSLcalib(forward, forward, maturity, tgt_vol, alpha, beta, rho,
                        lambda);

  if (guess) {
    sig0 = guess[0];
  }

  err = op_sabrSL_adi(forward, &strike, 1, maturity, 1.0, SRT_CALL, sig0, alpha,
                      beta, rho, lambda, floor_std, nt, nx, nz, &premium, 0,
                      NULL, 1, 0, premium_tgt, &newSig0, &premium2);

  if (err) {
    goto FREE_RETURN;
  }

  error = fabs(premium - premium_tgt);
  err = srt_f_optimpvol(premium, forward, strike, maturity, 1.0, SRT_CALL,
                        input_vol_type, &bsvol);

  errorv = fabs(bsvol - tgt_vol);

  if (errorv >= precision || err) {
    if (err || newSig0 < 0.0 || fabs(alpha < 0.001)) {
      if (err) {
        newSig0 = sig0 * 0.5;
      } else {
        newSig0 = sig0 * sqrt(premium_tgt * premium_tgt / premium / premium *
                              maturity);
      }
    }

    if (!err) {
      res_iter[0][0] = sig0;
      res_iter[0][1] = premium;
      nb_iter++;
    }
  } else {
    newSig0 = sig0;
  }

  /* Now do a Newton algorithm	*/

  k = 0;
  while ((k < nbIter) && (errorv >= precision)) {
    err = op_sabrSL_adi(forward, &strike, 1, maturity, 1.0, SRT_CALL, newSig0,
                        alpha, beta, rho, lambda, floor_std, nt, nx, nz,
                        &premium2, 0, NULL, 0, 0, 0, NULL, NULL);
    if (err) {
      goto FREE_RETURN;
    }

    err = srt_f_optimpvol(premium, forward, strike, maturity, 1.0, SRT_CALL,
                          input_vol_type, &bsvol);

    errorv = fabs(bsvol - tgt_vol);

    if (premium2 > 100 * premium_tgt || premium2 < 0.0 || err) {
      /* We have a problem */
      newSig0 *= 0.7;
    } else {
      /*
      newSig0 = sig0;
      */

      vega = (premium2 - premium) / (newSig0 - sig0);
      if (vega < 0) {
        /* problem !!! */
      } else {
        l = 0;
        while (l < nb_iter && res_iter[l][0] < newSig0) {
          l++;
        }

        for (j = nb_iter - 1; j >= l; j--) {
          res_iter[j + 1][0] = res_iter[j][0];
          res_iter[j + 1][1] = res_iter[j][1];
        }

        res_iter[l][0] = newSig0;
        res_iter[l][1] = premium2;
        nb_iter++;

        sig0 = newSig0;
        premium = premium2;

        newSig0 = solve_for_next_coef(res_iter, nb_iter, premium_tgt, 1);

        if (newSig0 < 0.0) {
          newSig0 = res_iter[0][0] * 0.25;
        }
      }
    }

    k++;
  }

  *res = newSig0;

  if (errorv <= precision) {
    smessage("Success at iteration %d", 1);
  } else {
    smessage("May have failed");
  }

  t2 = clock();
  smessage("Calibration time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  if (res_iter)
    free_dmatrix(res_iter, 0, nbIter, 0, 1);

  return err;
}

Err op_sabrSL_MC(double forward, double *strike, int nb_strike, double maturity,
                 double sigma_beta, double alpha, double beta, double rho,
                 double lambda, int npaths, int nsteps, double **res) {
  double a, b;
  double dt, sqdt;
  double rho2;
  long i, j;
  double S, U, U0, sig;
  double w1, w2;
  double const1, const2, const3;

  long seed = -123456789;

  dt = maturity / nsteps;
  sqdt = sqrt(dt);

  rho2 = sqrt(1.0 - rho * rho);
  a = beta * exp((beta - 1.0) * log(forward));
  b = (1.0 - beta) * exp(beta * log(forward));
  U0 = log(a * forward + b) / a;

  const1 = -0.5 * a * dt;
  const2 = lambda * dt;
  const3 = alpha * sqdt;

  for (j = 0; j < nb_strike + 2; j++) {
    res[j][0] = 0.0;
    res[j][1] = 0.0;
  }

  for (i = 0; i < npaths; i++) {
    U = U0;
    sig = sigma_beta;

    for (j = 0; j < nsteps; j++) {
      w1 = gauss_sample(&seed);
      w2 = rho * w1 + rho2 * gauss_sample(&seed);

      U += sig * (sig * const1 + sqdt * w1);
      sig += -const2 * (sig - sigma_beta) + const3 * sig * w2;
    }

    S = (exp(a * U) - b) / a;

    for (j = 0; j < nb_strike; j++) {
      if (S > strike[j]) {
        res[j][0] += (S - strike[j]) / npaths;
        res[j][1] += (S - strike[j]) * (S - strike[j]) / npaths;
      }
    }

    res[nb_strike][0] += S / npaths;
    res[nb_strike][1] += S * S / npaths;
    res[nb_strike + 1][0] += sig / npaths;
    res[nb_strike + 1][1] += sig * sig / npaths;
  }

  for (j = 0; j <= nb_strike + 1; j++) {
    res[j][1] = sqrt((res[j][1] - res[j][0] * res[j][0]) / npaths);
  }

  return NULL;
}

Err op_sabrSL_MC2(double forward, double *strike, int nb_strike,
                  double maturity, double sigma_beta, double alpha, double beta,
                  double rho, double lambda, int npaths, int nsteps,
                  int do_balsam, double **res) {
  double a, b;
  double dt, sqdt;
  double rho2;
  long i, j;
  double Z0, fwdZ, lnZ, lnZ0, sig;
  double w1;
  double const1, const2, const3, const4, const5;
  double lnsig, lnsig0;
  double price, var, sqvar, d1, d2;
  double intsig;
  double **matrix = NULL;
  double *matrixi;

  long seed = -123456789;
  Err err = NULL;

  dt = maturity / nsteps;
  sqdt = sqrt(dt);

  a = beta * exp((beta - 1.0) * log(forward));
  b = (1.0 - beta) * exp(beta * log(forward));

  Z0 = forward + b / a;
  lnZ0 = log(Z0);
  lnsig0 = log(sigma_beta);

  rho2 = a * a * (1.0 - rho * rho) * dt;

  const1 = -0.5 * a * a * rho * rho * dt;
  const2 = a * rho * sqdt;

  const3 = -0.5 * alpha * alpha * dt - lambda * dt;
  const4 = lambda * sigma_beta * dt;
  const5 = alpha * sqdt;

  for (j = 0; j < nb_strike; j++) {
    res[j][0] = 0.0;
    res[j][1] = 0.0;
    strike[j] += b / a;
  }

  /* fill the Brownian matrix */
  /* npaths has to be odd */
  npaths = 2 * ((long)(npaths / 2)) + 1;

  if (do_balsam) {
    matrix = dmatrix(0, npaths - 1, 0, nsteps - 1);
    if (!matrix) {
      return NULL;
    }

    err = balsam_generation(npaths, nsteps, matrix);
  }

  if (fabs(rho) < 1.0E-10) {
    for (i = 0; i < npaths; i++) {
      lnsig = lnsig0;
      var = 0.0;

      if (do_balsam) {
        matrixi = matrix[i];
      }

      for (j = 0; j < nsteps; j++) {
        if (do_balsam) {
          w1 = matrixi[j];
        } else {
          w1 = gauss_sample(&seed);
        }

        sig = exp(lnsig);
        var += sig * sig;
        lnsig += const3 + const4 / sig + const5 * w1;
      }

      var = var * rho2;
      sqvar = sqrt(var);
      var *= 0.5;

      for (j = 0; j < nb_strike; j++) {
        d1 = (log(Z0 / strike[j]) + var) / sqvar;
        d2 = d1 - sqvar;

        price = Z0 * norm(d1) - strike[j] * norm(d2);

        res[j][0] += price / npaths;
        res[j][1] += price * price / npaths;
      }
    }
  } else {
    for (i = 0; i < npaths; i++) {

      lnZ = lnZ0;
      lnsig = lnsig0;

      var = 0.0;
      intsig = 0.0;

      if (do_balsam) {
        matrixi = matrix[i];
      }

      for (j = 0; j < nsteps; j++) {
        if (do_balsam) {
          w1 = matrixi[j];
        } else {
          w1 = gauss_sample(&seed);
        }

        sig = exp(lnsig);
        var += sig * sig;
        intsig += sig * w1;
        lnsig += const3 + const4 / sig + const5 * w1;
      }

      lnZ += const1 * var + const2 * intsig;

      var = var * rho2;
      sqvar = sqrt(var);
      var *= 0.5;

      fwdZ = exp(lnZ);

      for (j = 0; j < nb_strike; j++) {
        d1 = (log(fwdZ / strike[j]) + var) / sqvar;
        d2 = d1 - sqvar;

        price = fwdZ * norm(d1) - strike[j] * norm(d2);

        res[j][0] += price / npaths;
        res[j][1] += price * price / npaths;
      }
    }
  }

  for (j = 0; j <= nb_strike + 1; j++) {
    res[j][1] = sqrt((res[j][1] - res[j][0] * res[j][0]) / npaths);
  }

  if (matrix)
    free_dmatrix(matrix, 0, npaths - 1, 0, nsteps - 1);

  return NULL;
}

double op_sabrSL(double F, double K, double T, double sigma, double alpha,
                 double beta, double rho, double lambda) {
  double Fb, gam1, gam2, Z, XZ;
  double coefa, coefb;
  double res;

  if (fabs(lambda) > 1.0E-10) {
    alpha *= sqrt((1.0 - exp(-2.0 * lambda * T)) / (2.0 * lambda * T));
  }

  coefa = beta * exp((beta - 1.0) * log(F));
  coefb = (1.0 - beta) * exp(beta * log(F));

  Fb = 0.5 * (F + K);
  gam1 = coefa / (coefa * Fb + coefb);
  gam2 = 0.0;
  Z = alpha / sigma * (F - K) / (coefa * Fb + coefb);
  XZ = log((sqrt(1.0 - 2.0 * rho * Z + Z * Z) + Z - rho) / (1.0 - rho));

  if (fabs(F - K) < 1.0E-10) {
    res = sigma * (coefa * F + coefb) / F;
  } else {
    res = sigma * log(F / K) * coefa /
          (log(coefa * F + coefb) - log(coefa * K + coefb)) * Z / XZ;
  }

  res *= (1.0 + ((2.0 * gam2 - gam1 * gam1 + 1.0 / Fb / Fb) * Fb * Fb / 24.0 +
                 alpha * Fb / sigma / (coefa * Fb + coefb) *
                     (0.25 * rho * gam1 * Fb + (2.0 - 3.0 * rho * rho) / 24.0 *
                                                   alpha * Fb / sigma /
                                                   (coefa * Fb + coefb))) *
                    sigma * sigma * (coefa * Fb + coefb) *
                    (coefa * Fb + coefb) / Fb / Fb * T);

  return res;
}

double op_sabrSLcalib(double F, double K, double T, double sigma, double alpha,
                      double beta, double rho, double lambda) {
  long k;
  double sigma_beta1, sigma_beta2;
  double vol1, vol2, vega;
  double error;

  sigma_beta1 = sigma * exp((1.0 - beta) * log(K));
  vol1 = op_sabrSL(F, K, T, sigma_beta1, alpha, beta, rho, lambda);

  error = fabs(vol1 - sigma);

  /* First run */
  if (vol1 < sigma) {
    sigma_beta2 = sigma_beta1 * 1.005;
  } else {
    sigma_beta2 = sigma_beta1 * 0.995;
  }

  vol2 = op_sabrSL(F, K, T, sigma_beta2, alpha, beta, rho, lambda);

  sigma_beta2 = sigma_beta1 +
                (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);
  vol2 = op_sabrSL(F, K, T, sigma_beta2, alpha, beta, rho, lambda);

  error = fabs(vol2 - sigma);
  k = 0;

  while (error > 0.000001 && k < 50) {
    vega = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
    sigma_beta1 = sigma_beta2;
    vol1 = vol2;

    sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;
    vol2 = op_sabrSL(F, K, T, sigma_beta2, alpha, beta, rho, lambda);

    error = fabs(vol2 - sigma);

    k++;
  }

  return sigma_beta2;
}

Err op_sabrQuad_MC(double forward, double *strike, int nb_strike,
                   double maturity, double sigma_beta, double alpha, double a,
                   double b, double c, double rho, double lambda, int npaths,
                   int nsteps, int do_balsam, double **res) {
  double dt, sqdt;
  double rho2;
  long i, j;
  double Z0, Z, sig, locvol;
  double w1, w2;
  double const1, const2, const3, const4, const5;
  double lnZ, lnZ0, lnsig, lnsig0;
  double price;
  double varZ, covar, fwdZ, varOpt;

  long seed = -123456789;
  Err err = NULL;

  dt = maturity / nsteps;
  sqdt = sqrt(dt);

  Z0 = forward;
  lnZ0 = log(Z0);
  lnsig0 = log(sigma_beta);

  rho2 = sqrt(1.0 - rho * rho);

  const1 = sigma_beta * sqdt;
  const2 = a * rho * sqdt;

  const3 = -0.5 * alpha * alpha * dt - lambda * dt;
  const4 = lambda * sigma_beta * dt;
  const5 = alpha * sqdt;

  for (j = 0; j < nb_strike + 3; j++) {
    res[j][0] = 0.0;
    res[j][1] = 0.0;
  }

  a *= sqdt;
  b *= sqdt;
  c *= sqdt;

  /* fill the Brownian matrix */
  /* npaths has to be odd */
  npaths = 2 * ((long)(npaths / 2)) + 1;

  for (i = 0; i < npaths; i++) {
    lnZ = lnZ0;
    Z = Z0;
    lnsig = lnsig0;

    for (j = 0; j < nsteps; j++) {
      w1 = gauss_sample(&seed);
      w2 = rho * w1 + rho2 * gauss_sample(&seed);

      sig = exp(lnsig);

      /*
      Z = exp(lnZ);
      locvol = sig * (a * Z + b + c / Z);

      lnZ += locvol * (-0.5 * locvol + w1);
      */

      locvol = sig * (Z * (a * Z + b) + c);

      Z += locvol * w1;

      if (Z < 1.0E-08) {
        Z = 1.0E-08;
      }

      lnsig += const3 + const4 / sig + const5 * w2;
    }

    /*
    Z = exp(lnZ);
    */

    sig = exp(lnsig);

    for (j = 0; j < nb_strike; j++) {
      price = max(Z - strike[j], 0);

      res[j][0] += price / npaths;
      res[j][1] += price * price / npaths;
      res[j][2] += price * Z / npaths;
    }

    res[nb_strike][0] += sig / npaths;
    res[nb_strike][1] += sig * sig / npaths;
    res[nb_strike + 1][0] += sig * sig / npaths;
    res[nb_strike + 1][1] += sig * sig * sig * sig / npaths;
    res[nb_strike + 2][0] += Z / npaths;
    res[nb_strike + 2][1] += Z * Z / npaths;
  }

  /* Compute the delta */
  fwdZ = res[nb_strike + 2][0];
  varZ = res[nb_strike + 2][1] - res[nb_strike + 2][0] * res[nb_strike + 2][0];

  for (j = 0; j < nb_strike; j++) {
    varOpt = res[j][1] - res[j][0] * res[j][0];
    covar = res[j][2] - res[j][0] * fwdZ;
    res[j][2] = covar / varZ;
    res[j][0] += res[j][2] * (forward - fwdZ);
    res[j][1] =
        sqrt((varOpt + res[j][2] * (res[j][2] * varZ - 2.0 * covar)) / npaths);
  }

  for (j = nb_strike; j < nb_strike + 3; j++) {
    res[j][1] = sqrt((res[j][1] - res[j][0] * res[j][0]) / npaths);
  }

  return NULL;
}