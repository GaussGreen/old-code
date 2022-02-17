
/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

static void local_Uniform_With_Antithetic_Noise(double **m, long pathindexstart,
                                                long pathindexend,
                                                long stepindexstart,
                                                long stepindexend, long *seed);

/*******************************************************************************
 *
 * FUNCTION      : srt_f_optmltspr(...)
 *
 * PURPOSE       : price an option with payoff
 *		 	( a X + bY + cZ - K )+
 *
 * DESCRIPTION   : uses a MonteCarlo type numerical integration with UWAN
 *		  technique to reduce the variance of the result
 *
 * CALLS         : local_Uniform_With_Antithetic_Noise
 *
 * PARAMETERS    : fwdx          - forward price of 1st underlying
 *               : fwdy          - forward price of 2nd underlying
 *               : fwdz          - forward price of 3rd underlying
 *               : sigx          - vol of 1st und
 *               : sigy          - vol of 2nd und
 *               : sigz          - vol of 3rd und
 *               : rhoxy         - correlation between x and y
 *               : rhoxz         - correlation between x and z
 *               : rhoyz         - correlation between y and z
 *				: a_x			- gearing of 1st und
 *               : b_y           - gearing of 2nd und
 *               : c_z           - gearing of 3rd und
 *               : k             - strike
 *               : mat           - maturity
 *               : disc          - discount factor
 *               : call_put      - option type
 *		: paths		- number of MC paths used for integration
 *               : greek         - greek wanted (only premium)
 *
 * RETURNS       : ??            - ??
 *
 *******************************************************************************/

/*
        Please note that the following conditions on correlations must
        be satisfied:

 (rhoxy*rhoxy+rhoyz*rhoyz-2*rhoxy*rhoxz*rhoyz) * ( 1 - rhoxz * rhoxz  )  <= 1
 (rhoyz*rhoyz+rhoxz*rhoxz-2*rhoxy*rhoxz*rhoyz) * ( 1 - rhoxy * rhoxy  )  <= 1
 (rhoxy*rhoxy+rhoxz*rhoxz-2*rhoxy*rhoxz*rhoyz) * ( 1 - rhoyz * rhoyz  )  <= 1

*/
double srt_f_optmltspd(double fwdx, double fwdy, double fwdz, double sigx,
                       double sigy, double sigz, double rhoxy, double rhoxz,
                       double rhoyz, double a_x, double b_y, double c_z,
                       double strike, double mat, double disc,
                       SrtCallPutType call_put, long paths,
                       SrtGreekType greek) {

  double **m, u = 0.0, v = 0.0, y, z;
  double fwd_xyz, strike_xyz, bs_xyz, d_xyz;
  double varx, vary, varz;
  double sqrt_varx, sqrt_vary, sqrt_varz;
  double alpha, beta, gamma;
  double sum, price;
  long path_numb = paths;
  long seed = RANDINIT;
  long i;
  double cp;

  cp = (call_put == SRT_CALL) ? 1 : -1;

  varx = sigx * sigx * mat;
  vary = sigy * sigy * mat;
  varz = sigz * sigz * mat;

  sqrt_varx = sqrt(varx);
  sqrt_vary = sqrt(vary);
  sqrt_varz = sqrt(varz);

  beta = (rhoxy - rhoxz * rhoyz) / (1 - rhoyz * rhoyz);
  gamma = (rhoxz - rhoyz * rhoxy) / (1 - rhoyz * rhoyz);
  alpha = sqrt(1 - beta * beta - gamma * gamma - 2 * beta * gamma * rhoyz);

  m = dmatrix(1, path_numb, 1, 2);

  local_Uniform_With_Antithetic_Noise(m, 1, path_numb, 1, 2, &seed);

  sum = 0.0;

  for (i = 1; i <= path_numb; i++) {
    u = inv_cumnorm_newton(0, m[i][1]);
    v = inv_cumnorm_newton(0, m[i][2]);

    y = sqrt(0.5 * (1 + rhoyz)) * u + sqrt(0.5 * (1 - rhoyz)) * v;
    z = sqrt(0.5 * (1 + rhoyz)) * u - sqrt(0.5 * (1 - rhoyz)) * v;

    fwd_xyz = a_x * fwdx * exp(0.5 * varx * (alpha * alpha - 1)) *
              exp(sqrt_varx * (beta * y + gamma * z));

    strike_xyz = strike - b_y * fwdy * exp(y * sqrt_vary - 0.5 * vary) -
                 c_z * fwdz * exp(z * sqrt_varz - 0.5 * varz);

    if (strike_xyz > 0.0) {
      d_xyz = log(fwd_xyz / strike_xyz) / (alpha * sqrt_varx);
      bs_xyz = cp * (fwd_xyz * norm(cp * (d_xyz + 0.5 * alpha * sqrt_varx)) -
                     strike_xyz * norm(cp * (d_xyz - 0.5 * alpha * sqrt_varx)));
    } else
      bs_xyz = fwd_xyz - strike_xyz;
    sum += bs_xyz;
  }
  price = sum / path_numb;
  price *= disc;

  free_dmatrix(m, 1, path_numb, 1, 2);
  return (price);
}

static void local_Uniform_With_Antithetic_Noise(double **m, long pathindexstart,
                                                long pathindexend,
                                                long stepindexstart,
                                                long stepindexend, long *seed) {

  long j, k;
  double unif;

  long step_numb, path_numb, index, N_2, length, left_paths;
  long newstart;

  step_numb = stepindexend - stepindexstart + 1;
  path_numb = pathindexend - pathindexstart + 1;

  N_2 = (long)(path_numb / 2);
  length = (long)pow((double)N_2, 1.0 / (double)step_numb);
  /* Closest 1/nth power that matches the number of paths */
  N_2 = (long)pow((double)length, (double)step_numb);

  left_paths = path_numb - 2 * N_2;

  for (j = 0; j < N_2; j++) {
    for (k = 0; k < step_numb; k++) {
      unif = uniform(seed);
      index = (long)((double)j / pow((double)length, k)) % length;
      m[pathindexstart + 2 * j][stepindexstart + k] =
          ((double)index + unif) / length;
      m[pathindexstart + 2 * j + 1][stepindexstart + k] =
          ((double)(index + 1) - unif) / length;
    }
  }
  if (left_paths > 0) {
    newstart = pathindexend - left_paths + 1;
    for (j = 0; j < left_paths; j++) {
      for (k = 0; k < step_numb; k++)
        m[newstart + j][stepindexstart + k] = uniform(seed);
    }
  }
}
