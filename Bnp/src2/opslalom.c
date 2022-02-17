/*******************************************************************************
**
**	OPSLALOM.C
**
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

double srt_f_optslalom(double fwd, double spot, double barrier, double vol,
                       double disc_fact, int today, int nb_dates,
                       int date[MAX_NUM_DATES], int up_do[MAX_NUM_DATES],
                       int step_num, SrtGreekType greek) {

  int i, j, n, opt_mat, c;
  double u, d, p_u, p_d, h, t_n, t_n_f;
  double a, bb, drift, df_n;
  double asset, disc_premium, intrinsic;
  double premium[MAX_STEP];

  opt_mat = date[nb_dates - 1];
  drift = log(fwd / spot) / opt_mat;

  h = (double)((opt_mat - today) / (double)step_num / 365.0);
  a = exp(drift * h);
  bb = a * a * (exp(vol * vol * h) - 1);
  u = ((a * a + bb + 1) +
       sqrt((a * a + bb + 1) * (a * a + bb + 1) - 4 * a * a)) /
      (2 * a);
  d = 1 / u;
  p_u = (a - d) / (u - d);
  p_d = 1 - p_u;

  asset = spot * pow(u, (double)(step_num));
  for (i = 0; i <= step_num; i++) {
    intrinsic = (asset - barrier) * (1 - 2 * up_do[nb_dates - 1]);

    if (intrinsic > 0)
      premium[i] = 100;
    else
      premium[i] = 0.0;

    asset /= u * u;
  }

  df_n = pow(disc_fact, 1 / (double)step_num);
  ;
  for (n = step_num - 1; n >= 0; n--) {
    t_n = (double)(n * h * 365 + today);
    t_n_f = (double)((n + 1) * h * 365 + today);

    c = 0;
    for (j = 0; j < nb_dates; j++) {
      if ((date[j] >= t_n) && (date[j] < t_n_f)) {
        c = 1;
        break;
      }
    };

    asset = spot * pow(u, (double)(n));
    for (i = 0; i <= n; i++) {
      disc_premium = (p_u * premium[i] + p_d * premium[i + 1]) * df_n;

      if (c == 1) {
        intrinsic = (asset - barrier) * (1 - 2 * up_do[j]);
        if (intrinsic > 0)
          premium[i] = disc_premium;
        else
          premium[i] = 0.0;
      } else {
        premium[i] = disc_premium;
      }
      asset /= u * u;
    }
  }
  return (premium[0]);
}

/* ========================================================================== */
