/* ====================================================================================
   FILENAME:  srt_h_ytat.h

   PURPOSE:   The parameters needed at a date t to compute a DF maturity T
              in an interest rate model
   ====================================================================================
 */

#ifndef SRT_H_YTAT_H
#define SRT_H_YTAT_H

/* ---------------------------------------------------------------------------
  TYPE            :YTatt_param
  DESCRIPTION     :parameter for computing the zero coupon yield to time T
                   having reached time t
  DEFINITION      :

  ----------------------------------------------------------------------------
*/

typedef struct YTatt_param {
  double fwd_zero_rate;
  double x_coeff[2];
  double phi_coeff[2][2];

  double vasicek_sr_coeff;
  double vasicek_mean_coeff;
  double vasicek_var_coeff;

  double A_t;
  double A_T;

  double h_t;
  double h_T;

  double beta;
  double eta;
  double lambda_t;
  double lambda_T;
  double zeta_t;
  double zeta_T;
  double **M;

} YTatt_param;

/* ----------------------------------------------------------------- */

Err Vasicek_Y_T_at_param(Ddate fix_date, Ddate *pay_dates, long dim_pay_date,
                         SrtUndPtr und, YTatt_param *YT_sam);

Err Y_T_at_t_param(Ddate t, Ddate *T_mat, long dim_Tmat, SrtUndPtr und,
                   YTatt_param *YT_sam);

Err Y_T_at_t_compute(long dim_Tmat, SrtSample *sam, YTatt_param *YT_sam,
                     double *YTt_yield, int index, SrtMdlDim mdl_dim,
                     SrtMdlType mdl_type);

Err Vasicek_Y_T_at_t_compute(long dim_Tmat, SrtSample *sam, YTatt_param *YT_sam,
                             double *YTt_yield, int index, SrtMdlDim mdl_dim,
                             SrtMdlType mdl_type);

#endif