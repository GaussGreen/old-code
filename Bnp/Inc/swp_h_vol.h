/* ===========================================================================

         FILENAME:     swp_h_vol.h

     PURPOSE:      The function that will have to be used everywhere to compute
                       a volatility for a given strike

   ===========================================================================
 */

#ifndef SWP_H_VOL_H
#define SWP_H_VOL_H

Err swp_f_vol(char *vol_id, Ddate start, Ddate end, double strike,
              double *volatility, double *power);

Err swp_f_SABRvol(char *vol_id, Ddate start, Ddate end, double strike,
                  double *volatility, double *power, int component);

Err swp_f_transformvol(char *yc_id, char *vol_id, char *freq_std,
                       char *basis_std, char *refRateCode, Ddate start,
                       Ddate end, double strike, char *freq, char *basis,
                       double *volatility, double *power);

Err swp_f_transformvol_func(
    char *yc_id, char *vol_id,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *freq_std, char *basis_std, char *refRateCode, Ddate start, Ddate end,
    double strike, char *freq, char *basis, double *volatility, double *power);

Err swp_f_IsSmileVol(char *vol_id,
                     double *yesno); // return 1 if Model is SABR or SABRAF.

#endif
