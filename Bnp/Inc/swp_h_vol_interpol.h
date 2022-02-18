/* ===========================================================================

   FILENAME:     swp_h_vol_interpol.h

   PURPOSE:      Provide functions to interpolate on a volatility matrix

  ============================================================================= */

#ifndef SWP_H_VOL_INTERPOL_H
#define SWP_H_VOL_INTERPOL_H

Err vol_interp_from_frd(
    double** array,
    int      nrow,
    int      ncol,
    double   date1,
    double   daysfromstart,
    int      method,
    double   today,
    double*  vol);

double vol_interpol_function(
    double   tgt_exp_date,
    double   tgt_und_mat,
    double*  exp_dates_vec,
    long     l_exp_dates,
    double*  und_mats_vec,
    long     l_und_mats,
    double** mkt_vol_2darray,
    long     m_exp_dates,
    long     m_und_mats);

#endif
