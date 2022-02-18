/*--------------------------------------------------------------
        FILE: srt_h_ra_calib.h
        PURPOSE: Range Accrual product SABR sigma-beta calibration
        AUTHOR: Dimitri Mayevski
        DATE: 01/08/2002
  --------------------------------------------------------------*/

#ifndef __SRT_H_RA_CALIB_H__
#define __SRT_H_RA_CALIB_H__

Err srt_f_calib_range_accrual(
    char*       yc_d,
    char*       yc_f,
    SrtProduct* product,
    int         nex,
    long*       ex_dates,
    double*     sig_d,
    double*     sig_f,
    double*     sig_fx,
    double      lam_d,
    double      lam_f,
    double      rho_df,
    double      rho_ffx,
    int         nx);

#endif /* #ifndef __SRT_H_RA_CALIB_H__ */