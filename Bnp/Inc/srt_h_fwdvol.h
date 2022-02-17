/*
        FILE NAME				SRT_F_FWDVOL.C
        AUTHOR					NONO SAVINE
        PURPOSE					COMPUTE FORWARD VOLATILITY  ,
   EITHER THROUGH AUTOCALIBRATION  , OR PARABOLAS
*/

#ifndef _SRT_H_FWDVOL_H
#define _SRT_H_FWDVOL_H

/*	Main function for autocal option */
Err srt_f_fwdvol_autocal(
    String model,   /* Model (LGM  , CHEY  , CHEY_BETA  , 1F or 2F) */
    long *tau_date, /* Tau Term Struct */
    double *tau_value, long num_tau, int fix_tau_flag,

    double alpha,                              /* 2F parameters */
    double beta, double rho, double chey_beta, /* Chey Beta parameter */
    long num_period, /* Number of required forward vols */
    long *vol_start, /* Vol Start Dates array */
    long *vol_end,   /* Vol End Dates array */

    int option_type, /* Option Type */

    long *und_start_act,                   /* Und start date array */
    char *und_tenor, SrtCompounding compd, /* Und compd */
    SrtBasisCode basis,                    /* Und basi */
    String ref_rate,                       /* Reference Rate Code */
    double *und_val,

    String yc_id,       /* YC id */
    String bs_vol_type, /* Lognormal or normal */
    Err(*GetVol)        /* GetVol function */
    (Ddate start, Ddate end, double strike, double dForward, double dSpread,
     double *bs_vol),
    SrtPriceType price_type, double *output_array, /* Output array */
    String fwd_vol_type);

Err srt_f_fwdvol_parabola(double level_slope_correl, /* Parameters */
                          double slope_vol,
                          long num_period, /* Number of required forward vols */
                          Ddate *vol_start,   /* Vol Start Dates array */
                          Ddate *vol_end,     /* Vol End Dates array */
                          String ref_rate,    /* Reference Rate Code */
                          String yc_id,       /* YC id */
                          String bs_vol_type, /* Lognormal or normal */
                          Err(*GetVol)        /* GetVol function */
                          (Ddate start, Ddate end, double strike,
                           double dForward, double dSpread, double *bs_vol),
                          double *fwd_vol, /* Output array */
                          long max_strike, long max_vol, double delta_strike,
                          double tol, long nvol);

Err srt_f_spfwdnormvol(
    Date today, Date vol_start, Date vol_end, SrtCompounding compd,
    double fwd_cms_start, double fwd_cms_end, double level_slope_correl,
    double slope_vol,
    Err (*GetVol)(Ddate, Ddate, double, double, double, double *),
    char *vol_type, double *fwd_vol, long max_strike, long max_vol,
    double delta_strike, double tol, long nvol);

#endif