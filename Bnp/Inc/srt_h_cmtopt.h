

Err srt_f_get_cmt_vol(
    String           swp_yc_name,
    SrtDiffusionType swp_vol_type,

    double yr_to_exp,
    double cms_rate,
    double cmt_rate,
    double cmt_strike,

    /* SABR CMS PARAMETERS */
    double atm_sigma_beta,
    double sabr_alpha,
    double sabr_beta,
    double sabr_rho,

    /* SPREAD DYNAMICS PARAMETERS */
    double      spot_spread,      /* THE SPREAD INITIAL CONDITION */
    double      spread_tau,       /* ONE OVER THE MEAN REVERSION */
    double      spread_fwd_level, /* THE MEAN REVERSION LEVEL */
    double      spread_loc_vol,   /* THE SPREAD VOLATILITY */
    double      rho,              /* CORRELATION BETWEEN CMS AND SPREAD */
    SRT_Boolean use_mkt_spread,

    double* cmt_vol);

Err srt_f_calibrate_CMT_CMS_spread_parameters(
    double yr_to_exp,
    String swp_yc_name,
    String swp_vol_type_str,
    double atm_sigma_beta,
    double sabr_alpha,
    double sabr_beta,
    double sabr_rho,

    double cmt_rate,
    double spot_spread,
    double fwd_spread,

    double spread_mean_reversion_guess,
    double spread_fwd_level_guess,
    double spread_loc_vol_guess,
    double rho_guess,

    long    num_strike,
    double* cmt_calibration_strikes,
    double* cmt_vols,

    double** param);