#ifndef __CMS_SPREAD_TIMESWAP_H
#define __CMS_SPREAD_TIMESWAP_H

#define LINEAR_METHOD 0
#define PIECEWISE_CONSTANT_METHOD 1

#define EXOTIC_LEG_PAY 0
#define EXOTIC_LEG_REC 1

Err cms_spread_digital_put(
    long today,
    long fixing_date,
    long start_date,
    long pay_date,
    int  exo_pay_rec,

    double fixing_time,
    double delay,

    double cms1_weight,
    long   cms1_end_date,
    char*  cms1_frequency,
    char*  cms1_basis,
    char*  cms1_refrate,
    char*  cms1_dom_for,

    char* cms1_mkt_ID,
    char* cms1_yc_name,
    char* cms1_vc_name,
    int   cms1_mkt_isSABR,

    double         cms1_nfp,
    SrtCompounding cms1_compounding,
    double         cms1_dRateConv,

    double*          cms1_strikes_in_vol,
    int              n_cms1_strikes_in_vol,
    SrtDiffusionType cms1_vol_type,
    int              cms1_cash_vol,

    double cms2_weight,
    long   cms2_end_date,
    char*  cms2_frequency,
    char*  cms2_basis,
    char*  cms2_refrate,
    char*  cms2_dom_for,

    char* cms2_mkt_ID,
    char* cms2_yc_name,
    char* cms2_vc_name,
    int   cms2_mkt_isSABR,

    double         cms2_nfp,
    SrtCompounding cms2_compounding,
    double         cms2_dRateConv,

    double*          cms2_strikes_in_vol,
    int              n_cms2_strikes_in_vol,
    SrtDiffusionType cms2_vol_type,
    int              cms2_cash_vol,

    double cms1_cms2_corr_bid,
    double cms1_cms2_corr_ask,

    double cms1_quanto_corr,
    double cms2_quanto_corr,
    double fwd_fx_vol,

    double       strike[],
    double       call_spread,
    int          use_cms_vol,
    long         copula_npts,
    long         copula_nsim,
    int          copula_degree,
    int          copula_nconv,
    SrtMCSamType copula_mctype,
    int          calculate_copula_price,
    double       copula_correlation_bump,

    double* lower_digital_lognormal_daily,
    double* lower_digital_normal_daily,
    double* lower_digital_copula_daily,

    double* upper_digital_lognormal_daily,
    double* upper_digital_normal_daily,
    double* upper_digital_copula_daily,

    double* atm_lognormal_daily,
    double* atm_normal_daily,
    double* atm_copula_daily);

Err cms_spread_timeswap_leg(
    char* dom_yc_name,
    long  today,
    int   exo_pay_rec,
    char* exo_basis,

    long*   exo_start_dates,
    long*   exo_end_dates,
    double* exo_coupons,
    double* exo_notionals,
    double* exo_lower_barriers,
    double* exo_upper_barriers,
    int     n_exo_coupons,

    long** index_fixing_dates,
    int*   n_index_fixing_dates,
    int    index_fixing_lag,

    double cms1_weight,
    char*  cms1_tenor,
    char*  cms1_frequency,
    char*  cms1_basis,
    char*  cms1_refrate,
    char*  cms1_dom_for,

    char* cms1_mkt_ID,
    char* cms1_yc_name,
    char* cms1_vc_name,
    int   cms1_mkt_isSABR,
    char* cms1_mkt_frequency,
    char* cms1_mkt_basis,
    char* cms1_mkt_refrate,

    double*          cms1_strikes_in_vol,
    int              n_cms1_strikes_in_vol,
    SrtDiffusionType cms1_vol_type,
    int              cms1_cash_vol,

    double cms2_weight,
    char*  cms2_tenor,
    char*  cms2_frequency,
    char*  cms2_basis,
    char*  cms2_refrate,
    char*  cms2_dom_for,

    char* cms2_mkt_ID,
    char* cms2_yc_name,
    char* cms2_vc_name,
    int   cms2_mkt_isSABR,
    char* cms2_mkt_frequency,
    char* cms2_mkt_basis,
    char* cms2_mkt_refrate,

    double*          cms2_strikes_in_vol,
    int              n_cms2_strikes_in_vol,
    SrtDiffusionType cms2_vol_type,
    int              cms2_cash_vol,

    double* cms1_cms2_corr_dates,
    double* cms1_cms2_corr_values_bid,
    double* cms1_cms2_corr_values_ask,
    int     n_cms1_cms2_corr,

    double* cms1_quanto_corr_dates,
    double* cms1_quanto_corr_values,
    int     n_cms1_quanto_corr,

    double* cms2_quanto_corr_dates,
    double* cms2_quanto_corr_values,
    int     n_cms2_quanto_corr,

    double* fwd_fx_vol_dates,
    double* fwd_fx_vol_values,
    int     n_fwd_fx_vol,

    double* lower_digital_lognormal,
    double* lower_digital_normal,
    double* lower_digital_copula,

    double* upper_digital_lognormal,
    double* upper_digital_normal,
    double* upper_digital_copula,

    double* CMSSpreadTimeSwapLeg_lognormal,
    double* CMSSpreadTimeSwapLeg_normal,
    double* CMSSpreadTimeSwapLeg_copula,

    double* atm_lognormal,
    double* atm_normal,
    double* atm_copula,

    double       call_spread,
    int          monitoring_frequency,
    int          use_cms_vol,
    long         copula_npts,
    long         copula_nsim,
    int          copula_degree,
    int          copula_nconv,
    SrtMCSamType copula_mctype,
    int          calculate_copula_price,
    double       copula_correlation_bump);

#endif