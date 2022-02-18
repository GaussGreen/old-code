/*
        FILE NAME				SRT_F_FWDVOLFNC.C
        AUTHOR					NONO SAVINE
        PURPOSE					UYTILITY FUNCTIONS USED IN SRT_F_FWDVOL.C
*/

#ifndef SRT_H_FWDVOLFNC_H
#define SRT_H_FWDVOLFNC_H

Err capfwd_vol_set_calibration_instruments(
    long           num_period,
    long*          vol_start,
    long*          vol_end,
    String         ref_rate,
    SrtCompounding compd,
    SrtBasisCode   basis,
    String         yc_id,
    String         bs_vol_type,
    Err (*GetVol)(
        Ddate start, Ddate end, double strike, double dForward, double dSpread, double* bs_vol),
    long*    cal_numinst,
    Date**   cal_start,
    Date**   cal_end,
    String** cal_freq,
    String** cal_basis,
    double** cal_str,
    double** cal_bndstr,
    String** cal_type,
    String** cal_recpay,
    String** cal_refrate,
    double** cal_price,
    double** cal_vega);

Err cmsfwd_vol_set_calibration_instruments(
    long  num_period,
    int   fix_tau_flag,
    long* vol_start,
    long* vol_end,

    char*          und_tenor,
    SrtCompounding swap_compd,
    SrtBasisCode   swap_basis,
    String         ref_rate,

    String yc_id,
    String bs_vol_type,
    Err(*GetVol) /* GetVol function */
    (Ddate start, Ddate end, double strike, double dForward, double dSpread, double* bs_vol),
    long*    cal_numinst,
    long**   cal_start,
    long**   cal_end,
    String** cal_freq,
    String** cal_basis,
    double** cal_str,
    double** cal_bndstr,
    String** cal_type,
    String** cal_recpay,
    String** cal_refrate,
    double** cal_price,
    double** cal_vega);

void fwd_vol_set_grfn_options(
    SrtMdlType mdl_type, String** grfn_param, String** grfn_value, long* num_grfn_param);

void fwd_vol_set_calib_options(
    SrtMdlType mdl_type,
    int        option_type,
    int        fix_tau_flag,
    String**   calib_param,
    String**   calib_value,
    long*      num_calib_param);

Err srt_f_fwd_resetcaplet(
    SrtUndPtr       und,
    SrtGrfnParam*   grfnparam,
    String          yc_id,
    Date            today,
    String          ref_rate,
    SrtMdlType      mdl_type,
    SrtReceiverType rec_pay,
    Ddate           str_fix,
    Ddate           spot_fix,
    Ddate           start_act,
    Ddate           end_th,
    String          compdStr,
    String          basisStr,
    double*         premium);

Err srt_f_fwd_resetcms(
    SrtUndPtr     irundptr,
    SrtGrfnParam* grfnparam,
    Date          strike_fix,
    Date          cms_fix,
    Date          start_act,
    Date          end_theo,
    String        basisStr,
    String        compStr,
    String        ref_rate,
    String        yc_id,
    Err (*GetVol)(Ddate, Ddate, double, double, double, double*),
    String  bs_vol_type,
    double* premium);

Err parabola(
    Date           today,
    Date           start,
    SrtCompounding compd,
    double         fwd_cms,
    SrtCallPutType call_put,
    Err (*GetVol)(Ddate, Ddate, double, double, double, double*),
    char*   vol_type,
    double* answer,
    long    max_strike,
    long    max_vol,
    double  delta_strike,
    double  tol,
    long    nvol);

#endif