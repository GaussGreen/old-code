/* SWP_H_CMT_INIT.H*/

#ifndef SWP_H_CMT_INIT_H
#define SWP_H_CMT_INIT_H

Err swp_f_strip_CMT(
    String  fwd_T_name,
    String* tenor_names,
    double* spread_rates,
    int     num_spreads,
    Date    cmt_today_long,
    String  yc_name,
    String  cmt_index_string,
    String  vc_name,
    String  mkt_name,
    Err (*GetVol)(long, long, double, double*),
    CMT_Param** cmt_param);

Err swp_f_store_fwd_spreads(
    String  fwd_spread_name,
    double* fwd_spread_dates,
    double* fwd_spread_values,
    int     num_fwd_spreads,
    String  yc_name,
    String  cmt_index_string,
    String  vc_name,
    String  mkt_name,
    Err (*GetVol)(long, long, double, double*),
    CMT_Param* cmt_param);

Err swp_f_strip_treas(
    String  fwd_T_name,
    String* tenor_names,
    double* spread_rates,
    int     num_spreads,
    Date    cmt_today_long,
    String  yc_name,
    String  cmt_index_string,
    String  vc_name,
    String  mkt_name,
    Err (*GetVol)(long, long, double, double*),
    CMT_Param** cmt_param);
#endif
