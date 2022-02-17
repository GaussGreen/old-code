
#ifndef __CTS_RESERVE_H
#define __CTS_RESERVE_H

void shift_LGMSV_model(long today, char *yc, LGMSV_MODEL model, double shift,
                       double log_norm);

Err cts_get_one_time_infos(CTS_MKT mkt, CTS_ADI_ARG adi_arg_for_midat, CTS cts,
                           double vol_shift, double *onetime_reserve,
                           int *one_time_calc, double *european_prices,
                           double *european_prices_bumped,
                           double *european_fees, int *nb_one_time_calc,
                           double *most_expensive_pv, CTS_IV fwd_iv_info);

Err cts_calibrate_european(CTS_UND und, CTS market_cts, CTS model_cts,
                           int for_fund, int calc_fwd_iv, int adj_fee,

                           CTS_ADI_ARG adi_arg, LGMSV_NUMERPARAMS numer_params,

                           int nb_european, double *european_values,

                           int nb_iter,

                           double *fitting_error);

Err cts_convert_to_equi_midat(CTS_MKT mkt, CTS mkt_cts, CTS mdl_cts);

Err cts_calibrate_equi_midat(CTS_MKT mkt, CTS_UND und, CTS cts, char *fund_freq,
                             char *fund_basis,

                             int notperiod, int *one_time_calc, double *prices,
                             double *fees,

                             double lambda,

                             double numer_tstar, double min_time);

Err cts_optimpvol(double premium, double fwd_price, double strike, double mat,
                  double disc, SrtCallPutType call_put,
                  SrtDiffusionType lognormal_normal, double *implied_vol);

Err cts_calc_reserve(/*	Initial results */
                     CTS cts_mkt, CTS cts_mdl, CTS_UND und, CTS_ADI_ARG adi_arg,

                     int for_fund,

                     double call, double *onetime_call,

                     /*	Model parameters */
                     int nsmilepar, double *smilepartime, double *alphaeps,
                     double *rhoeps, double *ldaeps, double *rho2eps,
                     double tstar,

                     /*	Calibration parameters */
                     LGMSV_NUMERPARAMS NumerParams, double numer_tstar,
                     char *cal_tenor, char *cal_ref, char *cal_freq,
                     char *cal_basis, int force_atm, double max_std_long,
                     double max_std_short, double vol_shift_long,
                     DIAGCALIB_VOLTYPE vol_type_long,
                     DIAGCALIB_SHIFTTYPE vol_shift_type_long,
                     double vol_shift_short, DIAGCALIB_VOLTYPE vol_type_short,
                     DIAGCALIB_SHIFTTYPE vol_shift_type_short,
                     double lambda_shift,
                     int calib_strategy, /*	-1: autocal  , 0: swaptions /
                                            cap  , 1: cap / swaptions */
                     int fix_lambda, char *short_tenor, char *short_refrate,
                     char *short_freq, char *short_basis, int fix_smile,
                     int smile_calib_months, /* 0: co-terminal swaption  ,
                                                otherwise underlyings with
                                                required nb months */
                     LGMSV_CalibParams *lgmsv_calib_params, double min_time,
                     int skip_last,   /*	If 1  , the last option is disregarded
                                         and the forward volatility is flat from
                                         option n-1 */
                     double min_fact, /*	Maximum down jump on variance */
                     double max_fact, /*	Maximum up jump on variance */
                     int use_jumps,   /*	1: we allow jumps on vol  , 0: we don't
                                       */
                     double prec, int maxiter, int keep_first,
                     int long_strike_flag,  /*	0: ATM 1: Coupon 2: Eq (PV/Lvl)
                                             */
                     int short_strike_flag, /*	0: ATM  , 1: implied digital
                                               caplet strike 2: same number of
                                               std */

                     int calc_fwd_iv, int adj_fee,

                     /*	Reserve parameters */
                     double lambda_reserve, double lgm_alpha, double lgm_gamma,
                     double lgm_rho, double onetime_vega, int vega_shift_type,

                     int reserve_method,    /*	0: old method  , 1: new method
                                             */
                     int lgm_reserve,       /*	1: calculate the reserve in the
                                               one factor */
                     int midat_reserve,     /*	1: reserve calculated using the
                                               midat-replication */
                     int recalib_european,  /*	If set to 1  , then european are
                                               recalibrated when changing lambda
                                             */
                     int euro_nb_iter,      /*	Number of Levenberg iteration
                                               for calibration of european */
                     int one_time_index,    /*	One Time callable used for
                                               reserve calculation */
                     int recalc_one_factor, /*	Recalculate MCEB in one factor
                                               for reserve */

                     int reserve_nstpt,   /*	number of time steps for reserve
                                             calculation */
                     int reserve_nstpx,   /*	number of space steps for
                                             reserve calculation */
                     int reserve_nstppsi, /*	number of time steps for reserve
                                             calculation */
                     int reserve_nstpvol, /*	number of space steps for
                                             reserve calculation */
                     double reserve_integ_mintime,

                     /*	Results */
                     double *switch_reserve, double *onetime_reserve,

                     CTS_IV fwd_iv_info, int save_extra_infos,
                     CTS_EXTRA_INFOS extra_infos);

#endif