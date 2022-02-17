#ifndef __DOUBLE_TRIGGER_SWAP_QUANTO_LEG_H
#define __DOUBLE_TRIGGER_SWAP_QUANTO_LEG_H

#define LINEAR_METHOD 0
#define PIECEWISE_CONSTANT_METHOD 1

#define EXOTIC_LEG_PAY 0
#define EXOTIC_LEG_REC 1

// Prices a double put with payoff E[max(strikeX - X_fixed_at_timeX,0) * max(strikeY - Y_fixed_at_timeY,0)
// where X and Y are both lognormals with fwds fwdX and fwdY
// and vols sigmaX and sigmaY and instantaneous correlation rhoXY
// X is fixed at timeX and Y is fixed at time Y.
Err   double_put_lognormal(double	x0,
						   double	B,
						   double	sigma_x,
						   double	T_x,
						   double	y0,
						   double	K,
						   double	sigma_y,
						   double	T_y,
						   double	correl,
						   double*	result);

//Calculates the PV of a double trigger swap quanto leg
Err   double_trigger_swap_quanto_leg(
									 char*				dom_yc_name,
									 char*				dom_vc_name,
									 int				dom_spot_lag,		   
									 double*			dom_strikes_in_vol,
									 int				n_dom_strikes_in_vol,
									 SrtDiffusionType	dom_vol_type,
									 int				dom_cash_vol,

									 char*				for_yc_name,
									 char*				for_vc_name,
									 int				for_spot_lag,
									 double*			for_strikes_in_vol,
									 int				n_for_strikes_in_vol,
									 SrtDiffusionType	for_vol_type,
									 int				for_cash_vol,

									 long				today,
									 
									 char*				cms1_tenor,
									 char*				cms1_frequency,
									 char*				cms1_basis,
									 char*				cms1_refrate,
									 int				cms1_dom_for,

									 char*				cms2_tenor,
									 char*				cms2_frequency,
									 char*				cms2_basis,
									 char*				cms2_refrate,
									 int				cms2_dom_for,

									 long*				exo_start_dates,
									 long*				exo_end_dates,
									 long*				exo_pay_dates,
									 char**				exo_basis,
									 double*			exo_notionals,
									 long*				exo_cms1_fixing_dates,
									 double*			exo_cms1_barriers,
									 long*				exo_cms2_fixing_dates,
									 double*			exo_cms2_barriers,
									 double*			exo_m12,
									 double*			exo_g2,
									 double*			exo_m2,
									 double*			exo_g3,
									 double*			exo_m3,
									 double*			exo_g4_X,
									 double*			exo_g4_Y,
									 double*			exo_m4,
									 double*			exo_cms1_past_fixings,
									 double*			exo_cms2_past_fixings,
									 int				n_exo_coupons,

									 double*			cms1_cms2_corr_times,
									 double*			cms1_cms2_corr_values_bid,
									 double*			cms1_cms2_corr_values_ask,
									 int				n_cms1_cms2_corr,

									 double*			cms1_quanto_corr_times,
									 double*			cms1_quanto_corr_values,
									 int				n_cms1_quanto_corr,

									 double*			cms2_quanto_corr_times,
									 double*			cms2_quanto_corr_values,
									 int				n_cms2_quanto_corr,

									 double*			fwd_fx_vols_times,
									 double*			fwd_fx_vols_values,
									 int				n_fwd_fx_vols,

									 int				pay_rec,
									 double				call_spread,
									 int				use_CMSRate_CMSTEC,
									 int				eod_fix_flag,
									 int				eod_pay_flag,

									 double*			cpn1,
									 double*			cpn2,
									 double*			cpn3,
									 double*			cpn4,
									 double*			DoubleTriggerSwapQuantoLeg);


/*
Err   frequency_and_basis_to_compounding_and_dRateConv(char* frequency,
													   char* basis,
												       SrtCompounding* compounding,
												       double* dRateConv);

Err   cms_get_tenor_num_and_unit_from_tenor(char* cms_tenor,
											int* tenor_num,
											char* tenor_unit);

Err	  constant_or_linear_interpolation(double* dates,
								 	   double* values,
									   int n_dates,
									   int method,
									   double date_wanted,
									   double* value_wanted);

Err   cms_get_end_date_from_tenor(long cms_start_date,
								  char* cms_tenor,
								  long* cms_end_date);

Err digital_quanto_float_leg_get_lognormal_vol(char* vc_name,
											   double fwd,
											   double fixing_time,
											   long start_date,
											   long end_date,
											   double strike,
											   double* result_volatility);

Err   digital_quanto_float_leg_cms_rate_quanto(double fwd,
											   double cms_fixing_time,
											   double cms_start_time,
											   double cms_nfp,

											   SrtCompounding cms_compounding,
											   double cms_delay,
											   double cms_dRateConv,

											   SrtDiffusionType cms_vol_type,
											   int cms_flat_vol_or_smile_approx,
											   long cms_start_date,
											   long cms_theo_end_date,
											   int cms_cash_vol,
											   double fwd_spread,
											   char* cms_vc_name,
											   int n_cms_strikes_in_vol,
											   double* cms_strikes_in_vol,

											   int cms_dom_for,

											   double* cms_quanto_corr_times,
											   double* cms_quanto_corr_values,
											   int n_cms_quanto_corr,

											   double* fwd_fx_vols_times,
											   double* fwd_fx_vols_values,
											   int n_fwd_fx_vols,

											   double* result_cms_rate_quanto);



  */
#endif
