 
#ifndef Fx5FCalib_h
#define Fx5FCalib_h

/*	Fx Implied volatility */
/*  The rates maturity dates have to be merged ! */
Err Fx5DtsImpliedVol(	double	opt_maturity,
						double	start_date,
						double	end_date,
						double	*maturity_rates,
						long	nbMat,
						double	*sig_curve_dom,
						double	lda_dom,
						double	alpha_dom,
						double	gamma_dom,
						double	*sig_curve_for,
						double	lda_for,
						double	alpha_for,
						double	gamma_for,
						double	*maturity_fx,
						double	*sig_curve_fx,
						long	nbrMat_fx, 
						double	*correl[],
						double	*fx_vol);

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx5DtsCalibration(	
					double	*exercise_opt,
					double	*maturity_opt,
					double	*vol_opt,
					long	nbrOpt,
					double	*maturity_rates,
					long	nbrMat,
					double	*sig_curve_dom,
					double	lda_dom,
					double	alpha_dom,
					double	gamma_dom,
					double	*sig_curve_for,
					double	lda_for,
					double	alpha_for,
					double	gamma_for,
					double	**correl,
					double	**fx_vol_curve);

/*	Implied vol direct from underlying */
Err Fx5DImpliedVol(
				   char		*fx_underlying,
				   double	**correl,
				   double	val_time,
				   double	start_time,
				   double	end_time,
				   double	*vol
				   );

Err Fill_5F_correl_matrix(
							double	*maturity_rates,
							long	nbrMat,
							double	*sig_curve_dom,
							double	lda_dom,
							double	alpha_dom,
							double	gamma_dom,
							double	rho_dom,
							double	*sig_curve_for,
							double	lda_for,
							double	alpha_for,
							double	gamma_for,
							double	rho_for,
							double	correl_dom_for,
							double	correl_dom_fx,
							double	correl_for_fx,
							double	**correl);

/*	Fx calibration direct from underlying */
Err Fx5DCalibration(
					char	*dom_underlying,
					char	*for_underlying,
					double	**correl,
					double	*exercise_opt,
					double	*maturity_opt,
					double	*vol_opt,
					long	nbropt,
					double	**fx_vol_curve);

Err Get_FX_StochRate_TermStructures5F(	char	*underlying,
										double	**sigma_date_dom,
										double	**sigma_dom,
										long	*sigma_n_dom,
										double	*fixed_tau_dom,
										double	*fixed_alpha_dom,
										double	*fixed_gamma_dom,
										double	*fixed_rho_dom,
										double	**sigma_date_for,
										double	**sigma_for,
										long	*sigma_n_for,
										double	*fixed_tau_for,
										double	*fixed_alpha_for,
										double	*fixed_gamma_for,
										double	*fixed_rho_for,
										double	**sigma_date_fx,
										double	**sigma_fx,
										long	*sigma_n_fx,
										double	*correl_dom_for,
										double	*correl_dom_fx,
										double	*correl_for_fx);

Err	calibrate_correl5D_from_histo(double	mat,		/* maturity of the long rate */
								double	lam_dom,
								double	alpha_dom,
								double	gamma_dom,
								double	rho_dom,
								int		calib_rho_dom,	/* if one we calibrate rho_dom							*/
								double	lam_for,
								double	alpha_for,
								double	gamma_for,
								double	rho_for,
								int		calib_rho_for,	/* if one we calibrate rho_for							*/
								double	**corr_histo,	/* input of the historical correl SD, LD, SF, LF, Fx	*/
								double	**corr_model);	/* ouput for correl between the brownian of the model	*/

Err	get_correlSL_from_model(	double	mat,		/* maturity of the long rate */
								double	lam_dom,
								double	alpha_dom,
								double	gamma_dom,
								double	lam_for,
								double	alpha_for,
								double	gamma_for,
								double	**corr_model,	/* input of the model correl  between brownian	*/
								double	**corr_sd);		/* ouput for correl between SD, LD, SF, LF, Fx	*/

#endif