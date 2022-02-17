#ifndef __OTCCALLER_H
#define __OTCCALLER_H


/*	Caller for one time callable power duals */
/*	---------------------------------------- */

Err otc_caller(
/*	Today's date */
long		today,
/*	The underlying */
int			use_calib,						/*	0: use fx3dund, 1: calibrate */
/*		if calib */	
double		fx_spot,						/*	2bd fwd */
long		fx_spot_date,
int			dom_calib,						/*	Calibrate domestic underlying */
char		*dom_und,						/*	If no, domestic underlying to be used */
char		*dom_yc,						/*	Domestic yc */
char		*dom_vc,						/*	Domestic vc (only if calib) */
char		*dom_ref,						/*	Domestic ref rate (only if calib) */
char		*dom_swap_freq,					/*	Domestic swap freq (only if calib) */
char		*dom_swap_basis,				/*	Domestic swap basis (only if calib) */
double		dom_lam,						/*	Domestic lambda */
int			for_calib,						/*	Same for foreign */
char		*for_und,
char		*for_yc,
char		*for_vc,
char		*for_ref,
char		*for_swap_freq,
char		*for_swap_basis,
double		for_lam,	

double		min_fact,						/*	Maximum down jump on variance */
double		max_fact,						/*	Maximum up jump on variance */
int			use_jumps,						/*	Allow vol term structure to jump */

double		*corr_times,
double		*correl_dom_for,				/*	Correlations */
double		*correl_dom_fx,
double		*correl_for_fx,
double		**cumulative_corr_matrix,		// Dom-For || Dom-Fx || For-Fx
long		corr_n_times,
CPDBETADLMPARAMS	cpd_dlm_params,
Err			(*get_ir_cash_vol)(				/*	Function to get IR cash vol from the markets */
				char	*vol_curve_name,	
				double	start_date, 
				double	end_date,
				double	cash_strike,
				int		zero,
				char	*ref_rate_name,
				double	*vol,
				double	*power),
/*	Fx vol from the market */
long		*fx_mkt_vol_date,
double		*fx_mkt_vol,
int			num_fx_mkt_vol,
/*	Fx SABR parameters from the market */
double		*fx_mkt_smile_alpha,
double		*fx_mkt_smile_beta,
double		*fx_mkt_smile_rho,
double		*fx_mkt_smile_pi,
int			use_sabr,				/*	0: no smile adj-t, 1: only for underlying, 2: fees adj-t for the call/KO */
int			smile_spec_type,		//	0: lognormal vol + SABR params, 1: sigma-beta + SABR params, 2: BMM (not yet supported)
/*		if no calilb */
char		*fx3dund,
/*	The structure */
long		start_date,				/*	Date at which initial notional exchange occurs */
/*		funding */
double		fund_not,				/*	Notional */
int			fund_ccy,				/*	0: domestic, 1: foreign 2: other */
char		*fund_ccy_yc,			/*	If different from domestic or foreign (fund_ccy = 2) */
double		fx_fund_dom,			/*	If different from domestic or foreign (fund_ccy = 2) 2 bd fwd */
long		fx_fund_dom_spot_date,
int			fund_ncpn,				/*	Number of coupons */
long		*fund_fix,				/*	Fixing dates */
long		*fund_start,			/*	Start dates */
long		*fund_pay,				/*	Pay dates */
char		**fund_basis,			/*	Basis */
double		*fund_spr,				/*	Forward spreads */
double		*fund_mrg,				/*	Margins */
double		*fund_fix_cpn,			/*	Past coupon fixing if relevant,
											includes spr, but not mrg, cvg and notional */
/*		pd */
double		pd_not,					/*	Notional */
int			pd_ncpn,				/*	Number of coupons */
long		*pd_fix,				/*	Fx fixing dates */
long		*pd_start,				/*	Start dates */
long		*pd_pay,				/*	Pay dates */
char		**pd_basis,				/*	Basis */
double		*pd_alpha,				/*	Coupon = alpha + beta * fx [capped, floored] */
double		*pd_beta,
int			*pd_floored,
double		*pd_floor,
int			*pd_capped,
double		*pd_cap,
/*		pd interp coupon specification */
int			*pd_nfxpts,				/*	Number of coupon interpolation points */
double		**pd_fxpts,				/*	Coupon interpolation points */
double		**pd_cpn_at_pts,		/*	Coupon values at interpolation points (interpolation is linear) */
int			*pd_lin_xtrpl_l,		/*	0 = flat extrapolation to the left, 1 = linear w first 2 pts slope */
int			*pd_lin_xtrpl_r,		/*	0 = flat extrapolation to the right, 1 = linear w last 2 pts slope */
double		*pd_fix_fx,				/*	Past Fx fixing if relevant */
/*		pd not refund */
long		*pd_not_ref_fix,			//	fx fixing dates FOR EACH CALL DATE + in the end if relevant
double		pd_not_ref_alpha,			/*	Final notional on PD leg */
double		pd_not_ref_beta,
int			pd_not_ref_floored,
double		pd_not_ref_floor,
int			pd_not_ref_capped,
double		pd_not_ref_cap,
/*		pd not interp specification */
int			*pd_not_ref_nfxpts,			/*	Number of notional interpolation points FOR EACH CALL DATE + in the end */
double		**pd_not_ref_fxpts,
double		**pd_not_ref_cpn_at_pts,
int			*pd_not_ref_lin_xtrpl_l,	/*	0 = flat extrapolation to the left, 1 = linear w first 2 pts slope */
int			*pd_not_ref_lin_xtrpl_r,	/*	0 = flat extrapolation to the right, 1 = linear w last 2 pts slope */
double		*pd_not_ref_fix_fx,			/*	fx fixings FOR EACH CALL DATE + in the end if relevant */
/*		calls */
int			*call_type,				/*	0: call, 1: KO */
int			ncall,					/*	Number of calls */
int			pay_rec,				/*	0: rec pd, 1: pay pd */
long		*ex_date,				/*	Call dates */
long		*set_date,				/*	Settlement dates */
double		*barrier,				/*	in case of a pure KO or a Callable KO */
int			*bar_type,				/*	0: up and in, 1: down and in */
double		*fees,					/*  fees if deal is called in domestic currency */
int			TARN_Do,				//	1: Prices a Target note Powerdual
int			use_GMA,				//	0: Does not output the multicallable price | 1: GMA1 | 1: GMA2 | 3: GMA1 & GMA2
/*	Numerical params */
long		req_stp,				/*	Number of time steps in the tree */
long		req_pth,				/*	Number of paths in the MC */
double		bar_smooth,				/*	Smoothing factor for barriers */
int			do_pecs,				/*	Do PECS in the MC */
int			forcetree,				/*	If equal to 1 then the valuation is done in a tree	*/
int			do_optim,				/*	If equal to 1 then the call are replaced by optimal KO	*/
int			force_optim,			/*	If equal to 1 then all call will be replaced by optimal KO	*/
int			fx_bound,				/*	If equal to 1 then optimisation on the Fx, on the IV otherwise	*/
int			use_bound,				/*	If equal to 1 then prices the call as UO on the Fx using a provided boundary	*/
int			do_infos,				/*	infos on callable right */
/*	EOD Fixing Flag */
int			eod_fix_flag,			/*	0: I, 1: E */
/*	EOD Payment Flag */
int			eod_pay_flag,			/*	0: I, 1: E */
/*	EOD Exercise Flag */
int			eod_ex_flag,			/*	0: I, 1: E */
/*	Exercised flag */
int			exercised,				/*	Flag */
long		ex_date_ex,				/*	Date when exercised */
long		ex_date_set,			/*	Corresponding settlement date */
/*  Parameters */
int			nSimul, 
int			Do_pecs, 
int			nPoints, 
int			nStd,
double		CummulPrecision,
int			CummulLinear,
double		std,
int			nIter,
int			PayoffFunction,
double		fwdSmileVisuNstd,
			/* otc */
int			fwdVolMethod,			/*	0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding Cvg Sbeta; 3=Cvg Cvg Sbeta */
int			smileOtc,				/*	use smile parameters for the OTC */
int			smileFee,				/*	use smile parameters for the OTC Fees */
int			otc,					/*	which call to keep */
int			smileModel,				/*	0=SSL; 1=Log Mix; 2=Beta Mix */
double		BMpi,					/*	if smileModel = 1 or 2 then BMpi is the probability of the state 1 */
/* Correl */
int			firstIndex,
int			secondIndex,
int			firstLong,
int			secondLong,
long		correlTstar,
char		*CPDsigma,
char		*CPDalpha,
char		*CPDbeta,
char		*CPDrho,
/* Change the Funding for speed */
int			FundingSpeedUp,
/* Change the calibration strategy */
int			LongShort,
/* Do not calculate all OTC*/
int			nStart,
int			oneOutOfN,
/*Fast MC*/
int			FMC_do,
double		FMC_precision,
int			FMC_min_paths,
/*	Results */
double		**Results);

Err FxFwdSmilePrice(	long	today,
						long	spot_date,
						double	spot_fx,
						long	forwardFixDate,
						double	forwardFixTime,
						long	forwardValDate,
						double	forwardValTime,
						long	fixDate,
						double	fixTime,
						long	valDate,
						double	valTime,
						char	*dom_yc,
						char	*for_yc,
						double	*sigma_time_dom, 
						int		sigma_n_dom,
						double	*sigma_dom,
						double	lda_dom,
						double	*sigma_time_for, 
						int		sigma_n_for,
						double	*sigma_for, 
						double	lda_for,
						double	*sigma_time_fx, 
						double	*sigma_fx, 
						int		sigma_n_fx,
						double	*corr_times, 
						double	*correl_dom_for, 
						double	*correl_dom_fx, 
						double	*correl_for_fx,
						int		corr_n_times,
						double	smile_vol,
						double	smile_alpha,
						double	smile_beta,
						double	smile_rho,
						double	SigmaFwd,
						double	AlphaFwd,
						double	BetaFwd,
						double	RhoFwd,
						int		nSimul,
						int		nPoints,
						int		nIter,
						int		nStd,
						double	std,
						int		smileModel,
						int		nStrikes,
						double	*Strikes,
						int		isCall,
						int		OutputVol,
						int		StochVolFudge,
						double	SwitchLevel,
						double	Proba,
						double	AlphaFudge,
						double	*Results);

Err otc_pricer(	CPD_STR				cpd,
				CPD_UND				und,
				otcpd_params		*OTC_params,
				otcpd_precalc		*OTC_precalc,
				SMILE_VOL_MARKET	smile_mkt,
				SMILE_PARAMETERS	smile_params,
				double				pd_not,
				int					mergeNtimes,
				double				*mergeTimes,
				double				*sigDom, 
				double				*sigFor, 
				double				*sigFx, 
				double				*corDF, 
				double				*corDFx, 
				double				*corFFx,
				long				start_date,				/*	Date at which initial notional exchange occurs */
				double				**cumulative_corr_matrix,
				int					erasing_call_done,		//	Don't pass to GMA all the OTC
				int					*erased_call_list,		//	List of erasedz
				/*	Results */
				double				**Results);

Err ko_pricer(	CPD_STR				cpd,
				CPD_UND				und,
				otcpd_params		*OTC_params,
				otcpd_precalc		*OTC_precalc,
				SMILE_VOL_MARKET	smile_mkt,
				SMILE_PARAMETERS	smile_params,
				double				pd_not,
				int					mergeNtimes,
				double				*mergeTimes,
				double				*sigDom, 
				double				*sigFor, 
				double				*sigFx, 
				double				*corDF, 
				double				*corDFx, 
				double				*corFFx,
				/*	The structure */
				long		start_date,				/*	Date at which initial notional exchange occurs */
				/*	Results */
				double		**Results);
#endif
