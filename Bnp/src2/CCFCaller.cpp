
#include "srt_h_all.h"
#include "srtaccess.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "opfnctns.h"
#include "math.h"
#include "RainbowOpt.h"

#define	NUM_HERMITE	10

static double dens(double x)
{
	return INV_SQRT_TWO_PI * exp(-x*x/2.0);
}

static double cms_rate(				 
						 double	delay,
						 double	nfp,
						 double	compd,
						 double	forward,
						 double	nvar)
{
static double	lvl,
				dlvl,
				cms_adj,
				del_adj,
				tot_adj;

	lvl = (1.0 - pow (1.0 + forward / compd, - nfp)) / (forward / compd);
	dlvl = 
		(nfp / compd) * pow (1.0 + forward / compd, - nfp - 1) * forward / compd
		- (1.0 - pow (1.0 + forward / compd, - nfp)) / compd;
	dlvl /= (forward / compd) * (forward / compd);

	/*	Adjustments	*/
	cms_adj = - dlvl/lvl;
	del_adj = - delay / (1.0 + forward * delay);
	tot_adj = cms_adj + del_adj;

	return (forward + tot_adj * nvar);
}

#define POS_VAL(X)							((X) > 0? (X): 0)
	
#define CALL_VAL_N(FWD, STRIKE, STD, D)			(((FWD) - (STRIKE))* norm ((D)) + (STD) * dens ((D)))

#define PUT_VAL_N(FWD, STRIKE, STD, D)			(((STRIKE) - (FWD))* norm (-(D)) + (STD) * dens ((D)))

#define OPT_VAL_MACRO_N(TYPE, FWD, STRIKE, STD) (													\
	(TYPE) == 0?																					\
		0.0:(																						\
			(TYPE) == 1?																			\
			POS_VAL ((FWD) - (STRIKE)):(															\
				(TYPE) == 2?																		\
					POS_VAL ((STRIKE) - (FWD)):(													\
					(TYPE) == 3?																	\
						CALL_VAL_N ((FWD), (STRIKE), (STD), ((FWD) - (STRIKE))/(STD)):			\
							PUT_VAL_N ((FWD), (STRIKE), (STD), ((FWD) - (STRIKE))/(STD))))))	

Err ccf_caller(
			/*	Today's date */
			long		today,
			/*	The underlying */
			int			use_calib,						/*	0: use lgm2fund, 1: calibrate */
			/*		if calib */
			char		*yc,							/*	yc */
			char		*vc,							/*	vc */
			char		*ref,							/*	ref rate (only if calib) */
			char		*swap_freq,						/*	swap freq (only if calib) */
			char		*swap_basis,					/*	swap basis (only if calib) */
			int			lam_ts,							/*	0: use unique lambda, 1: use ts */
			double		lambda,							/*	lambda if unique */
			int			tsnlam,							/*	number of lambdas if ts */
			double		*tslamtime,						/*	lambda times i.e. (date - today) / 365 if ts */
			double		*tslam,							/*	corresponding lambdas if ts */
			double		alpha,							/*	alpha */
			double		gamma,							/*	gamma */
			double		rho,							/*	rho */
			/*	End of calib params */
			Err			(*get_cash_vol)(				/*	function to get IR cash vol from the markets */
							char	*vol_curve_name,	
							double	start_date, 
							double	end_date,
							double	cash_strike,
							int		zero,
							char	*ref_rate_name,
							double	*vol,
							double	*power),
			/*		if no calilb */
			char		*lgm2dund,
			/*	The structure */
			long		start_date,				/*	Date at which initial notional exchange occurs */
			/*		funding */
			int			fund_ccy,				/*	0: domestic, 1: other */
			double		fund_not,				/*	If different from domestic or foreign (fund_ccy = 1) */
			char		*fund_ccy_yc,			/*	If different from domestic or foreign (fund_ccy = 1) */
			double		fx_fund_dom,			/*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
			long		fx_fund_dom_spot_date,
			int			fund_ncpn,
			long		*fund_fix,
			long		*fund_start,
			long		*fund_pay,
			char		**fund_basis,
			double		*fund_spr,
			double		*fund_mrg,
			double		*fund_fix_cpn,			/*	Past coupon fixing if relevant,
											includes spr, but not mrg, cvg and notional */
			/*		cf */
			double		cf_not,
			int			cf_ncpn,
			long		*cf_fix,
			long		*cf_start,
			long		*cf_pay,
			char		**cf_basis,
			char		**cf_cms_tenor1,
			char		*cf_cms_freq1,
			char		*cf_cms_basis1,
			double		*cf_cms_spread1,
			char		**cf_cms_tenor2,
			char		*cf_cms_freq2,
			char		*cf_cms_basis2,
			double		*cf_cms_spread2,
			long		spread_vol_n,
			double		*spread_vol_time,
			double		*spread_vol_floor,
			double		*spread_vol_cap,
			
			double		*spread_slbeta1,		/* Shifted log beta on the CMS1 */
			double		*spread_slbeta2,		/* Shifted log beta on the CMS2 */

			double		*cf_alpha,
			double		*cf_beta,
			double		*cf_gamma,
			int			*cf_floored,
			double		*cf_floor,
			int			*cf_capped,
			double		*cf_cap,
			double		*cf_fix_cpn,			/*	Past coupon fixing if relevant */

			int			cf_nopt,			/* Number of spread options */
			double		**cf_notopt,		/* Notional of the spread options */	
			double		**cf_strikeopt,		/* spread option strikes */
			int			**cf_typeopt,		/* spread option type 0 call 1 put */

			/*		calls */
			int			ncall,
			int			pay_rec,				/*	0: rec pd, 1: pay pd */
			long		*ex_date,
			long		*set_date,
			double		*fee,
			/*	Numerical params */
			int			req_stp,
			int			req_stpx,
			long			req_paths,
			/*	Calib params */
			int			force_atm,
			double		max_std_long,
			double		max_std_short,
			int			fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
																	to diagonal */
			int			cal_vol_shift_type,				/*	vol shift type for volatility */
			double		cal_vol_shift,					/*	vol shift */
			double		cal_lambda_shift,				/*	shift on lambda after calibration */
			int			one_f_equi,						/*	1F equivalent flag:
																	if set to 1, then 2F lambda will calibrate
																	to the cap priced within calibrated 1F
																	with the given lambda */
			int			skip_last,				/*	If 1, the last option is disregarded and the forward volatility is flat from option n-1 */
			double		long_prec,				/*	Precision on primary instruments */
			double		short_prec,				/*	Precision on secondary instruments */
			double		min_fact,				/*	Maximum down jump on variance */
			double		max_fact,				/*	Maximum up jump on variance */
			int			use_jumps,				/*	Allow vol term structure to jump */
			int			proba_weight,
			

			int			calc_fwdiv,				/*	output the model and market fwdiv  */
			int			adjust_fee,				/*	adjust the fee to price correctly the fwd iv */

			/*	EOD Flags */
			int			eod_fix_flag,			/*	0: I, 1: E */
			int			eod_pay_flag,			/*	0: I, 1: E */
			int			eod_ex_flag,			/*	0: I, 1: E */
			/*	CMS can be valued with smile */
			int			cms_adj,				/*	1: adjust for CMS effect
													0: don't */
			int			cms_for_iv_only,		/*	1: adjust for CMS effect only IV
													0: adjust both IV and call */
			int			use_cms_smile,			/*	1: value CMS with ATM
													0: value CMS with smile */
			int			cms_vol_adj,			/*	1: adjust for CMS vol effect
													0: don't */
			double		cms_beta1,				/*	How adjustment varies with rates */
			double		cms_beta2,
			int			num_strikes_in_vol,		/*	Array of strikes in vol matrix */
			double		*strikes_in_vol,
			SrtDiffusionType vol_type,			/*	Type of vol in matrix, SRT_NORMAL or SRT_LOGNORMAL */
			int			cash_vol,				/*	1: matrix is a cash vol
													0: matrix is a swap vol */
			int			cvg_sv,					/*	1: coverging model on spread vol/correlation
													0: sliding model on spread vol/correlation */
			int			is_corr,				/*	1: spread vols are correlations
													0: spread vols are vols */

			int			use_SL,					/*  1: use Shifted Log 
												    0: don't use				*/
			int			calib_SL,				/*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
													0: use the Shifted Log beta given */
			
			double		NbStdcalib_SL,			/*  Nb of std for the calibration of to the skew */

			int			calib_correl_SL,		/*	1: Calibrate the correlation between the two SL to get the same ATM
														normal spread vol 
													0: use the normal spread correl for the sl correl */
			
			int			use_cfoptions,			/*  1: Use the spread options and don't take into account 
														the floor and cap in the cf coupons
													0: take into account the floor and cap in the cf coupons */

			/*	Exercised flag */
			int			exercised,				/*	Flag */
			long		ex_date_ex,				/*	Date when exercised */
			long		ex_date_set,			/*	Corresponding settlement date */
			double		ex_fee,					/*	Corresponding fee */
			/*	Results */
			double		*fund_val,				/*	Value of the funding leg */
			double		*cf_val,				/*	Value of the Power Dual leg */
			double		*call_val,				/*	Value of the callable feature */			
			int			export_ts,				/*	1: Export TS, 0: don't */
			CCF_UND		und_exp)
{
	ccf_str			*ccf		= NULL;
	ccf_und			*und		= NULL;
	ccf_adi_arg		*adi_arg	= NULL;	

	CF_FUND_LEG		fund_leg;
	CF_FUND_CPN		fund_cpn;
	CF_EXO_LEG		exo_leg;
	CF_EXO_CPN		exo_cpn;
	CF_CMS_DESC		cms;

	int				call_feat;
	double			fund_leg_pv, exo_leg_pv;	
	SrtBasisCode	bas;	

	int				floortype, captype;
	double			df, stdfloor, stdcap, floor, cap, coupon;

	double			call;
	int				i, j;
	int				free_struct = 0;
	
	double			temp;

	double			vol, spread, forward, cmsrate[2];

	int				for_fund;
	long			fund_start_date, fin_not_date;
	double			eq_final_ex, eq_init_ex;
	
	double			tmp1, tmp2;

	CF_CALL			callt;
	double			fact;

	Err				err			= NULL;

	double			OtherCcyfundNot;
	char			Copyfund_ccy_yc[255];

	int				iNumOpt;
	double			spread_corr, spread_slcorr;
	double			optionvalue;

	
	/* Memory Allocation */
	ccf = calloc(1, sizeof(ccf_str));
	und = calloc(1, sizeof(ccf_und));
	adi_arg = calloc(1, sizeof(ccf_adi_arg));

	if (!ccf || !und || !adi_arg)
	{
		err = "memory allocation failure in ccf_caller";
		goto FREE_RETURN;
	}
		
	/*	If exercised */
	if (exercised)
	{
		/* Consider only the coupons which start before the ex date */
		i = 0;
		while (i < cf_ncpn && cf_start[i] < ex_date_ex)
		{
			i++;
		}
		cf_ncpn= i;

		i = 0;
		while (i < fund_ncpn && fund_start[i] < ex_date_ex)
		{
			i++;
		}
		fund_ncpn = i;

		/* Checks */
		/*if ((cf_ncpn > 0) && (ex_date_set<cf_pay[cf_ncpn-1]))
		{
			err = "Exercised settlement should be > than payment date of the last not called cf coupon ";
			goto FREE_RETURN;
		}

		if ((fund_ncpn > 0) && (ex_date_set<fund_pay[fund_ncpn-1]))
		{
			err = "Exercised settlement should be > than payment date of the last not called funding coupon ";
			goto FREE_RETURN;
		}
		*/

		/* Gestion of the particular case when the Ex settlement is past *
		 * or all the coupons (cf and funding) are called				 */
		if ((ex_date_set < today + eod_pay_flag)||((cf_ncpn == 0)&&(fund_ncpn == 0)))
		{	
			/* Fill the und */
			strcpy (und->name, "CALIB");
			und->today = today;
			strcpy (und->yc, yc);
			strcpy (und->vc, vc);
			strcpy (und->ref, ref);
			strcpy (und->swap_freq, swap_freq);
			strcpy (und->swap_basis, swap_basis);
			
			und->sigma_n = 0;
			und->lambda_n = 0;
			und->spread_vol_n = 0;

			/* Fill the output */
			*fund_val = *cf_val = *call_val = 0.0;

			/* Initial exchange of notional */
			if (start_date >= today + eod_pay_flag)
			{
				if (fund_ccy)
				{
					*fund_val -= swp_f_df (today, start_date, fund_ccy_yc) * fund_not 
						* fx_fund_dom * swp_f_df (today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (today, fx_fund_dom_spot_date, fund_ccy_yc);
					*cf_val -= cf_not * swp_f_df (today, start_date, yc);
				}
			}
			
			/* Final exchange of notional */
			if (ex_date_set >= today + eod_pay_flag)
			{
				if (fund_ccy)
				{
					/* Final exchange of notional */
					*fund_val += swp_f_df (today, ex_date_set, fund_ccy_yc) * fund_not 
						* fx_fund_dom * swp_f_df (today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (today, fx_fund_dom_spot_date, fund_ccy_yc);
					*cf_val += cf_not * swp_f_df (today, ex_date_set, yc);
					*call_val += -ex_fee * swp_f_df (today, ex_date_set, yc);
				}
			}

			/* export the ts if asked */
			if (export_ts)
			{
				ccf_copy_und (und, und_exp);
			}

			return NULL;
		}
		
		/*  those particular cases should not happen */
		if ((cf_ncpn > 0) && (fund_ncpn == 0))
		{
			err = "All the funding coupons are called but not all the cf coupons : Not allowed ";
			goto FREE_RETURN;
		}
		if ((cf_ncpn == 0) && (fund_ncpn > 0))
		{
			err = "All the cf coupons are called but not all the funding coupons : Not allowed ";
			goto FREE_RETURN;
		}
		
		/*  Case where (cf_ncpn > 0) && (fund_ncpn > 0) */
		ncall = 0;
	}

	/*	If cms adjustment only for IV */
	if (!exercised && cms_adj && cms_for_iv_only)
	{
		/* First, fo the non_callable */
		err = ccf_caller(
			today,
			use_calib,
			yc,
			vc,
			ref,
			swap_freq,
			swap_basis,
			lam_ts,
			lambda,
			tsnlam,
			tslamtime,
			tslam,
			alpha,
			gamma,
			rho,
			get_cash_vol,
			lgm2dund,
			start_date,
			fund_ccy,
			fund_not,
			fund_ccy_yc,
			fx_fund_dom,
			fx_fund_dom_spot_date,
			fund_ncpn,
			fund_fix,
			fund_start,
			fund_pay,
			fund_basis,
			fund_spr,
			fund_mrg,
			fund_fix_cpn,
			cf_not,
			cf_ncpn,
			cf_fix,
			cf_start,
			cf_pay,
			cf_basis,
			cf_cms_tenor1,
			cf_cms_freq1,
			cf_cms_basis1,
			cf_cms_spread1,
			cf_cms_tenor2,
			cf_cms_freq2,
			cf_cms_basis2,
			cf_cms_spread2,
			spread_vol_n,
			spread_vol_time,
			spread_vol_floor,
			spread_vol_cap,

			spread_slbeta1,
			spread_slbeta2,

			cf_alpha,
			cf_beta,
			cf_gamma,
			cf_floored,
			cf_floor,
			cf_capped,
			cf_cap,
			cf_fix_cpn,

			cf_nopt,
			cf_notopt,
			cf_strikeopt,
			cf_typeopt,

			0,
			pay_rec,
			ex_date,
			set_date,
			fee,
			req_stp,
			req_stpx,
			req_paths,
			force_atm,
			max_std_long,
			max_std_short,
			fix_lambda,				
			cal_vol_shift_type,
			cal_vol_shift,
			cal_lambda_shift,
			one_f_equi,
			skip_last,
			long_prec,
			short_prec,
			min_fact,
			max_fact,
			use_jumps,
			proba_weight,
			calc_fwdiv,
			adjust_fee,
			eod_fix_flag,
			eod_pay_flag,
			eod_ex_flag,
			1,
			0,
			use_cms_smile,
			cms_vol_adj,
			cms_beta1,
			cms_beta2,
			num_strikes_in_vol,
			strikes_in_vol,
			vol_type,
			cash_vol,
			cvg_sv,
			is_corr,

			use_SL,
			calib_SL,
			NbStdcalib_SL,
			calib_correl_SL,
			use_cfoptions,

			exercised,
			ex_date_ex,
			ex_date_set,
			ex_fee,
			fund_val,
			cf_val,
			call_val,
			export_ts,
			und_exp);
		if (err)
		{
			goto FREE_RETURN;
		}

		/*	Then value the call feature */
		err = ccf_caller(
			today,
			use_calib,
			yc,
			vc,
			ref,
			swap_freq,
			swap_basis,
			lam_ts,
			lambda,
			tsnlam,
			tslamtime,
			tslam,
			alpha,
			gamma,
			rho,
			get_cash_vol,
			lgm2dund,
			start_date,
			fund_ccy,
			fund_not,
			fund_ccy_yc,
			fx_fund_dom,
			fx_fund_dom_spot_date,
			fund_ncpn,
			fund_fix,
			fund_start,
			fund_pay,
			fund_basis,
			fund_spr,
			fund_mrg,
			fund_fix_cpn,
			cf_not,
			cf_ncpn,
			cf_fix,
			cf_start,
			cf_pay,
			cf_basis,
			cf_cms_tenor1,
			cf_cms_freq1,
			cf_cms_basis1,
			cf_cms_spread1,
			cf_cms_tenor2,
			cf_cms_freq2,
			cf_cms_basis2,
			cf_cms_spread2,
			spread_vol_n,
			spread_vol_time,
			spread_vol_floor,
			spread_vol_cap,

			spread_slbeta1,
			spread_slbeta2,

			cf_alpha,
			cf_beta,
			cf_gamma,
			cf_floored,
			cf_floor,
			cf_capped,
			cf_cap,
			cf_fix_cpn,

			cf_nopt,
			cf_notopt,
			cf_strikeopt,
			cf_typeopt,

			ncall,
			pay_rec,
			ex_date,
			set_date,
			fee,
			req_stp,
			req_stpx,
			req_paths,
			force_atm,
			max_std_long,
			max_std_short,
			fix_lambda,				
			cal_vol_shift_type,
			cal_vol_shift,
			cal_lambda_shift,
			one_f_equi,
			skip_last,
			long_prec,
			short_prec,
			min_fact,
			max_fact,
			use_jumps,
			proba_weight,
			calc_fwdiv,
			adjust_fee,
			eod_fix_flag,
			eod_pay_flag,
			eod_ex_flag,
			0,
			0,
			use_cms_smile,
			cms_vol_adj,
			cms_beta1,
			cms_beta2,
			num_strikes_in_vol,
			strikes_in_vol,
			vol_type,
			cash_vol,
			cvg_sv,
			is_corr,

			use_SL,
			calib_SL,
			NbStdcalib_SL,
			calib_correl_SL,
			use_cfoptions,

			exercised,
			ex_date_ex,
			ex_date_set,
			ex_fee,
			&tmp1,
			&tmp2,
			call_val,
			export_ts,
			und_exp);
			
		goto FREE_RETURN;
	}

	if (fund_ccy == 1)
	{
		/* save the initial fund not in the third ccy 
		   and the third Ccy yield curve name			*/
		OtherCcyfundNot = fund_not;
		strcpy (Copyfund_ccy_yc, fund_ccy_yc);

		/* Convert into domestic */
		fund_ccy = 0;
		for_fund = 1;
		err = convert_funding_to_domestic(
								today,
								start_date,
								eod_fix_flag,
								eod_pay_flag,
								fx_fund_dom,
								fx_fund_dom_spot_date,
								cf_not,
								yc,
								fund_ncpn,
								fund_fix,
								fund_start,
								fund_pay,
								fund_basis,
								fund_ccy_yc,		
								&fund_not,
								fund_spr,
								fund_mrg,
								fund_fix_cpn,
								&fund_start_date,
								&eq_final_ex,
								&eq_init_ex);
		if (err)
		{
			return err;
		}
	}
	else
	{
		for_fund = 0;
	}

	/*	Initialise structures */
	free_struct = 0;
	err = ccf_fill_check_all_struct(
		today,
		use_calib,
		yc,
		vc,
		ref,
		swap_freq,
		swap_basis,
		lam_ts,
		lambda,
		tsnlam,
		tslamtime,
		tslam,
		alpha,
		gamma,
		rho,
		get_cash_vol,
		lgm2dund,
		fund_not,		
		fund_ncpn,
		fund_fix,
		fund_start,
		fund_pay,
		fund_basis,
		fund_spr,
		fund_mrg,
		cf_not,
		cf_ncpn,
		cf_fix,
		cf_start,
		cf_pay,
		cf_basis,
		cf_cms_tenor1,
		cf_cms_freq1,
		cf_cms_basis1,
		cf_cms_spread1,
		cf_cms_tenor2,
		cf_cms_freq2,
		cf_cms_basis2,
		cf_cms_spread2,
		spread_vol_n,
		spread_vol_time,
		spread_vol_floor,
		spread_vol_cap,

		spread_slbeta1,		/* Shifted log beta on the CMS1 */
		spread_slbeta2,		/* Shifted log beta on the CMS2 */

		cvg_sv,
		is_corr,
		cf_alpha,
		cf_beta,
		cf_gamma,
		cf_floored,
		cf_floor,
		cf_capped,
		cf_cap,	
		
		cf_nopt,		/* Number of spread options */
		cf_notopt,		/* Notional of the spread options */	
		cf_strikeopt,	/* spread option strikes */
		cf_typeopt,		/* spread option type 0 call 1 put */
		use_SL,					/*  1: use Shifted Log 
									0: don't use				*/
		calib_SL,				/*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
									0: use the Shifted Log beta given */
		NbStdcalib_SL,			/*  Nb of std for the calibration of to the skew */
		calib_correl_SL,		/*	1: Calibrate the correlation between the two SL to get the same ATM
									   normal spread vol 
									0: use the normal spread correl for the sl correl */
			
		use_cfoptions,			/*  1: Use the spread options and don't take into account 
														the floor and cap in the cf coupons
									0: take into account the floor and cap in the cf coupons */

		cms_adj,
		use_cms_smile,
		cms_vol_adj,
		cms_beta1,
		cms_beta2,
		num_strikes_in_vol,
		strikes_in_vol,
		vol_type,
		cash_vol,
		ncall,
		pay_rec,
		ex_date,
		set_date,
		fee,
		req_stp,
		req_stpx,	
		req_paths,
		force_atm,
		max_std_long,
		max_std_short,
		fix_lambda,				
		cal_vol_shift_type,
		cal_vol_shift,
		cal_lambda_shift,
		one_f_equi,
		skip_last,
		long_prec,
		short_prec,
		min_fact,
		max_fact,
		use_jumps,
		proba_weight,
		eod_fix_flag,
		eod_ex_flag,
		ccf,
		und,
		&call_feat,
		adi_arg);

	if (err)
	{
		goto FREE_RETURN;
	}
	free_struct = 1;

	/*  0) calculate the fwd iv in the model */

	if (calc_fwdiv && ccf->num_calls > 0)
	{
		und->has_fwd_iv = 1;
		und->nb_fwdiv = ccf->num_calls;
		
		und->exercise_date = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->market_fwdiv = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->model_fwdiv = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->extra_fees = (double *) calloc(und->nb_fwdiv, sizeof(double));

		if (!und->exercise_date || !und->market_fwdiv || !und->model_fwdiv || !und->extra_fees)
		{
			err = "Memory allocation faillure in ccf caller";
			goto FREE_RETURN;
		}

		err = ccf_calc_mdl_iv_fwd(	ccf,
							und,
							adi_arg,
							NUM_HERMITE,
							und->model_fwdiv);

		if (err)
		{
			goto FREE_RETURN;
		}

		if (pay_rec == 1)
		{
			fact = 1.0;
		}
		else
		{
			fact = -1.0;
		}

		for (i=0; i<ccf->num_calls; i++)
		{
			callt = ccf->call + i;
			und->exercise_date[i] = callt->ex_date;
		}
	}
		
	/*	1)	Value funding leg */			

	fund_leg = ccf->fund_leg;
	fund_leg_pv = 0.0;	


	/*	Cash libor */
	if (fund_leg->num_cpn > 0)
	{
		if (for_fund)
		{
			fund_leg_pv += swp_f_df (today, fund_leg->cpn[0].start_date, yc)
				* eq_final_ex;

			if (calc_fwdiv)
			{
				/* initialisation of the initial notional */
				for (i=0; i<ccf->num_calls; i++)
				{
					callt = ccf->call + i;
					und->market_fwdiv[i] = fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
					* eq_final_ex;
				}
			}
		}
		else
		{
			fund_leg_pv += swp_f_df (today, fund_leg->cpn[0].start_date, yc)
				* fund_leg->notional;

			if (calc_fwdiv)
			{
				/* initialisation of the initial notional */
				for (i=0; i<ccf->num_calls; i++)
				{
					callt = ccf->call + i;
					und->market_fwdiv[i] = fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
					* fund_leg->notional;
				}
			}
		}
	}
	else
	{
		fin_not_date = fund_pay[fund_ncpn-1];

		if (fin_not_date >= today + eod_pay_flag)
		{
			if (for_fund)
			{
				temp = swp_f_df (today, fin_not_date, yc) * eq_final_ex;
				fund_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += fact * temp;
					}
				}
			}
			else
			{
				/*fund_leg_pv += swp_f_df (today, fin_not_date, yc) * fund_leg->notional;*/
			}
		}	
	}

	/*	Coupons: spread + margin */
	for (i=0; i<fund_leg->num_cpn; i++)
	{
		fund_cpn = fund_leg->cpn + i;
		temp = swp_f_df (today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
		fund_leg_pv += temp;

		if (calc_fwdiv)
		{
			j = 0;
			while (j < ccf->num_calls && (ccf->call + j)->fund_idx <= i)
			{
				und->market_fwdiv[j] += fact * temp;
				j++;
			}
		}
	}

	/*	Notional exchange */
	if (for_fund)
	{
		if (start_date >= today + eod_pay_flag)
		{
			fund_leg_pv -= swp_f_df (today, start_date, yc) * eq_init_ex;
		}

		if (calc_fwdiv)
		{
			for (i=0; i<ccf->num_calls; i++)
			{
				callt = ccf->call + i;
				und->market_fwdiv[i] -= fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
				* eq_final_ex;
			}
		}
	}
	else
	{
		if (fund_leg->num_cpn > 0)
		{
			temp = swp_f_df (today, fund_leg->cpn[fund_leg->num_cpn-1].pay_date, yc) 
			* fund_leg->notional;
			fund_leg_pv -= temp;

			if (calc_fwdiv)
			{				
				for (i=0; i<ccf->num_calls; i++)
				{
					und->market_fwdiv[i] -= fact * temp;
				}
			}
		}
	}

	/*	PV of coupons fixed in the past and not yet paid */
	i = 0;
	while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
	{
		if (fund_pay[i] >= today + eod_pay_flag)
		{
			err = interp_basis (fund_basis[i], &bas);
			if (err)
			{
				goto FREE_RETURN;	
			}
			
			fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i])
				* coverage (fund_start[i], fund_pay[i], bas)
				* fund_not
				* swp_f_df (today, fund_pay[i], yc);
		}

		i++;
	}

	/* Exerciced Case */
	/* Force the Notional to be paid at ex_settle_date */
	if ((exercised) && (fund_pay[fund_ncpn-1] != ex_date_set))
	{
		/* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at ex_settle_date */
		if (for_fund)
		{
			/* Other currency case */
			fund_leg_pv += (swp_f_df (und->today, ex_date_set, Copyfund_ccy_yc) - 
						swp_f_df (und->today, fund_pay[fund_ncpn-1], Copyfund_ccy_yc)) *
						OtherCcyfundNot 
						* fx_fund_dom * swp_f_df (und->today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);
		}
	}
		
	/*	2)	Value cf leg */

	exo_leg = ccf->cf_leg;
	exo_leg_pv = 0.0;

	
	/*	Coupons */

	for (i=0; i<exo_leg->num_cpn; i++)
	{
		exo_cpn = exo_leg->cpn + i;

		/*	Discount */
		df = swp_f_df (today, exo_cpn->pay_date, yc);
		
		/*	CMS */
		spread = 0.0;
		for (j=0; j<exo_cpn->ncms; j++)
		{
			cms = &(exo_cpn->cms[j]);
			
			if (exo_cpn->cms_fix_date > today)
			{
				vol = cms->atmvar * exo_cpn->cms_fix_time;
			}
			else
			{
				vol = 0.0;
			}
			forward = cms->fwd + cms->swap_spread;

			if (cms_adj && vol > 1.0E-08)
			{
				if (use_cms_smile)
				{
					cmsrate[j] = cms->fwd_cms;					
				}
				else
				{
					cmsrate[j] = cms_rate(
						 cms->delay,
						 cms->nfp,
						 cms->cpnd,
						 forward,
						 vol);
				}
			}
			else
			{
				cmsrate[j] = forward;
			}

			spread += exo_cpn->alphabeta[j] * cmsrate[j];
		}

	
		/* Case no use of string of options */
		if (!use_cfoptions )
		{
			if (!use_SL)
			{
				coupon = spread + exo_cpn->gamma;

				/*	Spread Vol */
				
				switch (exo_cpn->type)
				{
				case 0:
				default:
					/*	Midat */
					break;

				case 1:
				case 2:
					/*	CIF/CMS */
					
					if (exo_cpn->cms_fix_time > 1.0E-08)
					{
						stdfloor = fabs (exo_cpn->alphabeta[0]) * sqrt (exo_cpn->cms[0].floorvar * exo_cpn->cms_fix_time);
						stdcap = fabs (exo_cpn->alphabeta[0]) * sqrt (exo_cpn->cms[0].capvar * exo_cpn->cms_fix_time);
					}
					else
					{
						stdfloor = stdcap = 0.0;
					}
					break;

				case 3:
					if (exo_cpn->cms_fix_time > 1.0E-08)
					{
						stdfloor = interp(
							spread_vol_time,
							spread_vol_floor,
							spread_vol_n,
							exo_cpn->cms_fix_time,
							0,
							&temp);
						stdcap = interp(
							spread_vol_time,
							spread_vol_cap,
							spread_vol_n,
							exo_cpn->cms_fix_time,
							0,
							&temp);
						
						if (is_corr)
						{
							stdfloor = sqrt(
								exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].floorvar
								+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].floorvar
			+ 2.0 * stdfloor * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] * sqrt (exo_cpn->cms[0].floorvar * exo_cpn->cms[1].floorvar));

							stdcap = sqrt(
								exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].capvar
								+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].capvar
			+ 2.0 * stdcap * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] * sqrt (exo_cpn->cms[0].capvar * exo_cpn->cms[1].capvar));
					
						}
							
						stdfloor *= sqrt (exo_cpn->cms_fix_time);
						stdcap *= sqrt (exo_cpn->cms_fix_time);
					}
					else
					{
						stdfloor = stdcap = 0.0;
					}
					break;
				}


				if (exo_cpn->floored)
				{						
					if (stdfloor > 1.0E-08)
					{
						floortype = 4;	/*	Put */
					}
					else
					{
						floortype = 2;	/*	Put IV */
					}						
				}
				else
				{
					floortype = 0;		/*	None */
				}

				if (exo_cpn->capped)
				{						
					if (stdcap > 1.0E-08)
					{
						captype = 3;	/*	Call */
					}
					else
					{
						captype = 1;	/*	Call IV */
					}						
				}
				else
				{
					captype = 0;		/*	None */
				}

				floor = OPT_VAL_MACRO_N(
					floortype, 
					spread, 
					exo_cpn->floor - exo_cpn->gamma, 
					stdfloor);

				cap = OPT_VAL_MACRO_N(
					captype, 
					spread, 
					exo_cpn->cap - exo_cpn->gamma, 
					stdcap);
				
				/*	Coupon pv */

				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
			else
			{
				/* Use SL Model wo string of options*/
				coupon = spread + exo_cpn->gamma;
				floor = cap = 0.0;

				if (exo_cpn->cms_fix_time > 1.0E-08)
				{						
					/* Extract the correlation */
					spread_slcorr = interp(
											spread_vol_time,
											spread_vol_floor,
											spread_vol_n,
											exo_cpn->cms_fix_time,
											0,
											&temp);		

					/* Evaluation of the floor */
					if (exo_cpn->floored)
					{
						err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->floor - exo_cpn->gamma
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									SRT_PUT,
									&floor);
					}

					/* Evaluation of the cap */
					if (exo_cpn->capped)
					{
						err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->cap - exo_cpn->gamma
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									SRT_CALL,
									&cap);
					}
				}
				
				/*	Coupon pv */

				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
		}
		else
		{
			/* use a string of options */

			if (!use_SL)
			{
				coupon = exo_cpn->gamma;
				floor = cap = 0.0;

				/* Put all the value of the options in the floor variable */
				if (exo_cpn->cms_fix_time > 1.0E-08)
				{
					/* Extract the correlation */
					spread_corr = interp(spread_vol_time,
										spread_vol_floor,
										spread_vol_n,
										exo_cpn->cms_fix_time,
										0,
										&temp);	

					for (iNumOpt=0; iNumOpt<exo_cpn->nopt; iNumOpt++)
					{
						if (exo_cpn->notopt[iNumOpt] != 0)
						{
							stdfloor = sqrt(exo_cpn->cms_fix_time * (
							exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].optvar[iNumOpt]
							+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].optvar[iNumOpt]
							+ 2.0 * spread_corr * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] 
								* sqrt (exo_cpn->cms[0].optvar[iNumOpt] * exo_cpn->cms[1].optvar[iNumOpt])));

							floor += exo_cpn->notopt[iNumOpt] * OPT_VAL_MACRO_N(
																			   (exo_cpn->typeopt[iNumOpt] == SRT_CALL ? (stdfloor == 0 ? 1 : 3) : (stdfloor == 0 ? 2 : 4)),
																				spread,
																				exo_cpn->strikeopt[iNumOpt],
																				stdfloor);
						}
					}
				}
				
				/*	Coupon pv */
				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
			else
			{
				/* Use a Shifted log model for options */
				coupon = exo_cpn->gamma;
				floor = cap = 0.0;

				/* Put all the value of the options in the floor variable */

				/*	Spread Vol */			
				if (exo_cpn->cms_fix_time > 1.0E-08)
				{
					/* Extract the correlation */
					spread_slcorr = interp(spread_vol_time,
										spread_vol_floor,
										spread_vol_n,
										exo_cpn->cms_fix_time,
										0,
										&temp);	

					for (iNumOpt=0; iNumOpt<exo_cpn->nopt; iNumOpt++)
					{
						if (exo_cpn->notopt[iNumOpt] != 0)
						{
								err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->strikeopt[iNumOpt]
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									exo_cpn->typeopt[iNumOpt],
									&optionvalue);

								floor += exo_cpn->notopt[iNumOpt] * optionvalue;
						}
					}
				}				
				
				/*	Coupon pv */
				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
		}
	

		if (calc_fwdiv)
		{
			j = 0;
			while (j < ccf->num_calls && (ccf->call + j)->cf_idx <= i)
			{
				und->market_fwdiv[j] += -fact * temp;
				j++;
			}
		}
	}

	/*	PV of coupons fixed in the past and not yet paid */
	i = 0;
	while (i < cf_ncpn && cf_fix[i] < today + eod_fix_flag)
	{
		if (cf_pay[i] >= today + eod_pay_flag)
		{
			err = interp_basis (cf_basis[i], &bas);
			if (err)
			{
				goto FREE_RETURN;	
			}
			
			exo_leg_pv += cf_fix_cpn[i] 
				* coverage (cf_start[i], cf_pay[i], bas)
				* cf_not
				* swp_f_df (today, cf_pay[i], yc);
		}

		i++;
	}

	/*	Initial and final exchange */
	if (for_fund)
	{
		/*	Final */
		
		if (exo_leg->num_cpn > 0)
		{	
			{
				temp = swp_f_df (today, fund_leg->cpn[fund_leg->num_cpn-1].pay_date, yc)
				* cf_not;
				exo_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += -fact * temp;
					}
				}
			}
		}
		else
		{
			fin_not_date = fund_pay[fund_ncpn-1];
			if (fin_not_date >= today + eod_pay_flag)
			{
				temp = swp_f_df (today, fin_not_date, yc) * cf_not;
				exo_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += -fact * temp;
					}
				}
			}
		}

		/*	Initial */
		if (start_date >= today + eod_pay_flag)
		{
			exo_leg_pv -= swp_f_df (today, start_date, yc) * cf_not;
		}

		if (calc_fwdiv)
		{
			for (i=0; i<ccf->num_calls; i++)
			{
				callt = ccf->call + i;
				und->market_fwdiv[i] -= -fact * swp_f_df (today, (ccf->cf_leg->cpn + callt->cf_idx)->start_date, yc)
				* cf_not;
			}
		}
	}

	if (calc_fwdiv)
	{
		for (i=0; i<ccf->num_calls; i++)
		{
			callt = ccf->call + i;
			und->extra_fees[i] = -(und->market_fwdiv[i] - und->model_fwdiv[i]) / swp_f_df (today, callt->set_date, yc);
			
			if (adjust_fee)
			{
				callt->fee += und->extra_fees[i];
			}
		}
	}

	/*	4)	If there is at least one call after today, value call feature */

	if (call_feat == 1)
	{
		smessage ("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
		
		err = ccf_launch_adi (ccf, und, adi_arg, &call);
		
		if (err)
		{
			goto FREE_RETURN;
		}
	}
	else	
	{
		call = 0.0;
	}

	/* Exerciced Case */
	/* Force the Notional to be paid at ex_settle_date */
	if ((exercised)&& (fund_pay[fund_ncpn-1] != ex_date_set))
	{
		/* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at ex_settle_date */
		if (for_fund)
		{
			exo_leg_pv += (swp_f_df (und->today, ex_date_set, yc) - 
							 swp_f_df (und->today, fund_pay[fund_ncpn-1], yc)) * cf_not;
		}

		/* Add fee */
		if (ex_date_set >= today + eod_pay_flag)
		{
			call += -ex_fee * swp_f_df (today, ex_date_set, yc);
		}
	}

	*fund_val = fund_leg_pv;
	*cf_val = exo_leg_pv;
	*call_val = call;	

	if (export_ts)
	{
		ccf_copy_und (und, und_exp);
	}

FREE_RETURN:

	if (free_struct)
	{
		ccf_free_all_struct (ccf, und, call_feat, adi_arg);
	}
	if (ccf) free (ccf);
	if (und) free (und);
	if (adi_arg) free (adi_arg);

	return err;
}


Err ccf_caller_LambdaShift(
			/*	Today's date */
			long		today,
			/*	The underlying */
			int			use_calib,						/*	0: use lgm2fund, 1: calibrate */
			/*		if calib */
			char		*yc,							/*	yc */
			char		*vc,							/*	vc */
			char		*ref,							/*	ref rate (only if calib) */
			char		*swap_freq,						/*	swap freq (only if calib) */
			char		*swap_basis,					/*	swap basis (only if calib) */
			int			lam_ts,							/*	0: use unique lambda, 1: use ts */
			double		lambda,							/*	lambda if unique */
			int			tsnlam,							/*	number of lambdas if ts */
			double		*tslamtime,						/*	lambda times i.e. (date - today) / 365 if ts */
			double		*tslam,							/*	corresponding lambdas if ts */
			double		alpha,							/*	alpha */
			double		gamma,							/*	gamma */
			double		rho,							/*	rho */
			/*	End of calib params */
			Err			(*get_cash_vol)(				/*	function to get IR cash vol from the markets */
							char	*vol_curve_name,	
							double	start_date, 
							double	end_date,
							double	cash_strike,
							int		zero,
							char	*ref_rate_name,
							double	*vol,
							double	*power),
			/*		if no calilb */
			char		*lgm2dund,
			/*	The structure */
			long		start_date,				/*	Date at which initial notional exchange occurs */
			/*		funding */
			int			fund_ccy,				/*	0: domestic, 1: other */
			double		fund_not,				/*	If different from domestic or foreign (fund_ccy = 1) */
			char		*fund_ccy_yc,			/*	If different from domestic or foreign (fund_ccy = 1) */
			double		fx_fund_dom,			/*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
			long		fx_fund_dom_spot_date,
			int			fund_ncpn,
			long		*fund_fix,
			long		*fund_start,
			long		*fund_pay,
			char		**fund_basis,
			double		*fund_spr,
			double		*fund_mrg,
			double		*fund_fix_cpn,			/*	Past coupon fixing if relevant,
											includes spr, but not mrg, cvg and notional */
			/*		cf */
			double		cf_not,
			int			cf_ncpn,
			long		*cf_fix,
			long		*cf_start,
			long		*cf_pay,
			char		**cf_basis,
			char		**cf_cms_tenor1,
			char		*cf_cms_freq1,
			char		*cf_cms_basis1,
			double		*cf_cms_spread1,
			char		**cf_cms_tenor2,
			char		*cf_cms_freq2,
			char		*cf_cms_basis2,
			double		*cf_cms_spread2,
			long		spread_vol_n,
			double		*spread_vol_time,
			double		*spread_vol_floor,
			double		*spread_vol_cap,
			
			double		*spread_slbeta1,		/* Shifted log beta on the CMS1 */
			double		*spread_slbeta2,		/* Shifted log beta on the CMS2 */

			double		*cf_alpha,
			double		*cf_beta,
			double		*cf_gamma,
			int			*cf_floored,
			double		*cf_floor,
			int			*cf_capped,
			double		*cf_cap,
			double		*cf_fix_cpn,			/*	Past coupon fixing if relevant */

			int			cf_nopt,			/* Number of spread options */
			double		**cf_notopt,		/* Notional of the spread options */	
			double		**cf_strikeopt,		/* spread option strikes */
			int			**cf_typeopt,		/* spread option type 0 call 1 put */

			/*		calls */
			int			ncall,
			int			pay_rec,				/*	0: rec pd, 1: pay pd */
			long		*ex_date,
			long		*set_date,
			double		*fee,
			/*	Numerical params */
			int			req_stp,
			int			req_stpx,
			long			req_paths,
			/*	Calib params */
			int			force_atm,
			double		max_std_long,
			double		max_std_short,
			int			fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
																	to diagonal */
			int			cal_vol_shift_type,				/*	vol shift type for volatility */
			double		cal_vol_shift,					/*	vol shift */
			double		cal_lambda_shift,				/*	shift on lambda after calibration */
			int			one_f_equi,						/*	1F equivalent flag:
																	if set to 1, then 2F lambda will calibrate
																	to the cap priced within calibrated 1F
																	with the given lambda */
			int			skip_last,				/*	If 1, the last option is disregarded and the forward volatility is flat from option n-1 */
			double		long_prec,				/*	Precision on primary instruments */
			double		short_prec,				/*	Precision on secondary instruments */
			double		min_fact,				/*	Maximum down jump on variance */
			double		max_fact,				/*	Maximum up jump on variance */
			int			use_jumps,				/*	Allow vol term structure to jump */
			int			proba_weight,
			

			int			calc_fwdiv,				/*	output the model and market fwdiv  */
			int			adjust_fee,				/*	adjust the fee to price correctly the fwd iv */

			/*	EOD Flags */
			int			eod_fix_flag,			/*	0: I, 1: E */
			int			eod_pay_flag,			/*	0: I, 1: E */
			int			eod_ex_flag,			/*	0: I, 1: E */
			/*	CMS can be valued with smile */
			int			cms_adj,				/*	1: adjust for CMS effect
													0: don't */
			int			cms_for_iv_only,		/*	1: adjust for CMS effect only IV
													0: adjust both IV and call */
			int			use_cms_smile,			/*	1: value CMS with ATM
													0: value CMS with smile */
			int			cms_vol_adj,			/*	1: adjust for CMS vol effect
													0: don't */
			double		cms_beta1,				/*	How adjustment varies with rates */
			double		cms_beta2,
			int			num_strikes_in_vol,		/*	Array of strikes in vol matrix */
			double		*strikes_in_vol,
			SrtDiffusionType vol_type,			/*	Type of vol in matrix, SRT_NORMAL or SRT_LOGNORMAL */
			int			cash_vol,				/*	1: matrix is a cash vol
													0: matrix is a swap vol */
			int			cvg_sv,					/*	1: coverging model on spread vol/correlation
													0: sliding model on spread vol/correlation */
			int			is_corr,				/*	1: spread vols are correlations
													0: spread vols are vols */

			int			use_SL,					/*  1: use Shifted Log 
												    0: don't use				*/
			int			calib_SL,				/*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
													0: use the Shifted Log beta given */
			int			calib_correl_SL,		/*	1: Calibrate the correlation between the two SL to get the same ATM
														normal spread vol 
													0: use the normal spread correl for the sl correl */
			
			int			use_cfoptions,			/*  1: Use the spread options and don't take into account 
														the floor and cap in the cf coupons
													0: take into account the floor and cap in the cf coupons */

			/*	Exercised flag */
			int			exercised,				/*	Flag */
			long		ex_date_ex,				/*	Date when exercised */
			long		ex_date_set,			/*	Corresponding settlement date */
			double		ex_fee,					/*	Corresponding fee */
			double		dLambdaShift,
			/*	Results */
			double		*fund_val,				/*	Value of the funding leg */
			double		*cf_val,				/*	Value of the Power Dual leg */
			double		*call_val,				/*	Value of the callable feature */			
			int			export_ts,				/*	1: Export TS, 0: don't */
			CCF_UND		und_exp)
{
	ccf_str			*ccf		= NULL;
	ccf_und			*und		= NULL;
	ccf_adi_arg		*adi_arg	= NULL;	

	CF_FUND_LEG		fund_leg;
	CF_FUND_CPN		fund_cpn;
	CF_EXO_LEG		exo_leg;
	CF_EXO_CPN		exo_cpn;
	CF_CMS_DESC		cms;

	int				call_feat;
	double			fund_leg_pv, exo_leg_pv;	
	SrtBasisCode	bas;	

	int				floortype, captype;
	double			df, stdfloor, stdcap, floor, cap, coupon;

	double			call;
	int				i, j;
	int				free_struct = 0;
	
	double			temp;

	double			vol, spread, forward, cmsrate[2];

	int				for_fund;
	long			fund_start_date, fin_not_date;
	double			eq_final_ex, eq_init_ex;
	
	double			tmp1, tmp2;

	CF_CALL			callt;
	double			fact;

	Err				err			= NULL;

	double			OtherCcyfundNot;
	char			Copyfund_ccy_yc[255];

	int nShiftLambda;

	int				iNumOpt;
	double			spread_corr, spread_slcorr;
	double			optionvalue;

	
	/* Memory Allocation */
	ccf = calloc(1, sizeof(ccf_str));
	und = calloc(1, sizeof(ccf_und));
	adi_arg = calloc(1, sizeof(ccf_adi_arg));

	if (!ccf || !und || !adi_arg)
	{
		err = "memory allocation failure in ccf_caller";
		goto FREE_RETURN;
	}
		
	/*	If exercised */
	if (exercised)
	{
		/* Consider only the coupons which start before the ex date */
		i = 0;
		while (i < cf_ncpn && cf_start[i] < ex_date_ex)
		{
			i++;
		}
		cf_ncpn= i;

		i = 0;
		while (i < fund_ncpn && fund_start[i] < ex_date_ex)
		{
			i++;
		}
		fund_ncpn = i;

		/* Checks */
		if ((cf_ncpn > 0) && (ex_date_set<cf_pay[cf_ncpn-1]))
		{
			err = "Exercised settlement should be > than payment date of the last not called cf coupon ";
			goto FREE_RETURN;
		}

		if ((fund_ncpn > 0) && (ex_date_set<fund_pay[fund_ncpn-1]))
		{
			err = "Exercised settlement should be > than payment date of the last not called funding coupon ";
			goto FREE_RETURN;
		}

		/* Gestion of the particular case when the Ex settlement is past *
		 * or all the coupons (cf and funding) are called				 */
		if ((ex_date_set < today + eod_pay_flag)||((cf_ncpn == 0)&&(fund_ncpn == 0)))
		{	
			/* Fill the und */
			strcpy (und->name, "CALIB");
			und->today = today;
			strcpy (und->yc, yc);
			strcpy (und->vc, vc);
			strcpy (und->ref, ref);
			strcpy (und->swap_freq, swap_freq);
			strcpy (und->swap_basis, swap_basis);
			
			und->sigma_n = 0;
			und->lambda_n = 0;
			und->spread_vol_n = 0;

			/* Fill the output */
			*fund_val = *cf_val = *call_val = 0.0;

			/* Initial exchange of notional */
			if (start_date >= today + eod_pay_flag)
			{
				if (fund_ccy)
				{
					*fund_val -= swp_f_df (today, start_date, fund_ccy_yc) * fund_not 
						* fx_fund_dom * swp_f_df (today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (today, fx_fund_dom_spot_date, fund_ccy_yc);
					*cf_val -= cf_not * swp_f_df (today, start_date, yc);
				}
			}
			
			/* Final exchange of notional */
			if (ex_date_set >= today + eod_pay_flag)
			{
				if (fund_ccy)
				{
					/* Final exchange of notional */
					*fund_val += swp_f_df (today, ex_date_set, fund_ccy_yc) * fund_not 
						* fx_fund_dom * swp_f_df (today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (today, fx_fund_dom_spot_date, fund_ccy_yc);
					*cf_val += cf_not * swp_f_df (today, ex_date_set, yc);
					*call_val += -ex_fee * swp_f_df (today, ex_date_set, yc);
				}
			}

			/* export the ts if asked */
			if (export_ts)
			{
				ccf_copy_und (und, und_exp);
			}

			return NULL;
		}
		
		/*  those particular cases should not happen */
		if ((cf_ncpn > 0) && (fund_ncpn == 0))
		{
			err = "All the funding coupons are called but not all the cf coupons : Not allowed ";
			goto FREE_RETURN;
		}
		if ((cf_ncpn == 0) && (fund_ncpn > 0))
		{
			err = "All the cf coupons are called but not all the funding coupons : Not allowed ";
			goto FREE_RETURN;
		}
		
		/*  Case where (cf_ncpn > 0) && (fund_ncpn > 0) */
		ncall = 0;
	}

	/*	If cms adjustment only for IV */
	if (!exercised && cms_adj && cms_for_iv_only)
	{
		/* First, fo the non_callable */
		err = ccf_caller(
			today,
			use_calib,
			yc,
			vc,
			ref,
			swap_freq,
			swap_basis,
			lam_ts,
			lambda,
			tsnlam,
			tslamtime,
			tslam,
			alpha,
			gamma,
			rho,
			get_cash_vol,
			lgm2dund,
			start_date,
			fund_ccy,
			fund_not,
			fund_ccy_yc,
			fx_fund_dom,
			fx_fund_dom_spot_date,
			fund_ncpn,
			fund_fix,
			fund_start,
			fund_pay,
			fund_basis,
			fund_spr,
			fund_mrg,
			fund_fix_cpn,
			cf_not,
			cf_ncpn,
			cf_fix,
			cf_start,
			cf_pay,
			cf_basis,
			cf_cms_tenor1,
			cf_cms_freq1,
			cf_cms_basis1,
			cf_cms_spread1,
			cf_cms_tenor2,
			cf_cms_freq2,
			cf_cms_basis2,
			cf_cms_spread2,
			spread_vol_n,
			spread_vol_time,
			spread_vol_floor,
			spread_vol_cap,

			spread_slbeta1,
			spread_slbeta2,

			cf_alpha,
			cf_beta,
			cf_gamma,
			cf_floored,
			cf_floor,
			cf_capped,
			cf_cap,
			cf_fix_cpn,

			cf_nopt,
			cf_notopt,
			cf_strikeopt,
			cf_typeopt,

			0,
			pay_rec,
			ex_date,
			set_date,
			fee,
			req_stp,
			req_stpx,
			req_paths,
			force_atm,
			max_std_long,
			max_std_short,
			fix_lambda,				
			cal_vol_shift_type,
			cal_vol_shift,
			cal_lambda_shift,
			one_f_equi,
			skip_last,
			long_prec,
			short_prec,
			min_fact,
			max_fact,
			use_jumps,
			proba_weight,
			calc_fwdiv,
			adjust_fee,
			eod_fix_flag,
			eod_pay_flag,
			eod_ex_flag,
			1,
			0,
			use_cms_smile,
			cms_vol_adj,
			cms_beta1,
			cms_beta2,
			num_strikes_in_vol,
			strikes_in_vol,
			vol_type,
			cash_vol,
			cvg_sv,
			is_corr,

			use_SL,
			calib_SL,
			0.5,
			calib_correl_SL,
			use_cfoptions,

			exercised,
			ex_date_ex,
			ex_date_set,
			ex_fee,
			fund_val,
			cf_val,
			call_val,
			export_ts,
			und_exp);
		if (err)
		{
			goto FREE_RETURN;
		}

		/*	Then value the call feature */
		err = ccf_caller(
			today,
			use_calib,
			yc,
			vc,
			ref,
			swap_freq,
			swap_basis,
			lam_ts,
			lambda,
			tsnlam,
			tslamtime,
			tslam,
			alpha,
			gamma,
			rho,
			get_cash_vol,
			lgm2dund,
			start_date,
			fund_ccy,
			fund_not,
			fund_ccy_yc,
			fx_fund_dom,
			fx_fund_dom_spot_date,
			fund_ncpn,
			fund_fix,
			fund_start,
			fund_pay,
			fund_basis,
			fund_spr,
			fund_mrg,
			fund_fix_cpn,
			cf_not,
			cf_ncpn,
			cf_fix,
			cf_start,
			cf_pay,
			cf_basis,
			cf_cms_tenor1,
			cf_cms_freq1,
			cf_cms_basis1,
			cf_cms_spread1,
			cf_cms_tenor2,
			cf_cms_freq2,
			cf_cms_basis2,
			cf_cms_spread2,
			spread_vol_n,
			spread_vol_time,
			spread_vol_floor,
			spread_vol_cap,

			spread_slbeta1,
			spread_slbeta2,

			cf_alpha,
			cf_beta,
			cf_gamma,
			cf_floored,
			cf_floor,
			cf_capped,
			cf_cap,
			cf_fix_cpn,

			cf_nopt,
			cf_notopt,
			cf_strikeopt,
			cf_typeopt,

			ncall,
			pay_rec,
			ex_date,
			set_date,
			fee,
			req_stp,
			req_stpx,
			req_paths,
			force_atm,
			max_std_long,
			max_std_short,
			fix_lambda,				
			cal_vol_shift_type,
			cal_vol_shift,
			cal_lambda_shift,
			one_f_equi,
			skip_last,
			long_prec,
			short_prec,
			min_fact,
			max_fact,
			use_jumps,
			proba_weight,
			calc_fwdiv,
			adjust_fee,
			eod_fix_flag,
			eod_pay_flag,
			eod_ex_flag,
			0,
			0,
			use_cms_smile,
			cms_vol_adj,
			cms_beta1,
			cms_beta2,
			num_strikes_in_vol,
			strikes_in_vol,
			vol_type,
			cash_vol,
			cvg_sv,
			is_corr,

			use_SL,
			calib_SL,
			0.5,
			calib_correl_SL,
			use_cfoptions,

			exercised,
			ex_date_ex,
			ex_date_set,
			ex_fee,
			&tmp1,
			&tmp2,
			call_val,
			export_ts,
			und_exp);
			
		goto FREE_RETURN;
	}

	if (fund_ccy == 1)
	{
		/* save the initial fund not in the third ccy 
		   and the third Ccy yield curve name			*/
		OtherCcyfundNot = fund_not;
		strcpy (Copyfund_ccy_yc, fund_ccy_yc);

		/* Convert into domestic */
		fund_ccy = 0;
		for_fund = 1;
		err = convert_funding_to_domestic(
								today,
								start_date,
								eod_fix_flag,
								eod_pay_flag,
								fx_fund_dom,
								fx_fund_dom_spot_date,
								cf_not,
								yc,
								fund_ncpn,
								fund_fix,
								fund_start,
								fund_pay,
								fund_basis,
								fund_ccy_yc,		
								&fund_not,
								fund_spr,
								fund_mrg,
								fund_fix_cpn,
								&fund_start_date,
								&eq_final_ex,
								&eq_init_ex);
		if (err)
		{
			return err;
		}
	}
	else
	{
		for_fund = 0;
	}

	/*	Initialise structures */
	free_struct = 0;
	err = ccf_fill_check_all_struct(
		today,
		use_calib,
		yc,
		vc,
		ref,
		swap_freq,
		swap_basis,
		lam_ts,
		lambda,
		tsnlam,
		tslamtime,
		tslam,
		alpha,
		gamma,
		rho,
		get_cash_vol,
		lgm2dund,
		fund_not,		
		fund_ncpn,
		fund_fix,
		fund_start,
		fund_pay,
		fund_basis,
		fund_spr,
		fund_mrg,
		cf_not,
		cf_ncpn,
		cf_fix,
		cf_start,
		cf_pay,
		cf_basis,
		cf_cms_tenor1,
		cf_cms_freq1,
		cf_cms_basis1,
		cf_cms_spread1,
		cf_cms_tenor2,
		cf_cms_freq2,
		cf_cms_basis2,
		cf_cms_spread2,
		spread_vol_n,
		spread_vol_time,
		spread_vol_floor,
		spread_vol_cap,

		spread_slbeta1,		/* Shifted log beta on the CMS1 */
		spread_slbeta2,		/* Shifted log beta on the CMS2 */

		cvg_sv,
		is_corr,
		cf_alpha,
		cf_beta,
		cf_gamma,
		cf_floored,
		cf_floor,
		cf_capped,
		cf_cap,	
		
		cf_nopt,		/* Number of spread options */
		cf_notopt,		/* Notional of the spread options */	
		cf_strikeopt,	/* spread option strikes */
		cf_typeopt,		/* spread option type 0 call 1 put */
		use_SL,					/*  1: use Shifted Log 
									0: don't use				*/
		calib_SL,				/*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
									0: use the Shifted Log beta given */
		0.5,
		calib_correl_SL,		/*	1: Calibrate the correlation between the two SL to get the same ATM
									   normal spread vol 
									0: use the normal spread correl for the sl correl */
			
		use_cfoptions,			/*  1: Use the spread options and don't take into account 
														the floor and cap in the cf coupons
									0: take into account the floor and cap in the cf coupons */

		cms_adj,
		use_cms_smile,
		cms_vol_adj,
		cms_beta1,
		cms_beta2,
		num_strikes_in_vol,
		strikes_in_vol,
		vol_type,
		cash_vol,
		ncall,
		pay_rec,
		ex_date,
		set_date,
		fee,
		req_stp,
		req_stpx,	
		req_paths,
		force_atm,
		max_std_long,
		max_std_short,
		fix_lambda,				
		cal_vol_shift_type,
		cal_vol_shift,
		cal_lambda_shift,
		one_f_equi,
		skip_last,
		long_prec,
		short_prec,
		min_fact,
		max_fact,
		use_jumps,
		proba_weight,
		eod_fix_flag,
		eod_ex_flag,
		ccf,
		und,
		&call_feat,
		adi_arg);

	for(nShiftLambda = 0; nShiftLambda < und->lambda_n; ++ nShiftLambda)
	{
		und->lambda[nShiftLambda] += (dLambdaShift/100.);
	}

	/// adding shift to lambda


	if (err)
	{
		goto FREE_RETURN;
	}
	free_struct = 1;

	/*  0) calculate the fwd iv in the model */

	if (calc_fwdiv)
	{
		und->has_fwd_iv = 1;
		und->nb_fwdiv = ccf->num_calls;
		
		und->exercise_date = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->market_fwdiv = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->model_fwdiv = (double *) calloc(und->nb_fwdiv, sizeof(double));
		und->extra_fees = (double *) calloc(und->nb_fwdiv, sizeof(double));

		if (!und->exercise_date || !und->market_fwdiv || !und->model_fwdiv || !und->extra_fees)
		{
			err = "Memory allocation faillure in ccf caller";
			goto FREE_RETURN;
		}

		err = ccf_calc_mdl_iv_fwd(	ccf,
							und,
							adi_arg,
							NUM_HERMITE,
							und->model_fwdiv);

		if (err)
		{
			goto FREE_RETURN;
		}

		if (pay_rec == 1)
		{
			fact = 1.0;
		}
		else
		{
			fact = -1.0;
		}

		for (i=0; i<ccf->num_calls; i++)
		{
			callt = ccf->call + i;
			und->exercise_date[i] = callt->ex_date;
		}
	}
		
	/*	1)	Value funding leg */			

	fund_leg = ccf->fund_leg;
	fund_leg_pv = 0.0;	


	/*	Cash libor */
	if (fund_leg->num_cpn > 0)
	{
		if (for_fund)
		{
			fund_leg_pv += swp_f_df (today, fund_leg->cpn[0].start_date, yc)
				* eq_final_ex;

			if (calc_fwdiv)
			{
				/* initialisation of the initial notional */
				for (i=0; i<ccf->num_calls; i++)
				{
					callt = ccf->call + i;
					und->market_fwdiv[i] = fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
					* eq_final_ex;
				}
			}
		}
		else
		{
			fund_leg_pv += swp_f_df (today, fund_leg->cpn[0].start_date, yc)
				* fund_leg->notional;

			if (calc_fwdiv)
			{
				/* initialisation of the initial notional */
				for (i=0; i<ccf->num_calls; i++)
				{
					callt = ccf->call + i;
					und->market_fwdiv[i] = fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
					* fund_leg->notional;
				}
			}
		}
	}
	else
	{
		fin_not_date = fund_pay[fund_ncpn-1];

		if (fin_not_date >= today + eod_pay_flag)
		{
			if (for_fund)
			{
				temp = swp_f_df (today, fin_not_date, yc) * eq_final_ex;
				fund_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += fact * temp;
					}
				}
			}
			else
			{
				/*fund_leg_pv += swp_f_df (today, fin_not_date, yc) * fund_leg->notional;*/
			}
		}	
	}

	/*	Coupons: spread + margin */
	for (i=0; i<fund_leg->num_cpn; i++)
	{
		fund_cpn = fund_leg->cpn + i;
		temp = swp_f_df (today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
		fund_leg_pv += temp;

		if (calc_fwdiv)
		{
			j = 0;
			while (j < ccf->num_calls && (ccf->call + j)->fund_idx <= i)
			{
				und->market_fwdiv[j] += fact * temp;
				j++;
			}
		}
	}

	/*	Notional exchange */
	if (for_fund)
	{
		if (start_date >= today + eod_pay_flag)
		{
			fund_leg_pv -= swp_f_df (today, start_date, yc) * eq_init_ex;
		}

		if (calc_fwdiv)
		{
			for (i=0; i<ccf->num_calls; i++)
			{
				callt = ccf->call + i;
				und->market_fwdiv[i] -= fact * swp_f_df (today, (ccf->fund_leg->cpn + callt->fund_idx)->start_date, yc)
				* eq_final_ex;
			}
		}
	}
	else
	{
		if (fund_leg->num_cpn > 0)
		{
			temp = swp_f_df (today, fund_leg->cpn[fund_leg->num_cpn-1].pay_date, yc) 
			* fund_leg->notional;
			fund_leg_pv -= temp;

			if (calc_fwdiv)
			{				
				for (i=0; i<ccf->num_calls; i++)
				{
					und->market_fwdiv[i] -= fact * temp;
				}
			}
		}
	}

	/*	PV of coupons fixed in the past and not yet paid */
	i = 0;
	while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
	{
		if (fund_pay[i] >= today + eod_pay_flag)
		{
			err = interp_basis (fund_basis[i], &bas);
			if (err)
			{
				goto FREE_RETURN;	
			}
			
			fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i])
				* coverage (fund_start[i], fund_pay[i], bas)
				* fund_not
				* swp_f_df (today, fund_pay[i], yc);
		}

		i++;
	}

	/* Exerciced Case */
	/* Force the Notional to be paid at ex_settle_date */
	if ((exercised) && (fund_pay[fund_ncpn-1] != ex_date_set))
	{
		/* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at ex_settle_date */
		if (for_fund)
		{
			/* Other currency case */
			fund_leg_pv += (swp_f_df (und->today, ex_date_set, Copyfund_ccy_yc) - 
						swp_f_df (und->today, fund_pay[fund_ncpn-1], Copyfund_ccy_yc)) *
						OtherCcyfundNot 
						* fx_fund_dom * swp_f_df (und->today, fx_fund_dom_spot_date, yc) 
						/ swp_f_df (und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);
		}
	}
		
	/*	2)	Value cf leg */

	exo_leg = ccf->cf_leg;
	exo_leg_pv = 0.0;

	
	/*	Coupons */

	for (i=0; i<exo_leg->num_cpn; i++)
	{
		exo_cpn = exo_leg->cpn + i;

		/*	Discount */
		df = swp_f_df (today, exo_cpn->pay_date, yc);
		
		/*	CMS */
		spread = 0.0;
		for (j=0; j<exo_cpn->ncms; j++)
		{
			cms = &(exo_cpn->cms[j]);
			
			if (exo_cpn->cms_fix_date > today)
			{
				vol = cms->atmvar * exo_cpn->cms_fix_time;
			}
			else
			{
				vol = 0.0;
			}
			forward = cms->fwd + cms->swap_spread;

			if (cms_adj && vol > 1.0E-08)
			{
				if (use_cms_smile)
				{
					cmsrate[j] = cms->fwd_cms;					
				}
				else
				{
					cmsrate[j] = cms_rate(
						 cms->delay,
						 cms->nfp,
						 cms->cpnd,
						 forward,
						 vol);
				}
			}
			else
			{
				cmsrate[j] = forward;
			}

			spread += exo_cpn->alphabeta[j] * cmsrate[j];
		}

	
		/* Case no use of string of options */
		if (!use_cfoptions )
		{
			if (!use_SL)
			{
				coupon = spread + exo_cpn->gamma;

				/*	Spread Vol */
				
				switch (exo_cpn->type)
				{
				case 0:
				default:
					/*	Midat */
					break;

				case 1:
				case 2:
					/*	CIF/CMS */
					
					if (exo_cpn->cms_fix_time > 1.0E-08)
					{
						stdfloor = fabs (exo_cpn->alphabeta[0]) * sqrt (exo_cpn->cms[0].floorvar * exo_cpn->cms_fix_time);
						stdcap = fabs (exo_cpn->alphabeta[0]) * sqrt (exo_cpn->cms[0].capvar * exo_cpn->cms_fix_time);
					}
					else
					{
						stdfloor = stdcap = 0.0;
					}
					break;

				case 3:
					if (exo_cpn->cms_fix_time > 1.0E-08)
					{
						stdfloor = interp(
							spread_vol_time,
							spread_vol_floor,
							spread_vol_n,
							exo_cpn->cms_fix_time,
							0,
							&temp);
						stdcap = interp(
							spread_vol_time,
							spread_vol_cap,
							spread_vol_n,
							exo_cpn->cms_fix_time,
							0,
							&temp);
						
						if (is_corr)
						{
							stdfloor = sqrt(
								exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].floorvar
								+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].floorvar
			+ 2.0 * stdfloor * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] * sqrt (exo_cpn->cms[0].floorvar * exo_cpn->cms[1].floorvar));

							stdcap = sqrt(
								exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].capvar
								+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].capvar
			+ 2.0 * stdcap * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] * sqrt (exo_cpn->cms[0].capvar * exo_cpn->cms[1].capvar));
					
						}
							
						stdfloor *= sqrt (exo_cpn->cms_fix_time);
						stdcap *= sqrt (exo_cpn->cms_fix_time);
					}
					else
					{
						stdfloor = stdcap = 0.0;
					}
					break;
				}


				if (exo_cpn->floored)
				{						
					if (stdfloor > 1.0E-08)
					{
						floortype = 4;	/*	Put */
					}
					else
					{
						floortype = 2;	/*	Put IV */
					}						
				}
				else
				{
					floortype = 0;		/*	None */
				}

				if (exo_cpn->capped)
				{						
					if (stdcap > 1.0E-08)
					{
						captype = 3;	/*	Call */
					}
					else
					{
						captype = 1;	/*	Call IV */
					}						
				}
				else
				{
					captype = 0;		/*	None */
				}

				floor = OPT_VAL_MACRO_N(
					floortype, 
					spread, 
					exo_cpn->floor - exo_cpn->gamma, 
					stdfloor);

				cap = OPT_VAL_MACRO_N(
					captype, 
					spread, 
					exo_cpn->cap - exo_cpn->gamma, 
					stdcap);
				
				/*	Coupon pv */

				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
			else
			{
				/* Use SL Model wo string of options*/
				coupon = spread + exo_cpn->gamma;
				floor = cap = 0.0;

				if (exo_cpn->cms_fix_time > 1.0E-08)
				{						
					/* Extract the correlation */
					spread_slcorr = interp(
											spread_vol_time,
											spread_vol_floor,
											spread_vol_n,
											exo_cpn->cms_fix_time,
											0,
											&temp);		

					/* Evaluation of the floor */
					if (exo_cpn->floored)
					{
						err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->floor - exo_cpn->gamma
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									SRT_PUT,
									&floor);
					}

					/* Evaluation of the cap */
					if (exo_cpn->capped)
					{
						err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->cap - exo_cpn->gamma
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									SRT_CALL,
									&cap);
					}
				}
				
				/*	Coupon pv */

				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
		}
		else
		{
			/* use a string of options */

			if (!use_SL)
			{
				coupon = exo_cpn->gamma;
				floor = cap = 0.0;

				/* Put all the value of the options in the floor variable */
				if (exo_cpn->cms_fix_time > 1.0E-08)
				{
					/* Extract the correlation */
					spread_corr = interp(spread_vol_time,
										spread_vol_floor,
										spread_vol_n,
										exo_cpn->cms_fix_time,
										0,
										&temp);	

					for (iNumOpt=0; iNumOpt<exo_cpn->nopt; iNumOpt++)
					{
						if (exo_cpn->notopt[iNumOpt] != 0)
						{
							stdfloor = sqrt(exo_cpn->cms_fix_time * (
							exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] * exo_cpn->cms[0].optvar[iNumOpt]
							+ exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] * exo_cpn->cms[1].optvar[iNumOpt]
							+ 2.0 * spread_corr * exo_cpn->alphabeta[0] * exo_cpn->alphabeta[1] 
								* sqrt (exo_cpn->cms[0].optvar[iNumOpt] * exo_cpn->cms[1].optvar[iNumOpt])));

							floor += exo_cpn->notopt[iNumOpt] * OPT_VAL_MACRO_N(
																			   (exo_cpn->typeopt[iNumOpt] == SRT_CALL ? (stdfloor == 0 ? 1 : 3) : (stdfloor == 0 ? 2 : 4)),
																				spread,
																				exo_cpn->strikeopt[iNumOpt],
																				stdfloor);
						}
					}
				}
				
				/*	Coupon pv */
				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
			else
			{
				/* Use a Shifted log model for options */
				coupon = exo_cpn->gamma;
				floor = cap = 0.0;

				/* Put all the value of the options in the floor variable */

				/*	Spread Vol */			
				if (exo_cpn->cms_fix_time > 1.0E-08)
				{
					/* Extract the correlation */
					spread_slcorr = interp(spread_vol_time,
										spread_vol_floor,
										spread_vol_n,
										exo_cpn->cms_fix_time,
										0,
										&temp);	

					for (iNumOpt=0; iNumOpt<exo_cpn->nopt; iNumOpt++)
					{
						if (exo_cpn->notopt[iNumOpt] != 0)
						{
								err = OptSpreadNew(
									cmsrate[0] + exo_cpn->cms[0].slshift,									/* fwd x */
									exo_cpn->alphabeta[0],													/*  nx,  */
									exo_cpn->cms[0].slvol,													/* sigx, */
									cmsrate[1] + exo_cpn->cms[1].slshift,									/* fwdy, */
									exo_cpn->alphabeta[1],													/* ny,   */
									exo_cpn->cms[1].slvol,													/* sigy, */
									exo_cpn->strikeopt[iNumOpt]
										+ exo_cpn->alphabeta[0] * exo_cpn->cms[0].slshift
										+ exo_cpn->alphabeta[1] * exo_cpn->cms[1].slshift,					/* K,	 */
									exo_cpn->cms_fix_time,													/* mat,  */
									spread_slcorr,															/* rho,  */
									exo_cpn->typeopt[iNumOpt],
									&optionvalue);

								floor += exo_cpn->notopt[iNumOpt] * optionvalue;
						}
					}
				}				
				
				/*	Coupon pv */
				temp = df * (coupon + floor - cap) * exo_cpn->cvg;
				exo_leg_pv += temp;
			}
		}
	

		if (calc_fwdiv)
		{
			j = 0;
			while (j < ccf->num_calls && (ccf->call + j)->cf_idx <= i)
			{
				und->market_fwdiv[j] += -fact * temp;
				j++;
			}
		}
	}

	/*	PV of coupons fixed in the past and not yet paid */
	i = 0;
	while (i < cf_ncpn && cf_fix[i] < today + eod_fix_flag)
	{
		if (cf_pay[i] >= today + eod_pay_flag)
		{
			err = interp_basis (cf_basis[i], &bas);
			if (err)
			{
				goto FREE_RETURN;	
			}
			
			exo_leg_pv += cf_fix_cpn[i] 
				* coverage (cf_start[i], cf_pay[i], bas)
				* cf_not
				* swp_f_df (today, cf_pay[i], yc);
		}

		i++;
	}

	/*	Initial and final exchange */
	if (for_fund)
	{
		/*	Final */
		
		if (exo_leg->num_cpn > 0)
		{	
			{
				temp = swp_f_df (today, fund_leg->cpn[fund_leg->num_cpn-1].pay_date, yc)
				* cf_not;
				exo_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += -fact * temp;
					}
				}
			}
		}
		else
		{
			fin_not_date = fund_pay[fund_ncpn-1];
			if (fin_not_date >= today + eod_pay_flag)
			{
				temp = swp_f_df (today, fin_not_date, yc) * cf_not;
				exo_leg_pv += temp;

				if (calc_fwdiv)
				{				
					for (i=0; i<ccf->num_calls; i++)
					{
						und->market_fwdiv[i] += -fact * temp;
					}
				}
			}
		}

		/*	Initial */
		if (start_date >= today + eod_pay_flag)
		{
			exo_leg_pv -= swp_f_df (today, start_date, yc) * cf_not;
		}

		if (calc_fwdiv)
		{
			for (i=0; i<ccf->num_calls; i++)
			{
				callt = ccf->call + i;
				und->market_fwdiv[i] -= -fact * swp_f_df (today, (ccf->cf_leg->cpn + callt->cf_idx)->start_date, yc)
				* cf_not;
			}
		}
	}

	if (calc_fwdiv)
	{
		for (i=0; i<ccf->num_calls; i++)
		{
			callt = ccf->call + i;
			und->extra_fees[i] = -(und->market_fwdiv[i] - und->model_fwdiv[i]) / swp_f_df (today, callt->set_date, yc);
			
			if (adjust_fee)
			{
				callt->fee += und->extra_fees[i];
			}
		}
	}

	/*	4)	If there is at least one call after today, value call feature */

	if (call_feat == 1)
	{
		smessage ("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
		
		err = ccf_launch_adi (ccf, und, adi_arg, &call);
		
		if (err)
		{
			goto FREE_RETURN;
		}
	}
	else	
	{
		call = 0.0;
	}

	/* Exerciced Case */
	/* Force the Notional to be paid at ex_settle_date */
	if ((exercised)&& (fund_pay[fund_ncpn-1] != ex_date_set))
	{
		/* Substract the notional received at fund_pay[fund_ncpn-1] and add notional paid at ex_settle_date */
		if (for_fund)
		{
			exo_leg_pv += (swp_f_df (und->today, ex_date_set, yc) - 
							 swp_f_df (und->today, fund_pay[fund_ncpn-1], yc)) * cf_not;
		}

		/* Add fee */
		if (ex_date_set >= today + eod_pay_flag)
		{
			call += -ex_fee * swp_f_df (today, ex_date_set, yc);
		}
	}

	*fund_val = fund_leg_pv;
	*cf_val = exo_leg_pv;
	*call_val = call;	

	if (export_ts)
	{
		ccf_copy_und (und, und_exp);
	}

FREE_RETURN:

	if (free_struct)
	{
		ccf_free_all_struct (ccf, und, call_feat, adi_arg);
	}
	if (ccf) free (ccf);
	if (und) free (und);
	if (adi_arg) free (adi_arg);

	return err;
}

