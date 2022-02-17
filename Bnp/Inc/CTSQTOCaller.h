
#ifndef __CTSQTO_CALLER_H
#define __CTSQTO_CALLER_H

#include "LGMQuantoUnd.h"


/*	-------------------------------------------- */
/*	Caller for callable time swap quanto options */
/*	-------------------------------------------- */


Err ctsqto_caller(
			/*	Today's date */
			long		today,
			/*	The underlying */
			int			use_calib,						/*	0: use lgm2fund, 1: calibrate */
			/*		if calib */
			char		*dom_yc,							/*	dom yc */
			char		*dom_vc,							/*	dom vc */
			char		*dom_ref,							/*	dom ref rate */
			char		*dom_swap_freq,						/*	dom swap freq */
			char		*dom_swap_basis,					/*	dom swap basis */
			double		dom_lambda,							/*	dom lambda if unique */

			char		*for_yc,							/*	for yc */
			char		*for_vc,							/*	for vc */
			char		*for_ref,							/*	for ref rate */
			char		*for_instr_freq,					/*	for instr freq */
			char		*for_instr_basis,					/*	for instr basis */
			double		for_lambda,							/*	for lambda if unique */

			int			forcalib,							/*	0 : RA und, 1 : Diag */

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
			char		*fx3dund,
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
			/*		ra */
			double		ra_not,
			int			ra_cpn_type,
			int			ra_ncpn,
			double		*ra_cpns,
			long		*ra_start,
			long		*ra_pay,
			char		*ra_refrate,
			char		*ra_basis,
			int			*ra_nfixings,
			long		**ra_fixingdates,
			double		**ra_fixings,			/*	Past coupon fixing if relevant */

			// RA floating coupons
			int			ra_float_refrate_is_dom_for,
			char		*ra_float_refrate,
			long		*ra_float_fixl,
			double		*ra_float_past_fixings,
			double		*ra_float_gearings,			
			
			double		*upper_barr,
			double		*lower_barr,
			double		c_spread,
			long		obs_freq_iv,
			long		obs_freq_op,
			double		rho_df,

			// Params for numeraire adjustment if a floating coupon is paid
			double		correl_start,
			double		correl_end,
			int			float_adj_strike,

			int			n_fxvol_dates,
			long		*fxvol_dates,
			double		*fxvol,
			double		*qtocorrel,

			int			typeVol,

			double		*corr_times,
			double		*correl_dom_for,
			double		*correl_dom_fx,
			double		*correl_for_fx,
			int			corr_n_times,

			/*		calls */
			int			ncall,
			int			pay_rec,				/*	0: rec pd, 1: pay pd */
			long		*ex_date,
			long		*set_date,
			double		*fee,
			/*	Numerical params */
			int			req_stp,
			int			req_stpx,
			/*	Calib params */
			int			dom_force_atm,
			int			for_force_atm,
			double		max_std_long,
			double		max_std_short,
			int			fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
																	to diagonal */
			int			one_f_equi,						/*	1F equivalent flag:
																	if set to 1, then 2F lambda will calibrate
																	to the cap priced within calibrated 1F
																	with the given lambda */
			int			skip_last,						/*	If 1, the last option is disregarded
																	and the forward volatility is flat from option
																	n-1 */

			int			calc_fwdiv,						/* output the model and market fwdiv  */
			int			adjust_fee,						/* adjust the fee to price correctly the fwd iv */

			/*	EOD Flags */
			int			eod_fix_flag,			/*	0: I, 1: E */
			int			eod_pay_flag,			/*	0: I, 1: E */
			int			eod_ex_flag,			/*	0: I, 1: E */

			/*	Exercised flag */
			int			exercised,				/*	Flag */
			long		ex_date_ex,				/*	Date when exercised */
			long		ex_date_set,			/*	Corresponding settlement date */
			double		ex_fee,					/*	Corresponding fee */
			/*	Results */
			double		*fund_val,				/*	Value of the funding leg */
			double		*ra_val,				/*	Value of the Range Accrual leg */
			double		*call_val,				/*	Value of the callable feature */
			int			export_ts,				/*	1: Export TS, 0: don't */
			LGMQTO_UND		und_exp);


#endif