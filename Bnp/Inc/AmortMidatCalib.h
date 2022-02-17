
#ifndef __AMORTMIDAT_CALIB_H
#define __AMORTMIDAT_CALIB_H

#include "CPDCalib.h"
#include "srt_h_und_struct.h"

Err AmortMidat_lgmprcapgivenlambda(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_fee[],							/*	Exercise Fee for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
//double			swaption_strikes[],					/*	Strikes for standard swaptions */
double			swaption_cpn[],					/*	Discounted Cash-Flows */
//double			cpn_cvg_standard[],					/*	standard cvg from i-1 to i */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
double			lambda,								/*	Lambda */
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last,							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
int				price_cap,							/*	0: just calibrate */
double			*ex_sprice);							/*	Cap price as output */


double amortMidat_lgmcapval1F_b(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	double		cpn[],								/*	Notional */
	double		df[],								/*	Df to cash flow dates */
	double		cvg[],								/*	cvg from i-1 to i */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	int			nex,								/*	Number of exercise dates */
	double		ex_sstrike[],						/*	Strikes for cap */
	double		floatcoupon[],						/*	Float Coupon*/
	int			ex_cpn[],							/*	For each exercise date, first coupon
															to be exercised */
	int			ex_ncpn[],							/*	For each exercise date, number of coupons
															to be exercised */
	double		ex_zeta1[],							/*	Z1 at exercise date */
	double		ex_G1[]);								/*	G1 at exercise date */


double amortMidat_lgmcapval2F_b(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	double		cpn[],								/*	Notional */
	double		df[],								/*	Df to cash flow dates */
	double		cvg[],								/*	cvg from i-1 to i */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		cpn_G2[],							/*	G2 at cash-flow dates */
	int			nex,								/*	Number of exercise dates */
	double		ex_sstrike[],						/*	Strikes for cap */
	double		floatcoupon[],						/*	Float Coupon*/
	int			ex_cpn[],							/*	For each exercise date, first coupon
															to be exercised */
	int			ex_ncpn[],							/*	For each exercise date, number of coupons
															to be exercised */
	double		ex_zeta1[],							/*	Z1 at exercise date */
	double		ex_zeta2[],							/*	Z2 at exercise date */
	double		ex_zeta12[],						/*	Z12 at exercise date */
	double		ex_G1[],							/*	G1 at exercise date */
	double		ex_G2[]);							/*	G2 at exercise date */


double amortMidat_lgmcapval2F(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	double		cpn[],								/*	Notional */
	double		df[],								/*	Df to cash flow dates */
	double		cvg[],								/*	cvg from i-1 to i */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		cpn_G2[],							/*	G2 at cash-flow dates */
	int			nex,								/*	Number of exercise dates */
	int			ex_cpn[],							/*	For each exercise date, first coupon
															to be exercised */
	int			ex_ncpn[],							/*	For each exercise date, number of coupons
															to be exercised */
	double		ex_zeta1[],							/*	Z1 at exercise date */
	double		ex_zeta2[],							/*	Z2 at exercise date */
	double		ex_zeta12[],						/*	Z12 at exercise date */
	double		ex_G1[],							/*	G1 at exercise date */
	double		ex_G2[]);							/*	G2 at exercise date */


Err AmortMidat_lgmcalibzeta2F_ts(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flow */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
double			cpn_G1[],							/*	G1 at cash-flow dates */
double			cpn_G2[],							/*	G2 at cash-flow dates */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
double			ex_G1[],							/*	G1 at exercise date */
double			ex_G2[],							/*	G2 at exercise date */
//double			strike[],							/*	Strikes */
double			mkt_price[],						/*	Market prices */
double			ex_fee[],							/*	Exercise Fees for Diagonal */
double			ex_zeta[],							/*	Output: zetas (1) */
double			floatnotional[],					/*	Float Notional	*/
/*	Lambda, Alpha, gamma, rho */
int nlambda,
double lambda_time[], 
double lambda[],
double			alpha,
double			gamma,
double			rho,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
/*	Calibrate zeta to amortized diagonal given G: 2F case */
Err zcAmortMidat_lgmcalibzeta2F(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flow */
double			zc_last_cpn[],						/*	Discounted Last Fix CF */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
double			cpn_G1[],							/*	G1 at cash-flow dates */
double			cpn_G2[],							/*	G2 at cash-flow dates */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
double			ex_G1[],							/*	G1 at exercise date */
double			ex_G2[],							/*	G2 at exercise date */
//double			strike[],							/*	Strikes */
double			mkt_price[],						/*	Market prices */
double			ex_zeta[],							/*	Output: zetas (1) */
double			floatnotional[],					/*	Float Notional	*/
/*	Lambda, Alpha, gamma, rho */
double			lambda,
double			alpha,
double			gamma,
double			rho,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */

Err zcAmortMidat_lgmcalibzeta1F(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flow */
double			zc_last_cpn[],						/*	Discounted Last Fix CF */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
double			cpn_G[],							/*	G at cash-flow dates */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
double			ex_G[],								/*	G at exercise date */
double			mkt_price[],						/*	Market prices */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
double			lambda,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
Err AmortMidat_lgmprcapgivenlambda_ts(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_fee[],							/*	Exercise Fee for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
//double			swaption_strikes[],					/*	Strikes for standard swaptions */
double			swaption_cpn[],					/*	Discounted Cash-Flows */
//double			cpn_cvg_standard[],					/*	standard cvg from i-1 to i */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/

int	nNumLambdas,
double *pdLambdaValue,
double *pdLambdaTime,				

int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last,							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
int				price_cap,							/*	0: just calibrate */
double			*ex_sprice);							/*	Cap price as output */

Err zcAmortMidat_lgmprcapgivenlambda(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			zc_last_cpn[],						/*	Discounted Last Fix CF */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
double			lambda,								/*	Lambda */
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last,							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
int				price_cap,							/*	0: just calibrate */
double			*ex_sprice);							/*	Cap price as output */


Err amortMidat_lgmcalibzetalambda(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_fee[],							/*	Exercise Fee for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
double			ex_sprice,							/*	Market price for cap */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
double			swaption_cpn[],						/*	Standard Swaption Coupons	*/
//double			standard_cvg[],						/*	Standard Coverages	*/
int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
int				cap_or_swaption,					/*	0: calib lambda to cap, 1: calib lambda to sum of diag swaptions */
double			*lambda,							/*	Lambda: may be changed in the process */
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */

Err amortMidat_lgmcalibzetalambda_ts(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_fee[],							/*	Exercise Fee for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
double			ex_sprice,							/*	Market price for cap */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
double			swaption_cpn[],						/*	Standard Swaption Coupons	*/
//double			standard_cvg[],						/*	Standard Coverages	*/
int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
int				cap_or_swaption,					/*	0: calib lambda to cap, 1: calib lambda to sum of diag swaptions */

int				nNumLambda,
double		*pdLambdaValue,
double		*pdLambdaTime,

int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
Err amortMidat_modify_dates2(
	int				Nfix,
	long			*fix_pay_datesIn,				/*	Pay date for the amortised swap */

	int				Nfloat,
	long			*float_pay_datesIn,
	
	long			*fix_pay_datesOut,				/*	Pay date for the amortised swap */

	long			*float_pay_datesOut);			/*	Pay date of the amortised swap */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err zcamortMidat_lgmcalibzetalambda(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn[],								/*	Discounted Cash-Flows */
double			zc_last_cpn[],						/*	Discounted Last Fix CF */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn_df[],							/*	Df to cash-flow dates */
double			cpn_cvg[],							/*	cvg from i-1 to i */
int				nex,								/*	Total number of exercise dates */
int				ex_bool[],							/*	Exercise or not */
double			ex_time[],							/*	Exercise times */
int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
int				ex_sncpn[],							/*	Number of coupons in each caplet */
double			ex_lprice[],						/*	Market prices for diagonal */
double			ex_sstrike[],						/*	Strikes for cap */
double			floatcoupon[],						/*	Float Coupon */
double			ex_sprice,							/*	Market price for cap */
double			ex_zeta[],							/*	Output: zetas */
double			floatnotional[],					/*	Float Notional	*/
int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
double			*lambda,							/*	Lambda: may be changed in the process */
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
int				skip_last);							/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */



double amortMidat_lgmsumsopval1F_b(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	int			nex,								/*	Number of exercise dates */
	int			ex_cpn[],							/*	Index of the first cash-flow to be exercised */
	double		df[],								/*	Df to cash flow dates */
	double		cpn_G[],							/*	G at cash-flow dates */
	double		ex_zeta[],							/*	Z at exercise dates */
	double		ex_G[],								/*	G at exercise dates */
	double		swaption_cpn[]);						/*	Standard Swaption Coupons */								


double amortMidat_lgmsumsopval2F_b(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	int			nex,								/*	Number of exercise dates */
	int			ex_cpn[],							/*	Index of the first cash-flow to be exercised */
	double		df[],								/*	Df to cash flow dates */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		cpn_G2[],							/*	G2 at cash-flow dates */
	double		ex_zeta1[],							/*	Z1 at exercise dates */
	double		ex_zeta2[],							/*	Z2 at exercise dates */
	double		ex_zeta12[],							/*	Z12 at exercise dates */
	double		ex_G1[],								/*	G1 at exercise dates */
	double		ex_G2[],								/*	G2 at exercise dates */
	double		swaption_cpn[]);							/*	Standard Swaption Coupons */								


/*	Calibrate lgm on amortised swaptions : main function */
Err amortMidat_cpd_calib_diagonal(
	int				notperiod,
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */
	char			*default_ref_rate_name,					/*	Name of the reference rate */
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	double			vol_shift,
	int				shift_type,						/*	0:	Additive
														1:	Multiplicative */
	int				*ex_bool,						/*	Exercise or not
															NULL = True */
	long			start_date,						/*	Start date of the amortised swap */
	long			end_date,						/*	End date for the amortised swap */
	double			*long_strike,					/*	Diagonal swaption strikes
															NULL = ATM */
	double			*short_strike,					/*	Short swaption strikes
															NULL = ATM */
	int				strike_type,					/*	0: ATM
														1: CASH
														2: SWAP
														3: STD */
	double			*diag_prices,					/* Diagonal Prices */
	double			*ex_fee,						/* Exercise Fees of Diagonal Swaptions */
	double			*fixNotional,						/* Amortised Notional
													NULL = No amortisation */
	double			*floatNotional,						/* Amortised Notional
													NULL = No amortisation */
	double			*margin,						

	double			max_std_short,
	char			*ref_rate_name,					
	char			*swaption_freq,					/*	Frequency, basis and ref. rate of calibrated swaptions */
	char			*swaption_basis,
	int				fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
	int				one_f_equi,						/*	1F equivalent flag:
															if set to 1, then 2F lambda will calibrate
															to the cap priced within calibrated 1F
															with the given lambda */
	int				skip_last,						/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
	double			max_var_jump,					/*	Maximum multiplicative variation in vol */

	double			*lambda,						/*	Lambda: may be changed in the process */
	int				one2F,							/*	Number of factors */
	/*	Alpha, Gamma, Rho (2F only) */
	double			alpha,
	double			gamma,
	double			rho,
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data);					/*	NULL = don't save calibration instrument data */

Err amortMidat_cpd_calib_diagonal_new(
	int				notperiod,
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */

	char			*default_ref_rate_name,			/*	Name of the default reference rate */
	char			*swaption_basis,				/*	 */
	char			*swaption_freq,					/*	 */
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	double			vol_shift,
	int				shift_type,						/*	0:	Additive
														1:	Multiplicative */
	
	int				*ex_bool,						/*	Exercise or not : NULL = True */
	double			*diag_prices,					/* Diagonal Prices */
	double			*ex_fee,						/* Exercise Fees of Diagonal Swaptions */

	char			*fix_freq,
	char			*fix_basis,
	int				Nfix,
	long			*fix_start_dates,				/*	Start date of the amortised swap */
	long			*fix_end_dates,					/*	End date for the amortised swap */
	long			*fix_pay_dates,					/*	Pay date for the amortised swap */
	double			*fix_rates,						/*	Diagonal swaption strikes */
	double			*fix_notionals,					/* Amortised Notional */

	char			*float_freq,
	char			*float_basis,
	int				Nfloat,
	long			*float_start_dates,				/*	Start date of the amortised swap */
	long			*float_end_dates,				/*	End date for the amortised swap */
	long			*float_pay_dates,				/*	Pay date for the amortised swap */
	double			*float_margins,					/*	Margins */
	double			*float_spreads,					/*	Spreads */
	double			*float_notionals,				/*	Amortised Notional */

	/* Calibration Parameters*/
	double			*short_strike,					/*	For Calibration : Short swaption strikes NULL = ATM */
	int				strike_type,					/*	0: ATM
														1: CASH
														2: SWAP
														3: STD
														4: IRR Flat
														5: IRR with vol
														6: SWAPTION */

	double			max_std_short,

	int				fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
	int				one_f_equi,						/*	1F equivalent flag:
															if set to 1, then 2F lambda will calibrate
															to the cap priced within calibrated 1F
															with the given lambda */
	int				skip_last,						/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
	int				use_jump,						/*	Allow Jump (and vol = 0) in vol TS */
	double			max_var_jump,					/*	Maximum multiplicative variation in variance */

	double			*lambda,						/*	Lambda: may be changed in the process */
	int				one2F,							/*	Number of factors */
	/*	Alpha, Gamma, Rho (2F only) */
	double			alpha,
	double			gamma,
	double			rho,

	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig);


//	Calibrate lgm: main function 
Err amortMidat_cpd_calib_diagonal_new_ts(
	int				notperiod,
	char			*yc_name,						
	char			*vol_curve_name,				

	char			*default_ref_rate_name,			
	char			*swaption_basis,				
	char			*swaption_freq,					
	Err				(*get_cash_vol)(				
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	double			vol_shift,
	int				shift_type,						
													
	
	int				*ex_bool,						
	double			*diag_prices,					
	double			*ex_fee,						

	char			*fix_freq,
	char			*fix_basis,
	int				Nfix,
	long			*fix_start_dates,				
	long			*fix_end_dates,					
	long			*fixpaydates,					

	double			*fix_rates,						
	double			*fix_notionals,					

	char			*float_freq,
	char			*float_basis,
	int				Nfloat,
	long			*float_start_dates,				
	long			*float_end_dates,				
	long			*floatpaydates,				

	double			*float_margins,					
	double			*float_spreads,					
	double			*float_notionals,				

	// Calibration Parameters
	double			*short_strike,					
	int				strike_type,					
	double			max_std_short,

	int				fix_lambda,						
	int				one_f_equi,						
	int				skip_last,						
													
	int				use_jump,						
	double			max_var_jump,					

	int				nNum_lambda,						
	double		 *pdLambdaValue,
	double		 *pdLambdaDate,
	
	int				one2F,							
	//	Alpha, Gamma, Rho (2F only) 
	double			alpha,
	double			gamma,
	double			rho,

	int				*num_sig,						
	double			**sig_time,
	double			**sig);


Err ZCamortMidat_cpd_calib_diagonal(
	int				notperiod,
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */
	char			*default_ref_rate_name,					/*	Name of the reference rate */
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	double			vol_shift,
	int				shift_type,						/*	0:	Additive
														1:	Multiplicative */
	int				*ex_bool,						/*	Exercise or not
															NULL = True */
	long			start_date,						/*	Start date of the amortised swap */
	long			end_date,						/*	End date for the amortised swap */
	double			*long_strike,					/*	Diagonal swaption strikes
															NULL = ATM */
	double			*short_strike,					/*	Short swaption strikes
															NULL = ATM */
	int				strike_type,					/*	0: ATM
														1: CASH
														2: SWAP
														3: STD
														4: IRR Flat
														5: IRR with Vol*/
	double			*diag_prices,					/* Diagonal Prices */
	double			*fixNotional,						/* Amortised Notional
													NULL = No amortisation */
	double			*floatNotional,						/* Amortised Notional
													NULL = No amortisation */
	double			*margin,

	double			max_std_short,
	char			*ref_rate_name,					
	char			*swaption_freq,					/*	Frequency, basis and ref. rate of calibrated swaptions */
	char			*swaption_basis,
	int				fix_lambda,						/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
	int				one_f_equi,						/*	1F equivalent flag:
															if set to 1, then 2F lambda will calibrate
															to the cap priced within calibrated 1F
															with the given lambda */
	int				skip_last,						/*	If 1, the last option is disregarded
															and the forward volatility is flat from option
															n-1 */
	double			max_var_jump,					/*	Maximum multiplicative variation in variance */

	double			*lambda,						/*	Lambda: may be changed in the process */
	int				one2F,							/*	Number of factors */
	/*	Alpha, Gamma, Rho (2F only) */
	double			alpha,
	double			gamma,
	double			rho,
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data);					/*	NULL = don't save calibration instrument data */

Err europSwaption_clsdfrm(
							   char			*yc_name,
							   SrtUndPtr	und,
							   int			notperiod,
							   long			lStartDate,
							   long			lEndDate,
							   char			*fixFreq,
							   char			*fixBasis,
							   char			*refRateName,
							   double		exer_fee,
							   double		strike,
							   double		margin,
							   char			*recPayStr, 
							   double		*price);

Err europAmortSwaption_clsdfrm_new(
							   char			*yc_name,
							   SrtUndPtr	und,
							   int			notperiod,

							   char			*fixbasis,
							   int			NFix,
							   long			*FixStartDates,
							   long			*FixEndDates,
							   double		*FixRates,
							   double		*FixNotionals,
							   
							   char			*refrate,
							   int			NFloat,
							   long			*FloatStartDates,
							   long			*FloatEndDates,
							   double		*FloatMargins,
							   double		*FloatSpreads,
							   double		*FloatNotionals,

							   double		exer_fee,
							   char			*recPayStr, 
							   double		*price);

Err europAmortSwaption_clsdfrm(
							   char *yc_name,
							   SrtUndPtr und,
							   int notperiod,
							   long lStartDate,
							   long lEndDate,
							   char *fixFreq,
							   char *fixBasis,
							   char *refRateName,
							   double exer_fee,
							   int Nfix,
							   double *fixNotional,
							   double *fixRates,
							   int Nfloat,
							   double *floatNotional,
							   double *margin,
							   char	*recPayStr, 
							   double *price);

Err europZCSwaption_clsdfrm(
							   char *yc_name,
							   SrtUndPtr und,
							   int notperiod,
							   long lStartDate,
							   long lEndDate,
							   char *fixFreq,
							   char *fixBasis,
							   char *refRateName,
							   int Nfix,
							   double *fixNotional,
							   double *fixRates,
							   int Nfloat,
							   double *floatNotional,
							   double *margin,
							   double *price);

Err amortMidat_modify_dates(
	int				Nfix,
	long			*fix_start_datesIn,				/*	Start date of the amortised swap */
	long			*fix_end_datesIn,					/*	End date for the amortised swap */

	int				Nfloat,
	long			*float_start_datesIn,			/*	Start date of the amortised swap */
	long			*float_end_datesIn,
	
	long			*fix_start_datesOut,				/*	Start date of the amortised swap */
	long			*fix_end_datesOut,				/*	End date for the amortised swap */

	long			*float_start_datesOut,			/*	Start date of the amortised swap */
	long			*float_end_datesOut	);				/*	End date for the amortised swap */


double amortMidat_lgmamortsopval1F(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	double		cpn[],								/*	Notional */
	double		df[],								/*	Df to cash flow dates */
	double		cvg[],								/*	cvg from i-1 to i */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		ex_zeta1,							/*	Z1 at exercise date */
	double		ex_G1);								/*	G1 at exercise date */

double amortMidat_lgmamortsopval2F(
	int			ncpn,								/*	Number of cash-flow dates, including
															start and end date */
	double		cpn[],								/*	Notional */
	double		df[],								/*	Df to cash flow dates */
	double		cvg[],								/*	cvg from i-1 to i */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		cpn_G2[],							/*	G2 at cash-flow dates */
	double		ex_zeta1,							/*	Z1 at exercise date */
	double		ex_zeta2,							/*	Z2 at exercise date */
	double		ex_zeta12,							/*	Z12 at exercise date */
	double		ex_G1,								/*	G1 at exercise date */
	double		ex_G2);								/*	G2 at exercise date */

Err Get_LGM_TermStructureOneOrTwoFact(char *underlying, double **sigma_date, double **sigma, long *sigma_n,
						  double **tau_date, double **tau, long *tau_n,
						  double *alpha, double *gamma, double *rho);

#endif