#ifndef LGMSVCALIBH
#define	LGMSVCALIBH

#include "srt_h_all.h"
#include "cpdcalib.h"
#include "lgmsvpde.h"

#define MAX_NCPN	512

/*	Underlying term structs */
typedef struct
{
	int		ncpn;
	int		start_index;
	double	cpn_time[MAX_NCPN];
	double	cpn_df[MAX_NCPN];
	double	cpn_cvg[MAX_NCPN];
	double	gam[MAX_NCPN];
	double	gam2[MAX_NCPN];
	double	strike;

} lgmsv_swap_param, *LGMSV_SWAP_PARAM;


Err lgmSVprcapgivenlambda_num(
			int				ncpn,								/*	Total number of cash-flow dates */
			double			cpn_time[],							/*	Cash-Flow times */
			double			cpn_df[],							/*	Df to cash-flow dates */
			double			cpn_cvg[],							/*	cvg from i-1 to i */
			int				nex,								/*	Total number of exercise dates */
			double			ex_time[],							/*	Exercise times */
			int				ex_cpn[],							/*	Index of the first cash-flow to be exercised */
			int				ex_endcpn[],						/*	Index of the last cash-flow to be exercised */
			int				ex_sncpn[],							/*	Number of coupons in each caplet */
			double			ex_lstrike[],						/*	Strikes for diagonal */
			double			ex_lprice[],						/*	Market prices for diagonal */
			double			ex_sstrike[],						/*	Strikes for cap */
			double			vol[],								/*	Output: volatility structure */
			double			lambda,								/*	Lambda */

			/*	Market */
			long			today,
			char			*cYieldCurve,
			
			/*	Lambda Alpha and Rho */
			double			alpha,
			double			lambdaEps,
			double			rho,
						
			int				skip_last,							/*	If 1, the last option is disregarded
																		and the forward volatility is flat from option
																		n-1 */
			int				price_cap,							/*	0: just calibrate */
			double			*ex_sprice,							/*	Cap price as output */
			
			/*	Space Discretisation	*/
			int				nbT,
			int				nbPhi,
			int				nbX,
			int				nbEps,
			
			/*	Newton Parameters		*/
			int				nbIter,
			double			precision,
			
			/*	Model Params */
			LGMSVParam		ModelParams
			);

/*	Calibrate lgm: main function */
Err lgmsv_calib_diagonal(
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */
	char			*ref_rate_name,					/*	Name of the reference rate */
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
						/*	If ex_date is NULL,
						exercise dates will be generated 2bd before start */
	int				num_ex_dates,					/*	Exercise dates, 
															all supposed to be on or after today */
	long			*ex_date,						/*	Supposed to be sorted 
															NULL = 2bd before each coupon */
	char			**end_tenor,					/*	Tenors of the underlying instruments
														or "DIAG" */
	long			end_date,						/*	End date for diagonal */
	double			*long_strike,					/*	Diagonal swaption strikes
															NULL = ATM */
	double			*short_strike,					/*	Short swaption strikes
															NULL = ATM */
	int				strike_type,					/*	0: ATM
														1: CASH
														2: SWAP
														3: STD */
	double			max_std_long,
	double			max_std_short,
	char			*swaption_freq,					/*	Frequency and basis of underlying swaptions */
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
	double			*lambda,						/*	Lambda: may be changed in the process */

	/*	Alpha, Gamma, Rho (2F only) */
	double			alpha,
	double			lambdaEps,
	double			rho,
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,

	/*	Space Discretisation	*/
	int				nbT,
	int				nbPhi,
	int				nbX,
	int				nbEps,

	/*	Newton Parameters		*/
	int				nbIterMax,
	double			precision,

	/*	Model Params */
	LGMSVParam		ModelParams,

	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data);					/*	NULL = don't save calibration instrument data */

#endif