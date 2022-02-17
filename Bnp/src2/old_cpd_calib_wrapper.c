// Old cpd_calib_diagonal wrapper -> gate to cpd_calib_diagonal_dlm
// Only works in the case fix_lambda = 1 for now

#include "srt_h_all.h"
#include "CPDCalib.h"
#include "DiagCalibDLM.h"
#include "old_cpd_calib_wrapper.h"

Err old_cpd_calib_diagonal_wrapper(
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
	double			min_fact,						/*	Maximum down jump on variance */
	double			max_fact,						/*	Maximum up jump on variance */
	int				use_jumps,						/*	Allow vol term structure to jump */
	int				proba_weight,					/*	Proba weighting for caplet */
	double			*proba,

	double			*lambda,						/*	Lambda: may NOT be changed in the process */
	int				one2F,							/*	Number of factors */
	/*	Alpha, Gamma, Rho (2F only) */
	double			alpha,
	double			gamma,
	double			rho,
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data)					/*	NULL = don't save calibration instrument data */
{
	Err						err = NULL;
	int						i, *cal_datel = NULL;
	char					**end_tenorl = NULL;
	double					*strike = NULL;
	long					*used_ex_date = NULL;
	double					lam_time = 1.0, lam_shift = 0.0;
	cpd_diag_calib_param	calib_param;
	diag_calib_lm_params	lm_params;

	SrtCurvePtr		yc_ptr;
	SrtCompounding	ifreq;
	SrtBasisCode	ibasis;
	long			today, theo_date, act_date;
	int				ncpn;

	if (fix_lambda == 0)
		return serror("old_cpd_calib_diagonal_wrapper can only be used with fix_lambda = 1 for now");

	/* Get the exercise dates */
	if (ex_date)
	{
		/* just recopy */
		used_ex_date = calloc(num_ex_dates, sizeof(long));

		if (!used_ex_date)
		{
			err = "Memory allocation faillure in old_cpd_calib_diagonal_wrapper";
			return err;
		}

		memcpy(used_ex_date, ex_date, num_ex_dates * sizeof(long));
	}
	else
	{
		/* need to create the schedule ! */

		yc_ptr = lookup_curve (yc_name);

		if (!yc_ptr)
		{
			err = "Yield Curve not found";
			return err;
		}

		today = get_today_from_curve (yc_ptr);

		err = interp_compounding (swaption_freq, &ifreq);
		if (err) return err;

		err = interp_basis (swaption_basis, &ibasis);
		if (err) return err;

		theo_date = end_date;
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);

		ncpn = 1;

		while (act_date > today)
		{
			theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
			act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
			ncpn++;
		}

		ncpn--;

		if (ncpn < 2)
		{
			err = "Not enough coupons in cpd_calib_diagonal";
			return err;		
		}

		num_ex_dates = ncpn - 1;

		used_ex_date = calloc(num_ex_dates, sizeof(long));

		if (!used_ex_date)
		{
			err = "Memory allocation faillure in old_cpd_calib_diagonal_wrapper";
			return err;
		}

		theo_date = end_date;
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		i = ncpn - 1;		

		while (i >= 0)
		{
			if (i < num_ex_dates)
			{
				used_ex_date[i] = add_unit (act_date, - 2, SRT_BDAY, MODIFIED_SUCCEEDING);		
			}

			theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
			act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);

			i--;
		}		
	}

	cal_datel = (int *) calloc(num_ex_dates, sizeof(int));
	end_tenorl = svector_size(0, num_ex_dates-1, 5);
	strike = calloc(num_ex_dates, sizeof(double));

	for (i=0; i < num_ex_dates; i++)
	{
		cal_datel[i] = 1;
		strcpy(end_tenorl[i], "DIAG");
		strike[i] = 0.0;
	}

	cpd_calib_set_default_param(&calib_param);

	calib_param.vol_shift = vol_shift;
	calib_param.shift_type = shift_type;
	calib_param.strike_type = strike_type;
	calib_param.max_std = max_std_long;
	calib_param.skip_last = skip_last;
	calib_param.min_fact = min_fact;
	calib_param.max_fact = max_fact;
	calib_param.use_jumps = use_jumps;
	calib_param.precision = 0.00001;

	diag_calib_lm_params_set_default_param(&lm_params);

	err = cpd_calib_diagonal_dlm(
								yc_name,
								vol_curve_name,
								get_cash_vol,
								ref_rate_name,
								swaption_freq,
								swaption_basis,
								ref_rate_name,
								num_ex_dates,
								used_ex_date,
								cal_datel,
								end_tenorl,
								end_date,
								(long_strike ? long_strike : strike),
								&calib_param,
								NULL,
								NULL,
								NULL,
								0,
								NULL,
								NULL,
								NULL,
								0,
								NULL,
								NULL,
								NULL,
								1,
								one_f_equi,
								1,
								&lam_time,
								lambda,
								&lam_shift,
								one2F,
								alpha,
								gamma,
								rho,

								0,
								NULL,
								NULL,
								0,
								NULL,
								NULL,

								num_sig,
								sig_time,
								sig,
								&lm_params,
								inst_data );

	free(used_ex_date);
	free(cal_datel);
	free_svector_size(end_tenorl, 0, num_ex_dates-1, 5);
	free(strike);

	return err;
}
