
#include "srt_h_all.h"
#include "opfnctns.h"
#include "math.h"


#include "LGMSVCalib.h"
#include "LGMSVClosedForm.h"
#include "LGMSVGrfn.h"
#include "LGMSVMC.h"
#include "LGMSVPDE.h"
#include "LGMSVUtil.h"
#include "LGMSVClosedFormApprox.h"
#include "LGMSVCalibApprox.h"

#include "CTSProdStruct.h"

#define	LGMSV_MINVOL		0.0001
#define	LGMSV_MAXVOL		0.1
#define	MAX_EURO			200
#define	BUMP_VOL			0.0005

static	CTS_UND				stat_und;
static	CTS					stat_market_cts,
							stat_model_cts;
							
static	int					stat_for_fund,
							stat_calc_fwd_iv,
							stat_adj_fee,
							stat_nb_euro;
							
static	CTS_ADI_ARG			stat_adi_arg;

static	LGMSV_NUMERPARAMS	stat_num_params;

static	double				stat_price[MAX_EURO],
							*stat_target,
							stat_sensi[MAX_EURO][MAX_EURO];
							
static	double				stat_call,
							stat_euro[MAX_EURO];

#define OPT_EPS				1.0e-15
#define	VOL_MULT			5.0
#define	VOL_MULT_ITER		5

Err	static_cts_price_european(double	index,
							  double	vol[],
							  double	*price,
							  double	*gradient,
							  int		nvol)
{
	static int		i, j, ind;
	static double	sum2;
	Err		err	= NULL;

	ind = (int) (index);	

	if (ind == 0)
	{
		/* We calculate everything */

		/* Initialise the modell */

		for (i=0; i<nvol; i++)
		{
			stat_und->model.dSigma[i] = vol[i+1];

			if (stat_und->model.dSigma[i] < LGMSV_MINVOL)
			{
				stat_und->model.dSigma[i] = LGMSV_MINVOL;
			}
			if (stat_und->model.dSigma[i] > LGMSV_MAXVOL)
			{
				stat_und->model.dSigma[i] = LGMSV_MAXVOL;
			}
		}
				
		/* Calculate Sensi */
		for (j=0; j<nvol; j++)
		{
			stat_und->model.dSigma[j] += BUMP_VOL;			
			
			/* prices */
			err = cts_calc_model_fwd_iv(stat_und,
										stat_model_cts,
										stat_calc_fwd_iv,
										stat_adj_fee,
										stat_num_params);

			if (err) goto FREE_RETURN;

			err = cts_adjust_model_fwd_iv(	stat_und,
											stat_market_cts,
											stat_model_cts,
											stat_for_fund,
											stat_calc_fwd_iv,
											stat_adj_fee,
											0,
											0,
											NULL);

			if (err) goto FREE_RETURN;

			err = cts_launch_adi (stat_model_cts, stat_und, stat_adi_arg, &stat_call, &(stat_euro[0]));

			if (err) goto FREE_RETURN;
			
			for (i=0; i<stat_nb_euro; i++)
			{				
				stat_sensi[i][j] = stat_euro[i];
			}

			stat_und->model.dSigma[j] -= BUMP_VOL;
		}
		
		err = cts_calc_model_fwd_iv(stat_und,
									stat_model_cts,
									stat_calc_fwd_iv,
									stat_adj_fee,
									stat_num_params);

		if (err) goto FREE_RETURN;

		err = cts_adjust_model_fwd_iv(	stat_und,
										stat_market_cts,
										stat_model_cts,
										stat_for_fund,
										stat_calc_fwd_iv,
										stat_adj_fee,
										0,
										0,
										NULL);

		if (err) goto FREE_RETURN;

		err = cts_launch_adi (stat_model_cts, stat_und, stat_adi_arg, &stat_call, &(stat_euro[0]));

		if (err) goto FREE_RETURN;
		
		/* Calculate Prices */

		sum2 = 0.0;

		for (i=0; i<stat_nb_euro; i++)
		{						
			for (j=0; j<nvol; j++)
			{
				stat_sensi[i][j] = (stat_sensi[i][j] - stat_euro[i]) / BUMP_VOL;
			}

			stat_price[i] = stat_euro[i] - stat_target[i];

			sum2 += stat_price[i] * stat_price[i];
		}

		sum2 = sqrt(sum2);
	}
		
	/* Return the asked result */	

	(*price) = stat_price[ind];

	for (j=0; j<nvol; j++)
	{
		gradient[j+1] = stat_sensi[ind][j];
	}

FREE_RETURN:

	return err;
}

Err	cts_calibrate_european(	CTS_UND				und,
							CTS					market_cts,
							CTS					model_cts,
							int					for_fund,
							int					calc_fwd_iv,
							int					adj_fee,

							CTS_ADI_ARG			adi_arg,
							LGMSV_NUMERPARAMS	numer_params,

							int					nb_european,
							double				*european_values,
												
							int					nb_iter,
												
							double				*fitting_error)
{
Err	err = NULL;

double	data[MAX_CPN],
		weight[MAX_CPN],
		target[MAX_CPN],
		vol[MAX_CPN];

int		i;
	
	for (i=0; i<nb_european; i++)
	{		
		data[i+1] = i;
		weight[i+1] = 1.0 / european_values[i];
		target[i+1] = 0.0;
	}

	for (i=0; i<und->model.iNbPWTime; i++)
	{
		vol[i+1] = und->model.dSigma[i];
	}

	stat_und = und;
	stat_market_cts = market_cts;
	stat_model_cts = model_cts;
	stat_for_fund = for_fund;
	stat_calc_fwd_iv = calc_fwd_iv;
	stat_adj_fee = adj_fee;
	stat_adi_arg = adi_arg;	
	stat_num_params = numer_params;

	stat_nb_euro = nb_european;
	stat_target = european_values;

	err  = levenberg_marquardt(
								data,
								target,
								weight,
								nb_european,
								vol,
								und->model.iNbPWTime,
								nb_iter,
								static_cts_price_european,
								fitting_error);

	return err;
}

static double optval(double fwd,double strike,double vol,double mat,double disc_fact,
					SrtCallPutType call_put, SrtGreekType greek, SrtDiffusionType log_or_norm)
{
double premium;

	if (log_or_norm == SRT_LOGNORMAL)  
	{    
		premium = srt_f_optblksch(fwd,strike,vol,mat,disc_fact,call_put,greek);  
	}  
	else    
	{    
		premium = srt_f_optblknrm(fwd,strike,vol,mat,disc_fact,call_put,greek);  
	}

	return premium;
}

Err 	cts_optimpvol(
		double 	           premium, 
		double 	           fwd_price, 
		double 	           strike,
		double 	           mat, 
		double 	           disc, 
		SrtCallPutType 	   call_put,
		SrtDiffusionType   lognormal_normal,
		double             *implied_vol)
{ 	
int 	  i;
double    vol_up, vol_down, vol_middle, vol_temp, vol_acc, vol_shift;
double	  vol_diffold, vol_diff;
double	  store_fwd_price;
double 	  deriv;
double	  prem_new, prem_up, prem_down;
double    intrinsic;

	/* value by default */
	*implied_vol = 0.0;

	/* Compute option intrinsic value */   
	intrinsic  = optval(	fwd_price,
					strike,
					NULL_VOL,
					mat, 
					disc, 
					call_put, 
					PREMIUM,
					lognormal_normal);  


	/* Checks target premium is above intrinsic value */   
	if (intrinsic > premium + OPT_EPS)
	{
		return serror ("Intrinsic higher than option premium");
	}
	
	/* Renormalise the option price by the spot value and factorise df */
	if ( fwd_price != 0.00)
	{
		store_fwd_price = fabs(fwd_price);
		strike /= store_fwd_price;
		premium /= store_fwd_price;
		fwd_price = fwd_price / store_fwd_price;

		premium /= disc;
		disc = 1.0;
	}

	/* Sets initial guess for vol */
	vol_up = 20;
	vol_down =  0.000000001;
	vol_shift = 0.0000000001;
	vol_acc = 0.000001;

	prem_up =  optval(fwd_price,
						strike,
						vol_up,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal);

	prem_down = optval(fwd_price,
						strike,
						vol_down,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal);

	if ((fabs(prem_up-premium) < DBL_EPSILON) && (fabs(prem_down-premium) < DBL_EPSILON))
	{
		return serror("Too many solution in range");	
	}

	if (fabs(prem_down-premium) < DBL_EPSILON)
	{
		*implied_vol = vol_down;
		if (lognormal_normal == SRT_NORMAL)
			*implied_vol *= store_fwd_price;
		return NULL;
	}

	if (fabs(prem_up-premium) < DBL_EPSILON)
	{
		*implied_vol = vol_up;
		if (lognormal_normal == SRT_NORMAL)
			*implied_vol *= store_fwd_price;
		return NULL;
	}

	/* same sign for the interval */
	if ((prem_down-premium) * (prem_up-premium) > 0.0)
	{		
		/* modification for very low forward in JPY */
		if (prem_up < premium)
		{
			/* one more try with vol_up * 5.0 */
			i = 0;
			while (i < VOL_MULT_ITER && prem_up < premium)
			{
				vol_up *= VOL_MULT;
				prem_up =  optval(fwd_price,
						strike,
						vol_up,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal);
				i++;
			}

			if (prem_up < premium)
			{
				return serror("No solution in range");
			}
		}
		else
		{			
			return serror("No solution in range");
		}
	}

	/* orientation of the search */
	if (prem_up < 0)
	{
		vol_middle = vol_up;
		vol_up = vol_down;
		vol_down = vol_middle;
	}

	vol_middle = 0.5 * (vol_up + vol_down);
	vol_diffold = fabs(vol_up - vol_down);
	vol_diff = vol_diffold;

	prem_new = optval(fwd_price,
						strike,
						vol_middle,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal);
  
	deriv = (optval(fwd_price,
						strike,
						vol_middle + vol_shift,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal) - prem_new) / vol_shift;


	for (i = 0; i <= MAX_ITER; i++)
	{
		if ( (((vol_middle - vol_up) * deriv - (prem_new - premium)) * ((vol_middle - vol_down) * deriv - (prem_new - premium)) >= 0.0) 
			|| (fabs(2.0*(prem_new-premium)) > fabs(vol_diffold * deriv)))
		{
			/* bissection if Newton is out of range, or not decreasing fast enough */
			vol_diffold = vol_diff;
			vol_diff = 0.5 * (vol_up - vol_down);
			vol_middle = vol_down + vol_diff;
			if (vol_down == vol_middle) /* The change is negligible */
			{
				*implied_vol = vol_middle;
				if (lognormal_normal == SRT_NORMAL)
					*implied_vol *= store_fwd_price;
				return NULL;
			}
		}
		else
		{
			/* the change in newton is acceptable, take it */
			vol_diffold = vol_diff;
			vol_diff = (prem_new - premium) / deriv;
			vol_temp = vol_middle;
			vol_middle -= vol_diff;
			if (vol_temp == vol_middle)
			{
				*implied_vol = vol_middle;
				if (lognormal_normal == SRT_NORMAL)
					*implied_vol *= store_fwd_price;
				return NULL;
			}
		}

		if (fabs(vol_diff) < vol_acc)
		{
			*implied_vol = vol_middle;
			if (lognormal_normal == SRT_NORMAL)
				*implied_vol *= store_fwd_price;
			return NULL;
		}

		prem_new = optval(fwd_price,
						strike,
						vol_middle,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal);  

		deriv = (optval(fwd_price,
						strike,
						vol_middle + vol_shift,
						mat, 
						disc, 
						call_put, 
						PREMIUM,
						lognormal_normal) - prem_new) / vol_shift;

		/* maintain the bracket on the root */
		if ((prem_new-premium) < 0.0)
		{
			vol_down = vol_middle;
		}
		else
		{
			vol_up = vol_middle;
		}
	}

	*implied_vol = 0.0;

	return NULL;

} /* srt_f_optimpvol() */