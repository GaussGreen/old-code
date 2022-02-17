
#include "srt_h_all.h"
#include "srt_h_lgmtypes.h"
#include "srt_h_lgmprotos.h"
#include "opfnctns.h"
#include "srt_h_allFx3F.h"
#include "swp_h_swap_pricing.h"
#include "diagCalibDLM.h"
#include "math.h"


#define			BUMP_LAM_LM		5.0E-04
#define			BUMP_LAM_NEWTON	1.0E-02
#define			LIMIT_LAM_SV	-0.30

static int		NB_MOM_LM,
				INDEX_MOM[MAX_INST],
				BREAK_INST,
				FREQ_SHORT_LM,
				SHIFT_FREQ_LM,
				NEXS_LM,
				MAX_ITER_NEWTON;

static double	*EX_PRICES_LM;

static double	PRICE_LM[MAX_CPN],
				SENSI_LM[MAX_CPN][MAX_CPN],
				PRECISION_LM;


static Err		(*CALIBRATION_FUNCTION)(void *AllPrimaryInst,
										void *Model,
										void *Params),

				(*PRICING_FUNCTION)(void *AllSecondaryInst,
									void *Model,
									void *Params,
									int		iIndex,
									double	*dPrice),

				(*BUMPINGMODEL_FUNCTION)(void	*Model,
										int		iIndex,
										double	NewValue);

void			*ALLPRIMINST,
				*ALLSECINST,
				*MODEL_LM,
				*PARAMS;

				
void init_static_moment_dlm(void	*AllPrimaryInst,
							void	*AllSecondaryInst,
							void	*Model,
							void	*Params,

							Err		(*CalibrationFunction)(	void *AllPrimaryInst,
															void *Model,
															void *Params),

							Err		(*PricingFunction)(	void *AllSecondaryInst,
														void *Model,
														void *Params,
														int		iIndex,
														double	*dPrice),

							Err		(*BumpingModelFunction)(	void	*Model,
																int		iIndex,
																double	NewValue),

							int							nb_exe_dates,
							double						*target_prices,
															
							DIAG_CALIB_LM_PARAMS		ln_params)							
{
int	i, nbm;
	
	NEXS_LM = nb_exe_dates;
	EX_PRICES_LM = target_prices;
	
	ALLPRIMINST = AllPrimaryInst;
	ALLSECINST = AllSecondaryInst;
	MODEL_LM = Model;
	PARAMS = Params;

	CALIBRATION_FUNCTION = CalibrationFunction;
	PRICING_FUNCTION = PricingFunction;
	BUMPINGMODEL_FUNCTION = BumpingModelFunction;

	NB_MOM_LM = ln_params->nb_moment;
	FREQ_SHORT_LM = ln_params->freq_short;
	SHIFT_FREQ_LM = ln_params->shift_freq;
	PRECISION_LM = ln_params->precision;
	MAX_ITER_NEWTON = ln_params->nb_iter;

	BREAK_INST = ln_params->break_moment;

	if (BREAK_INST)
	{
		nbm = (int) (nb_exe_dates / NB_MOM_LM + 0.5);		

		INDEX_MOM[0] = 0;

		for (i=1; i<NB_MOM_LM; i++)
		{
			INDEX_MOM[i] = INDEX_MOM[i-1] + nbm;
		}

		INDEX_MOM[NB_MOM_LM] = nb_exe_dates - 1;
	}
}

/* This method try to match the momentum */
Err	static_pricer_momentum_dlm(double	index,
										  double	lam[],
										  double	*price,
										  double	*gradient,
										  int		nlam)
{
	static int	i, j, k, ind, pos;
	Err		err;
	double	*lambda = &(lam[1]);
	static double temp, temp_pow[MAX_CPN];

	ind = (int) (index);

	if (ind == 0)
	{		
		/* Calculate Sensi */
		for (j=0; j<nlam; j++)
		{
			lambda[j] += BUMP_LAM_LM;

			err = BUMPINGMODEL_FUNCTION(MODEL_LM,
										j,
										lambda[j]);

			if (err)
			{
				return err;
			}

			err = CALIBRATION_FUNCTION(	ALLPRIMINST,
										MODEL_LM,
										PARAMS);

			if (err)
			{
				return err;
			}
			
			for (i=0; i<NEXS_LM; i++)
			{
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					err = PRICING_FUNCTION(	ALLSECINST, MODEL_LM, PARAMS, i, &SENSI_LM[i][j]);

					if (err)
					{
						return err;
					}
											
					SENSI_LM[i][j] = (SENSI_LM[i][j] - EX_PRICES_LM[i]);
				}
			}

			lambda[j] -= BUMP_LAM_LM;

			err = BUMPINGMODEL_FUNCTION(MODEL_LM,
										j,
										lambda[j]);

			if (err)
			{
				return err;
			}
		}

		err = CALIBRATION_FUNCTION(	ALLPRIMINST,
										MODEL_LM,
										PARAMS);

		if (err)
		{
			return err;
		}		
		
		/* Calculate Prices */
		for (i=0; i<NEXS_LM; i++)
		{
			if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
			{
				err = PRICING_FUNCTION(	ALLSECINST, MODEL_LM, PARAMS, i, &PRICE_LM[i]);

				if (err)
				{
					return err;
				}

				PRICE_LM[i] = (PRICE_LM[i] - EX_PRICES_LM[i]);
			}
		}		

		/* Computes the momentum */
		if (BREAK_INST)
		{
			pos = -1;

			for (k=0; k<NB_MOM_LM; k++)
			{
				temp = 0.0;

				for (i=pos+1; i<=INDEX_MOM[k+1]; i++)
				{
					temp += PRICE_LM[i];
				}

				PRICE_LM[NEXS_LM + k] = temp;
				
				pos = i-1;
			}

			for (j=0; j<nlam; j++)
			{
				pos = -1;

				for (k=0; k<NB_MOM_LM; k++)
				{
					temp = 0.0;

					for (i=pos+1; i<=INDEX_MOM[k+1]; i++)
					{
						temp += SENSI_LM[i][j];
					}

					SENSI_LM[NEXS_LM + k][j] = (temp - PRICE_LM[NEXS_LM + k]) / BUMP_LAM_LM;
					
					pos = i-1;
				}				
			}
		}
		else
		{
			for (k=0; k<NB_MOM_LM; k++)
			{
				temp = 0.0;

				for (i=0; i<NEXS_LM; i++)
				{
					if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
					{
						temp += pow(PRICE_LM[i], k + 1);
					}
				}
				
				PRICE_LM[NEXS_LM + k] = temp / fabs(temp) * exp(log(fabs(temp)) / (k+1));
			}

			for (j=0; j<nlam; j++)
			{
				for (k=0; k<NB_MOM_LM; k++)
				{
					temp = 0.0;				

					for (i=0; i<NEXS_LM; i++)
					{
						if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
						{
							temp += pow(SENSI_LM[i][j], k + 1);
						}
					}				

					SENSI_LM[NEXS_LM + k][j] = (temp / fabs(temp) * exp(log(fabs(temp)) / (k+1)) - PRICE_LM[NEXS_LM + k]) / BUMP_LAM_LM;
				}
			}
		}
	}

	/* Return the asked result */

	(*price) = PRICE_LM[NEXS_LM + ind];

	for (j=0; j<nlam; j++)
	{
		gradient[j+1] = SENSI_LM[NEXS_LM + ind][j];
	}

	return NULL;
}

static Err	calib_lambda_simpleNewton_dlm(	double	*lam,
											double	*price)
{
	int		i, j, k, l;	
	double	temp, price1, price2;
	double	lam1, lam2;
	double	**res_iter = NULL;
	Err		err;
	
	res_iter = dmatrix(0, MAX_ITER_NEWTON, 0, 1);	

	if (!res_iter)
	{
		err = "Memory allocation faillure in lgmprcapgiventauts_momentum_dlm";
		goto FREE_RETURN;
	}

	lam1 = *lam;
	
	if (lam1 < LIMIT_LAM_SV)
	{
		lam1 = LIMIT_LAM_SV;
	}
	
	/* First guess */
	err = BUMPINGMODEL_FUNCTION(MODEL_LM,
								0,
								lam1);

	err = CALIBRATION_FUNCTION(	ALLPRIMINST,
								MODEL_LM,
								PARAMS);

	if (err)
	{
		return err;
	}		
	
	/* Calculate Prices */
	
	price1 = 0.0;

	for (i=0; i<NEXS_LM; i++)
	{
		if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
		{
			err = PRICING_FUNCTION(	ALLSECINST, MODEL_LM, PARAMS, i, &temp);

			if (err)
			{
				return err;
			}

			price1 += (temp - EX_PRICES_LM[i]);
		}
	}

	if (fabs(price1) < PRECISION_LM)
	{
		*lam = lam1;
		*price = price1;
		goto FREE_RETURN;
	}

	res_iter[0][0] = lam1;
	res_iter[0][1] = price1;

	/* second guess */
	lam2 = lam1 + BUMP_LAM_NEWTON;				

	k = 1;

	while (fabs(price1) > PRECISION_LM && k < MAX_ITER_NEWTON)
	{
		k++;

		err = BUMPINGMODEL_FUNCTION(MODEL_LM,
									0,
									lam2);

		err = CALIBRATION_FUNCTION(	ALLPRIMINST,
									MODEL_LM,
									PARAMS);

		if (err)
		{
			return err;
		}		
		
		/* Calculate Prices */
		
		price2 = 0.0;

		for (i=0; i<NEXS_LM; i++)
		{
			if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
			{
				err = PRICING_FUNCTION(	ALLSECINST, MODEL_LM, PARAMS, i, &temp);

				if (err)
				{
					return err;
				}

				price2 += (temp - EX_PRICES_LM[i]);
			}
		}

		/* Save Res */
		l = 0;
		while (l < k-1 && res_iter[l][1] < price2)
		{
			l++;
		}

		if (l < k-1)
		{
			for (j=k-2; j>=l; j--)
			{
				res_iter[j+1][0] = res_iter[j][0];
				res_iter[j+1][1] = res_iter[j][1];
			}
		}

		res_iter[l][0] = lam2;
		res_iter[l][1] = price2;

		lam1 = lam2;
		price1 = price2;

		lam2 = solve_for_next_coef_dlm(	res_iter,
										k,
										0.0,
										2);

		if (lam2 < LIMIT_LAM_SV)
		{
			smessage("DiagCalibSVDLM: Limit lambda has been reached, calibration stopped");
			lam2 = LIMIT_LAM_SV;
			break;
		}
	}
	
	/* Return the asked result */

	*lam = lam2;
	*price = price2;

FREE_RETURN:

	if (res_iter) free_dmatrix(res_iter, 0, MAX_ITER_NEWTON, 0, 1);

	return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err LGMSVcalibration_dlm(

int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
int				nlam,								/*	Lambda TS: may NOT be changed in the process */
double			lam_time[],
double			lam[],

void	*AllPrimaryInst,
void	*AllSecondaryInst,
void	*Model,
void	*Params,

Err		(*CalibrationFunction)(	void *AllPrimaryInst,
								void *Model,
								void *Params),

Err		(*PricingFunction)(	void *AllSecondaryInst,
							void *Model,
							void *Params,
							int		iIndex,
							double	*dPrice),

Err		(*BumpingModelFunction)(	void	*Model,
									int		iIndex,
									double	NewValue),

int							nb_exe_dates,
double						*target_prices,
double						*target_vegas,
								
DIAG_CALIB_LM_PARAMS		lm_params)
{
	int	i;	
	Err err = NULL;
	double	data[MAX_CPN],
			weight[MAX_CPN],
			target[MAX_CPN],
			lambda[MAX_CPN];
	
	double	fitting_error;

	if (!fix_lambda && nb_exe_dates < 1)
	{
		return serror ("Cannot calibrate lambda and zeta with less than 2 exercises - choose fix lambda");
	}

	if (fix_lambda)
	{
		err = CalibrationFunction(AllPrimaryInst, Model, Params);
		
		return err;
	}
	else
	{
		for (i=1; i<=nlam; i++)
		{			
			lambda[i] = lam[i-1];
		}

		lm_params->shift_freq = (int) (fmod(lm_params->shift_freq, nb_exe_dates) + 0.5);

		if (lm_params->use_moment)
		{
			/* Momentum method */
			lambda[0] = lm_params->nb_moment;

			for (i=1; i<=lm_params->nb_moment; i++)
			{					
				data[i] = i-1;
				if  (lm_params->vega_weight)
				{
					weight[i] = (i * 1.0);
				}
				else
				{
					weight[i] = 1.0;
				}

				target[i] = 0.0;
			}

			init_static_moment_dlm(	AllPrimaryInst,
									AllSecondaryInst,
									Model,
									Params,
									CalibrationFunction,
									PricingFunction,
									BumpingModelFunction,
									nb_exe_dates,
									target_prices,
									lm_params);

			if (nlam == 1)
			{
				err = calib_lambda_simpleNewton_dlm(&lambda[1], &fitting_error);

				if (err) return err;

				smessage ("Calibration of lambda, error in bp: %.00f", 10000 * fitting_error);
			}
			else
			{
				err  = levenberg_marquardt(
									data,
									target,
									weight,
									lm_params->nb_moment,
									lambda,
									nlam,
									lm_params->nb_iter,
									static_pricer_momentum_dlm,
									&fitting_error);

				if (err) return err;
			}
		}
		else
		{
			for (i=1; i<=nb_exe_dates; i++)
			{					
				data[i] = (int) (lm_params->shift_freq + (i - 1) * lm_params->freq_short + 0.5);
				if  (lm_params->vega_weight)
				{
					weight[i] = target_vegas[(int)data[i]];
				}
				else
				{
					weight[i] = 1.0;
				}
				target[i] = target_prices[(int)data[i]];
			}			
		}
		
		for (i=1; i<=nlam; i++)
		{
			lam[i-1] = lambda[i];		
		}
	}

	return NULL;
}

