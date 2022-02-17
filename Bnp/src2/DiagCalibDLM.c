
#include "srt_h_all.h"
#include "srt_h_lgmtypes.h"
#include "srt_h_lgmprotos.h"
#include "opfnctns.h"
#include "srt_h_allFx3F.h"
#include "swp_h_swap_pricing.h"
#include "swp_h_vol.h"
#include "diagCalibDLM.h"
#include "DiagCalibGen.h"
#include "math.h"

#define	MIN_CALIB_TIME	0.03

#define NUM_HERMITE		6
#define	CALPRES			1.0e-08
#define	NITER			10
#define	BUMP_LAM_LM		5.0E-04
#define	BUMP_LAM_NEWTON	1.0E-02

#define EPS_MAT_TO_LONG_DATE_CONVERSION  1.0e-8	

#define ONE_MONTH 0.083333333
double x[NUM_HERMITE+1], w[NUM_HERMITE+1];

//Interpolate linearly or piecewise constant
void constant_or_linear_interpolation_dlm(	double	*dates,
											double	*values,
											int		n_dates,
											double	date_wanted,
											double	*value_wanted)
{
	int i = 0;
	double result;
	
	//Flat intrapolation if before
	if(date_wanted <= dates[0]) result = values[0];
	else
	{	
		//Checks if the index is below the limit before comparing the dates
		//to avoid release compilation problems
		while ( i < n_dates && date_wanted > dates[i]) i++;
		i--;

		if (i < n_dates - 1)
		{
			//Interpolate linearly
			result = values[i] + (values[i+1] - values[i])/(dates[i+1] - dates[i]) * (date_wanted - dates[i]);
		}
		else
		{
			//Flat extrapolation
			result = values[n_dates - 1];
		}
	}

	*value_wanted = result;
}

static Err get_end_date(
				long	ex_date,
				long	struct_end_date,
				char	*tenor,
				int		theo_act,
				long	*end_date)
{
	long			start_date;
	Err				err;

	start_date = add_unit (ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	strupper (tenor);
	strip_white_space (tenor);
	if (*tenor == 'D')
	{
		*end_date = struct_end_date;
	}
	else
	{
		err = add_tenor (start_date, tenor, NO_BUSDAY_CONVENTION, end_date);
		if (err)
		{
			return err;
		}
	}

	if (theo_act)
	{
		*end_date = bus_date_method (*end_date, MODIFIED_SUCCEEDING);
	}

	return NULL;
}

static CPD_DIAG_CALIB_PARAM PARAML_LM,
							PARAMS_LM;

static int		NCPNL_LM,
				NEXL_LM,
				*EX_CPNL_LM,
				*EX_ENDCPNL_LM,

				NCPNS_LM,
				NEXS_LM,
				*EX_CPNS_LM,
				*EX_ENDCPNS_LM,
				
				ONE2F_LM,
				NB_MOM_LM,
				INDEX_MOM[MAX_INST],
				BREAK_INST,
				FREQ_SHORT_LM,
				SHIFT_FREQ_LM,				

				NLAM_LM,

				USE_JUMPS;

static double	*CPN_TIMEL_LM,
				*CPN_DFL_LM,
				*CPN_CVGL_LM,
				*EX_TIMEL_LM,
				*EX_STRIKEL_LM,
				*EX_PRICEL_LM,
				*EX_VEGAL_LM,
				CPN_GL_LM[MAX_CPN],
				CPN_G2L_LM[MAX_CPN], 
				EX_GL_LM[MAX_CPN],
				EX_G2L_LM[MAX_CPN],
				*EX_ZETAL_LM,
				
				*CPN_TIMES_LM,
				*CPN_DFS_LM,
				*CPN_CVGS_LM,
				*EX_TIMES_LM,
				*EX_STRIKES_LM,
				*EX_WEIGHTS_LM,
				*EX_PRICES_LM,
				CPN_GS_LM[MAX_CPN],
				CPN_G2S_LM[MAX_CPN], 
				EX_GS_LM[MAX_CPN],
				EX_G2S_LM[MAX_CPN],
				*EX_ZETAS_LM,				
				
				ZETA2_LM[MAX_CPN],
				ZETA12_LM[MAX_CPN],

				*LAM_TIME_LM,
				*LAM_LM_INIT,				
				LAM_LM[MAX_CPN],
				ALPHA_LM,
				GAMMA_LM,
				RHO_LM,
				
				PREC_LM,
				
				PRICE_LM[MAX_CPN],
				SENSI_LM[MAX_CPN][MAX_CPN];		

double		LIMIT_LAM_LM_DOWN;
double		LIMIT_LAM_LM_UP;
				
void init_static_lgmprcapgiventauts_dlm(

int				ncpnl,								/*	Total number of cash-flow dates */
double			cpn_timel[],							/*	Cash-Flow times */
double			cpn_dfl[],							/*	Df to cash-flow dates */
double			cpn_cvgl[],							/*	cvg from i-1 to i */
int				nexl,								/*	Total number of exercise dates */
double			ex_timel[],							/*	Exercise times */
int				ex_cpnl[],							/*	Index of the first cash-flow to be exercised */
int				ex_endcpnl[],
double			ex_strikel[],						/*	Strikes */
double			ex_pricel[],						/*	Market prices */
double			ex_vegal[],						/*	Market prices */
double			ex_zetal[],							/*	Output: zetas */

int				ncpns,								/*	Total number of cash-flow dates */
double			cpn_times[],							/*	Cash-Flow times */
double			cpn_dfs[],							/*	Df to cash-flow dates */
double			cpn_cvgs[],							/*	cvg from i-1 to i */
int				nexs,								/*	Total number of exercise dates */
double			ex_times[],							/*	Exercise times */
int				ex_cpns[],							/*	Index of the first cash-flow to be exercised */
int				ex_endcpns[],
double			ex_strikes[],						/*	Strikes */
double			ex_weights[],						/*	Weights on secondary instruments */
double			ex_prices[],						/*	Market prices */
double			ex_zetas[],							/*	Output: zetas */

int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
int				nlam,								/*	Lambda TS: may NOT be changed in the process */
double			lam_time[],
double			lam[],
double			alpha,
double			gamma,
double			rho,
CPD_DIAG_CALIB_PARAM paraml,
CPD_DIAG_CALIB_PARAM params,
DIAG_CALIB_LM_PARAMS ln_params)
{
int	i, nbm;

	NCPNL_LM = ncpnl;	
	NEXL_LM= nexl;
	EX_CPNL_LM = ex_cpnl;
	EX_ENDCPNL_LM = ex_endcpnl;
	CPN_TIMEL_LM= cpn_timel;
	CPN_DFL_LM = cpn_dfl;
	CPN_CVGL_LM = cpn_cvgl;
	EX_TIMEL_LM = ex_timel;
	EX_STRIKEL_LM = ex_strikel;
	EX_PRICEL_LM = ex_pricel;
	EX_VEGAL_LM = ex_vegal;
	EX_ZETAL_LM = ex_zetal;
	
	NCPNS_LM = ncpns;	
	NEXS_LM= nexs;
	EX_CPNS_LM = ex_cpns;
	EX_ENDCPNS_LM = ex_endcpns;
	CPN_TIMES_LM= cpn_times;
	CPN_DFS_LM = cpn_dfs;
	CPN_CVGS_LM = cpn_cvgs;
	EX_TIMES_LM = ex_times;
	EX_STRIKES_LM = ex_strikes;
	EX_WEIGHTS_LM = ex_weights;
	EX_PRICES_LM = ex_prices;
	EX_ZETAS_LM = ex_zetas;

	ONE2F_LM = one2F;	
	LAM_TIME_LM = lam_time;
	NLAM_LM = nlam;	
	LAM_LM_INIT = lam;

	ALPHA_LM = alpha;
	GAMMA_LM = gamma;
	RHO_LM = rho;

	PARAML_LM = paraml;
	PARAMS_LM = params;

	NB_MOM_LM = ln_params->nb_moment;
	FREQ_SHORT_LM = ln_params->freq_short;
	SHIFT_FREQ_LM = ln_params->shift_freq;

	BREAK_INST = ln_params->break_moment;

	if (BREAK_INST)
	{
		nbm = (int) (nexs / NB_MOM_LM + 0.5);

		INDEX_MOM[0] = 0;

		for (i=1; i<NB_MOM_LM; i++)
		{
			INDEX_MOM[i] = INDEX_MOM[i-1] + nbm;
		}

		INDEX_MOM[NB_MOM_LM] = nexs - 1;
	}	
}

/* This method try to match the momentum */
Err	static_lgmprcapgiventauts_momentum_dlm(double	index,
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
		/* We calculate everything */
		if (ONE2F_LM == 1)
		{			
			/* Calculate Sensi */
			for (j=0; j<nlam; j++)
			{
				lambda[j] += BUMP_LAM_LM;
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

				err = lgmcalibzeta1F_tauts2(
									NCPNL_LM,
									CPN_TIMEL_LM,
									CPN_DFL_LM,
									CPN_CVGL_LM,
									CPN_GL_LM,
									NEXL_LM,
									EX_TIMEL_LM,
									EX_CPNL_LM,
									EX_ENDCPNL_LM,
									EX_GL_LM,
									EX_STRIKEL_LM,
									EX_PRICEL_LM,
									EX_VEGAL_LM,
									EX_ZETAL_LM,
									nlam,
									LAM_TIME_LM,
									lambda,
									0,
									PARAML_LM->precision,
									PARAML_LM->vega_prec,
									PARAML_LM->min_fact,
									PARAML_LM->max_fact,
									PARAML_LM->use_jumps);

				if (err)
				{
					return err;
				}

				static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

				for (i=0; i<NEXS_LM; i++)
				{
					if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
					{
						SENSI_LM[i][j] = (lgmsopval1F(
											EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
											&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
											&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
											&(CPN_GS_LM[EX_CPNS_LM[i]]), 
											EX_ZETAS_LM[i], 
											EX_GS_LM[i],
											EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
					}
				}

				lambda[j] -= BUMP_LAM_LM;
			}

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

			err = lgmcalibzeta1F_tauts2(
									NCPNL_LM,
									CPN_TIMEL_LM,
									CPN_DFL_LM,
									CPN_CVGL_LM,
									CPN_GL_LM,
									NEXL_LM,
									EX_TIMEL_LM,
									EX_CPNL_LM,
									EX_ENDCPNL_LM,
									EX_GL_LM,
									EX_STRIKEL_LM,
									EX_PRICEL_LM,
									EX_VEGAL_LM,
									EX_ZETAL_LM,
									nlam,
									LAM_TIME_LM,
									lambda,
									0,
									PARAML_LM->precision,
									PARAML_LM->vega_prec,
									PARAML_LM->min_fact,
									PARAML_LM->max_fact,
									PARAML_LM->use_jumps);

			if (err)
			{
				return err;
			}

			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);
			
			/* Calculate Prices */
			for (i=0; i<NEXS_LM; i++)
			{
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					PRICE_LM[i] = (lgmsopval1F(
											EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
											&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
											&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
											&(CPN_GS_LM[EX_CPNS_LM[i]]), 
											EX_ZETAS_LM[i], 
											EX_GS_LM[i],
											EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
				}
			}
		}
		else
		{
			/* Calculate Sensi */
			for (j=0; j<nlam; j++)
			{
				lambda[j] += BUMP_LAM_LM;	
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
				static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
				static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

				err = lgmcalibzeta2F_tauts2(
										NCPNL_LM,
										CPN_TIMEL_LM,
										CPN_DFL_LM,
										CPN_CVGL_LM,
										CPN_GL_LM,
										CPN_G2L_LM,
										NEXL_LM,
										EX_TIMEL_LM,
										EX_CPNL_LM,
										EX_ENDCPNL_LM,
										EX_GL_LM,
										EX_G2L_LM,
										EX_STRIKEL_LM,
										EX_PRICEL_LM,
										EX_VEGAL_LM,
										EX_ZETAL_LM,
										nlam,
										LAM_TIME_LM,
										lambda,
										ALPHA_LM,
										GAMMA_LM,
										RHO_LM,
										0,
										PARAML_LM->precision,
										PARAML_LM->vega_prec,
										PARAML_LM->min_fact,
										PARAML_LM->max_fact,
										PARAML_LM->use_jumps);

				if (err)
				{
					return err;
				}

				static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

				static_lgmcalczeta2zeta12_tauts(
					NEXS_LM,
					EX_TIMES_LM,
					EX_ZETAS_LM,
					nlam,
					LAM_TIME_LM,
					lambda,
					ALPHA_LM,
					GAMMA_LM,
					RHO_LM,
					ZETA2_LM,
					ZETA12_LM);

				for (i=0; i<NEXS_LM; i++)
				{	
					if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
					{
						SENSI_LM[i][j] = (lgmsopval2F (
												EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
												&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
												&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
												&(CPN_GS_LM[EX_CPNS_LM[i]]), 
												&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
												EX_ZETAS_LM[i], 
												ZETA2_LM[i], 
												ZETA12_LM[i], 
												EX_GS_LM[i],
												EX_G2S_LM[i],
												EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
					}
				}

				lambda[j] -= BUMP_LAM_LM;
			}

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
			static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

			err = lgmcalibzeta2F_tauts2(
										NCPNL_LM,
										CPN_TIMEL_LM,
										CPN_DFL_LM,
										CPN_CVGL_LM,
										CPN_GL_LM,
										CPN_G2L_LM,
										NEXL_LM,
										EX_TIMEL_LM,
										EX_CPNL_LM,
										EX_ENDCPNL_LM,
										EX_GL_LM,
										EX_G2L_LM,
										EX_STRIKEL_LM,
										EX_PRICEL_LM,
										EX_VEGAL_LM,
										EX_ZETAL_LM,
										nlam,
										LAM_TIME_LM,
										lambda,
										ALPHA_LM,
										GAMMA_LM,
										RHO_LM,
										0,
										PARAML_LM->precision,
										PARAML_LM->vega_prec,
										PARAML_LM->min_fact,
										PARAML_LM->max_fact,
										PARAML_LM->use_jumps);

			if (err)
			{
				return err;
			}

			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

			static_lgmcalczeta2zeta12_tauts(
				NEXS_LM,
				EX_TIMES_LM,
				EX_ZETAS_LM,
				nlam,
				LAM_TIME_LM,
				lambda,
				ALPHA_LM,
				GAMMA_LM,
				RHO_LM,
				ZETA2_LM,
				ZETA12_LM);

			for (i=0; i<NEXS_LM; i++)
			{	
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					PRICE_LM[i] = (lgmsopval2F (
										EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
										&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
										&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
										&(CPN_GS_LM[EX_CPNS_LM[i]]), 
										&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
										EX_ZETAS_LM[i], 
										ZETA2_LM[i], 
										ZETA12_LM[i], 
										EX_GS_LM[i],
										EX_G2S_LM[i],
										EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
				}
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

Err	static_lgmprcapgiventauts_dlm(double	index,
							  double	lam[],
							  double	*price,
							  double	*gradient,
							  int		nlam)
{
	static int	i, j, ind;
	Err		err;
	double	*lambda = &(lam[1]);	

	ind = (int) (index);

	if (ind - SHIFT_FREQ_LM == 0)
	{
		/* We calculate everything */
		if (ONE2F_LM == 1)
		{			
			/* Calculate Sensi */
			for (j=0; j<nlam; j++)
			{
				lambda[j] += BUMP_LAM_LM;
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

				err = lgmcalibzeta1F_tauts2(
									NCPNL_LM,
									CPN_TIMEL_LM,
									CPN_DFL_LM,
									CPN_CVGL_LM,
									CPN_GL_LM,
									NEXL_LM,
									EX_TIMEL_LM,
									EX_CPNL_LM,
									EX_ENDCPNL_LM,
									EX_GL_LM,
									EX_STRIKEL_LM,
									EX_PRICEL_LM,
									EX_VEGAL_LM,
									EX_ZETAL_LM,
									nlam,
									LAM_TIME_LM,
									lambda,
									0,
									PARAML_LM->precision,
									PARAML_LM->vega_prec,
									PARAML_LM->min_fact,
									PARAML_LM->max_fact,
									PARAML_LM->use_jumps);

				if (err)
				{
					return err;
				}

				static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

				for (i=0; i<NEXS_LM; i++)
				{
					if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
					{
						SENSI_LM[i][j] = lgmsopval1F(
											EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
											&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
											&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
											&(CPN_GS_LM[EX_CPNS_LM[i]]), 
											EX_ZETAS_LM[i], 
											EX_GS_LM[i],
											EX_STRIKES_LM[i]);
					}
				}

				lambda[j] -= BUMP_LAM_LM;
			}

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

			err = lgmcalibzeta1F_tauts2(
									NCPNL_LM,
									CPN_TIMEL_LM,
									CPN_DFL_LM,
									CPN_CVGL_LM,
									CPN_GL_LM,
									NEXL_LM,
									EX_TIMEL_LM,
									EX_CPNL_LM,
									EX_ENDCPNL_LM,
									EX_GL_LM,
									EX_STRIKEL_LM,
									EX_PRICEL_LM,
									EX_VEGAL_LM,
									EX_ZETAL_LM,
									nlam,
									LAM_TIME_LM,
									lambda,
									0,
									PARAML_LM->precision,
									PARAML_LM->vega_prec,
									PARAML_LM->min_fact,
									PARAML_LM->max_fact,
									PARAML_LM->use_jumps);
			if (err)
			{
				return err;
			}

			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);
			
			/* Calculate Prices */
			for (i=0; i<NEXS_LM; i++)
			{
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					PRICE_LM[i] = lgmsopval1F(
											EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
											&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
											&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
											&(CPN_GS_LM[EX_CPNS_LM[i]]), 
											EX_ZETAS_LM[i], 
											EX_GS_LM[i],
											EX_STRIKES_LM[i]);

					for (j=0; j<nlam; j++)
					{
						SENSI_LM[i][j] = (SENSI_LM[i][j] - PRICE_LM[i]) / BUMP_LAM_LM;
					}
				}
			}
		}
		else
		{
			/* Calculate Sensi */
			for (j=0; j<nlam; j++)
			{
				lambda[j] += BUMP_LAM_LM;	
				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
				static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

				static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
				static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

				err = lgmcalibzeta2F_tauts2(
										NCPNL_LM,
										CPN_TIMEL_LM,
										CPN_DFL_LM,
										CPN_CVGL_LM,
										CPN_GL_LM,
										CPN_G2L_LM,
										NEXL_LM,
										EX_TIMEL_LM,
										EX_CPNL_LM,
										EX_ENDCPNL_LM,
										EX_GL_LM,
										EX_G2L_LM,
										EX_STRIKEL_LM,
										EX_PRICEL_LM,
										EX_VEGAL_LM,
										EX_ZETAL_LM,
										nlam,
										LAM_TIME_LM,
										lambda,
										ALPHA_LM,
										GAMMA_LM,
										RHO_LM,
										0,
										PARAML_LM->precision,
										PARAML_LM->vega_prec,
										PARAML_LM->min_fact,
										PARAML_LM->max_fact,
										PARAML_LM->use_jumps);

				if (err)
				{
					return err;
				}

				static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

				static_lgmcalczeta2zeta12_tauts(
					NEXS_LM,
					EX_TIMES_LM,
					EX_ZETAS_LM,
					nlam,
					LAM_TIME_LM,
					lambda,
					ALPHA_LM,
					GAMMA_LM,
					RHO_LM,
					ZETA2_LM,
					ZETA12_LM);

				for (i=0; i<NEXS_LM; i++)
				{	
					if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
					{
						SENSI_LM[i][j] = lgmsopval2F (
												EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
												&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
												&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
												&(CPN_GS_LM[EX_CPNS_LM[i]]), 
												&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
												EX_ZETAS_LM[i], 
												ZETA2_LM[i], 
												ZETA12_LM[i], 
												EX_GS_LM[i],
												EX_G2S_LM[i],
												EX_STRIKES_LM[i]);
					}
				}

				lambda[j] -= BUMP_LAM_LM;
			}

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

			static_lgmsetupG_tauts (nlam, LAM_TIME_LM, lambda, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
			static_lgmsetupG2_tauts (nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

			err = lgmcalibzeta2F_tauts2(
										NCPNL_LM,
										CPN_TIMEL_LM,
										CPN_DFL_LM,
										CPN_CVGL_LM,
										CPN_GL_LM,
										CPN_G2L_LM,
										NEXL_LM,
										EX_TIMEL_LM,
										EX_CPNL_LM,
										EX_ENDCPNL_LM,
										EX_GL_LM,
										EX_G2L_LM,
										EX_STRIKEL_LM,
										EX_PRICEL_LM,
										EX_VEGAL_LM,
										EX_ZETAL_LM,
										nlam,
										LAM_TIME_LM,
										lambda,
										ALPHA_LM,
										GAMMA_LM,
										RHO_LM,
										0,
										PARAML_LM->precision,
										PARAML_LM->vega_prec,
										PARAML_LM->min_fact,
										PARAML_LM->max_fact,
										PARAML_LM->use_jumps);

			if (err)
			{
				return err;
			}

			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, nlam, LAM_TIME_LM, lambda, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

			static_lgmcalczeta2zeta12_tauts(
				NEXS_LM,
				EX_TIMES_LM,
				EX_ZETAS_LM,
				nlam,
				LAM_TIME_LM,
				lambda,
				ALPHA_LM,
				GAMMA_LM,
				RHO_LM,
				ZETA2_LM,
				ZETA12_LM);

			for (i=0; i<NEXS_LM; i++)
			{	
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					PRICE_LM[i] = lgmsopval2F (
										EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
										&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
										&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
										&(CPN_GS_LM[EX_CPNS_LM[i]]), 
										&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
										EX_ZETAS_LM[i], 
										ZETA2_LM[i], 
										ZETA12_LM[i], 
										EX_GS_LM[i],
										EX_G2S_LM[i],
										EX_STRIKES_LM[i]);

					for (j=0; j<nlam; j++)
					{
						SENSI_LM[i][j] = (SENSI_LM[i][j] - PRICE_LM[i]) / BUMP_LAM_LM;
					}
				}
			}		
		}
	}

	/* Return the asked result */	

	(*price) = PRICE_LM[ind];
	for (j=0; j<nlam; j++)
	{
		gradient[j+1] = SENSI_LM[ind][j];
	}

	return NULL;
}

double solve_for_next_coef_dlm(	double	**res_iter,
								int		nb_iter,
								double	premium_tgt,
								int		method)			/* 0: linear 1: quadratic */
{
static double	vega, coef, a, b, c, delta, res1, res2;
static int		i, type;

	if (nb_iter == 2 || method == 0)
	{
		/* linear interpolation */

		i = 0;

		while (i < nb_iter && res_iter[i][1] < premium_tgt)
		{
			i++;
		}

		if (i > 0)
		{
			if (i == nb_iter)
			{
				i = i - 2;
			}
			else
			{
				i = i - 1;
			}
		}

		/* interpolation between i and i+1 */

		vega = (res_iter[i+1][1] - res_iter[i][1]) / (res_iter[i+1][0] - res_iter[i][0]);
		coef = res_iter[i][0] + (premium_tgt - res_iter[i][1]) / vega;
		
	}
	else if (nb_iter > 2)
	{
		/* quadratic interpolation */

		if (nb_iter == 3)
		{
			i = 1;
		}
		else
		{
			i = 0;

			while (i < nb_iter && res_iter[i][1] < premium_tgt)
			{
				i++;
			}

			if (i <= 1)
			{
				i = 1;
			}
			else if (i >= nb_iter - 1)
			{
				i = nb_iter - 2;
			}
			else
			{
				if ((premium_tgt - res_iter[i][1]) < 0.5 * (res_iter[i+1][1] - res_iter[i][1]))
				{
					i--;
				}
			}
		}

		/* interpolation between i-1, i and i+1 */

		/* check that numbers are different */
		if (fabs(res_iter[i][0] - res_iter[i-1][0]) < 1.0E-12 || fabs(res_iter[i][0] - res_iter[i+1][0]) < 1.0E-12)
		{
			coef = res_iter[i][0];
			return coef;
		}

		/* check that the function is really monotonic */
		if (fabs(res_iter[i][1] - res_iter[i-1][1]) < 1.0E-08)
		{
			if ((premium_tgt - res_iter[i][1]) * (premium_tgt - res_iter[i-1][1]) < 0.0)
			{
				coef = 0.5 * (res_iter[i][0] + res_iter[i-1][0]);
			}
			else
			{
				coef = 0.5 * (res_iter[i][0] + res_iter[i+1][0]);
			}

			return coef;
		}

		if (fabs(res_iter[i][1] - res_iter[i+1][1]) < 1.0E-08)
		{
			if ((premium_tgt - res_iter[i][1]) * (premium_tgt - res_iter[i+1][1]) < 0.0)
			{
				coef = 0.5 * (res_iter[i][0] + res_iter[i+1][0]);
			}
			else
			{
				coef = 0.5 * (res_iter[i][0] + res_iter[i-1][0]);
			}

			return coef;
		}

		/* solve quadratic equation */

		a =		res_iter[i-1][1] / (res_iter[i-1][0] - res_iter[i][0]) / (res_iter[i-1][0] - res_iter[i+1][0])
			+	res_iter[i][1] / (res_iter[i][0] - res_iter[i-1][0]) / (res_iter[i][0] - res_iter[i+1][0])
			+	res_iter[i+1][1] / (res_iter[i+1][0] - res_iter[i-1][0]) / (res_iter[i+1][0] - res_iter[i][0]);

		if (fabs(a) < 1.0E-12)
		{
			/* we are linear */
			coef = solve_for_next_coef_dlm(res_iter, nb_iter, premium_tgt, 0);
			return coef;
		}

		b =		 res_iter[i-1][1] / (res_iter[i-1][0] - res_iter[i][0]) / (res_iter[i-1][0] - res_iter[i+1][0])
			*	(res_iter[i][0] + res_iter[i+1][0]) + res_iter[i][1] / (res_iter[i][0] - res_iter[i-1][0]) / (res_iter[i][0] - res_iter[i+1][0])
			*	(res_iter[i-1][0] + res_iter[i+1][0]) + res_iter[i+1][1] / (res_iter[i+1][0] - res_iter[i-1][0]) / (res_iter[i+1][0] - res_iter[i][0]) * (res_iter[i][0] + res_iter[i-1][0]);
		b = -b;

		c = res_iter[i-1][1] / (res_iter[i-1][0] - res_iter[i][0]) / (res_iter[i-1][0] - res_iter[i+1][0]) 
			* res_iter[i][0] * res_iter[i+1][0] + res_iter[i][1] / (res_iter[i][0] - res_iter[i-1][0]) / (res_iter[i][0] - res_iter[i+1][0])
			* res_iter[i-1][0] * res_iter[i+1][0] + res_iter[i+1][1] / (res_iter[i+1][0] - res_iter[i-1][0]) / (res_iter[i+1][0] - res_iter[i][0]) * res_iter[i-1][0] * res_iter[i][0] 
			- premium_tgt;

		delta = b * b - 4 * a * c;

		if (delta > 0)
		{
			delta = sqrt(delta);
			res1 = (-b + delta) / (2 * a);
			res2 = (-b - delta) / (2 * a);

			/* choose the one */
			if (fabs(res2 - res_iter[i][0]) < fabs(res1 - res_iter[i][0]))
			{
				coef = res2;
			}
			else
			{
				coef = res1;
			}
		}
		else
		{
			coef = solve_for_next_coef_dlm(	res_iter,
											nb_iter,
											premium_tgt,
											0);
		}
	}
	else
	{
		/* only one point is available */
		if (res_iter[0][1] > premium_tgt)
		{
			coef = res_iter[0][0] * 0.9;
		}
		else
		{
			coef = res_iter[0][0] * 1.1;
		}
	}

	return coef;
}

static Err	lgmprcapgiventauts_simpleNewton_dlm(	double	*lam,
													double	*price,
													double	*lam_sens)
{
	int		i, j, k, l;	
	double	price1, price2;
	double	lam1, lam2;
	double	**res_iter = NULL;
	int		already_tried;
	Err		err;
	
	res_iter = dmatrix(0, PARAMS_LM->nb_iter_max, 0, 1);	

	if (!res_iter)
	{
		err = "Memory allocation faillure in lgmprcapgiventauts_momentum_dlm";
		goto FREE_RETURN;
	}

	*lam_sens = 0.0;

	lam1 = 0.0;

	/* We calculate everything */
	if (ONE2F_LM == 1)
	{
		/* First guess */		
		for (i=0; i<NLAM_LM; i++)
		{
			LAM_LM[i] = LAM_LM_INIT[i] + lam1;
		}		

		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

		err = lgmcalibzeta1F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM,
									NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_STRIKEL_LM,
									EX_PRICEL_LM, EX_VEGAL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, 0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

		if (err) goto FREE_RETURN;
		
		static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

		price1 = 0.0;

		for (i=0; i<NEXS_LM; i++)
		{
			if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
			{
				price1 += EX_WEIGHTS_LM[i] * (lgmsopval1F(
									EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i],
									&(CPN_DFS_LM[EX_CPNS_LM[i]]),
									&(CPN_CVGS_LM[EX_CPNS_LM[i]]),
									&(CPN_GS_LM[EX_CPNS_LM[i]]),
									EX_ZETAS_LM[i],
									EX_GS_LM[i],
									EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
			}
		}

		if (fabs(price1) < PARAMS_LM->precision)
		{			
			for (i=0; i<NLAM_LM; i++)
			{
				lam[i] = LAM_LM_INIT[i] + lam1;
			}

			*price = price1;
			goto FREE_RETURN;
		}

		res_iter[0][0] = lam1;
		res_iter[0][1] = price1;

		/* second guess */
		lam2 = lam1 + BUMP_LAM_NEWTON;		
				
		for (i=0; i<NLAM_LM; i++)
		{
			LAM_LM[i] = LAM_LM_INIT[i] + lam2;
		}		

		k = 1;

		while (fabs(price1) > PARAMS_LM->precision && k < PARAMS_LM->nb_iter_max)
		{
			k++;

			static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

			err = lgmcalibzeta1F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM,
										NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_STRIKEL_LM,
										EX_PRICEL_LM, EX_VEGAL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

			if (err) goto FREE_RETURN;
			
			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

			price2 = 0.0;

			for (i=0; i<NEXS_LM; i++)
			{
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					price2 += EX_WEIGHTS_LM[i] * (lgmsopval1F(
										EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i],
										&(CPN_DFS_LM[EX_CPNS_LM[i]]),
										&(CPN_CVGS_LM[EX_CPNS_LM[i]]),
										&(CPN_GS_LM[EX_CPNS_LM[i]]),
										EX_ZETAS_LM[i],
										EX_GS_LM[i],
										EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
				}
			}

			/* Save Res */
			l = 0;
			while (l < k-1 && res_iter[l][0] < lam2)
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

			lam2 = solve_for_next_coef_gen(	res_iter,
											k,
											0.0,
											PARAMS_LM->precision,
											2);
			
			for (i=0; i<NLAM_LM; i++)
			{
				LAM_LM[i] = LAM_LM_INIT[i] + lam2;
			}			

			if (lam2 < LIMIT_LAM_LM_DOWN || lam2 > LIMIT_LAM_LM_UP)
			{
				if (lam2 < LIMIT_LAM_LM_DOWN)
				{
					lam2 = LIMIT_LAM_LM_DOWN;															
				}
				else
				{
					lam2 = LIMIT_LAM_LM_UP;															
				}

				already_tried = 0;

				for (i=0; i<k; i++)
				{
					if (fabs(lam2 - res_iter[i][0]) < 1.0E-08)
					{
						already_tried = 1;
						*price = res_iter[i][1];
						break;
					}
				}

				for (i=0; i<NLAM_LM; i++)
				{
					LAM_LM[i] = LAM_LM_INIT[i] + lam2;
				}

				if (already_tried)
				{					
					/* recalibrate */
					static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
					static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

					err = lgmcalibzeta1F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM,
												NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_STRIKEL_LM,
												EX_PRICEL_LM, EX_VEGAL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

					/* Return the asked result */
					for (i=0; i<NLAM_LM; i++)
					{
						lam[i] = LAM_LM[i];
					}

					goto FREE_RETURN;
				}
			}				
		}

		/* Last recalibration */
		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);

		err = lgmcalibzeta1F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM,
									NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_STRIKEL_LM,
									EX_PRICEL_LM, EX_VEGAL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

		if (err) goto FREE_RETURN;
	}
	else
	{			
		/* First guess */
		for (i=0; i<NLAM_LM; i++)
		{
			LAM_LM[i] = LAM_LM_INIT[i] + lam1;
		}

		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
		static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
		static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

		err = lgmcalibzeta2F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM, CPN_G2L_LM,
									NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_G2L_LM, EX_STRIKEL_LM, EX_PRICEL_LM, EX_VEGAL_LM, 
									EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

		if (err) goto FREE_RETURN;
		
		static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

		static_lgmcalczeta2zeta12_tauts(NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM, ZETA2_LM, ZETA12_LM);

		price1 = 0.0;

		for (i=0; i<NEXS_LM; i++)
		{	
			if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
			{
				price1 += EX_WEIGHTS_LM[i] * (lgmsopval2F (
									EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
									&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
									&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
									&(CPN_GS_LM[EX_CPNS_LM[i]]), 
									&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
									EX_ZETAS_LM[i], 
									ZETA2_LM[i], 
									ZETA12_LM[i], 
									EX_GS_LM[i],
									EX_G2S_LM[i],
									EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
			}
		}

		if (fabs(price1) < PARAMS_LM->precision)
		{
			for (i=0; i<NLAM_LM; i++)
			{
				lam[i] = LAM_LM_INIT[i] + lam1;
			}
			
			*price = price1;
			goto FREE_RETURN;
		}

		res_iter[0][0] = lam1;
		res_iter[0][1] = price1;

		/* second guess */
		lam2 = lam1 + BUMP_LAM_NEWTON;	
		
		for (i=0; i<NLAM_LM; i++)
		{
			LAM_LM[i] = LAM_LM_INIT[i] + lam2;
		}

		k = 1;

		while (fabs(price1) > PARAMS_LM->precision && k < PARAMS_LM->nb_iter_max)
		{
			k++;

			static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
			static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

			static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
			static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

			err = lgmcalibzeta2F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM, CPN_G2L_LM,
										NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_G2L_LM, EX_STRIKEL_LM, EX_PRICEL_LM, EX_VEGAL_LM, 
										EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

			if (err) goto FREE_RETURN;
			
			static_interpolate_zeta(NEXL_LM, EX_TIMEL_LM, EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM);

			static_lgmcalczeta2zeta12_tauts(NEXS_LM, EX_TIMES_LM, EX_ZETAS_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM, ZETA2_LM, ZETA12_LM);

			price2 = 0.0;

			for (i=0; i<NEXS_LM; i++)
			{	
				if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM))
				{
					price2 += EX_WEIGHTS_LM[i] * (lgmsopval2F (
										EX_ENDCPNS_LM[i] + 1 - EX_CPNS_LM[i], 
										&(CPN_DFS_LM[EX_CPNS_LM[i]]), 
										&(CPN_CVGS_LM[EX_CPNS_LM[i]]), 
										&(CPN_GS_LM[EX_CPNS_LM[i]]), 
										&(CPN_G2S_LM[EX_CPNS_LM[i]]), 
										EX_ZETAS_LM[i], 
										ZETA2_LM[i], 
										ZETA12_LM[i], 
										EX_GS_LM[i],
										EX_G2S_LM[i],
										EX_STRIKES_LM[i]) - EX_PRICES_LM[i]);
				}
			}

			/* Save Res */
			l = 0;
			while (l < k-1 && res_iter[l][0] < lam2)
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

			lam2 = solve_for_next_coef_gen(	res_iter,
											k,
											0.0,
											PARAMS_LM->precision,
											2);

			for (i=0; i<NLAM_LM; i++)
			{
				LAM_LM[i] = LAM_LM_INIT[i] + lam2;
			}			

			if (lam2 < LIMIT_LAM_LM_DOWN || lam2 > LIMIT_LAM_LM_UP)
			{
				if (lam2 < LIMIT_LAM_LM_DOWN)
				{
					lam2 = LIMIT_LAM_LM_DOWN;															
				}
				else
				{
					lam2 = LIMIT_LAM_LM_UP;															
				}

				already_tried = 0;

				for (i=0; i<k; i++)
				{
					if (fabs(lam2 - res_iter[i][0]) < 1.0E-08)
					{
						already_tried = 1;
						*price = res_iter[i][1];
						break;
					}
				}

				for (i=0; i<NLAM_LM; i++)
				{
					LAM_LM[i] = LAM_LM_INIT[i] + lam2;
				}

				if (already_tried)
				{
					/* recalibrate */
					static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
					static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

					static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
					static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

					err = lgmcalibzeta2F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM, CPN_G2L_LM,
												NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_G2L_LM, EX_STRIKEL_LM, EX_PRICEL_LM, EX_VEGAL_LM, 
												EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

					/* Return the asked result */
					for (i=0; i<NLAM_LM; i++)
					{
						lam[i] = LAM_LM[i];
					}

					goto FREE_RETURN;
				}
			}
		}

		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_GL_LM, NEXL_LM, EX_TIMEL_LM, EX_GL_LM);
		static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNL_LM, CPN_TIMEL_LM, CPN_G2L_LM, NEXL_LM, EX_TIMEL_LM, EX_G2L_LM);

		static_lgmsetupG_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, NCPNS_LM, CPN_TIMES_LM, CPN_GS_LM, NEXS_LM, EX_TIMES_LM, EX_GS_LM);
		static_lgmsetupG2_tauts (NLAM_LM, LAM_TIME_LM, LAM_LM, GAMMA_LM, NCPNS_LM, CPN_TIMES_LM, CPN_G2S_LM, NEXS_LM, EX_TIMES_LM, EX_G2S_LM);

		err = lgmcalibzeta2F_tauts2(NCPNL_LM, CPN_TIMEL_LM, CPN_DFL_LM, CPN_CVGL_LM, CPN_GL_LM, CPN_G2L_LM,
									NEXL_LM, EX_TIMEL_LM, EX_CPNL_LM, EX_ENDCPNL_LM, EX_GL_LM, EX_G2L_LM, EX_STRIKEL_LM, EX_PRICEL_LM, EX_VEGAL_LM, 
									EX_ZETAL_LM, NLAM_LM, LAM_TIME_LM, LAM_LM, ALPHA_LM, GAMMA_LM, RHO_LM,  0, PARAML_LM->precision, PARAML_LM->vega_prec, PARAML_LM->min_fact, PARAML_LM->max_fact, PARAML_LM->use_jumps);

		if (err) goto FREE_RETURN;
	}

	/* Return the asked result */
	for (i=0; i<NLAM_LM; i++)
	{
		lam[i] = LAM_LM_INIT[i] + lam2;
	}

	*price = price2;

	solve_for_next_sensi(	res_iter,
							k,
							lam2,
							0.0005,
							lam_sens);

FREE_RETURN:

	if (res_iter) free_dmatrix(res_iter, 0, PARAMS_LM->nb_iter_max, 0, 1);

	return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err lgmcalibzetalambda_tauts_dlm(
int				ncpnl,								/*	Total number of cash-flow dates */
double			cpn_timel[],							/*	Cash-Flow times */
double			cpn_dfl[],							/*	Df to cash-flow dates */
double			cpn_cvgl[],							/*	cvg from i-1 to i */
int				nexl,								/*	Total number of exercise dates */
double			ex_timel[],							/*	Exercise times */
int				ex_cpnl[],							/*	Index of the first cash-flow to be exercised */
int				ex_ncpnl[],							/*	Number of coupons in each caplet */
double			ex_strikel[],						/*	Strikes for diagonal */
double			ex_pricel[],						/*	Market prices for diagonal */
double			ex_vegal[],							/*	Market vegas for diagonal */

int				ncpns,								/*	Total number of cash-flow dates */
double			cpn_times[],							/*	Cash-Flow times */
double			cpn_dfs[],							/*	Df to cash-flow dates */
double			cpn_cvgs[],							/*	cvg from i-1 to i */
int				nexs,								/*	Total number of exercise dates */
double			ex_times[],							/*	Exercise times */
int				ex_cpns[],							/*	Index of the first cash-flow to be exercised */
int				ex_ncpns[],							/*	Number of coupons in each caplet */
double			ex_strikes[],						/*	Strikes for diagonal */
double			ex_weights[],						/*	Weights for short */
double			ex_prices[],						/*	Market prices for diagonal */
double			ex_vegas[],							/*	Market vega for short */

double			ex_zeta[],							/*	Output: zetas */
int				fix_lambda,							/*	0: calib lambda to cap, 1: fix lambda calib
															to diagonal */
int				nlam,								/*	Lambda TS: may NOT be changed in the process */
double			lam_time[],
double			lam[],
int				one2F,								/*	Number of factors */
/*	Alpha, Gamma, Rho (2F only) */
double			alpha,
double			gamma,
double			rho,
CPD_DIAG_CALIB_PARAM paraml,
CPD_DIAG_CALIB_PARAM params,
DIAG_CALIB_LM_PARAMS lm_params,
double			*lam_sens)							
{
	int	i;
	double	ex_spriceres[MAX_CPN];
	Err err = NULL;
	double	data[MAX_CPN],
			weight[MAX_CPN],
			target[MAX_CPN],
			lambda[MAX_CPN],
			ex_zetas[MAX_INST];

	double	fitting_error;
	double	sum_weight;

	if (!fix_lambda && nexs < 1)
	{
		return serror ("Cannot calibrate lambda and zeta with less than 2 exercises - choose fix lambda");
	}

	if (fix_lambda)
	{
		err = lgmprcapgivenlambda_tauts2(
					ncpnl, cpn_timel, cpn_dfl, cpn_cvgl,
					nexl, ex_timel, ex_cpnl, ex_ncpnl, ex_strikel, ex_pricel, ex_vegal, ex_ncpns, ex_strikes, ex_zeta,
					one2F, nlam, lam_time, lam, alpha, gamma, rho, 0, paraml->precision, paraml->vega_prec, paraml->min_fact, paraml->max_fact, paraml->use_jumps, 0, ex_spriceres);
		
		return err;
	}
	else
	{
		for (i=1; i<=nlam; i++)
		{			
			lambda[i] = lam[i-1];
		}

		lm_params->shift_freq = (int) (fmod(lm_params->shift_freq, nexs) + 0.5);

		if (lm_params->use_moment)
		{
			/* Momentum method */
			lambda[0] = lm_params->nb_moment;
			sum_weight = 0.0;

			for (i=0; i<nexs; i++)
			{				
				if (ex_weights)
				{
					weight[i] = ex_weights[i];
				}
				else
				{
					if (lm_params->vega_weight)
					{
						weight[i] = ex_vegas[i];
					}
					else
					{
						weight[i] = 1.0;
					}
				}
				
				sum_weight += weight[i];
				target[i] = 0.0;
			}
			
			sum_weight /= (nexs * 1.0);

			/* Rescaling */
			/*
			for (i=0; i<nexs; i++)
			{
				weight[i] /= sum_weight;
			}
			*/

			init_static_lgmprcapgiventauts_dlm(ncpnl, cpn_timel, cpn_dfl, cpn_cvgl, nexl, ex_timel, ex_cpnl, ex_ncpnl, ex_strikel, ex_pricel, ex_vegal, ex_zeta,
				ncpns, cpn_times, cpn_dfs, cpn_cvgs, nexs, ex_times, ex_cpns, ex_ncpns, ex_strikes, weight, ex_prices, ex_zetas,
				one2F, nlam, lam_time, lam, alpha, gamma, rho, paraml, params, lm_params);

			/* find the limit lam */
			LIMIT_LAM_LM_DOWN = lam[0];
			LIMIT_LAM_LM_UP = lam[0];

			for (i=1; i<nlam; i++)
			{
				if (lam[i] < LIMIT_LAM_LM_DOWN)
				{
					LIMIT_LAM_LM_DOWN = lam[i];
				}

				if (lam[i] > LIMIT_LAM_LM_UP)
				{
					LIMIT_LAM_LM_UP = lam[i];
				}
			}

			if (nlam > 1)
			{				
				LIMIT_LAM_LM_DOWN = -0.60 - LIMIT_LAM_LM_DOWN;
				LIMIT_LAM_LM_UP = 0.60 - LIMIT_LAM_LM_UP;
			}
			else
			{
				LIMIT_LAM_LM_DOWN = params->lambda_min - LIMIT_LAM_LM_DOWN;
				LIMIT_LAM_LM_UP = params->lambda_max - LIMIT_LAM_LM_UP;
			}

			err = lgmprcapgiventauts_simpleNewton_dlm(lam, &fitting_error, lam_sens);

			if (err)
			{
				return err;
			}		

			smessage ("Calibration of lambda, error in bp: %.00f", 10000 * fitting_error);

			/*
			if (nlam == 1)
			{				
				err = lgmprcapgiventauts_simpleNewton_dlm(lambda, &fitting_error);

				if (err)
				{
					return err;
				}

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
									static_lgmprcapgiventauts_momentum_dlm,
									&fitting_error);

				if (err)
				{
					return err;
				}
			}
			*/
		}
		else
		{
			for (i=1; i<=nexs; i++)
			{					
				data[i] = (int) (lm_params->shift_freq + (i - 1) * lm_params->freq_short + 0.5);

				if (ex_weights)
				{
					weight[i] = ex_weights[(int)data[i]];
				}
				else
				{
					if  (lm_params->vega_weight)
					{
						weight[i] = ex_vegas[(int)data[i]];
					}
					else
					{
						weight[i] = 1.0;
					}
				}

				target[i] = ex_prices[(int)data[i]];
			}

			init_static_lgmprcapgiventauts_dlm(ncpnl, cpn_timel, cpn_dfl, cpn_cvgl, nexl, ex_timel, ex_cpnl, ex_ncpnl, ex_strikel, ex_pricel, ex_vegal, ex_zeta,
				ncpns, cpn_times, cpn_dfs, cpn_cvgs, nexs, ex_times, ex_cpns, ex_ncpns, ex_strikes, weight, ex_prices,
				ex_zetas, one2F, nlam, lam_time, lambda, alpha, gamma, rho, paraml, params, lm_params);

			err  = levenberg_marquardt(
								data,
								target,
								weight,
								nexs / lm_params->freq_short,
								lambda,
								nlam,
								lm_params->nb_iter,
								static_lgmprcapgiventauts_dlm,
								&fitting_error);

			for (i=1; i<=nlam; i++)
			{
				lam[i-1] = lambda[i];
			}
		}				
	}

	return NULL;
}

Err	Construct_MultiSchedule(
							char	*yc_name,
							long	today,
							char	*instr_freq,
							char	*instr_basis,
							int		num_ex_dates,					/*	Exercise dates */
							long	*ex_date,						/*	Supposed to be sorted */							
							char	**end_tenor,					/*	Tenors of the underlying instruments */
							long	end_date,

							long	*ncpn,
							long	*cpn_date,
							double	*cpn_time,
							double	*cpn_cvg,
							double	*cpn_df,
							long	*theo_end_dates,
							long	*act_end_dates)
{
int				i;
SrtCompounding	ifreq;
SrtBasisCode	ibasis;
long			theo_date, act_date, temp_date, temp_date2;
Err				err = NULL;

	err = interp_compounding (instr_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (instr_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Find the end date as the longest total maturity */
	theo_date = end_date;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	for (i=0; i<num_ex_dates; i++)
	{
		err = get_end_date (ex_date[i], end_date, end_tenor[i], 0, &(theo_end_dates[i]));
		if (err)
		{
			goto FREE_RETURN;
		}
		act_end_dates[i] = bus_date_method (theo_end_dates[i], MODIFIED_SUCCEEDING);	
	}

	for (i=0; i<num_ex_dates; i++)
	{
		if (theo_end_dates[i] > theo_date || act_end_dates[i] > act_date)
		{
			theo_date = theo_end_dates[i];
			act_date = act_end_dates[i];
		}
	}

	*ncpn = 1;
	temp_date = theo_date;
	temp_date2 = act_date;

	while (act_date > today)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		(*ncpn)++;
	}
	
	(*ncpn)--;

	if (*ncpn < 2)
	{
		err = "Not enough coupons in cpd_calib_diagonal";
		goto FREE_RETURN;		
	}

	theo_date = temp_date;
	act_date = temp_date2;
	i = *ncpn - 1;

	while (i >= 0)
	{
		cpn_time[i] = (act_date - today) * YEARS_IN_DAY; 
		cpn_date[i] = act_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		cpn_cvg[i] = coverage (temp_date, act_date, ibasis);
		cpn_df[i] = swp_f_df (today, act_date, yc_name);
		act_date = temp_date;

		i--;
	}

	cpn_cvg[0] = 0.0;

FREE_RETURN:

	return err;
}


Err	Reduct_ExerciseDates(	
							long	today,

							int		old_num_ex_date,					/*	Exercise dates */
							long	*old_ex_date,						/*	Supposed to be sorted */
							double	*old_strike,
							double	*old_weight,
							long	*old_theo_end_dates,
							long	*old_act_end_dates,
							int		ncpn,
							long	*cpn_date,

							int		*cal_date,
							CPD_DIAG_CALIB_PARAM param,

							int		*num_ex_date,
							long	**ex_date,
							double	**strike,
							double	**weight,
							long	**theo_end_dates,
							long	**act_end_dates)
{
int		i, j, k, l;
long	min_cal_date;
Err		err = NULL;

	/* first copy old into new */

	memcpy((*ex_date), old_ex_date, old_num_ex_date * sizeof(long));
	memcpy((*strike), old_strike, old_num_ex_date * sizeof(double));
	if (old_weight && weight) memcpy((*weight), old_weight, old_num_ex_date * sizeof(double));	
	memcpy((*theo_end_dates), old_theo_end_dates, old_num_ex_date * sizeof(long));
	memcpy((*act_end_dates), old_act_end_dates, old_num_ex_date * sizeof(long));

	(*num_ex_date) = old_num_ex_date;

	/*	Remove non-calibration dates */
	for (i=(*num_ex_date)-1; i>=0; i--)
	{		
		if (cal_date[i] == 0)
		{
			for (k=i-1; k>=0; k--)
			{
				(*ex_date)[k+1] = (*ex_date)[k];
				(*theo_end_dates)[k+1] = (*theo_end_dates)[k];
				(*act_end_dates)[k+1] = (*act_end_dates)[k];
				(*strike)[k+1] = (*strike)[k];
				if (old_weight && weight) (*weight)[k+1] = (*weight)[k];				
			}

			(*ex_date)++;
			(*theo_end_dates)++;
			(*act_end_dates)++;
			(*strike) ++;
			if (old_weight && weight) (*weight) ++;

			
			(*num_ex_date)--;			
			if ((*num_ex_date) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
	}

	/* remove short maturities */
	min_cal_date = (long) (today + param->min_calib_time * DAYS_IN_YEAR + 0.5);

	i = 0;
	while ((*ex_date)[i] < min_cal_date && i <(*num_ex_date)-1)
	{
		i++;

		(*ex_date)++;
		(*theo_end_dates)++;
		(*act_end_dates)++;
		(*strike) ++;
		if (old_weight && weight) (*weight) ++;

		
		(*num_ex_date)--;			
		if ((*num_ex_date) < 1)
		{
			err = "All exercise dates are past in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}

	
	/*	Remove redundant dates */
	j = ncpn - 1;
	l = ncpn + 1;
	for (i=(*num_ex_date)-1; i>=0; i--)
	{
		while (j > 0 && cpn_date[j] > (*ex_date)[i])
		{
			j--;
		}
		if (cpn_date[j] < (*ex_date)[i] && j < ncpn-1)
		{
			j++;
		}
		
		if (j >= ncpn - 1 || j == l)
		{
			for (k=i-1; k>=0; k--)
			{
				(*ex_date)[k+1] = (*ex_date)[k];
				(*theo_end_dates)[k+1] = (*theo_end_dates)[k];
				(*act_end_dates)[k+1] = (*act_end_dates)[k];
				(*strike)[k+1] = (*strike)[k];
				if (old_weight && weight) (*weight)[k+1] = (*weight)[k];
			}

			(*ex_date)++;
			(*theo_end_dates)++;
			(*act_end_dates)++;
			(*strike)++;
			if (old_weight && weight) (*weight) ++;
			
			(*num_ex_date)--;			
			if ((*num_ex_date) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			l = j;
		}	
	}

	/*	Remove close dates */
	j = (*num_ex_date) - 1;
	for (i=(*num_ex_date)-2; i>=param->keep_first; i--)
	{
		if (((*ex_date)[j] - (*ex_date)[i]) * YEARS_IN_DAY < param->min_time - ONE_MONTH)
		{
			for (k=i-1; k>=0; k--)
			{
				(*ex_date)[k+1] = (*ex_date)[k];
				(*theo_end_dates)[k+1] = (*theo_end_dates)[k];
				(*act_end_dates)[k+1] = (*act_end_dates)[k];
				(*strike)[k+1] = (*strike)[k];
				if (old_weight && weight) (*weight)[k+1] = (*weight)[k];
			}

			(*ex_date)++;
			(*theo_end_dates)++;
			(*act_end_dates)++;
			(*strike)++;
			if (old_weight && weight) (*weight) ++;
			
			(*num_ex_date)--;	
			j--;
			if ((*num_ex_date) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			j = i;
		}
	}

	/*	Remove last? */
	if (param->skip_last && (*num_ex_date) > 1)
	{
		(*num_ex_date)--;
	}

	/*	Remove past dates */
	while ((*ex_date)[0] <= today)
	{
		(*ex_date)++;
		(*theo_end_dates)++;
		(*act_end_dates)++;
		(*strike)++;
		if (old_weight && weight) (*weight) ++;
		
		(*num_ex_date)--;	
		j--;
		if ((*num_ex_date) < 1)
		{
			err = "All exercise dates are past in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}

FREE_RETURN:

	return err;
}

Err	shift_the_vol(double				dFwd,
				  double				dStrike,
				  double				dMaturity,
				  double				dInputVol,
				  double				dInputVolType,
									
				  double				dVolShift,
				  DIAGCALIB_VOLTYPE		eShiftVolType,
				  DIAGCALIB_SHIFTTYPE	eShiftType,
				  double				*dOutputVol,
				  double				*dOutputVolType)
{
double	dPrice;
Err		err = NULL;
	
	if (fabs(dVolShift) < 1.0E-08 || eShiftVolType == LGM_VOL)
	{
		*dOutputVol = dInputVol;
		*dOutputVolType = dInputVolType;
	}
	else
	{
		if (dInputVolType < 0.5 && eShiftVolType == LOGNORMAL_VOL)
		{
			dPrice = srt_f_optblknrm(dFwd, dStrike, dInputVol, dMaturity, 1.0, SRT_PUT, PREMIUM);

			err = srt_f_optimpvol(dPrice, dFwd, dStrike, dMaturity, 1.0, SRT_PUT, SRT_LOGNORMAL, &dInputVol);

			if (err) return err;						
		}
		else if (dInputVolType > 0.5 && eShiftVolType == NORMAL_VOL)
		{
			dPrice = srt_f_optblksch(dFwd, dStrike, dInputVol, dMaturity, 1.0, SRT_PUT, PREMIUM);

			err = srt_f_optimpvol(dPrice, dFwd, dStrike, dMaturity, 1.0, SRT_PUT, SRT_NORMAL, &dInputVol);

			if (err) return err;			
		}			

		if (eShiftType == MULTIPLICATIVE)
		{
			*dOutputVol = dInputVol * (1.0 + dVolShift);
		}
		else
		{
			*dOutputVol = dInputVol + dVolShift;
		}

		if (eShiftVolType == LOGNORMAL_VOL)
		{
			*dOutputVolType = 1.0;
		}
		else
		{
			*dOutputVolType = 0.0;
		}
	}

	return err;
}

Err	Calculate_Underlying(	long	today,
							char	*yc_name,						/*	Name of the yield curve */
							char	*vol_curve_name,				/*	Name of the market vol curve */	
							Err		(*get_cash_vol)(				/*	Function to get cash vol from the market */
												char	*vol_curve_name,	
												double	start_date, 
												double	end_date,
												double	cash_strike,
												int		zero,
												char	*ref_rate_name,
												double	*vol,
												double	*power),
							char	*vol_ref_rate_name,
							char	*instr_freq,
							char	*instr_basis,
							char	*ref_rate_name,
							
							CPD_DIAG_CALIB_PARAM param,

							long	ncpn,
							long	*cpn_date,
							double	*cpn_time,
							double	*cpn_cvg,
							double	*cpn_df,
							int		num_ex_dates,
							long	*ex_date,
							long	*theo_end_dates,
							long	*act_end_dates,
							double	*strike,
							double	*weights,
							int		vega_weight,

							/* shift informations */
							int		num_vol_shift,
							double	*vol_shift_time,
							double	*vol_shift,

							/* information in the case of the short instruments */
							int		num_ex_primdates,
							double	*ex_primtime,
							double	*ex_primnstd,

							double	*ex_time,
							long	*ex_startcpn,
							long	*ex_endcpn,
							double	*ex_lvl,
							double	*ex_fwd,
							double	*ex_strike,
							double	*ex_weights,
							double	*ex_price,
							double	*ex_atm_price,
							double	*ex_vega,
							double	*ex_vol,
							double	*ex_nstd)
{
int		i, j, k, l;
double	lvl, dfi, dff, std, power, swp_rte, spr;
double	shift;
SrtCcyParam *ptrCcyParam = 0;
SrtCrvPtr crvPtr;
char* szStdBasis = 0;
char* szStdFreq = 0;
Err	err = NULL;

/* if we are transforming the vol, then get the CCY defaults */
	if ( param->transform_vol )
	{
		crvPtr = lookup_curve( yc_name );
		ptrCcyParam = new_CcyParam();
		if ( err = swp_f_get_CcyParam_from_CcyStr( crvPtr->curve_ccy, &ptrCcyParam) )
			goto FREE_RETURN;
		if ( err = translate_basis(&szStdBasis, ptrCcyParam->swap_basis_code) ) 
			goto FREE_RETURN;
		if ( translate_compounding(&szStdFreq, ptrCcyParam->compd) )
			goto FREE_RETURN;
	}

	num_ex_dates = num_ex_dates;
	j = 0;		
	for (i=0; i<num_ex_dates; i++)
	{
		while (cpn_date[j] < ex_date[i])
		{
			j++;
		}

		ex_startcpn[i] = j;
		ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;

		k = j;
		while (cpn_date[k] < act_end_dates[i])
		{
			k++;
		}
		if (k > 0 && cpn_date[k] - act_end_dates[i] > act_end_dates[i] - cpn_date[k-1])
		{
			k--;
		}

		if (k <= j)
		{
			k = j + 1;
		}
		ex_endcpn[i] = k;
		
		if (j >= ncpn || k >= ncpn)
		{
			err = "Coupon date bug in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}

	for (i=0; i<num_ex_dates; i++)
	{
		j = ex_startcpn[i];
		l = ex_endcpn[i];

		lvl = 0.0;
		for (k=j+1; k<=l; k++)
		{
			lvl += cpn_cvg[k] * cpn_df[k];
		}
		dfi = swp_f_df (today, cpn_date[j], yc_name);
		dff = swp_f_df (today, cpn_date[l], yc_name);

		ex_lvl[i] = lvl;
		ex_fwd[i] = (dfi - dff) / lvl;

		/*	ATM std */
		if ( param->transform_vol )
			err = swp_f_transformvol_func(	yc_name, 
											vol_curve_name,
											get_cash_vol,
											szStdFreq, 
											szStdBasis,
											vol_ref_rate_name,
											(double) add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
											(double) theo_end_dates[i], 
											ex_fwd[i], 
											instr_freq,
											instr_basis,
											&std, 
											&power);
		else
			err = get_cash_vol (
				vol_curve_name, 
				add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
				theo_end_dates[i], 
				ex_fwd[i], 
				0, 
				vol_ref_rate_name, 
				&std, 
				&power);

		if (err)
		{
			goto FREE_RETURN;
		}

		/* Apply the shift if relevant */
		if (vol_shift_time && vol_shift)
		{
			/* Linear Interpolation on the shift */
			constant_or_linear_interpolation_dlm(	vol_shift_time,
													vol_shift,
													num_vol_shift,
													ex_time[i],
													&shift);			
		}
		else
		{
			shift = 0.0;
		}

		err = shift_the_vol(ex_fwd[i], 
							ex_fwd[i],
							ex_time[i], 
							std,
							power,
							param->vol_shift + shift,
							param->vol_type,
							param->shift_type,
							&std,
							&power);

		if (err) goto FREE_RETURN;				

		if (power > 0.5)
		{
			power = srt_f_optblksch(
				ex_fwd[i], 
				ex_fwd[i], 
				std, 
				ex_time[i], 
				1.0, 
				SRT_CALL, 
				PREMIUM);

			err = srt_f_optimpvol(
				power,
				ex_fwd[i],
				ex_fwd[i],
				ex_time[i],
				1.0,
				SRT_CALL,
				SRT_NORMAL,
				&std);
		}

		if (ex_atm_price)
		{
			ex_atm_price[i] = srt_f_optblknrm(
											ex_fwd[i],
											ex_fwd[i],
											std,
											ex_time[i],
											ex_lvl[i],
											SRT_PUT,
											PREMIUM);
		}

		std *= sqrt (ex_time[i]);
		
		/*	Strike */
		if (param->strike_type == 0 || param->strike_type == 5 || (strike[i] < 1.0e-04 && param->strike_type != 4))
		{
			ex_strike[i] = ex_fwd[i];
		}
		else
		if (param->strike_type == 1)
		{
			ex_strike[i] = strike[i];
		}
		else
		if (param->strike_type == 2)
		{
			if (err = swp_f_ForwardRate(
				cpn_date[j],
				theo_end_dates[i],
				instr_freq,
				instr_basis,
				yc_name,
				ref_rate_name,
				&swp_rte))
			{
				goto FREE_RETURN;
			}

			spr = swp_rte - ex_fwd[i];
			
			ex_strike[i] = strike[i] - spr;
		}
		else
		if (param->strike_type == 3)
		{
			ex_strike[i] = ex_fwd[i] + strike[i] * std;
		}
		else
		if (param->strike_type == 4)
		{
			ex_strike[i] = ex_fwd[i] + ex_primnstd[Get_Index(ex_time[i], ex_primtime, num_ex_primdates)] * std;
		}

		ex_nstd[i] = (ex_strike[i] - ex_fwd[i]) / std;

		/*	Apply max std */
		if (ex_strike[i] > ex_fwd[i] + param->max_std * std)
		{
			ex_strike[i] = ex_fwd[i] + param->max_std * std;
		}
		else if (ex_strike[i] < ex_fwd[i] - param->max_std * std)
		{
			ex_strike[i] = ex_fwd[i] - param->max_std * std;
		}

		/*	Make sure strikes are positive (actually more than 1bp)
				otherwise use ATM	*/
		if (ex_strike[i] < 1.0e-04)
		{
			ex_strike[i] = ex_fwd[i];
		}

		if ( param->transform_vol )
			err = swp_f_transformvol_func(	yc_name, 
											vol_curve_name,
											get_cash_vol,
											szStdFreq, 
											szStdBasis,
											vol_ref_rate_name,
											(double) add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
											(double) theo_end_dates[i], 
											ex_strike[i], 
											instr_freq,
											instr_basis,
											&(ex_vol[i]), 
											&power);
		else
			err = get_cash_vol (
				vol_curve_name, 
				add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
				theo_end_dates[i], 
				ex_strike[i], 
				0, 
				vol_ref_rate_name, 
				&(ex_vol[i]), 
				&power);

		if (err)
		{
			goto FREE_RETURN;
		}

		/* Apply the shift if relevant */
		if (vol_shift_time && vol_shift)
		{
			/* Linear Interpolation on the shift */
			constant_or_linear_interpolation_dlm(	vol_shift_time,
												vol_shift,
												num_vol_shift,
												ex_time[i],
												&shift);			
		}
		else
		{
			shift = 0.0;
		}

		err = shift_the_vol(ex_fwd[i], 
							ex_strike[i],
							ex_time[i], 
							ex_vol[i],
							power,
							param->vol_shift + shift,
							param->vol_type,
							param->shift_type,
							&(ex_vol[i]),
							&power);

		if (err) goto FREE_RETURN;		

		if (power > 0.5)
		{
			ex_price[i] = srt_f_optblksch(
				ex_fwd[i],
				ex_strike[i],
				ex_vol[i],
				ex_time[i],
				ex_lvl[i],
				SRT_PUT,
				PREMIUM);

			ex_vega[i] = srt_f_optblksch(
				ex_fwd[i],
				ex_strike[i],
				ex_vol[i],
				ex_time[i],
				ex_lvl[i],
				SRT_PUT,
				VEGA);
		}
		else
		{
			ex_price[i] = srt_f_optblknrm(
				ex_fwd[i],
				ex_strike[i],
				ex_vol[i],
				ex_time[i],
				ex_lvl[i],
				SRT_PUT,
				PREMIUM);

			ex_vega[i] = srt_f_optblknrm(
				ex_fwd[i],
				ex_strike[i],
				ex_vol[i],
				ex_time[i],
				ex_lvl[i],
				SRT_PUT,
				VEGA);
		}

		if (weights && ex_weights)
		{
			ex_weights[i] = weights[i];
		}
		else
		if (ex_weights)
		{
			if (vega_weight)
			{
				ex_weights[i] = ex_vega[i];
			}
			else
			{
				ex_weights[i] = 1.0;
			}
		}
	}

FREE_RETURN:
	if ( ptrCcyParam ) free_CcyParam( ptrCcyParam );
	return err;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Calibrate lgm: main function											   //
//	New version: calibrates not necessarily to diagonal						   //
//		with lambda calibration												   //
//	Now save market instruments prices in the inst_data structure PMc 17Nov03  //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
Err cpd_calib_diagonal_dlm(
	/*	Market */
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */	
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	
	/* Get vol ref */
	char			*vol_ref_rate_name,
	
	/* Long Instruments */
	char			*instr_long_freq,					/*	Frequency and basis of instruments */
	char			*instr_long_basis,
	char			*long_ref_rate_name,				/*	Name of the reference rate */

	int				num_ex_datesl,						/*	Long Exercise dates */
	long			*ex_datel_,							/*	Supposed to be sorted */
	int				*cal_datel,							/*	1: use ex_date as calibration date, 0: don't */
	char			**end_tenorl,						/*	Tenors of the underlying instruments */
	long			end_datel,							/*	End date for diagonal */
	double			*strikel_,							/*	Strikes */		

	CPD_DIAG_CALIB_PARAM paraml,

	/* Short Instruments */
	char			*instr_short_freq,					/*	Frequency and basis of instruments */
	char			*instr_short_basis,
	char			*short_ref_rate_name,				/*	Name of the reference rate */

	int				num_ex_datess,						/*	Short Exercise dates */
	long			*ex_dates_,							/*	Supposed to be sorted */
	int				*cal_dates,							/*	1: use ex_date as calibration date, 0: don't */
	char			**end_tenors,						/*	Tenors of the underlying instruments */
	long			end_dates,							/*	End date for diagonal */
	double			*strikes_,							/*	Strikes */
	double			*weights_,

	CPD_DIAG_CALIB_PARAM params,
						
	/*	Model */
	int				fix_lambda,
	int				one_factor_equi,				/*	Calibrate 2 Factor to 1 Factor price */
	int				nlam,							/*	Lambda TS: may NOT be changed in the process */
	double			lam_time[],
	double			lam[],
	double			lam_shift[],
	int				one2F,							/*	Number of factors */
	double			alpha,							/*	Alpha, Gamma, Rho (2F only) */
	double			gamma,
	double			rho,	

	/*	Shift Parameters */
	int				num_long_vol_shift,
	double			*long_vol_shift_time,
	double			*long_vol_shift,

	int				num_short_vol_shift,
	double			*short_vol_shift_time,
	double			*short_vol_shift,

	/*	Output */
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	
	/*	Parameters */	
	DIAG_CALIB_LM_PARAMS lm_params,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data)					/*	NULL = don't save calibration instrument data */
{
	int				ncpnl, ncpns;	

	double			ex_timel[MAX_INST],
					ex_fwdl[MAX_INST],
					ex_lvll[MAX_INST],
					ex_strikel[MAX_INST],
					ex_voll[MAX_INST],
					ex_pricel[MAX_INST],
					ex_atm_pricel[MAX_INST],
					ex_vegal[MAX_INST],
					ex_atm_vegal[MAX_INST],
					ex_nstdl[MAX_INST],

					ex_times[MAX_INST],
					ex_fwds[MAX_INST],
					ex_lvls[MAX_INST],
					ex_strikes[MAX_INST],
					ex_weights[MAX_INST],
					ex_prices[MAX_INST],
					ex_atm_prices[MAX_INST],
					ex_vols[MAX_INST],
					ex_vegas[MAX_INST],
					ex_nstds[MAX_INST],

					ex_zeta[MAX_INST],
					ex_zeta_temp[MAX_INST],
					ex_G_temp[MAX_INST];

	int				ex_startcpnl[MAX_INST],
					ex_endcpnl[MAX_INST],
					ex_startcpns[MAX_INST],
					ex_endcpns[MAX_CPN];

	long			cpn_datel[MAX_CPN],
					cpn_dates[MAX_CPN];

	double			cpn_timel[MAX_CPN],
					cpn_cvgl[MAX_CPN],
					cpn_dfl[MAX_CPN],
					cpn_times[MAX_CPN],
					cpn_cvgs[MAX_CPN],
					cpn_dfs[MAX_CPN],
					cpn_G_temp[MAX_CPN];

	long			tmplngl1[MAX_CPN],
					tmplngl2[MAX_CPN];
	long			tmplngs1[MAX_CPN],
					tmplngs2[MAX_CPN];
	long			*theo_end_datesl,
					*act_end_datesl,
					*theo_end_datess,
					*act_end_datess;

	int				shift_lam;

	double			*used_strikel,
					*used_pricel,
					*used_vegal,
					*used_strikes,
					*used_prices;
	
	long			today; 
	
	double			exp_fact, lam_sens;
	
	SrtCurvePtr		yc_ptr;

	int				i, j;

	Err				err				= NULL;

	long			ex_datel__[MAX_INST], *ex_datel;
	double			strikel__[MAX_INST], *strikel;
	long			ex_dates__[MAX_INST], *ex_dates;
	double			strikes__[MAX_INST], *strikes;
	double			weights__[MAX_INST], *weights;

	/*	Copy data so as not to change the original */
	ex_datel = &(ex_datel__[0]);
	strikel = &(strikel__[0]);	
	ex_dates = &(ex_dates__[0]);
	strikes = &(strikes__[0]);
	weights = &(weights__[0]);

	theo_end_datesl = &(tmplngl1[0]);
	act_end_datesl = &(tmplngl2[0]);

	theo_end_datess = &(tmplngs1[0]);
	act_end_datess = &(tmplngs2[0]);

	*sig_time = NULL;
	*sig = NULL;

	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}
	today = get_today_from_curve (yc_ptr);

	if (one2F == 2)
	{
		HermiteStandard (x, w, NUM_HERMITE);
	}

	/* Long instruments */
	err = Construct_MultiSchedule(yc_name,
								today,
								instr_long_freq,
								instr_long_basis,
								num_ex_datesl,
								ex_datel_,
								end_tenorl,
								end_datel,
								&ncpnl,
								cpn_datel,
								cpn_timel,
								cpn_cvgl,
								cpn_dfl,
								theo_end_datesl,
								act_end_datesl);

	if (err)
	{
		goto FREE_RETURN;
	}

	err = Reduct_ExerciseDates(	today,
								num_ex_datesl,
								ex_datel_,
								strikel_,
								NULL,
								theo_end_datesl,
								act_end_datesl,
								ncpnl,
								cpn_datel,
								cal_datel,
								paraml,
								&num_ex_datesl,
								&ex_datel,
								&strikel,
								NULL,
								&theo_end_datesl,
								&act_end_datesl);

	if (err)
	{
		goto FREE_RETURN;
	}
	
	err = Calculate_Underlying(	today,
								yc_name,
								vol_curve_name,
								get_cash_vol,
								vol_ref_rate_name,
								instr_long_freq,
								instr_long_basis,
								long_ref_rate_name,
								paraml,
								ncpnl,
								cpn_datel,
								cpn_timel,
								cpn_cvgl,
								cpn_dfl,
								num_ex_datesl,
								ex_datel,
								theo_end_datesl,
								act_end_datesl,
								strikel,
								NULL,
								0,

								num_long_vol_shift,
								long_vol_shift_time,
								long_vol_shift,

								num_ex_datesl,
								ex_timel,
								ex_nstdl,

								ex_timel,
								ex_startcpnl,
								ex_endcpnl,
								ex_lvll,
								ex_fwdl,
								ex_strikel,
								NULL,
								ex_pricel,
								ex_atm_pricel,
								ex_atm_vegal,
								ex_voll,
								ex_nstdl);

	if (err) goto FREE_RETURN;

	/* Short instruments */

	if (!fix_lambda || (one_factor_equi && one2F == 2))
	{
		err = Construct_MultiSchedule(yc_name,
									today,
									instr_short_freq,
									instr_short_basis,
									num_ex_datess,
									ex_dates_,
									end_tenors,
									end_dates,
									&ncpns,
									cpn_dates,
									cpn_times,
									cpn_cvgs,
									cpn_dfs,
									theo_end_datess,
									act_end_datess);

		if (err)
		{
			goto FREE_RETURN;
		}

		err = Reduct_ExerciseDates(	today,
									num_ex_datess,
									ex_dates_,
									strikes_,
									weights_,
									theo_end_datess,
									act_end_datess,
									ncpns,
									cpn_dates,
									cal_dates,
									params,
									&num_ex_datess,
									&ex_dates,
									&strikes,
									&weights,
									&theo_end_datess,
									&act_end_datess);

		if (err)
		{
			goto FREE_RETURN;
		}

		if (!weights_) weights = NULL;

		err = Calculate_Underlying(	today,
									yc_name,
									vol_curve_name,
									get_cash_vol,
									vol_ref_rate_name,
									instr_short_freq,
									instr_short_basis,
									short_ref_rate_name,
									params,
									ncpns,
									cpn_dates,
									cpn_times,
									cpn_cvgs,
									cpn_dfs,
									num_ex_datess,
									ex_dates,
									theo_end_datess,
									act_end_datess,
									strikes,
									weights,
									lm_params->vega_weight,

									num_short_vol_shift,
									short_vol_shift_time,
									short_vol_shift,

									num_ex_datesl,
									ex_timel,
									ex_nstdl,

									ex_times,
									ex_startcpns,
									ex_endcpns,
									ex_lvls,
									ex_fwds,
									ex_strikes,
									ex_weights,
									ex_prices,
									ex_atm_prices,
									ex_vegas,
									ex_vols,
									ex_nstds);

		if (err)
		{
			goto FREE_RETURN;
		}		
	}

	if (one_factor_equi && one2F == 2 && fix_lambda)
	{
		/* First Calibrate the one factor model */
		err = lgmcalibzetalambda_tauts_dlm(	ncpnl,
											cpn_timel,
											cpn_dfl,
											cpn_cvgl,
											num_ex_datesl,
											ex_timel,
											ex_startcpnl,
											ex_endcpnl,
											ex_strikel,
											ex_pricel,
											ex_vegal,

											ncpns,
											cpn_times,
											cpn_dfs,
											cpn_cvgs,
											num_ex_datess,
											ex_times,
											ex_startcpns,
											ex_endcpns,
											ex_strikes,
											ex_weights,
											ex_prices,
											ex_vegas,

											ex_zeta,
											fix_lambda,
											nlam,
											lam_time,
											lam,
											1,
											alpha,
											gamma,
											rho,
											paraml,
											params,
											lm_params,
											&lam_sens);

		if (err) goto FREE_RETURN;

		/* then reprice in the one factor each of the secondary instruments */
		
		static_interpolate_zeta(num_ex_datesl, ex_timel, ex_zeta, nlam, lam_time, lam, num_ex_datess, ex_times, ex_zeta_temp);
		static_lgmsetupG_tauts (nlam, lam_time, lam, ncpns, cpn_times, cpn_G_temp, num_ex_datess, ex_times, ex_G_temp);

		for (i=0; i<num_ex_datess; i++)
		{
			ex_prices[i] = lgmsopval1F(	ex_endcpns[i] - ex_startcpns[i] + 1,
										&(cpn_dfs[ex_startcpns[i]]),
										&(cpn_cvgs[ex_startcpns[i]]),
										&(cpn_G_temp[ex_startcpns[i]]),
										ex_zeta_temp[i],
										ex_G_temp[i],
										ex_strikes[i]);
		}

		fix_lambda = 0;
	}

	/*	2.)	Calibrate zeta */

	/* First seperate lambda calibration if needed */

	shift_lam = 0;

	for (i=0; i<nlam; i++)
	{
		if ((lam_shift && fabs(lam_shift[i]) > 1.0E-08) || (params && fabs(params->lambda_shift) > 1.0E-08))
		{
			shift_lam = 1;
			break;
		}
	}

	if (!fix_lambda && (params->strike_type == 5 || shift_lam))
	{		
		if (params->strike_type == 5)
		{
			/* lambda is first calibrated from ATM prices */
			used_strikel = ex_fwdl;
			used_pricel = ex_atm_pricel;
			used_vegal = ex_atm_vegal;
			used_strikes = ex_fwds;
			used_prices = ex_atm_prices;
		}
		else
		{
			/* lambda will be shifted */
			used_strikel = ex_strikel;
			used_pricel = ex_pricel;
			used_vegal = ex_vegal;
			used_strikes = ex_strikes;
			used_prices = ex_prices;
		}

		/* call the calibration function */
		err = lgmcalibzetalambda_tauts_dlm(	ncpnl,
											cpn_timel,
											cpn_dfl,
											cpn_cvgl,
											num_ex_datesl,
											ex_timel,
											ex_startcpnl,
											ex_endcpnl,
											used_strikel,
											used_pricel,
											used_vegal,

											ncpns,
											cpn_times,
											cpn_dfs,
											cpn_cvgs,
											num_ex_datess,
											ex_times,
											ex_startcpns,
											ex_endcpns,
											used_strikes,
											ex_weights,
											used_prices,
											ex_vegas,

											ex_zeta,
											fix_lambda,
											nlam,
											lam_time,
											lam,
											one2F,
											alpha,
											gamma,
											rho,											
											paraml,
											params,
											lm_params,
											&lam_sens);

		if (err) goto FREE_RETURN;		

		/* Fix the lambda */
		fix_lambda = 1;
	}

	if (fix_lambda && shift_lam)
	{
		/* Apply the shift before calibration */
		if (lam_shift)
		{
			lam[i] += lam_shift[i] + params->lambda_shift;			
		}
		else
		{
			lam[i] += params->lambda_shift;			
		}

		shift_lam = 0;
	}

	/* Recalibrate */
	err = lgmcalibzetalambda_tauts_dlm(
		
		ncpnl,
		cpn_timel,
		cpn_dfl,
		cpn_cvgl,
		num_ex_datesl,
		ex_timel,
		ex_startcpnl,
		ex_endcpnl,
		ex_strikel,
		ex_pricel,
		ex_vegal,

		ncpns,
		cpn_times,
		cpn_dfs,
		cpn_cvgs,
		num_ex_datess,
		ex_times,
		ex_startcpns,
		ex_endcpns,
		ex_strikes,
		ex_weights,
		ex_prices,
		ex_vegas,

		ex_zeta,
		fix_lambda,
		nlam,
		lam_time,
		lam,
		one2F,
		alpha,
		gamma,
		rho,
		paraml,
		params,
		lm_params,
		&lam_sens);

	if (err) goto FREE_RETURN;

	if (shift_lam)
	{
		/* Apply the shift */				
		for (i=0; i<nlam; i++)
		{
			if (lam_shift)
			{
				lam[i] += lam_shift[i] + params->lambda_shift;			
			}
			else
			{
				lam[i] += params->lambda_shift;			
			}
		}

		shift_lam = 0;
		fix_lambda = 1;

		err = lgmcalibzetalambda_tauts_dlm(		
											ncpnl,
											cpn_timel,
											cpn_dfl,
											cpn_cvgl,
											num_ex_datesl,
											ex_timel,
											ex_startcpnl,
											ex_endcpnl,
											ex_strikel,
											ex_pricel,
											ex_vegal,

											ncpns,
											cpn_times,
											cpn_dfs,
											cpn_cvgs,
											num_ex_datess,
											ex_times,
											ex_startcpns,
											ex_endcpns,
											ex_strikes,
											ex_weights,
											ex_prices,
											ex_vegas,

											ex_zeta,
											fix_lambda,
											nlam,
											lam_time,
											lam,
											one2F,
											alpha,
											gamma,
											rho,
											paraml,
											params,
											lm_params,
											&lam_sens);

		if (err) goto FREE_RETURN;
	}
	
	/*	3.)	Transform into sigma */

	*num_sig = num_ex_datesl;
	*sig_time = (double*) calloc (num_ex_datesl, sizeof (double));
	*sig = (double*) calloc (num_ex_datesl, sizeof (double));
	
	if (!sig_time || !sig)
	{
		err = "Allocation error (3) in cpd_calib_diagonal";
		goto FREE_RETURN;
	}

	(*sig_time)[0] = ex_timel[0];
	exp_fact = static_lgmcalcexpfact_tauts (0.0, ex_timel[0], nlam, lam_time, lam);
	(*sig)[0] = sqrt (ex_zeta[0] / exp_fact);

	for (i=1; i<num_ex_datesl; i++)
	{
		(*sig_time)[i] = ex_timel[i];
		if (ex_zeta[i] > ex_zeta[i-1])
		{
			exp_fact = static_lgmcalcexpfact_tauts (ex_timel[i-1], ex_timel[i], nlam, lam_time, lam);
			(*sig)[i] = sqrt ((ex_zeta[i] - ex_zeta[i-1]) / exp_fact);
		}
		else
		{
			smessage ("Diagonal calibration failed at exercise year %.2f - Calibration stopped", ex_timel[i]);
			for (j=i; j<num_ex_datesl; j++)
			{
				(*sig)[j] = (*sig)[i-1];
			}
			i = num_ex_datesl;
		}
	}

	/* 3.5) Shift the Vol if required */
	if (fabs(paraml->vol_shift) > 1.0E-08 && paraml->vol_type == LGM_VOL)
	{
		for (i=0; i<num_ex_datesl; i++)
		{
			if (paraml->shift_type == MULTIPLICATIVE)
			{
				(*sig)[j] += (*sig)[j] * paraml->vol_shift;
			}
			else
			{
				(*sig)[j] += paraml->vol_shift;
			}
		}
	}			

	/*	4.)	Save instrument data if required */
	if (inst_data)
	{				
		inst_data->num_inst = num_ex_datesl;
		inst_data->exer_dates_long = (long*) calloc (num_ex_datesl, sizeof (long));
		inst_data->start_dates = (long*) calloc (num_ex_datesl, sizeof (long));
		inst_data->end_dates = (long*) calloc (num_ex_datesl, sizeof (long));		
		inst_data->long_strikes = (double*) calloc (num_ex_datesl, sizeof (double));
		inst_data->market_prices_long = (double*) calloc (num_ex_datesl, sizeof (double));

		if (!fix_lambda)
		{
			inst_data->num_insts = num_ex_datess;
			inst_data->exer_dates_short = (long*) calloc (num_ex_datess, sizeof (long));
			inst_data->start_datess = (long*) calloc (num_ex_datess, sizeof (long));
			inst_data->end_datess = (long*) calloc (num_ex_datess, sizeof (long));		
			inst_data->short_strikes = (double*) calloc (num_ex_datess, sizeof (double));
			inst_data->short_weights = (double*) calloc (num_ex_datess, sizeof (double));
			inst_data->market_prices_short = (double*) calloc (num_ex_datess, sizeof (double));

			if (!inst_data->exer_dates_long
				|| !inst_data->start_dates 
				|| !inst_data->end_dates
				|| !inst_data->long_strikes
				|| !inst_data->market_prices_long
				|| !inst_data->exer_dates_short
				|| !inst_data->start_datess
				|| !inst_data->end_datess
				|| !inst_data->short_strikes
				|| !inst_data->short_weights
				|| !inst_data->market_prices_short)
			{
				err = "Allocation error (4) in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			inst_data->num_insts = 0;

			if (!inst_data->exer_dates_long
				|| !inst_data->start_dates 
				|| !inst_data->end_dates
				|| !inst_data->long_strikes
				|| !inst_data->market_prices_long)
			{
				err = "Allocation error (4) in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}


		for (i=0; i<num_ex_datesl; i++)
		{
			inst_data->exer_dates_long[i] = (long) (ex_timel[i] * DAYS_IN_YEAR + EPS_MAT_TO_LONG_DATE_CONVERSION) + today;
			inst_data->start_dates[i] = cpn_datel[ex_startcpnl[i]];
			inst_data->end_dates[i] = cpn_datel[ex_endcpnl[i]];			
			inst_data->long_strikes[i] = ex_strikel[i];
			inst_data->market_prices_long[i] = ex_pricel[i];
		}
			
		if (!fix_lambda)
		{
			for (i=0; i<num_ex_datess; i++)
			{
				inst_data->exer_dates_short[i] = (long) (ex_times[i] * DAYS_IN_YEAR + EPS_MAT_TO_LONG_DATE_CONVERSION) + today;
				inst_data->start_datess[i] = cpn_dates[ex_startcpns[i]];
				inst_data->end_datess[i] = cpn_dates[ex_endcpns[i]];			
				inst_data->short_strikes[i] = ex_strikes[i];
				inst_data->market_prices_short[i] = ex_prices[i];
			}
		}
	}

FREE_RETURN:

	if (err)
	{
		if (*sig_time) free (*sig_time);
		*sig_time = NULL;

		if (*sig) free (*sig);
		*sig = NULL;

		if (inst_data)
		{
			cpd_free_calib_inst_data (inst_data);
		}
	}

	return err;
}

Err	AllocateCalibInst(	int					iNbStrike,
						CALIBINSTRUMENTDLM	CalibInst)
{
Err	err = NULL;

	CalibInst->iNbStrike = iNbStrike;

	CalibInst->dStrike = NULL;
	CalibInst->dStrike = calloc(iNbStrike, sizeof(double));

	CalibInst->dPrice = NULL;
	CalibInst->dPrice = calloc(iNbStrike, sizeof(double));

	CalibInst->dNormalVol = NULL;
	CalibInst->dNormalVol = calloc(iNbStrike, sizeof(double));

	CalibInst->dVega = NULL;
	CalibInst->dVega = calloc(iNbStrike, sizeof(double));

	CalibInst->dNbStd = NULL;
	CalibInst->dNbStd = calloc(iNbStrike, sizeof(double));

	CalibInst->dModelPrice = NULL;
	CalibInst->dModelPrice = calloc(iNbStrike, sizeof(double));

	if (!CalibInst->dStrike || !CalibInst->dPrice || !CalibInst->dNormalVol
		|| !CalibInst->dVega || !CalibInst->dNbStd || !CalibInst->dModelPrice)
	{
		err = "Memory allocation failure in AllocateCalibInst";
		return err;
	}

	return err;
}

void FreeCalibInst(	CALIBINSTRUMENTDLM	CalibInst)
{
	if (CalibInst)
	{
		if (CalibInst->dStrike) free(CalibInst->dStrike);
		CalibInst->dStrike = NULL;

		if (CalibInst->dPrice) free(CalibInst->dPrice);
		CalibInst->dPrice = NULL;

		if (CalibInst->dNormalVol) free(CalibInst->dNormalVol);
		CalibInst->dNormalVol = NULL;

		if (CalibInst->dVega) free(CalibInst->dVega);
		CalibInst->dVega = NULL;

		if (CalibInst->dNbStd) free(CalibInst->dNbStd);
		CalibInst->dNbStd = NULL;	

		if (CalibInst->dModelPrice) free(CalibInst->dModelPrice);
		CalibInst->dModelPrice = NULL;	
	}
}

Err	AllocateAllCalibInst(int					iNbInst,
						 int					iNbStrike,
						 ALLCALININSTRUMENTSDLM	AllCalibInst)
{
int	i;
Err	err = NULL;

	AllCalibInst->sCalibInst = calloc(iNbInst, sizeof(CalibInstrumentDLM));

	if (!AllCalibInst)
	{
		err = "Memory Allocation faillure in AllocateAllCalibInst";
		return err;
	}

	AllCalibInst->iNbInst = iNbInst;

	for (i=0; i<iNbInst; i++)
	{
		err = AllocateCalibInst(iNbStrike,
								&(AllCalibInst->sCalibInst[i]));
		if (err)
		{
			return err;
		}
	}

	return err;
}

void FreeAllCalibInst(ALLCALININSTRUMENTSDLM	AllCalibInst)
{
int	i;

	if (AllCalibInst)
	{
		if (AllCalibInst->sCalibInst)
		{
			for (i=0; i<AllCalibInst->iNbInst; i++)
			{
				FreeCalibInst(&(AllCalibInst->sCalibInst[i]));
			}

			free(AllCalibInst->sCalibInst);
		}
	}
}

void	AllocateCalibExeSchedule(	long				*lExeDates,
									double				*dExeTimes,
									long				*lStartDates,
									long				*lTheoEndDates,
									long				*lActEndDates,
									double				*dStrikes,
									double				*dStrikesS1,
									double				*dStrikesS2,
									double				*dWeights,
									CALIBEXESCHEDULEDLM CalibExeSchedule)
{
	CalibExeSchedule->lExeDates = lExeDates;
	CalibExeSchedule->dExeTimes = dExeTimes;
	CalibExeSchedule->lStartDate = lStartDates;
	CalibExeSchedule->lTheoEndDates = lTheoEndDates;
	CalibExeSchedule->lActEndDates = lActEndDates;
	CalibExeSchedule->dStrikes = dStrikes;		
	CalibExeSchedule->dStrikesS1 = dStrikesS1;		
	CalibExeSchedule->dStrikesS2 = dStrikesS2;
	CalibExeSchedule->dWeights = dWeights;	
}

Err	Construct_CalibSchedule(
							char				*yc_name,
							long				today,
							char				*instr_freq,
							char				*instr_basis,
							int					num_ex_dates,					/*	Exercise dates */
							long				*ex_date,						/*	Supposed to be sorted */							
							int					*cal_date,						/*	1: use ex_date as calibration date, 0: don't */
							char				**end_tenor,					/*	Tenors of the underlying instruments */
							long				end_date,
							double				*strikes,
							double				*strikesS1,
							double				*strikesS2,
							double				*weights,

							CALIBCPNSCHEDULEDLM	CalibCpnSchedule,
							CALIBEXESCHEDULEDLM	CalibExeSchedule)
{
int				i,j;
SrtCompounding	ifreq;
SrtBasisCode	ibasis;
long			theo_date, act_date, temp_date, temp_date2;
Err				err = NULL;

	err = interp_compounding (instr_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (instr_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	CalibCpnSchedule->lToday = today;
	CalibExeSchedule->iNExe = num_ex_dates;
	
	for (i=0; i<num_ex_dates; i++)
	{
		err = get_end_date (ex_date[i], end_date, end_tenor[i], 0, &(CalibExeSchedule->lTheoEndDates[i]));
		if (err)
		{
			goto FREE_RETURN;
		}
		CalibExeSchedule->lActEndDates[i] = bus_date_method (CalibExeSchedule->lTheoEndDates[i], MODIFIED_SUCCEEDING);	
		CalibExeSchedule->lExeDates[i] = ex_date[i];
		CalibExeSchedule->dExeTimes[i] = (ex_date[i] - today) * YEARS_IN_DAY;

		if (strikes)
		{
			CalibExeSchedule->dStrikes[i] = strikes[i];
		}
				
		if (strikesS1 && CalibExeSchedule->dStrikesS1)
		{
			CalibExeSchedule->dStrikesS1[i] = strikesS1[i];
		}

		if (strikesS2 && CalibExeSchedule->dStrikesS2)
		{
			CalibExeSchedule->dStrikesS2[i] = strikesS2[i];
		}

		if (CalibExeSchedule->dWeights)
		{
			if (weights)
			{
				CalibExeSchedule->dWeights[i] = weights[i];
			}
			else
			{
				CalibExeSchedule->dWeights[i] = 1.0;
			}
		}
	}

	/*	Find the end date as the longest total maturity and cal_date == 1*/
	j = num_ex_dates-1;

	if (cal_date)
	{		
		while((cal_date[j] == 0) && (j>0))
		{
			j--;
		}
	}

	theo_date = CalibExeSchedule->lTheoEndDates[j];
	act_date = CalibExeSchedule->lActEndDates[j];

	for (i=0; i<j ; i++)
	{
		if (cal_date && cal_date[i])
		{
			if (CalibExeSchedule->lTheoEndDates[i] > theo_date || CalibExeSchedule->lActEndDates[i] > act_date)
			{
				theo_date = CalibExeSchedule->lTheoEndDates[i];
				act_date = CalibExeSchedule->lActEndDates[i];
			}
		}
	}	

	CalibCpnSchedule->iNCpn = 1;
	temp_date = theo_date;
	temp_date2 = act_date;

	while (act_date > today)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		(CalibCpnSchedule->iNCpn)++;
	}
	
	(CalibCpnSchedule->iNCpn)--;

	if (CalibCpnSchedule->iNCpn < 2)
	{
		err = "Not enough coupons in schedule: check your frequency";
		goto FREE_RETURN;		
	}

	theo_date = temp_date;
	act_date = temp_date2;
	i = CalibCpnSchedule->iNCpn - 1;

	while (i >= 0)
	{
		CalibCpnSchedule->dCpnTime[i] = (act_date - today) * YEARS_IN_DAY; 
		CalibCpnSchedule->lCpnDate[i] = act_date;
		CalibCpnSchedule->lCpnTheoDate[i] = theo_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		CalibCpnSchedule->dCpnCvg[i] = coverage (temp_date, act_date, ibasis);
		CalibCpnSchedule->dCpnDf[i] = swp_f_df (today, act_date, yc_name);
		act_date = temp_date;

		i--;
	}

	CalibCpnSchedule->dCpnCvg[0] = 0.0;

FREE_RETURN:

	return err;
}

Err	Reduct_ExeSchedule(	
						CALIBCPNSCHEDULEDLM		CalibCpnSchedule,
						CALIBEXESCHEDULEDLM		CalibExeSchedule,
						int						*cal_date,
						CPD_DIAG_CALIB_PARAM	param)
{
int		i, j, k, l;
Err		err = NULL;	

	/*	Remove non-calibration dates */
	for (i=(CalibExeSchedule->iNExe)-1; i>=0; i--)
	{		
		if (cal_date[i] == 0)
		{
			for (k=i-1; k>=0; k--)
			{
				(CalibExeSchedule->lExeDates)[k+1] = (CalibExeSchedule->lExeDates)[k];
				(CalibExeSchedule->dExeTimes)[k+1] = (CalibExeSchedule->dExeTimes)[k];
				(CalibExeSchedule->lTheoEndDates)[k+1] = (CalibExeSchedule->lTheoEndDates)[k];
				(CalibExeSchedule->lActEndDates)[k+1] = (CalibExeSchedule->lActEndDates)[k];
				(CalibExeSchedule->dStrikes)[k+1] = (CalibExeSchedule->dStrikes)[k];

				if (CalibExeSchedule->dStrikesS1)
				{
					CalibExeSchedule->dStrikesS1[k+1] = CalibExeSchedule->dStrikesS1[k];
				}

				if (CalibExeSchedule->dStrikesS2)
				{
					CalibExeSchedule->dStrikesS2[k+1] = CalibExeSchedule->dStrikesS2[k];
				}

				if (CalibExeSchedule->dWeights)
				{
					CalibExeSchedule->dWeights[k+1] = CalibExeSchedule->dWeights[k];
				}
			}

			(CalibExeSchedule->lExeDates)++;
			(CalibExeSchedule->dExeTimes)++;
			(CalibExeSchedule->lTheoEndDates)++;
			(CalibExeSchedule->lActEndDates)++;
			(CalibExeSchedule->dStrikes) ++;

			if (CalibExeSchedule->dStrikesS1) (CalibExeSchedule->dStrikesS1)++;
			if (CalibExeSchedule->dStrikesS2) (CalibExeSchedule->dStrikesS2)++;
			if (CalibExeSchedule->dWeights) (CalibExeSchedule->dWeights)++;

			(CalibExeSchedule->iNExe)--;
			if ((CalibExeSchedule->iNExe) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
	}

	/*	Remove short maturities */	
	while (CalibExeSchedule->iNExe > 0 && CalibExeSchedule->dExeTimes[0] < param->min_calib_time)
	{
		if (CalibExeSchedule->iNExe > 1 || CalibExeSchedule->dExeTimes[0] <= 0.0)
		{
			(CalibExeSchedule->lExeDates)++;
			(CalibExeSchedule->dExeTimes)++;
			(CalibExeSchedule->lTheoEndDates)++;
			(CalibExeSchedule->lActEndDates)++;
			(CalibExeSchedule->dStrikes) ++;
			
			if (CalibExeSchedule->dStrikesS1) (CalibExeSchedule->dStrikesS1)++;
			if (CalibExeSchedule->dStrikesS2) (CalibExeSchedule->dStrikesS2)++;
			if (CalibExeSchedule->dWeights) (CalibExeSchedule->dWeights)++;

			(CalibExeSchedule->iNExe)--;
		}
		else
			break;

		if ((CalibExeSchedule->iNExe) < 1)
		{
			err = "All exercise dates are past in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}
	
	/*	Remove redundant dates */
	j = CalibCpnSchedule->iNCpn - 1;
	l = CalibCpnSchedule->iNCpn + 1;
	for (i=(CalibExeSchedule->iNExe)-1; i>=0; i--)
	{
		while (j > 0 && CalibCpnSchedule->lCpnDate[j] > (CalibExeSchedule->lExeDates)[i])
		{
			j--;
		}
		if (CalibCpnSchedule->lCpnDate[j] < (CalibExeSchedule->lExeDates)[i] && j < CalibCpnSchedule->iNCpn - 1)
		{
			j++;
		}
		
		if (j >= CalibCpnSchedule->iNCpn - 1 || j == l)
		{
			for (k=i-1; k>=0; k--)
			{
				(CalibExeSchedule->lExeDates)[k+1] = (CalibExeSchedule->lExeDates)[k];
				(CalibExeSchedule->dExeTimes)[k+1] = (CalibExeSchedule->dExeTimes)[k];
				(CalibExeSchedule->lTheoEndDates)[k+1] = (CalibExeSchedule->lTheoEndDates)[k];
				(CalibExeSchedule->lActEndDates)[k+1] = (CalibExeSchedule->lActEndDates)[k];
				(CalibExeSchedule->dStrikes)[k+1] = (CalibExeSchedule->dStrikes)[k];

				if (CalibExeSchedule->dStrikesS1)
				{
					CalibExeSchedule->dStrikesS1[k+1] = CalibExeSchedule->dStrikesS1[k];
				}

				if (CalibExeSchedule->dStrikesS2)
				{
					CalibExeSchedule->dStrikesS2[k+1] = CalibExeSchedule->dStrikesS2[k];
				}

				if (CalibExeSchedule->dWeights)
				{
					CalibExeSchedule->dWeights[k+1] = CalibExeSchedule->dWeights[k];
				}
			}

			(CalibExeSchedule->lExeDates)++;
			(CalibExeSchedule->dExeTimes)++;
			(CalibExeSchedule->lTheoEndDates)++;
			(CalibExeSchedule->lActEndDates)++;
			(CalibExeSchedule->dStrikes)++;

			if (CalibExeSchedule->dStrikesS1) (CalibExeSchedule->dStrikesS1)++;
			if (CalibExeSchedule->dStrikesS2) (CalibExeSchedule->dStrikesS2)++;
			if (CalibExeSchedule->dWeights) (CalibExeSchedule->dWeights)++;
			
			(CalibExeSchedule->iNExe)--;			
			if ((CalibExeSchedule->iNExe) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			l = j;
		}	
	}

	/*	Remove close dates */
	j = (CalibExeSchedule->iNExe) - 1;
	for (i=(CalibExeSchedule->iNExe)-2; i>=param->keep_first; i--)
	{
		if (((CalibExeSchedule->lExeDates)[j] - (CalibExeSchedule->lExeDates)[i]) * YEARS_IN_DAY < param->min_time - ONE_MONTH)
		{
			for (k=i-1; k>=0; k--)
			{
				(CalibExeSchedule->lExeDates)[k+1] = (CalibExeSchedule->lExeDates)[k];
				(CalibExeSchedule->dExeTimes)[k+1] = (CalibExeSchedule->dExeTimes)[k];
				(CalibExeSchedule->lTheoEndDates)[k+1] = (CalibExeSchedule->lTheoEndDates)[k];
				(CalibExeSchedule->lActEndDates)[k+1] = (CalibExeSchedule->lActEndDates)[k];
				(CalibExeSchedule->dStrikes)[k+1] = (CalibExeSchedule->dStrikes)[k];

				if (CalibExeSchedule->dStrikesS1)
				{
					CalibExeSchedule->dStrikesS1[k+1] = CalibExeSchedule->dStrikesS1[k];
				}

				if (CalibExeSchedule->dStrikesS2)
				{
					CalibExeSchedule->dStrikesS2[k+1] = CalibExeSchedule->dStrikesS2[k];
				}

				if (CalibExeSchedule->dWeights)
				{
					CalibExeSchedule->dWeights[k+1] = CalibExeSchedule->dWeights[k];
				}
			}

			(CalibExeSchedule->lExeDates)++;
			(CalibExeSchedule->dExeTimes)++;
			(CalibExeSchedule->lTheoEndDates)++;
			(CalibExeSchedule->lActEndDates)++;
			(CalibExeSchedule->dStrikes)++;

			if (CalibExeSchedule->dStrikesS1) (CalibExeSchedule->dStrikesS1)++;
			if (CalibExeSchedule->dStrikesS2) (CalibExeSchedule->dStrikesS2)++;
			if (CalibExeSchedule->dWeights) (CalibExeSchedule->dWeights)++;
			
			(CalibExeSchedule->iNExe)--;	
			j--;
			if ((CalibExeSchedule->iNExe) < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			j = i;
		}
	}

	/*	Remove last? */
	if (param->skip_last && (CalibExeSchedule->iNExe) > 1)
	{
		(CalibExeSchedule->iNExe)--;
	}	

FREE_RETURN:

	return err;
}

Err	Calculate_CalibInst(long					today,
						char					*yc_name,						/*	Name of the yield curve */
						char					*vol_curve_name,				/*	Name of the market vol curve */	
						Err						(*get_cash_vol)(				/*	Function to get cash vol from the market */
															char	*vol_curve_name,	
															double	start_date, 
															double	end_date,
															double	cash_strike,
															int		zero,
															char	*ref_rate_name,
															double	*vol,
															double	*power),
						char					*vol_ref_rate_name,
						char					*instr_freq,
						char					*instr_basis,
						char					*ref_rate_name,
						
						int						is_for_calibration,
						
						CPD_DIAG_CALIB_PARAM	param,

						CALIBCPNSCHEDULEDLM		CalibCpnSchedule,
						CALIBEXESCHEDULEDLM		CalibExeSchedule,

						ALLCALININSTRUMENTSDLM	AllCalibInst,

						/* information in the case of the short instruments */
						ALLCALININSTRUMENTSDLM	AllCalibShortInst)
{
int		i, j, k, l, m, iLoopCoupon;
double	lvl, sumdf, dfi, dff, std, power, swp_rte, spr;
long	start_date;
SrtBasisCode		ibasis;
CALIBINSTRUMENTDLM	CalibInst;
double	*dStrikes;
double	vol_shift;

DIAGCALIB_VOLTYPE	vol_type;
DIAGCALIB_SHIFTTYPE shift_type;

int		strike_type;
double	dMainStrike;

Err	err = NULL;

	err = interp_basis (instr_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	j = 0;

	AllCalibInst->sCalibCpnSchedule = CalibCpnSchedule;
	AllCalibInst->sCalibExeSchedule = CalibExeSchedule;

	for (i=0; i<AllCalibInst->iNbInst; i++)
	{
		CalibInst = &(AllCalibInst->sCalibInst[i]);

		if (CalibExeSchedule->lStartDate[i] == 0)
		{
			start_date = add_unit(CalibExeSchedule->lExeDates[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);			
		}
		else
		{
			start_date = CalibExeSchedule->lStartDate[i];
		}

		/* we call the next coupon */
		while (j < CalibCpnSchedule->iNCpn && CalibCpnSchedule->lCpnDate[j] < start_date - 5)
		{
			j++;
		}
		if (j > 0 && fabs(start_date - CalibCpnSchedule->lCpnDate[j]) >
					fabs(start_date - CalibCpnSchedule->lCpnDate[j-1])
			&& CalibExeSchedule->lExeDates[i] < CalibCpnSchedule->lCpnDate[j-1])
		{
			j--;
		}
		
		CalibInst->iStartCpn = j;
		CalibExeSchedule->iStartCpn[i] = j;		

		k = j;
		while (k < CalibCpnSchedule->iNCpn && CalibCpnSchedule->lCpnDate[k] < CalibExeSchedule->lActEndDates[i] - 5)
		{
			k++;
		}
		if (k > 0 && fabs(CalibExeSchedule->lActEndDates[i] - CalibCpnSchedule->lCpnDate[k]) >
				fabs(CalibExeSchedule->lActEndDates[i] - CalibCpnSchedule->lCpnDate[k-1]))
		{
			k--;
		}

		if (k <= j)
		{
			k = j + 1;
		}

		CalibInst->iEndCpn = k;
		CalibExeSchedule->iEndCpn[i] = k;
		
		if (j >= CalibCpnSchedule->iNCpn || k >= CalibCpnSchedule->iNCpn)
		{
			err = "Coupon date bug in cpd_calib_diagonal";
			goto FREE_RETURN;
		}

		/* we calibrate to cash instruments */
		CalibInst->dSpread = 0.0;
		CalibInst->iNbCoupon = CalibInst->iEndCpn - CalibInst->iStartCpn + 1;
		CalibInst->dExeTime = CalibExeSchedule->dExeTimes[i];
		CalibInst->lExeDate = CalibExeSchedule->lExeDates[i];
		CalibInst->iIsCMS = 0;
		CalibInst->dPayTime = 0.0;
		CalibInst->lPayDate = 0;
		CalibInst->iIsLiborOption = 0;
	}

	for (i=0; i<AllCalibInst->iNbInst; i++)
	{
		CalibInst = &(AllCalibInst->sCalibInst[i]);

		j = CalibInst->iStartCpn;
		l = CalibInst->iEndCpn;		

		lvl = 0.0;
		sumdf = 0.0;

		for (k=j+1; k<=l; k++)
		{
			lvl += CalibCpnSchedule->dCpnCvg[k] * CalibCpnSchedule->dCpnDf[k];
			sumdf += CalibCpnSchedule->dCpnDf[k];
		}

		dfi = CalibCpnSchedule->dCpnDf[j];
		dff = CalibCpnSchedule->dCpnDf[l];

		CalibInst->dLevel = lvl;
		CalibInst->dSumDf = sumdf;
		CalibInst->dFwdCash = (dfi - dff) / lvl;
		CalibExeSchedule->dATMStrike[i] = CalibInst->dFwdCash;

		if (is_for_calibration)
		{
			/*	ATM std */
			err = get_cash_vol (
				vol_curve_name,
				CalibCpnSchedule->lCpnDate[j],
				CalibCpnSchedule->lCpnTheoDate[l],
				CalibInst->dFwdCash,
				0, 
				vol_ref_rate_name, 
				&std,
				&power);

			if (err) goto FREE_RETURN;

			err = shift_the_vol(CalibInst->dFwdCash, 
								CalibInst->dFwdCash, 								
								CalibInst->dExeTime, 
								std,
								power,
								param->vol_shift,
								param->vol_type,
								param->shift_type,
								&std,
								&power);

			if (err) goto FREE_RETURN;

			if (power > 0.5)
			{
				power = srt_f_optblksch(
					CalibInst->dFwdCash, 
					CalibInst->dFwdCash, 
					std, 
					CalibInst->dExeTime, 
					1.0, 
					SRT_PUT, 
					PREMIUM);

				err = srt_f_optimpvol(
					power,
					CalibInst->dFwdCash,
					CalibInst->dFwdCash,
					CalibInst->dExeTime,
					1.0,
					SRT_PUT,
					SRT_NORMAL,
					&std);
			}

			CalibInst->dATMVol = std;

			std *= sqrt (CalibInst->dExeTime);

			for (m=0; m<CalibInst->iNbStrike; m++)
			{
				switch (m)
				{
					case 0:
						dStrikes = CalibExeSchedule->dStrikes;
						vol_shift = param->vol_shift;
						vol_type = param->vol_type;
						shift_type = param->shift_type;
						strike_type = param->strike_type;
						break;
					case 1:
						dStrikes = CalibExeSchedule->dStrikesS1;
						vol_shift = param->smile_vol_shift;
						vol_type = param->smile_vol_type;
						shift_type = param->smile_shift_type;
						strike_type = param->smile_strike_type;
						break;
					case 2:
						dStrikes = CalibExeSchedule->dStrikesS2;
						vol_shift = param->smile_vol_shift;
						vol_type = param->smile_vol_type;
						shift_type = param->smile_shift_type;
						strike_type = param->smile_strike_type;
						break;
				}

				/*	Strike */
				if (strike_type == 0 || (fabs(dStrikes[i]) < 1.0e-04 && strike_type != 4))
				{
					CalibInst->dStrike[m] = CalibInst->dFwdCash;
				}
				else
				if (strike_type == 1)
				{
					CalibInst->dStrike[m] = dStrikes[i];
				}
				else
				if (strike_type == 2)
				{
					if (err = swp_f_ForwardRate(
						CalibCpnSchedule->lCpnDate[j],
						CalibCpnSchedule->lCpnTheoDate[l],
						instr_freq,
						instr_basis,
						yc_name,
						ref_rate_name,
						&swp_rte))
					{
						goto FREE_RETURN;
					}

					spr = swp_rte - CalibInst->dFwdCash;
					
					CalibInst->dStrike[m] = dStrikes[i] - spr;
				}
				else
				if (strike_type == 3)
				{
					CalibInst->dStrike[m] = CalibInst->dFwdCash + dStrikes[i] * std;
				}
				else
				if (strike_type == 4)
				{
					CalibInst->dStrike[m] = CalibInst->dFwdCash + (AllCalibShortInst->sCalibInst[Get_Index(CalibInst->dExeTime, AllCalibShortInst->sCalibExeSchedule->dExeTimes, AllCalibShortInst->sCalibExeSchedule->iNExe)]).dNbStd[m] * std;
				}
				else
				if (strike_type == 5)
				{
					if (m > 0)
					{
						CalibInst->dStrike[m] = CalibInst->dStrike[0] + dStrikes[i] * std;
					}
					else
					{
						CalibInst->dStrike[m] = CalibInst->dFwdCash + dStrikes[i] * std;
					}
				}

				CalibInst->dNbStd[m] = (CalibInst->dStrike[m] - CalibInst->dFwdCash) / std;

				if (m == 0)
				{
					/* Save Main Strike */
					dMainStrike = CalibInst->dStrike[0];
				}

				/*	Apply max std for main strike only */
				if (m == 0 || strike_type < 3)
				{
					if (CalibInst->dStrike[m] > CalibInst->dFwdCash + param->max_std * std)
					{
						if (m == 0)
						{
							CalibInst->dStrike[m] = CalibInst->dFwdCash + param->max_std * std;
						}
						else
						{
							/* Apply the original number of std */
							CalibInst->dStrike[m] = CalibInst->dFwdCash + (CalibInst->dStrike[m] - dMainStrike);
						}
					}
					else if (CalibInst->dStrike[m] < CalibInst->dFwdCash - param->max_std * std)
					{
						if (m == 0)
						{
							CalibInst->dStrike[m] = CalibInst->dFwdCash - param->max_std * std;
						}
						else
						{
							/* Apply the original number of std */
							CalibInst->dStrike[m] = CalibInst->dFwdCash + (CalibInst->dStrike[m] - dMainStrike);
						}
					}
				}

				/*	Make sure strikes are positive (actually more than 1bp)
						otherwise use ATM	*/
				if (CalibInst->dStrike[m] < 1.0e-04)
				{
					CalibInst->dStrike[m] = CalibInst->dFwdCash;
				}

				err = get_cash_vol (
					vol_curve_name, 
					CalibCpnSchedule->lCpnDate[j],
					CalibCpnSchedule->lCpnTheoDate[l],
					CalibInst->dStrike[m], 
					0, 
					vol_ref_rate_name, 
					&(CalibInst->dNormalVol[m]), 
					&power);

				if (err)
				{
					goto FREE_RETURN;
				}

				err = shift_the_vol(CalibInst->dFwdCash, 
									CalibInst->dStrike[m], 								
									CalibInst->dExeTime, 
									CalibInst->dNormalVol[m],
									power,
									vol_shift,
									vol_type,
									shift_type,
									&(CalibInst->dNormalVol[m]),
									&power);

				if (err) goto FREE_RETURN;

				if (power > 0.5)
				{
					CalibInst->dPrice[m] = srt_f_optblksch(
						CalibInst->dFwdCash,
						CalibInst->dStrike[m],
						CalibInst->dNormalVol[m],
						CalibInst->dExeTime,
						CalibInst->dLevel,
						SRT_PUT,
						PREMIUM);

					err = srt_f_optimpvol(	CalibInst->dPrice[m],
											CalibInst->dFwdCash,
											CalibInst->dStrike[m],
											CalibInst->dExeTime,
											CalibInst->dLevel,
											SRT_PUT,
											SRT_NORMAL,
											&(CalibInst->dNormalVol[m]));
				}
				else
				{
					CalibInst->dPrice[m] = srt_f_optblknrm(
						CalibInst->dFwdCash,
						CalibInst->dStrike[m],
						CalibInst->dNormalVol[m],
						CalibInst->dExeTime,
						CalibInst->dLevel,
						SRT_PUT,
						PREMIUM);					
				}

				CalibInst->dVega[m] = srt_f_optblknrm(
						CalibInst->dFwdCash,
						CalibInst->dStrike[m],
						CalibInst->dNormalVol[m],
						CalibInst->dExeTime,
						CalibInst->dLevel,
						SRT_PUT,
						VEGA) * 100;
			}

			CalibExeSchedule->dStrikes[i] = CalibInst->dStrike[0];
			CalibExeSchedule->dPrices[i] = CalibInst->dPrice[0];

			if (param->vega_prec && CalibExeSchedule->dWeights)
			{
				CalibExeSchedule->dWeights[i] = 1.0 / CalibInst->dVega[0];
			}

			/* Fill the Coupons */
			CalibInst->dCpnValues[0] = -CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn];

			for (iLoopCoupon=1; iLoopCoupon<CalibInst->iNbCoupon; iLoopCoupon++)
			{
				CalibInst->dCpnValues[iLoopCoupon] = CalibInst->dStrike[0] * CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn+iLoopCoupon] * CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn+iLoopCoupon];
			}

			CalibInst->dCpnValues[CalibInst->iNbCoupon-1] += CalibCpnSchedule->dCpnDf[CalibInst->iEndCpn];
		}
		else
		{
			/* Calculate the spread */
			if (err = swp_f_ForwardRate(
						CalibCpnSchedule->lCpnDate[j],
						CalibCpnSchedule->lCpnTheoDate[l],
						instr_freq,
						instr_basis,
						yc_name,
						ref_rate_name,
						&swp_rte))
			{
				goto FREE_RETURN;
			}

			CalibInst->dSpread = swp_rte - CalibInst->dFwdCash;
		}
	}

FREE_RETURN:

	return err;
}

/*	Function to get and convert the volatility */
Err	GetCashVolAndConvert(
				GET_CASH_VOL_FUNC	get_cash_vol,
				char				*vol_curve_name,	
				double				start_date, 
				double				end_date,
				double				mat,
				double				cash_fwd,
				double				cash_strike,
				char				*ref_rate_name,
				SrtDiffusionType	vol_type,
				double				*vol)
{
Err		err = NULL;
double	power, opt;

	if (cash_strike < 1.0E-08)
	{
		*vol = 1.0E-08;
		return NULL;
	}
	
	/* first get the cash vol */
	if (mat > 1.0E-08)
	{
		err = get_cash_vol(	vol_curve_name,
							start_date, 
							end_date,
							cash_strike,
							0,
							ref_rate_name,
							vol,
							&power);

		if (err) return err;

		if (power > 0.5 && vol_type == SRT_NORMAL)
		{
			/* we need to convert from lognormal to normal */
			opt =  srt_f_optblksch(	cash_fwd,
									cash_strike,
									*vol,
									mat,
									1.0,
									SRT_CALL,
									SRT_PREMIUM);

			if (mat >= YEARS_IN_DAY && opt > 1.0E-08)
			{
				err = srt_f_optimpvol(	opt,
										cash_fwd,
										cash_strike,
										mat,
										1.0,
										SRT_CALL,
										SRT_NORMAL,
										vol);

				if (err) return err;
			}
			else
			{
				*vol = 1.0E-08;
			}
		}
		else if (power < 0.5 && vol_type == SRT_LOGNORMAL)
		{
			/* we need to convert from normal to lognormal */
			opt =  srt_f_optblknrm(	cash_fwd,
									cash_strike,
									*vol,
									mat,
									1.0,
									SRT_CALL,
									SRT_PREMIUM);

			if (mat >= YEARS_IN_DAY && opt > 1.0E-08)
			{
				err = srt_f_optimpvol(	opt,
										cash_fwd,
										cash_strike,
										mat,
										1.0,
										SRT_CALL,
										SRT_LOGNORMAL,
										vol);

				if (err) return err;
			}
			else
			{
				*vol = 1.0E-08;
			}
		}
	}
	else
	{
		*vol = 1.0E-08;
	}

	return err;
}