/* ==========================================================================
   FILE_NAME:	Fx3FMc.c

   PURPOSE:		Monte Carlo Three Factors

   DATE:		08/28/00
   ========================================================================== */

#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "MCEBOptimisation.h"
#include "math.h"

/*	Main function */
/*	------------- */
Err	 mc_main_5dfx(	
					/*	Time data */
					long		nb_paths,
					int			nb_col,
					double		*time,
					double		*date,
					long		nb_dates,
					double		*dom_ifr,
					double		*dom_fwd1,
					double		*dom_fwd2,
					double		*dom_exp1,
					double		*dom_exp2,
					double		*dom_phi1,
					double		*dom_phi2,
					double		*dom_phi12,
					double		*dom_gam1_fwd,
					double		*dom_gam2_fwd,
					double		*dom_bond_pay,
					double		*dom_gam1_pay,
					double		*dom_gam2_pay,
					double		*for_ifr,
					double		*for_fwd1,
					double		*for_fwd2,
					double		*for_exp1,
					double		*for_exp2,			
					double		*for_phi1,
					double		*for_phi2,
					double		*for_phi12,
					double		*for_gam1_fwd,
					double		*for_gam2_fwd,
					double		*fx_fwd,
					double		***covar,
					/*	Product data */
					void		**func_parm_tab, 
					/*	Model data */
					double		dom_lam, 					
					double		dom_alpha,
					double		dom_gamma,					
					double		for_lam,
					double		for_alpha,
					double		for_gamma,					
					double		**correl,
					/*	Market data */
					double		spot_fx,
					char		*dom_yc,
					char		*for_yc,
					/* do PECS adjustment */
					int			do_pecs,					
					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					/*	Initialisation function to be called at the beggining of each path 
							or NULL if none */
					void (*init_func)(),
					/*	Payoff function */
					Err (*payoff_func)(
										/* Event */
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										/* Market data */
										double	spot_fx,
										double	R1D,
										double	R2D,
										double	R1F,
										double	R2F,
										double	Z,
										/* Results */
										int		nb_col,
										double	*res,
										int		*stop_path),
					/*	Results */
					double		**res)
{
int		stop_path;
long	i, j, k, l, m;
double	df;
double	R1D, R2D, R1F, R2F, Z;

double	*temp				= NULL,
		*res_evt			= NULL,
		*sum_price			= NULL,
		*sum_2price			= NULL,
		**matrix			= NULL,
		**chol				= NULL,
		***save_values		= NULL;

double	*matrixi;

double	tempc[5];
double	sum, logspot;
int		optim_today;
clock_t	time1, time2;

Err		err					= NULL;
	
	time1 = clock();
	
	logspot = log(spot_fx);

	/*	nb_paths has to be odd */
	nb_paths = 2 * ((long) (nb_paths / 2)) + 1;

	temp = dvector (0, nb_col - 1);
	res_evt = dvector (0, nb_col - 1);
	sum_price = dvector (0, nb_col - 1);
	sum_2price = dvector (0, nb_col - 1);
	matrix = dmatrix (0, nb_paths - 1, 0, 5 * (nb_dates - 1) - 1);		
	chol = dmatrix (0, 4, 0, 4);

	if (!temp || !res_evt || !sum_price || !sum_2price || !matrix || !covar)
	{
		err = "Memory allocation failure in mc_main_5dfx";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		err = mceb_allocate_savevalues_for_GRFN(nb_paths,
												nb_dates,
												params,
												&save_values);

		if (err) goto FREE_RETURN;
	}

	/* Initialisation */
	memset (sum_price, 0, nb_col * sizeof (double));
	memset (sum_2price, 0, nb_col * sizeof (double));	

	/* fill the Brownian matrix */
	err = balsam_generation (nb_paths, 5 * (nb_dates - 1), matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 1 -BalSam generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	time1 = clock();
	if (do_pecs)
	{
		adjust_emp_covar (matrix, nb_paths, 5 * (nb_dates - 1));
	}
	time2 = clock();
	smessage ("Phase 2 -PECS adjustment, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	time1 = clock();

	/* Adjust the  matrix to have the right covariance matrix */
	m = 0;
	for (j=1; j<nb_dates; j++)
	{	
		/* Find the Cholesky of the covar */
		nr_choldc (5, covar[j], chol);		

		for (i=0; i<nb_paths; i++)
		{
			for (l=0; l<5; l++)
			{
				sum = 0.0;
				for (k=0; k<5; k++)
				{
					sum += chol[l][k] * matrix[i][m + k];
				}

				tempc[l] = sum;
			}

			for (l=0; l<5; l++)
			{
				matrix[i][m + l] = tempc[l];
			}
		}

		m += 5;
	}

	time2 = clock();
	smessage ("Phase 3 -correlation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	time1 = clock();
	/* Do the simulation */
	for (i=0; i<nb_paths; i++)
	{
		R1D = R2D = R1F = R2F = 0.0;
		Z = logspot;
		matrixi = matrix[i];
		
		memset (temp, 0, nb_col * sizeof (double));
		stop_path = 0;
		
		if (init_func)
		{
			init_func();
		}

		/*	If today is an event date */
		if (func_parm_tab[0])
		{

			err = payoff_func(	date[0],
								time[0],
								func_parm_tab[0],
								spot_fx,																
								R1D,
								R2D,
								R1F,
								R2F,
								exp(Z),
								nb_col,
								res_evt,
								&stop_path);
			
			if (err)
			{
				goto FREE_RETURN;
			}

			df = dom_bond_pay[0];
			
			for (k=0; k<nb_col; k++)
			{
				temp[k] += res_evt[k] / df;
			}

			if (do_optimisation)
			{
				mceb_fill_savevalues_from_GRFN(	save_values[0],
												res_evt,
												i,
												df,
												params);
			}
		}
		
		m = 0;

		for (j=1; stop_path == 0 && j<nb_dates; j++)
		{			
			/* reconstruction formulae to calculate Fwd Bond */
			Z += fx_fwd[j] 
				- for_gam1_fwd[j-1] * R1F - for_gam2_fwd[j-1] * R2F
				+ dom_gam1_fwd[j-1] * R1D + dom_gam2_fwd[j-1] * R2D
				+ matrixi[m+4];

			R1D = R1D * dom_exp1[j] + dom_fwd1[j] + matrixi[m];
			R2D = R2D * dom_exp2[j] + dom_fwd2[j] + matrixi[m+1];

			R1F = R1F * for_exp1[j] + for_fwd1[j] + matrixi[m+2];
			R2F = R2F * for_exp2[j] + for_fwd2[j] + matrixi[m+3];

			df = dom_bond_pay[j] * exp(- dom_gam1_pay[j] * R1D - dom_gam2_pay[j] * R2D);
			
			err = payoff_func(	date[j],
								time[j],
								func_parm_tab[j],
								spot_fx,
								R1D,
								R2D,
								R1F,
								R2F,
								exp(Z),
								nb_col,
								res_evt,
								&stop_path);

			
			if (err)
			{
				goto FREE_RETURN;
			}

			/* Discount the payoff */
			for (k=0; k<nb_col; k++)
			{
				temp[k] += res_evt[k] / df;
			}
			
			if (do_optimisation)
			{
				mceb_fill_savevalues_from_GRFN(	save_values[j],
												res_evt,
												i,
												df,
												params);
			}

			m += 5;
		}

		for (k=0; k<nb_col; k++)
		{
			sum_price[k] += temp[k] / nb_paths;
			sum_2price[k] += temp[k] * temp[k] / nb_paths;
		}

		if (do_optimisation && params->iKnockInCol)
		{
			/* we recopy in the col pay the pv of the column */
			for (j=0; j<nb_dates; j++)
			{
				if (optimise[j])
				{
					save_values[j][params->iNbIndex][i] = temp[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
				}
			}
		}
	}
	
	/* Save results */
	for (k=0; k<nb_col; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = (sum_2price[k] - sum_price[k] * sum_price[k]) / nb_paths;

		if (fabs(res[k][1]) < 1.0E-10)
		{
			res[k][1] = 0.0;
		}
		else
		{
			res[k][1] = sqrt(res[k][1]);
		}
	}

	time2 = clock();
	smessage ("Phase 4 -evaluation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		time1 = clock();	

		/* free the useless memory */
		if (matrix)
		{
			free_dmatrix (matrix, 0, nb_paths - 1, 0, 3 * (nb_dates-1) - 1);
			matrix = NULL;
		}
		
		if (chol)
		{
			free_dmatrix(chol, 0, 4, 0, 4);
			chol = NULL;
		}
		
		if (time[0] < 1.0E-08 && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}

		err = find_and_optimise_boundary(	save_values,
											nb_dates,
											nb_paths,
											optimise,
											params,
											&(res[nb_col][0]),
											&(res[nb_col][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[nb_col][0] < 0.0)
				{
					res[nb_col][0] = 0.0;
					res[nb_col][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[nb_col][0])
				{
					res[nb_col][0] = save_values[0][params->iNbIndex][0];
					res[nb_col][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}		

		time2 = clock();
		smessage ("Phase 5 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}

FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, nb_col - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, nb_col - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, nb_col - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, nb_col - 1);
	}

	if (matrix)
	{
		free_dmatrix (matrix, 0, nb_paths - 1, 0, 5 * (nb_dates-1) - 1);
	}	

	if (chol)
	{
		free_dmatrix(chol, 0, 4, 0, 4);
	}	

	mceb_free_savevalues_for_GRFN(save_values, nb_paths, nb_dates, params);

	return err;
}

/*	Main Test function */
/*	------------- */
Err	 mc_main_5dfx_test(	
					/*	Time data */
					long		nb_paths,
					int			nb_col,
					double		*time,
					double		*date,
					long		nb_dates,
					double		*dom_ifr,		/*	Distributions */
					double		*dom_fwd,
					double		*dom_std,
					double		*dom_phi1,
					double		*dom_phi2,
					double		*dom_phi12,
					double		*dom_beta,
					double		*dom_bond_pay,
					double		*dom_beta_pay,
					double		*for_ifr,
					double		*for_fwd,
					double		*for_std,
					double		*for_phi1,
					double		*for_phi2,
					double		*for_phi12,
					double		*for_beta,
					double		*fx_fwd,
					double		*fx_std,
					double		*dom_for_cov,
					double		*dom_fx_cov,
					double		*for_fx_cov,
					/*	Product data */
					void		**func_parm_tab, 
					/*	Model data */
					double		dom_lam, 					
					double		dom_alpha,
					double		dom_gamma,
					double		dom_rho,
					double		for_lam,
					double		for_alpha,
					double		for_gamma,
					double		for_rho,
					double		**correl,
					/*	Market data */
					double		spot_fx,
					char		*dom_yc,
					char		*for_yc,
					/* do PECS adjustment */
					int			do_pecs,
					/*	Initialisation function to be called at the beggining of each path 
							or NULL if none */
					void (*init_func)(),
					/*	Payoff function */
					Err (*payoff_func)(
										/* Event */
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										/* Market data */
										double	spot_fx,
										double	R1D,
										double	R2D,
										double	R1F,
										double	R2F,
										double	Z,
										/* Results */
										int		nb_col,
										double	*res,
										int		*stop_path),
					/*	Results */
					double		**res)
{
int		stop_path;
long	i, j, k, l;
double	t1, t2;
double	df;
double	R1D, R2D, R1F, R2F, RD, RF, Z;
clock_t	time1, time2;
double	*temp				= NULL;
double	*res_evt			= NULL;
double	*sum_price			= NULL;
double	*sum_2price			= NULL;
double	**matrix			= NULL,
		**covar				= NULL,
		**chol				= NULL;

double	tempc[5];
double	sum, dt, logspot;
double	dom_lam2, for_lam2;

Err		err					= NULL;
	
	time1 = clock();
	
	dom_lam2 = dom_lam + dom_gamma;
	for_lam2 = for_lam + for_gamma;
	logspot = log(spot_fx);

	/* nb_paths has to be odd */
	nb_paths = 2 * ((long) (nb_paths / 2)) + 1;

	temp = dvector (0, nb_col - 1);
	res_evt = dvector (0, nb_col - 1);
	sum_price = dvector (0, nb_col - 1);
	sum_2price = dvector (0, nb_col - 1);
	matrix = dmatrix (0, nb_paths - 1, 0, 5 * (nb_dates - 1) - 1);	
	covar = dmatrix(0, 4, 0, 4);
	chol = dmatrix (0, 4, 0, 4);

	if (!temp || !res_evt || !sum_price || !sum_2price || !matrix || !covar)
	{
		err = "Memory allocation failure in mc_main_5dfx";
		goto FREE_RETURN;
	}	

	for (k=0; k<nb_col; k++)
	{
		sum_price[k] = sum_2price[k] = 0.0;
	}

	/* fill the Brownian matrix */
	err = balsam_generation (nb_paths, 5 * (nb_dates - 1), matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 1 -BalSam generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	if (do_pecs)
	{
		adjust_emp_covar (matrix, nb_paths, 5 * (nb_dates - 1));
	}
	time2 = clock();
	smessage ("Phase 2 -PECS adjustment, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	/* for each date */
	for (j=1; j<nb_dates; j++)
	{
		/* fill the covar matrix */
		covar[0][0] = dom_std[j] * dom_std[j];
		covar[0][1] = covar[1][0] = correl[0][1] * dom_std[j] * dom_alpha * dom_std[j];
		covar[0][2] = covar[2][0] = correl[0][2] * dom_std[j] * for_std[j];
		covar[0][3] = covar[3][0] = correl[0][3] * dom_std[j] * for_alpha * for_std[j];
		covar[0][4] = covar[4][0] = correl[0][4] * dom_std[j] * fx_std[j];
		covar[1][1] = dom_alpha * dom_std[j] * dom_alpha * dom_std[j];
		covar[1][2] = covar[2][1] = correl[1][2] * dom_alpha * dom_std[j] * for_std[j];
		covar[1][3] = covar[3][1] = correl[1][3] * dom_alpha * dom_std[j] * for_alpha * for_std[j];
		covar[1][4] = covar[4][1] = correl[1][4] * dom_alpha * dom_std[j] * fx_std[j];
		covar[2][2] = for_std[j] * for_std[j];
		covar[2][3] = covar[3][2] = correl[2][3] * for_std[j] * for_alpha * for_std[j];
		covar[2][4] = covar[4][2] = correl[2][4] * for_std[j] * fx_std[j];
		covar[3][3] = for_alpha * for_std[j] * for_alpha * for_std[j];
		covar[3][4] = covar[4][3] = correl[3][4] * for_alpha * for_std[j] * fx_std[j];
		covar[4][4] = fx_std[j] * fx_std[j];

		nr_choldc (5, covar, chol);

		for (i=0; i<nb_paths; i++)
		{
			for (l=0; l<5; l++)
			{
				sum = 0.0;
				for (k=0; k<5; k++)
				{
					sum += chol[l][k] * matrix[i][(j - 1) * 5 + k];
				}

				tempc[l] = sum;
			}

			for (l=0; l<5; l++)
			{
				matrix[i][(j - 1) * 5 + l] = tempc[l];
			}
		}
	}

	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 3 -correlation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	for (i=0; i<nb_paths; i++)
	{
		t1 = 0;
		R1D = R2D = R1F = R2F = 0.0;
		Z = logspot;
		
		memset (temp, 0, nb_col * sizeof (double));
		stop_path = 0;
		
		if (init_func)
		{
			init_func();
		}

		df = 0.0;

		/*	If today is an event date */
		if (func_parm_tab[0])
		{			
			
			err = payoff_func(	date[0],
								time[0],
								func_parm_tab[0],
								spot_fx,																
								R1D,
								R2D,
								R1F,
								R2F,
								exp(Z),
								nb_col,
								res_evt,
								&stop_path);
			
			if (err)
			{
				goto FREE_RETURN;
			}		
			
			for (k=0; k<nb_col; k++)
			{
				temp[k] += res_evt[k] * exp(-df);
			}
		}

		l = 0;

		for (j=1; stop_path == 0 && j<nb_dates; j++)
		{
			t2 = time[j];
			dt = t2 - t1;

			RD = dom_ifr[j-1] + R1D + R2D;
			RF = for_ifr[j-1] + R1F + R2F;

			/* reconstruction formulae to calculate Fwd Bond */
			Z += (RD - RF) * dt + fx_fwd[j] + matrix[i][l+4];

			R1D += (dom_phi1[j-1] + dom_phi12[j-1] - dom_lam * R1D) * dt + matrix[i][l];
			R2D += (dom_phi2[j-1] + dom_phi12[j-1] - dom_lam2 * R2D) * dt + matrix[i][l+1];

			R1F +=	(for_phi1[j-1] + for_phi12[j-1] - for_lam * R1F ) * dt
					- correl[2][4] * for_std[j] * fx_std[j]
					+ matrix[i][l+2];

			R2F +=	(for_phi2[j-1] + for_phi12[j-1] - for_lam2 * R2F) * dt
					- correl[3][4] * for_alpha * for_std[j] * fx_std[j]
					+ matrix[i][l+3];
			
			df += RD * dt;

			if (func_parm_tab[j])
			{				
				err = payoff_func(	date[j],
									t2,
									func_parm_tab[j],
									spot_fx,
									R1D,
									R2D,
									R1F,
									R2F,
									exp(Z),
									nb_col,
									res_evt,
									&stop_path);
			
				if (err)
				{
					goto FREE_RETURN;
				}

				for (k=0; k<nb_col; k++)
				{
					temp[k] += res_evt[k] * exp(-df);
				}
			}

			t1 = t2;
			l += 5;
		}

		for (k=0; k<nb_col; k++)
		{
			sum_price[k] += temp[k] / nb_paths;
			sum_2price[k] += temp[k] * temp[k] / nb_paths;
		}
	}
	
	for (k=0; k<nb_col; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = sqrt ((sum_2price[k] - sum_price[k] * sum_price[k]) / nb_paths);
	}

	time2 = clock();
	smessage ("Phase 4 -evaluation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, nb_col - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, nb_col - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, nb_col - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, nb_col - 1);
	}

	if (matrix)
	{
		free_dmatrix (matrix, 0, nb_paths - 1, 0, 5 * (nb_dates-1) - 1);
	}	

	if (covar)
	{
		free_dmatrix(covar, 0, 4, 0, 4);
	}

	if (chol)
	{
		free_dmatrix(chol, 0, 4, 0, 4);
	}	

	return err;
}

