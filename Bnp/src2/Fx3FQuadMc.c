/* ==========================================================================
   FILE_NAME:	Fx3FQuadMc.c

   PURPOSE:		Monte Carlo Three Factors

   DATE:		11/20/00
   ========================================================================== */

#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "math.h"

/*	Main function */
/*	------------- */
Err	 mcQuad_main_3dfx(	
					/*	Time data */
					long		npaths,
					long		nsteps,
					int			num_col,
					double		*time,
					double		*date,
					double		*dom_ifr,		
					double		*dom_std,
					double		*dom_phi,
					double		*for_ifr,
					double		*for_std,
					double		*for_phi,
					double		*fx_fwd,
					double		*fx_std,
					double		alpha,
					double		beta,
					double		gamma,
					double		sig0,

					/*	Product data */
					void		**func_parm_tab,
					int			*has_evt,
					/*	Model data */
					double		dom_lam, 
					double		for_lam, 
					double		*corr_dom_for,
					double		*corr_dom_fx,
					double		*corr_for_fx,
					/*	Market data */
					double		spot_fx,
					char		*dom_yc,
					char		*for_yc,
					/* do PECS adjustment */
					int			do_pecs,
					/*	Payoff function */
					Err (*payoff_func)(
										/* Event */
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										/* Market data */
										double	spot_fx,
										void	*dom_yc,
										double	dom_lam,
										double	dom_phi,
										void	*for_yc,
										double	for_lam,
										double	for_phi,
										double	Xdom,
										double	Yfor,
										double	Zfx,
										/* Results */
										int		num_col,
										double	*res,
										int		*stop_path),
					/*	Results */
					double		**res)
{
int		stop_path, num_evt;
long	i, j, k;
double	t1, t2;
double	df, dt, zc;
double	X, Y, Z, Z2, logSpot;
double	stdev1;

clock_t	time1, time2;
double	*temp				= NULL;
double	*res_evt			= NULL;
double	*sum_price			= NULL;
double	*sum_2price			= NULL;
double	**matrix			= NULL;
Err		err					= NULL;
	
	time1 = clock();
	logSpot = log(spot_fx);
	
	/* npaths has to be odd */
	npaths = 2 * ((long) (npaths / 2)) + 1;

	temp = dvector (0, num_col - 1);
	res_evt = dvector (0, num_col - 1);
	sum_price = dvector (0, num_col - 1);
	sum_2price = dvector (0, num_col - 1);
	matrix = dmatrix (0, npaths - 1, 0, 3 * (nsteps - 1) - 1);

	if (!temp || !res_evt || !sum_price || !sum_2price || !matrix)
	{
		err = "Memory allocation failure in mcBeta_main_3dfx";
		goto FREE_RETURN;
	}

	for (k=0; k<num_col; k++)
	{
		sum_price[k] = sum_2price[k] = 0.0;
	}

	/* fill the Brownian matrix */
	err = balsam_generation (npaths, 3 * (nsteps - 1), matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 2 -BalSam generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	if (do_pecs)
	{
		adjust_emp_covar (matrix, npaths, 3 * (nsteps - 1));
	}
	time2 = clock();
	smessage ("Phase 3 -PECS adjustment, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	err = correlate_variable_corr (&(dom_std[1]), &(for_std[1]), &(fx_std[1]),
							 &(corr_dom_for[1]), &(corr_dom_fx[1]), &(corr_for_fx[1]),
							 npaths, nsteps-1,
							 matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 4 -correlation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	time1 = clock();

	for (i=0; i<npaths; i++)
	{
		t1 = 0;
		X = Y = 0;
		Z = logSpot;
		
		memset (temp, 0, num_col * sizeof (double));
				
		stop_path = 0;
		num_evt = 0;

		/*	If today is an event date */
		if (func_parm_tab[0] && has_evt[0])
		{		
			
			err = payoff_func(	date[0],
								time[0],
								func_parm_tab[0],
								spot_fx,
								dom_yc,
								dom_lam,
								dom_phi[0],
								for_yc,
								for_lam,
								for_phi[0],
								X,
								Y,
								Z,
								num_col,
								res_evt,
								&stop_path);
			num_evt += 1;
			
			if (err)
			{
				goto FREE_RETURN;
			}		
			
			for (k=0; k<num_col; k++)
			{
				temp[k] += res_evt[k];
			}
		}

		zc = 0.0;

		for (j=1; stop_path == 0 && j<nsteps; j++)
		{
			t2 = time[j];
			stdev1 = sig0 * sqrt(t1);
			dt = t2 - t1;

			/* random X ,Y and Z */

			zc += (X + dom_ifr[j-1]) * dt; 			

			Z2 = fx_fwd[j-1] / exp(Z);
			
			if (Z2 > 1.0E-08)
			{
				if (j > 1)
				{
					Z2 = fx_std[j] * (	alpha * stdev1 * Z2 
										+ (1.0 - Z2) * (beta
													  + gamma / stdev1 / Z2 * (1.0 - beta)));
				}
				else
				{
					Z2 = alpha * stdev1;
				}

				Z +=	(X + dom_ifr[j-1] - Y - for_ifr[j-1]) * dt
					-0.5 * Z2 * Z2
					+ Z2 / fx_std[j] * matrix[i][3*(j-1) + 2];
			}
							
			X += (dom_phi[j-1] - dom_lam * X) * dt + matrix[i][3*(j-1)];

			Y += (for_phi[j-1] - for_lam * Y) * dt 
				- corr_for_fx[j] * for_std[j] * Z2 + matrix[i][3*(j-1) + 1];			
  
			if (has_evt[j])
			{
				err = payoff_func(	date[j],
									t2,
									func_parm_tab[num_evt],
									spot_fx,
									dom_yc,
									dom_lam,
									dom_phi[j],
									for_yc,
									for_lam,
									for_phi[j],
									X,
									Y,
									Z,
									num_col,
									res_evt,
									&stop_path);
				num_evt += 1;

				if (err)
				{
					goto FREE_RETURN;
				}

				df = exp(-zc);
				for (k=0; k<num_col; k++)
				{
					temp[k] += res_evt[k] * df;
				}

			}
			
			

			t1 = t2;
		}

		for (k=0; k<num_col; k++)
		{
			sum_price[k] += temp[k] / npaths;
			sum_2price[k] += temp[k] * temp[k] / npaths;
		}
	}
	
	for (k=0; k<num_col; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = sqrt ((sum_2price[k] - sum_price[k] * sum_price[k]) / npaths);
	}

	time2 = clock();
	smessage ("Phase 5 -evaluation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, num_col - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, num_col - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, num_col - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, num_col - 1);
	}

	if (matrix)
	{
		free_dmatrix (matrix, 0, npaths - 1, 0, 3 * (nsteps-1) - 1);
	}
	
	return err;
}