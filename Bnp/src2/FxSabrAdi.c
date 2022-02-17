/* ==========================================================================
   FILE_NAME:	LGM2FSabr.c

   PURPOSE:		ADI implementation of the Sabr Fx model.
				Discretisation is ADI.

   DATE:		03/16/01
   
   AUTHOR:		L.C.
   ========================================================================== */

#include "FxSabrAdi.h"
#include "opfnctns.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "FxSabrGrfn.h"
#include "math.h"

#define	NSTD_FX		5.0
#define	NSTD_VOL	5.0

#define MINNODE		25
#define	MAXNODE2	50


Err FxSabrPrecalculations(
							/*	Time data		*/
							int			nstp,
							double		*time,
							
							/*	Model data		*/
							double		sig0,
							double		*drift,
							double		alpha,
							double		beta,
							double		rho,
							double		lambda,

							/*	Market data */
							double		spot_fx,

							/*	Output			*/
							double		*exp_h,
							double		*expect_z,
							double		*drift_z,
							double		*std1)
{
double	spotBeta, z0, const1, const2, const3, const4, const5, x;
double	tmid, temp, exp_z, varxf, sig, dt;
int		i;

	spotBeta = pow (spot_fx, beta - 1.0);
	z0 = rho * alpha / (1.0 - beta) * (1.0 / spotBeta - 1) - sig0;

	const1 = -0.5 * alpha * beta * rho * sig0 * sig0 * spotBeta;
	const2 = -alpha * rho * (1.0 / spotBeta - 1.0) / (1.0 - beta);

	const5 = 2.0 * lambda - alpha * alpha;
	const3 = const1 * 2.0 * lambda / const5;
	const4 = -const1 * alpha * alpha / const5;

	/* Initialisation */
	exp_h[0] = 0.0;
	expect_z[0] = z0;

	temp = 0.0;
	exp_z = z0;
	varxf = 0.0;
	sig = 1.0;

	for (i=0; i<nstp-1; i++)
	{
		dt = (time[i+1] - time[i]);

		/* First update varxf */
		if (fabs(drift[i]) > 1.0E-10)
		{
			varxf += sig * sig * (exp(2.0 * drift[i] * dt) - 1.0) / (2.0 * drift[i]);
		}
		else
		{
			varxf += sig * sig * dt;
		}

		/* From ti to tmid*/
		dt /= 2.0;
		tmid = time[i] + dt;

		/* Calculation of the expectation at middle point */
		
		if (fabs(drift[i]) > 1.0E-10)
		{
			exp_z += const3 * sig * (exp(drift[i] * dt) - 1.0) / drift[i];

			x = drift[i] - const5;

			if (fabs(x) > 1.0E-10)
			{
				exp_z += const4 * sig * exp(-time[i] * drift[i]) 
								* (exp(x * tmid) - exp(x * time[i])) / x;
			}
			else
			{
				exp_z += const4 * sig * exp(-time[i] * drift[i]) * dt;
			}

			exp_z += const2 / sig * (1.0 - exp(-drift[i] * dt));
		}
		else
		{
			exp_z += const3 * sig * dt;

			x = -const5;

			if (fabs(x) > 1.0E-10)
			{
				exp_z += const4 * sig * (exp(x * tmid) - exp(x * time[i])) / x;
			}
			else
			{
				exp_z += const4 * sig * dt;
			}	
		}
		
		temp += drift[i] * dt;
		sig = exp(temp);

		expect_z[i] = exp_z;
		drift_z[i] = const3 * sig + const4 * sig * exp(-const5 * (time[i] + dt)) + const2 / sig * drift[i];
		exp_h[i] = sig;

		/* Update exp_z, temp and sig at time[i+1] */	
		
		if (fabs(drift[i]) > 1.0E-10)
		{
			exp_z += const3 * sig * (exp(drift[i] * dt) - 1.0) / drift[i];

			x = drift[i] - const5;

			if (fabs(x) > 1.0E-10)
			{
				exp_z += const4 * sig * exp(-tmid * drift[i]) 
								* (exp(x * time[i+1]) - exp(x * tmid)) / x;
			}
			else
			{
				exp_z += const4 * sig * exp(-tmid * drift[i]) * dt;
			}

			exp_z += const2 / sig * (1.0 - exp(-drift[i] * dt));
		}
		else
		{
			exp_z += const3 * sig * dt;

			x = -const5;

			if (fabs(x) > 1.0E-10)
			{
				exp_z += const4 * sig * (exp(x * time[i+1]) - exp(x * tmid)) / x;
			}
			else
			{
				exp_z += const4 * sig * dt;
			}	
		}
		
		temp += drift[i] * dt;
		sig = exp(temp);
	}

	*std1 = sig0 * sqrt(varxf) * spotBeta;	

	return NULL;
}

Err FxSabrPrecalculations2(
							/*	Time data		*/
							int			nstp,
							double		*time,
							
							/*	Model data		*/
							double		sig0,
							double		*drift,
							double		alpha,
							double		beta,
							double		rho,
							double		lambda,

							/*	Market data */
							double		spot_fx,

							/*	Output			*/
							double		*exp_h,
							double		*expect_z,
							double		*drift_z,
							double		*std1)
{
double	spotBeta, z0, const1;
double	temp, varxf, sig, dt;
int		i;

	spotBeta = pow (spot_fx, beta - 1.0);
	z0 = rho * alpha / (1.0 - beta) * (1.0 / spotBeta - 1) - sig0;

	const1 = rho * alpha * (1.0 / spotBeta - 1.0) / (1.0 - beta);

	/* Initialisation */
	exp_h[0] = 0.0;
	expect_z[0] = z0;

	varxf = 0.0;
	sig = 1.0;

	temp = 0.0;

	for (i=0; i<nstp-1; i++)
	{
		dt = (time[i+1] - time[i]);

		/* First update varxf */
		if (fabs(drift[i]) > 1.0E-10)
		{
			varxf += sig * sig * (exp(2.0 * drift[i] * dt) - 1.0) / (2.0 * drift[i]);
		}
		else
		{
			varxf += sig * sig * dt;
		}

		/* Then the moment */

		/*
		dt /= 2.0;		
		*/

		temp += drift[i] * dt;
		sig = exp(temp);

		exp_h[i] = sig;
		expect_z[i] = const1 / sig - sig0;
		drift_z[i] = -const1 * drift[i] / sig * (2.0 * dt);

		/*
		temp += drift[i] * dt;
		sig = exp(temp);
		*/
	}

	*std1 = sig0 * sqrt(varxf) * spotBeta;	

	return NULL;
}
						  
Err FxSabr_adi(	
				/*	Time data		*/
				int			nstp,					
				double		*time,
				double		*date,

				/*	Discretisation	*/
				int			nstepfx,
				int			nstepvol,
									
				/*	Model data		*/
				double		sig0,
				double		*drift,
				double		alpha,
				double		beta,
				double		rho,
				double		lambda,

				double		floorstd,

				/*	Product data */
				void		**func_parm_tab, 
				int			*eval_evt,
				double		*bar_lvl,
				int			*bar_col,
				int			*is_bar,
				
				/*	Market data */
				double		spot_fx,	/*	The cash one */
				char		*dom_yc,
				char		*for_yc,
				
				/*	Payoff function */
				Err (*payoff_func)(/* Event */
									double	evt_date,
									double	evt_time,
									void	*func_parm, 
									
									/* Market data	*/										
									long	today,
									double	spot_fx,	/*	The cash one */
									void	*dom_yc,
									void	*for_yc,
																			
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*x,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
				/*	Result */
				int			nprod, 
				double		*res)
{
Err					err = NULL;

long				today;				
int					i, j, k, step;
int					index_x, index_z, nstepx, nstepz;
double				fux, flx, dt, dft;
double				alpha2, alpharho, alpha2rho2, alpharhobeta5;
double				temp;
double				std1, std3, spotBeta;
double				z0;
					
double				*exp_h			= NULL,					
					*expect_z		= NULL,
					*expect_z2		= NULL,
					*drift_z		= NULL,
					*x				= NULL,
					*z				= NULL,
					***values		= NULL,
					***values_p1	= NULL,
					***values_temp	= NULL,
					**mux			= NULL,
					**muz			= NULL,
					**varx			= NULL,					
					**varz			= NULL,
					*varxinit		= NULL,
					*muzinit1		= NULL,
					*muzinit2		= NULL,
					*sigBeta		= NULL,
					**r				= NULL;
					
double				*varxi,
					*varzi,
					*muzi,
					**valuesi,
					**valuesip,
					**valuesim,
					*valuesij;

double				const_sig, sigBetaij;
double				const_sigi, const_varx, const_muz1, const_muz2, const_varz, const_sigi1, const_muz21;
double				adj_lambda, adj_lambda1, adj_lambda2;
					
int					lx, ux, lz, uz;

double				bar_pres, bar_c1, bar_c2;
long				bar_index;
					
clock_t				t1, t2;

CNPDE_TEMP_2D_ADI	pdestr, *pde = NULL;

double				logx, logmax, logmin;

FILE				*stream = fopen("C:\\toto.txt", "w+");


	/*	For stability reasons */

	if (alpha < 1.0e-05)
	{
		alpha = 1.0E-05;
	}

	if (fabs(beta - 1.0) < 1.0E-08)
	{
		beta = 1.0 - 1.0E-08;
	}

	t1 = clock();

	/* Constant	calculations */
	today = (long) (date[0] + 1.0e-06);		

	spotBeta = pow (spot_fx, beta - 1.0);

	alpha2 = alpha * alpha;	
	alpharho = alpha * rho;
	alpharhobeta5 = 0.5 * alpharho * beta;
	alpha2rho2 = alpha * alpha * (1.0 - rho * rho);
	const_sig = alpharho / (1.0 - beta);

	z0 = rho * alpha / (1.0 - beta) * (1.0 / spotBeta - 1) - sig0;
				
	/*	nstep has to be a odd nuber			*/
	nstepx = ((int) (nstepfx / 2)) * 2 + 1;
	nstepz = ((int) (nstepvol / 2)) * 2 + 1;

	/*	we want at least three points in each directions */
	if (nstepx < 3)
	{
		nstepx = 3;
	}
	if (nstepz < 3)
	{
		nstepz = 3;
	}
		
	/*	Memory allocations */

	x = calloc (nstepx, sizeof(double));
	z = dvector (0, nstepz-1);

	exp_h = dvector (0, nstp-1);
	expect_z = dvector (0, nstp-1);	
	drift_z = dvector (0, nstp-1);
	
	if (!exp_h || !expect_z || !drift_z || !x || !z)
	{
		err = "Memory allocation error (1) in FxSabr_adi";
		goto FREE_RETURN;
	}

	err = FxSabrPrecalculations2(nstp,
								time,
								sig0,
								drift,
								alpha,
								beta,
								rho,
								lambda,
								spot_fx,
								exp_h,
								expect_z,
								drift_z,
								&std1);
	
	/*	Calculate the standard deviation of X and the expectaion of Z	*/	

	/*
		flx = find_beta_lim(spot_fx,
							std1 / spotBeta,
							beta,
							0.995);

		fux = find_beta_lim(spot_fx,
							std1 / spotBeta,
							beta,
							0.005);
	*/

	

	std3 = sig0 * sqrt(alpha2rho2 * time[nstp-1]);
													
	/*	Then discretise space in the orthogonal system x / z */
	
	/*	Lognormal / Normal bounds of X */

	flx = log(spot_fx) -std1 * std1 / 2.0 - NSTD_FX * std1;
	fux = flx  + 2.0 * NSTD_FX * std1;
	flx = max(exp(flx), spot_fx / 1E06);

	if (beta >= 0.5)
	{
		/* lognormal upperbound */
		if (fux > 1.0E4)
		{
			fux = spot_fx * 10000;
		}
		else
		{			
			fux = exp(fux);

			if (fux < spot_fx)
			{
				fux = spot_fx * 1.005;
			}
		}
	}
	else
	{
		/* normal upperbound */
		fux = spot_fx * (1.0 + NSTD_FX * std1); 
	}
	
	/*	Bound for the calculation of max and min drift and var */

	if (floorstd == 5.0)
	{
		logmax = log(fux);
		logmin = log(flx);
	}
	else
	{
		logmax = log(spot_fx) -std1 * std1 / 2.0 + floorstd * std1;
		logmin = log(spot_fx) -std1 * std1 / 2.0 - floorstd * std1;	
	}
	
	/*	Discretisation of x and z */

	if (bar_lvl)
	{
		disc_linleft_logright_bar(&x, &nstepx, spot_fx, flx, fux, std1, NSTD_FX, bar_lvl, nstp, &bar_pres, &index_x);
	}
	else
	{
		if (beta >= 0.5)
		{
			disc_linleft_logright_center(x, nstepx, spot_fx, flx, fux, std1, NSTD_FX, &index_x);
		}
		else
		{
			disc_linleft_linright_center(x, nstepx, spot_fx, flx, fux, std1, NSTD_FX, &index_x);
		}
	}

	/*
	disc_normal_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3, NSTD_VOL, &index_z);
	*/

	disc_linleft_linright_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3, NSTD_VOL, &index_z);

	/*	Precalculate variances and expectations				*/	

	/* first allocate memory */	
	values = f3tensor (0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor (0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	mux = dmatrix (0, nstepx-1, 0, nstepz-1);
	muz = dmatrix (0, nstepx-1, 0, nstepz-1);
	muzinit1 = dvector (0, nstepx-1);
	muzinit2 = dvector (0, nstepx-1);
	varz = dmatrix (0, nstepx-1, 0, nstepz-1);
	r = dmatrix (0, nstepx-1, 0, nstepz-1);

	sigBeta = dvector (0, nstepx-1);
	varx = dmatrix (0, nstepx-1, 0, nstepz-1);
	varxinit = dvector (0, nstepx-1);

	if (!values || !values_p1 || !mux || !muz || !muzinit1 || !muzinit2 || !sigBeta
		|| !varx || !varxinit || !varz  || !r)
	{
		err = "Memory allocation error (2) in FxSabr_adi";
		goto FREE_RETURN;
	}

	for (i=0; i<nstepx; i++)
	{		
		logx = min(max(log(x[i]), logmin), logmax);
		temp = exp((beta - 1.0) * logx);
		
		sigBeta[i] = alpharho * (1.0 / temp - 1.0) / (1.0 - beta);
		varxinit[i] = exp(2.0 * beta * logx);							/* z^(2.0 * beta) */
		muzinit1[i] = alpharhobeta5 * temp;
		muzinit2[i] = alpharho * (1.0 / temp - 1.0 / spotBeta) / (1.0 - beta);
	}

	/*	Final payoff valuation					*/
	if (!eval_evt[nstp-1])
	{
		err = "No event at last step in FxSabr_Adi";
		goto FREE_RETURN;
	}

	/*	Eval payoff */
	err = payoff_func (	date[nstp-1],
						time[nstp-1],
						func_parm_tab[nstp-1],
						today,
						spot_fx,
						dom_yc,
						for_yc,
						0,
						nstepx-1,
						0,
						nstepz-1,
						x,						
						nprod,
						values_p1
						);


	/*
	if (nprod > 2)
	{
		for (j=0; j<nstepz; j++)
		{
			temp = z[j];

			for (i=0; i<nstepx; i++)
			{
				values_p1[i][j][2] = temp;
			}
		}
	}
	*/
	
	if (err)
	{
		goto FREE_RETURN;
	}	

	if (bar_lvl && is_bar[nstp-1])
	{
		bar_index = 0;

		for (i=0; i<nstepx; i++)
		{
			valuesi = values[i];

			if (fabs(x[i] - bar_lvl[nstp-1]) < bar_pres)
			{
				bar_index = i;
				i = nstepx + 1;
			}
		}

		if (bar_index > 0 && bar_index < nstepx - 1)
		{
			valuesi = values_p1[bar_index];
			valuesip = values_p1[bar_index + 1];
			valuesim = values_p1[bar_index - 1];

			bar_c1 = (x[bar_index] - x[bar_index-1]) / (x[bar_index+1] - x[bar_index-1]);
			bar_c2 = 1.0 - bar_c1;

			for (j=0; j<nstepz; j++)
			{
				valuesi[j][bar_col[nstp-1]] = bar_c1 * valuesip[j][bar_col[nstp-1]] + bar_c2 * valuesim[j][bar_col[nstp-1]];
			}
		}
	}

	/*	Initialize the CNPDE_TEMP_2D		*/

	pde = &pdestr;

	num_f_pde_init_2d_adi(	pde,
							nstepx,
							nstepz,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (2) in FxSabr_Adi";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;
	
	/*	now do the backward pde					*/

	for (step=nstp-2; step>=0; step--)
	{

		dt = time[step+1] - time[step];

		adj_lambda1 = lambda / exp_h[step] * dt;
		adj_lambda2 = lambda * sig0 * dt;

		const_sigi1 = expect_z[step] * exp_h[step];
		const_muz21 = 1.0 / exp_h[step] * drift[step] * dt;
		const_varz = alpha2rho2 / exp_h[step] / exp_h[step];

		for (i=lx; i<=ux; i++)
		{			
			varxi = varx[i];
			varzi = varz[i];
			muzi = muz[i];

			const_sigi = sigBeta[i] - const_sigi1;
			const_varx = varxinit[i];
			const_muz1 = -muzinit1[i] / exp_h[step];
			const_muz2 = -muzinit2[i] * const_muz21;

			for (j=lz; j<=uz; j++)
			{
				/* reconstruction of the vol */
				sigBetaij = max(const_sigi - z[j] * exp_h[step], 1.0E-08);

				adj_lambda = adj_lambda1 * sigBetaij - adj_lambda2;
				
				sigBetaij *= sigBetaij * dt;
				
				varxi[j] = const_varx * sigBetaij;
				muzi[j] = const_muz1 * sigBetaij + const_muz2 + adj_lambda;
				varzi[j] = const_varz * sigBetaij;
			}
		}

		if (step == nstp-2)
		{
			for (i=lx; i<=ux; i++)
			{						
				for (j=lz; j<=uz; j++)
				{
					fprintf(stream, "%d	%d	%f	%f	%f	%f	%f\n", i, j, x[i], z[j], muz[i][j], varx[i][j], varz[i][j]);
				}
			}
		}

		fclose(stream);


		/*	convolve							*/

		num_f_pde_one_step_backward_2f_adi(	pde,
											nstepx,
											x,
											nstepz,
											z,
											0,
											nprod-1,
											values_p1,
											mux,
											muz,
											varx,
											varz,
											r,
											values,
											lx,
											ux,
											lz,
											uz);

		/*  Apply discounting	*/
		
		dft = swp_f_df (date[step], date[step+1], dom_yc);
				
		if (bar_lvl && is_bar[step])
		{
			bar_index = 0;

			for (i=lx; i<=ux; i++)
			{
				valuesi = values[i];

				if (!bar_index && fabs(x[i] - bar_lvl[step]) < bar_pres)
				{
					bar_index = i;
				}

				for (j=lz; j<=uz; j++)
				{
					valuesij = valuesi[j];

					for (k=0; k<nprod; k++)
					{
						valuesij[k] *= dft;		
					}				
				}
			}
		}
		else
		{
			for (i=lx; i<=ux; i++)
			{
				valuesi = values[i];

				for (j=lz; j<=uz; j++)
				{
					valuesij = valuesi[j];

					for (k=0; k<nprod; k++)
					{
						valuesij[k] *= dft;		
					}				
				}
			}
		}
		
		/*	Eval payoff */
		if (eval_evt[step])
		{
			err = payoff_func (	date[step],
								time[step],
								func_parm_tab[step],						
								today,
								spot_fx,
								dom_yc,
								for_yc,
								lx,
								ux,
								lz,
								uz,
								x,
								nprod,
								values
								);
			if (err)
			{
				goto FREE_RETURN;
			}
		}

		if (bar_lvl && is_bar[step] && bar_index > lx && bar_index < ux)
		{
			valuesi = values[bar_index];
			valuesip = values[bar_index + 1];
			valuesim = values[bar_index - 1];

			bar_c1 = (x[bar_index] - x[bar_index-1]) / (x[bar_index+1] - x[bar_index-1]);
			bar_c2 = 1.0 - bar_c1;

			for (j=lz; j<=uz; j++)
			{
				valuesi[j][bar_col[step]] = bar_c1 * valuesip[j][bar_col[step]] + bar_c2 * valuesim[j][bar_col[step]];
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;		
	}

	/* copy the result					*/
	for (k=0; k<nprod; k++)
	{
		res[k] = values_p1[index_x][index_z][k];
	}
	
	t2 = clock();
	
	smessage ("Convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:


	if (exp_h) free_dvector (exp_h, 0, nstp-1);
	if (expect_z) free_dvector (expect_z, 0, nstp-1);
	if (drift_z) free_dvector (drift_z, 0, nstp-1);

	if (pde) num_f_pde_free_2d_adi (pde, nstepx, nstepz, nprod);

	if (x) free(x);
	if (z) free_dvector (z, 0, nstepz-1);
	if (values) free_f3tensor (values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor (values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (mux) free_dmatrix (mux, 0, nstepx-1, 0, nstepz-1);
	if (muz) free_dmatrix (muz, 0, nstepx-1, 0, nstepz-1);
	if (varx) free_dmatrix (varx, 0, nstepx-1, 0, nstepz-1);
	if (varz) free_dmatrix (varz, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix (r, 0, nstepx-1, 0, nstepz-1);

	if (varxinit) free_dvector (varxinit, 0, nstepx-1);
	if (muzinit1) free_dvector (muzinit1, 0, nstepx-1);
	if (muzinit2) free_dvector (muzinit2, 0, nstepx-1);
	if (sigBeta) free_dvector (sigBeta, 0, nstepx-1);	
	
	
			
	return err;
}


Err	FxSabrCalibration(
						/*	Underlying	*/
						char	*dom_yc,
						char	*for_yc,
						long	today,
						double	spot_fx,

						/*	Model Parameters			*/						
						double	alpha,
						double	beta,
						double	rho,
						double	lambda,

						double	floormu,

						/*	Options Parameters			*/
						double	*exercise,
						double	*maturity,
						double	*volatility,
						int		nbOpt,

						
						
						/*	Discretisation Parameters	*/
						long	nstpt_tgt,
						int		nstpfx,
						int		nstpvol,
						int		nbIter,
						double	precision,

						/*	Result						*/
						double	**sigma						
						)
{

Err		err = NULL;

int		i, j, k, l, nb_iter;
double	dt, df, bsvol, coef, newCoef, var, vega; 
double	premium, premium2, premium_tgt;
double	newVar, newVar2;
double	error, errorv;
double	sig0, newSig0;
double	vol1, vol2, var_diff;

double	a, b, c, delta;

long	spot_date;

double	*time		= NULL,
		*date		= NULL,
		*drift_res	= NULL,
		*drift_adi	= NULL,
		**res_iter	= NULL;

int		*eval_evt	= NULL;

int		maxstep, nstpt, old_index;

double	**func_parm	= NULL;

double	newCoef2;

time_t	t1, t2;

	t1 = clock();
	
	/*	For stability reasons */
	if (alpha < 1.0E-05)
	{
		alpha = 1.0E-05;
	}

	if (fabs(beta - 1.0) < 1.0e-08)
	{
		beta = 1.0 - 1.0E-08;
	}
	
	/* look for the spot date */
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	spot_fx = spot_fx * swp_f_df (today, spot_date, dom_yc) / swp_f_df (today, spot_date, for_yc);
	
	maxstep = nstpt_tgt + nbOpt;
	
	date = dvector(0, maxstep-1);
	drift_res = dvector(0, nbOpt-1);
	drift_adi = dvector(0, maxstep-1);
	eval_evt = lvector(0, maxstep-1);
	func_parm = dmatrix(0, maxstep-1, 0, 0);
	res_iter = dmatrix(0, nbIter, 0, 1);
	
	(*sigma) = (double *) calloc(nbOpt, sizeof(double));
		
	if (!date || !drift_res || !drift_adi || !eval_evt || !func_parm || !sigma || !res_iter)
	{
		err = "Memory allocation error (1) in FxSabrCalibration";
		goto FREE_RETURN;
	}
				
	date[0] = today;
	var = 0;

	smessage("Starting Calibration of Fx Sabr");

	/*	First Option to be calibrated through SABRADI formula 
		Not possible anymore because of the lambda
		
	err = op_sabr_calib_adi(spot_fx,
							spot_fx,
							exercise[0],
							volatility[0],
							SRT_LOGNORMAL,
							alpha,
							beta,
							rho,
							nstpt_tgt,
							nstpfx, 
							nstpvol,
							nbIter,
							precision,
							&sig0);


	if (err)
	{
		goto FREE_RETURN;
	}
	
	 */

	/* Calibration of sig0 to the first option */
	time = calloc(2, sizeof(double));
	if (!time)
	{
		err = "Memory allocation error (2) in FxSabrCalibration";
		goto FREE_RETURN;
	}

	time[0] = 0.0;
	time[1] = exercise[0];

	nstpt = 2;
	err = fill_time_vector (&time, 
							&nstpt,
							0,
							NULL,
							0,
							NULL, 
							nstpt_tgt);

	drift_adi[0] = 0.0;

	for (j=1; j<nstpt; j++)
	{			
		date[j] = today + time[j] * 365.0;
		drift_adi[j] = 0.0;
	}

	df = swp_f_df(today, (long) (today + exercise[0] * 365.0 + 1.0E-08), dom_yc);
	sig0 = volatility[0] * exp((1.0 - beta) * log(spot_fx));
	newVar = sig0 * sig0 * exercise[0];

	premium_tgt = srt_f_optblksch(	spot_fx,
									spot_fx,
									volatility[0],
									exercise[0],
									df,
									SRT_CALL,
									PREMIUM
									);

	eval_evt[nstpt-1] = 1;

	err = FxSabr_adi(	nstpt,
						time,
						date,
						nstpfx,
						nstpvol,
						sig0,
						drift_adi,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						func_parm,
						eval_evt,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_ATM_opt,
						1,
						&premium
						);

	if (err)
	{
		goto FREE_RETURN;
	}

	error = fabs(premium - premium_tgt);
	err = srt_f_optimpvol(		
							premium,
							spot_fx,
							spot_fx,
							exercise[0],								
							df,
							SRT_CALL,
							SRT_LOGNORMAL,
							&bsvol);

	errorv = fabs(bsvol - volatility[0]);

	if (errorv >= precision || err)
	{
		/* first adjustment assuming price is linar in stdev */

		newVar2 = premium_tgt * premium_tgt / premium / premium * newVar;
		newSig0 = sqrt(newVar2 / exercise[0]);
	}
	else
	{
		newSig0 = sig0;
	}

	/* Now do a Newton algorithm	*/

	k = 0;
	while ((k < nbIter) && (errorv >= precision))
	{
		err = FxSabr_adi(
						nstpt,
						time,
						date,
						nstpfx,
						nstpvol,
						newSig0,
						drift_adi,
						alpha,
						beta,
						rho,
						lambda,
						floormu,
						func_parm,
						eval_evt,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_ATM_opt,
						1,
						&premium2
						);
		if (err)
		{
			goto FREE_RETURN;
		}

		err = srt_f_optimpvol(		
							premium2,
							spot_fx,
							spot_fx,
							exercise[0],								
							df,
							SRT_CALL,
							SRT_LOGNORMAL,
							&bsvol);

		errorv = fabs(bsvol - volatility[0]);

		if (premium2 > 100 * premium_tgt || premium2 < 0.0 || err)
		{
			/* We have a problem */
			newSig0 *= 0.7;				
		}
		else
		{			
			error = fabs(premium2 - premium_tgt);

			vega = (premium2 - premium) / (newSig0 * newSig0 - sig0 * sig0);
			if (vega < 0)
			{
				/* convergence problem... 
				nstpfx = (int) (nstpfx * 1.5);
				nstpvol = (int) (nstpvol * 1.5);
				*/

				newSig0 *= 0.7;
			}
			else
			{
				sig0 = newSig0;
				newSig0 = sig0 * sig0 - (premium2 - premium_tgt) / vega;

				if (newSig0 < 0.0)
				{
					/* We have a problem */
					newSig0 *= 0.7;
				}
				else
				{
					newSig0 = sqrt(newSig0);
					premium = premium2;
				}
			}
		}

		k++;
	}

	sig0 = newSig0;
	eval_evt[nstpt-1] = 0;

	if (errorv <= precision)
	{
		smessage("Option %d: Success at iteration %d", 1, k);
	}
	else
	{
		smessage("Option %d: May have failed", 1);
	}

	if (time)
	{
		free (time);
		time = NULL;
	}

	(*sigma)[0] = sig0;
	drift_res[0] = 0.0;

	var = sig0 * sig0 * exercise[0];	

	/*	Next options to be calibrated */
	for (i=1; i<nbOpt; i++)
	{
		/* Initialisation	*/

		for (j=0; j<=nbIter; j++)
		{
			res_iter[j][0] = 0.0;
			res_iter[j][1] = 0.0;
		}

		nb_iter = 0;

		/* First construct the time steps */
		
		nstpt = i + 2;
		time = (double *) malloc(nstpt * sizeof(double));
		if (!time)
		{
			err = "Memory allocation error (2) in FxSabrCalibration";
			goto FREE_RETURN;
		}

		/* Add the exercise dates */
		time[0] = 0.0;
		for (j=0; j<=i; j++)
		{
			time[j+1] = exercise[j];
		}
				
		/*	Fill the time vector */
		err = fill_time_vector (&time, 
								&nstpt,
								0,
								NULL,
								0,
								NULL, 
								nstpt_tgt);

		for (j=1; j<nstpt; j++)
		{			
			date[j] = today + time[j] * DAYS_IN_YEAR;
		}

		/* Fill the Drift Vector */
		k = 0;
		for (j=0; j<i; j++)
		{
			while (time[k] < exercise[j] - 1.0E-08)
			{
				drift_adi[k] = drift_res[j];
				k++;
			}
		}

		old_index = k;

		/* Calculation of the first guess */
		vol1 = volatility[i-1] * exp((1.0 - beta) * log(spot_fx));
		vol2 = volatility[i] * exp((1.0 - beta) * log(spot_fx));
		var_diff = (vol2 * vol2 * exercise[i] - vol1 * vol1 * exercise[i-1]);
		
		if (var_diff < 0.0)
		{
			coef = 0.0;
		}
		else
		{			
			dt = exercise[i] - exercise[i-1];
			a = 2.0 / 3.0 * dt * dt;
			b = dt;
			c = 1.0 - var_diff / ((*sigma)[i-1] * (*sigma)[i-1] * dt);
			delta = b * b - 4.0 * a * c;

			if (delta >= 0.0)
			{
				coef = (-b + sqrt(delta)) / (2.0 * a);
			}
			else
			{
				coef = (var_diff / ((*sigma)[i-1] * (*sigma)[i-1]) - dt) / (dt * dt);
			}
		}
		
		/*	Fill term Structure	*/
		for (j=old_index; j<nstpt; j++)
		{
			drift_adi[j] = coef;
		}
		
		newVar = var + (*sigma)[i-1] * (*sigma)[i-1] / (2.0 * coef) * (exp(2.0 * coef * dt) - 1.0);
				
		
		/*	Adjustment in the PDE	*/
		
		df = swp_f_df(today, (long) (today + exercise[i] * 365.0 + 1.0E-08), dom_yc);
		
		premium_tgt = srt_f_optblksch(	spot_fx,
										spot_fx,
										volatility[i],
										exercise[i],
										df,
										SRT_CALL,
										PREMIUM
										);

		/*	Fill the param for the payoff funtion	*/
		
		eval_evt[nstpt-1] = 1;	
		j = 0;

		err = FxSabr_adi(	nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							sig0,
							drift_adi,
							alpha,
							beta,
							rho,
							lambda,
							floormu,
							func_parm,
							eval_evt,
							NULL,
							NULL,
							NULL,
							spot_fx,
							dom_yc,
							for_yc,
							payoff_fx_sabr_adi_ATM_opt,
							1,
							&premium
							);
				
		if (err)
		{
			goto FREE_RETURN;
		}

		error = fabs(premium - premium_tgt);
		err = srt_f_optimpvol(		
								premium,
								spot_fx,
								spot_fx,
								exercise[i],								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bsvol);

		errorv = fabs(bsvol - volatility[i]);

		if (errorv >= precision || err)
		{
			if (err)
			{
				newCoef = coef / 2;
			}
			else
			{
				/* Save result */						
				res_iter[0][0] = coef;
				res_iter[0][1] = premium;
				nb_iter++;

				/* We adjust the coef assuming that the price is more
					and less linear in the stdev						*/

				newVar2 = premium_tgt * premium_tgt / premium / premium * newVar;
				

				dt = exercise[i] - exercise[i-1];
				a = 2.0 / 3.0 * dt * dt;
				b = dt;
				c = 1.0 - (newVar2 - var) / ((*sigma)[i-1] * (*sigma)[i-1] * dt);
				delta = b * b - 4.0 * a * c;

				if (delta >= 0.0)
				{
					newCoef = (-b + sqrt(delta)) / (2.0 * a);
				}
				else
				{
					newCoef = ((newVar2 - var) / ((*sigma)[i-1] * (*sigma)[i-1]) - dt) / (dt * dt);
				}
			}			

			/*	Fill term Structure	*/
			for (j=old_index; j<nstpt; j++)
			{
				drift_adi[j] = newCoef;
			}
		}
		else
		{
			newCoef = coef;
		}

		/* Now do a Newton algorithm	*/

		k = 0;
		while ((k < nbIter) && (errorv >= precision))
		{
			err = FxSabr_adi(
							nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							sig0,
							drift_adi,
							alpha,
							beta,
							rho,
							lambda,
							floormu,
							func_parm,
							eval_evt,
							NULL,
							NULL,
							NULL,
							spot_fx,
							dom_yc,
							for_yc,
							payoff_fx_sabr_adi_ATM_opt,
							1,
							&premium2
							);

			if (err)
			{
				goto FREE_RETURN;
			}

			err = srt_f_optimpvol(		
								premium2,
								spot_fx,
								spot_fx,
								exercise[i],								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bsvol);

			errorv = fabs(bsvol - volatility[i]);

			if (premium2 > 100 * premium_tgt || premium2 < 0.0 || err)
			{
				/* We have a problem */

				newCoef = (coef + newCoef) / 2.0;

				for (j=old_index; j<nstpt; j++)
				{
					drift_adi[j] = newCoef;
				}

				err = NULL;
			}
			else
			{
				l = 0;
				while (l < nb_iter && res_iter[l][0] < newCoef)
				{
					l++;
				}

				for (j=nb_iter-1; j>=l; j--)
				{
					res_iter[j+1][0] = res_iter[j][0];
					res_iter[j+1][1] = res_iter[j][1];
				}

				res_iter[l][0] = newCoef;
				res_iter[l][1] = premium2;
				nb_iter++;


				newCoef2 = solve_for_next_coef(	res_iter,
												nb_iter,
												premium_tgt,
												0);

	
				error = fabs(premium2 - premium_tgt);

				vega = (premium2 - premium) / (newCoef - coef);
				
				if (vega < 0)
				{
					/* convergence problem... 
					nstpfx = (int) (nstpfx * 1.5);
					nstpvol = (int) (nstpvol * 1.5);
					*/

					newCoef = (coef + newCoef) / 2.0;
					for (j=old_index; j<nstpt; j++)
					{
						drift_adi[j] = newCoef;
					}
				}
				else
				{
					coef = newCoef;
					newCoef = coef - (premium2 - premium_tgt) / vega;

					newCoef = newCoef2;

					/*	Fill term Structure	*/
					for (j=old_index; j<nstpt; j++)
					{
						drift_adi[j] = newCoef;
					}

					premium = premium2;
				}
			}
			
			k++;
		}

		drift_res[i] = newCoef;
		

		(*sigma)[i] = (*sigma)[i-1] * exp(newCoef* (exercise[i] - exercise[i-1]));
		
		var += (*sigma)[i-1] * (*sigma)[i-1] / (2.0 * newCoef) * (exp(2.0 * newCoef * (exercise[i] - exercise[i-1])) - 1.0);
		
		eval_evt[nstpt-1] = 0;

		if (errorv <= precision)
		{
			smessage("Option %d: Success at iteration %d", i+1, k);
		}
		else
		{
			smessage("Option %d: May have failed", i+1);
		}

		if (time)
		{
			free (time);
			time = NULL;
		}
	}

	t2 = clock();
	smessage ("Calibration time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:
	
	if (time) free(time);
	if (date) free_dvector(date, 0, maxstep-1);
	if (eval_evt) free_lvector(eval_evt, 0, maxstep-1);
	if (func_parm)	free_dmatrix (func_parm, 0, maxstep-1, 0, 0);

	if (drift_adi) free_dvector(drift_adi, 0, maxstep-1);
	if (drift_res) free_dvector(drift_res, 0, nbOpt-1);
	if (res_iter) free_dmatrix(res_iter, 0, nbIter, 0, 1);

	return err;
}

void Interpolate_Payoff(	double	***values,
							double	*x,
							int		ux,
							int		nstepx,
							int		lz,
							int		uz,
							int		nprod,
							int		bar_col,
							double	barrier,
							int		interp_bar)
{
static	double	d1, d2, coef1, coef2, coef3;
static	int		j, k;
static  double	**valuesi1, 
				**valuesi2,
				**valuesi3;

	if (ux + 1 < nstepx && ux > 0)
	{
		valuesi1 = values[ux+1];
		valuesi2 = values[ux];
		valuesi3 = values[ux-1];

		/* adjust payoff */							
		
		if (nprod != 1)
		{
			d1 = x[ux+1] - x[ux];
			d2 = x[ux] - x[ux-1];

			coef1 = (barrier * (barrier - x[ux] - x[ux-1]) + x[ux] * x[ux-1])
					/ (d1 * (d1 + d2));
			coef2 = (barrier * (barrier - x[ux+1] - x[ux-1]) + x[ux+1] * x[ux-1])
					/ (-d1 * d2);			

			coef3 = 1.0 - coef1 - coef2;
		}

		for (j=lz; j<=uz; j++)
		{
			/* knock out payoff */
			if (interp_bar)
			{
				for (k=0; k<nprod; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
			}
			else
			{
				valuesi2[j][bar_col] = 0.0;
				
				/* other columns, quadratic interpolation */
				for (k=0; k<bar_col; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
				for (k=bar_col+1; k<nprod; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
			}
		}
	}
	else if (ux == nstepx - 1)
	{
		coef2 = (barrier - x[ux-1]) / (x[ux] - x[ux-1]);
		coef3 = 1.0 - coef2;

		valuesi2 = values[ux];
		valuesi3 = values[ux-1];

		for (j=lz; j<=uz; j++)
		{
			if (interp_bar)
			{			
				for (k=0; k<nprod; k++)
				{
					valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
			}
			else
			{
				/* knock out payoff */
				valuesi2[j][bar_col] = 0.0;
				
				/* other columns, quadratic interpolation */
				for (k=0; k<bar_col; k++)
				{
					valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
				for (k=bar_col+1; k<nprod; k++)
				{
					valuesi2[j][k] = coef2 * valuesi2[j][k] + coef3 * valuesi3[j][k];
				}
			}
		}
	}
	else
	{
		coef1 = (barrier - x[ux]) / (x[ux+1] - x[ux]);
		coef2 = 1.0 - coef1;

		valuesi1 = values[ux+1];
		valuesi2 = values[ux];

		for (j=lz; j<=uz; j++)
		{
			if (interp_bar)
			{
				for (k=0; k<nprod; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
				}
			}
			else
			{
				/* knock out payoff */
				valuesi2[j][bar_col] = 0.0;
				
				/* other columns, quadratic interpolation */
				for (k=0; k<bar_col; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
				}
				for (k=bar_col+1; k<nprod; k++)
				{
					valuesi2[j][k] = coef1 * valuesi1[j][k] + coef2 * valuesi2[j][k];
				}
			}
		}
	}
}

void	disc_normal_center(		double	*x,
								long	nstp,
								double	fwd,
								double	xmin,
								double	xmax,
								double	stdn,
								double	nstd,
								long	*index)
{
double	dmin, dmax;
double	meshx;
long	i;

	(*index) = (nstp - 1) / 2;
	
	dmin = norm((xmin - fwd) / stdn);
	dmax = norm((xmax - fwd) / stdn);

	meshx = (0.5 - dmin) / (nstp - 1) * 2;

	x[0] = xmin;

	for (i=1; i<=(*index); i++)
	{
		x[i] = fwd + stdn * inv_cumnorm_fast(dmin + i * meshx);
	}

	meshx = (dmax - 0.5) / (nstp - 1) * 2;

	for (i=(*index)+1; i<nstp-1; i++)
	{
		x[i] = fwd + stdn * inv_cumnorm_fast(0.5 + (i - (*index)) * meshx);
	}

	x[nstp-1] = xmax;
}

void	disc_linleft_logright(	double	*x,
								long	nstp,
								double	fwd,
								double	xmin,
								double	xmax,
								double	stdln,
								double	nstd,
								long	*index)
{
double	meshx;
long	i;

	meshx = log(xmax / xmin) / (nstp - 1);


	(*index) = (int) (log(fwd / xmin) / meshx + 0.5);
	
				
	x[0] = xmin;

	meshx = (fwd - xmin) / (*index);
	
	for (i=1; i<=(*index); i++)
	{
		x[i] = x[i-1] + meshx;
	}

	meshx = exp(log(xmax / fwd) / (nstp-1 - (*index)));

	for (i=(*index)+1; i<nstp; i++)
	{
		x[i] = x[i-1] * meshx;
	}
}

void	disc_linleft_logright_center(	double	*x,
										long	nstp,
										double	fwd,
										double	xmin,
										double	xmax,
										double	stdln,
										double	nstd,
										long	*index)
{
double	meshx;
long	i;


	(*index) = (int) (nstp / 2.0);
	
	x[0] = xmin;

	meshx = (fwd - xmin) / (*index);
	
	for (i=1; i<=(*index); i++)
	{
		x[i] = x[i-1] + meshx;
	}

	meshx = exp(log(xmax / fwd) / (nstp-1 - (*index)));

	for (i=(*index)+1; i<nstp; i++)
	{
		x[i] = x[i-1] * meshx;
	}
}

void	disc_linleft_linright_center(	double	*x,
										long	nstp,
										double	fwd,
										double	xmin,
										double	xmax,
										double	stdln,
										double	nstd,
										long	*index)
{
double	meshx;
long	i;


	(*index) = (int) (nstp / 2.0);
	
	x[0] = xmin;

	meshx = (fwd - xmin) / (*index);
	
	for (i=1; i<=(*index); i++)
	{
		x[i] = x[i-1] + meshx;
	}

	meshx = (xmax - fwd) / (nstp-1 - (*index));

	for (i=(*index)+1; i<nstp; i++)
	{
		x[i] = x[i-1] + meshx;
	}
}

void	disc_linleft_logright_bar(	double	**x,
									long	*nstp,
									double	fwd,
									double	xmin,
									double	xmax,
									double	stdln,
									double	nstd,
									double	*barrier,
									long	nb_bar,
									double	*bar_pres,
									long	*index)
{
double	meshx, meshx2;
long	i, j;
double	*bar;
long	nb_bar2;
int		nb_dec;
long	temp;
int		exit;

	bar = malloc(nb_bar * sizeof(double));	

	if (fwd >= 10)
	{
		nb_dec = 10000;
	}
	else if (fwd >= 1)
	{
		nb_dec = 10000;
	}
	else
	{
		nb_dec = 10000;
	}
		
	for (i=0; i<nb_bar; i++)
	{		
		temp = (int) (barrier[i] * nb_dec);
		bar[i] = (temp * 1.0) / nb_dec;
	}

	num_f_add_number(&nb_bar, &bar, fwd);
	num_f_add_number(&nb_bar, &bar, xmax);
	

	nb_bar2 = nb_bar;
	num_f_sort_vector(nb_bar2, bar);
	num_f_unique_vector(&nb_bar2, bar);

	*x = realloc(*x, (*nstp + nb_bar2) * sizeof(double));

	meshx = log(xmax / xmin) / (*nstp - 1);
	(*index) = (int) (log(fwd / xmin) / meshx + 0.5);
	meshx = (fwd - xmin) / (*index);
	meshx2 = exp(log(xmax / fwd) / (*nstp-1 - (*index)));
				
	(*x)[0] = xmin;

	/* first Barrier */ 

	j = 0;
	while (j < nb_bar2 && bar[j] < xmin)
	{
		j++;
	}

	i = 1;

	exit = 1;
	while (exit)
	{
		/* discretisation of x */
		

		if (j <nb_bar2 && bar[j] < (*x)[i-1] + meshx + 1.0E-08)
		{
			/* we add a barrier point */

			if (bar[j] - (*x)[i-1] > meshx / 2 || (j > 0 && fabs((*x)[i-1] - bar[j-1]) < 1.0E-06))
			{
				(*x)[i] = bar[j];							
			}
			else
			{
				i--;
				(*x)[i] = bar[j];				
			}
				
			if (fabs(bar[j] - fwd) < 1.0E-08)
			{		
				*index = i;
				exit = 0;
			}
			else
			{
				i++;
				(*x)[i] = 2.0 * bar[j] - (*x)[i-2];
				if (j < nb_bar2-1 && (*x)[i] > bar[j+1])
				{
					i--;
				}
			}			

			i++;
			j++;
		}
		else
		{
			/* we add a regular point */
			(*x)[i] = (*x)[i-1] + meshx;
			i++;
		}
	}

	

	exit = 1;	
	while (exit)
	{
		/* discretisation of x */
		if (j <nb_bar2 && bar[j] < (*x)[i-1] * meshx2 + 1.0E-08)
		{
			/* we add a barrier point */

			if (bar[j] > (*x)[i-1] * (1 + meshx2) / 2.0 || (j > 0 && fabs((*x)[i-1] - bar[j-1]) < 1.0E-06))
			{
				(*x)[i] = bar[j];
			}
			else
			{
				i--;
				(*x)[i] = bar[j];
			}

			if (fabs(bar[j] - xmax) < 1.0E-08)
			{
				*nstp = i + 1;
				exit = 0;
			}
			else
			{
				i++;
				(*x)[i] = 2.0 * bar[j] - (*x)[i-2];
				if (j < nb_bar2-1 && (*x)[i] > bar[j+1])
				{
					i--;
				}
			}

			i++;
			j++;
		}
		else
		{
			/* we add a regular point */
			(*x)[i] = (*x)[i-1] * meshx2;
			i++;
		}
	}

	*bar_pres = 1.0 / nb_dec + 1.0E-08;

	free (bar);

}

void	disc_linleft_logright_strike(	double	*x,
										long	nstp,
										double	fwd,
										double	xmin,
										double	xmax,
										double	stdln,
										double	nstd,
										double	strike,
										long	*index)
{
double	meshx;
long	i, index_s;

	meshx = log(xmax / xmin) / (nstp - 1);


	(*index) = (int) (log(fwd / xmin) / meshx + 0.5);
	
				
	x[0] = xmin;

	meshx = (fwd - xmin) / (*index);
	index_s = (int) ((strike - xmin) / meshx + 0.5);

	if (strike > fwd || strike < xmin || index_s == (*index) || index_s < 1)
	{				
		for (i=1; i<=(*index); i++)
		{
			x[i] = x[i-1] + meshx;
		}
	}
	else
	{		
		meshx = (strike - xmin) / index_s;
		
		for (i=1; i<=index_s; i++)
		{
			x[i] = x[i-1] + meshx;
		}

		meshx = (fwd - strike) / ((*index) - index_s);		
		
		for (i=index_s+1; i<=(*index); i++)
		{
			x[i] = x[i-1] + meshx;
		}
	}

	meshx = exp(log(xmax / fwd) / (nstp-1 - (*index)));

	index_s = (int) (log(strike / fwd) / log(meshx) + 0.5) + (*index);

	if (strike <fwd || strike > xmax || index_s == (*index) || index_s > nstp-2)
	{
		for (i=(*index)+1; i<nstp; i++)
		{
			x[i] = x[i-1] * meshx;
		}
	}
	else
	{
		meshx = exp(log(strike / fwd) / (index_s - (*index)));

		for (i=(*index)+1; i<=index_s; i++)
		{
			x[i] = x[i-1] * meshx;
		}

		meshx = exp(log(xmax / strike) / (nstp-1 - index_s));

		for (i=index_s+1; i<nstp; i++)
		{
			x[i] = x[i-1] * meshx;
		}
	}
}

void	disc_linleft_logright_center_strike(	double	*x,
												long	nstp,
												double	fwd,
												double	xmin,
												double	xmax,
												double	stdln,
												double	nstd,
												double	strike,
												long	*index)
{
double	meshx;
long	i, index_s;

	(*index) = (int) (nstp / 2.0);
				
	x[0] = xmin;
	meshx = (fwd - xmin) / (*index);

	if (strike < fwd * 0.9999 && strike > xmin)
	{
		index_s = (int) ((strike - xmin) / meshx + 0.5);

		if (index_s == (*index))
		{
			index_s = 0;
		}
		else
		{
			meshx = (strike - xmin) / index_s;
		}
		
		for (i=1; i<=index_s; i++)
		{
			x[i] = x[i-1] + meshx;
		}

		meshx = (fwd - strike) / ((*index) - index_s);		
		
		for (i=index_s+1; i<(*index); i++)
		{
			x[i] = x[i-1] + meshx;
		}
	}
	else
	{
		for (i=1; i<(*index); i++)
		{
			x[i] = x[i-1] + meshx;
		}
	}

	x[*index] = fwd;
	i++;

	meshx = exp(log(xmax / fwd) / (nstp-1 - (*index)));

	if (strike > fwd * 1.00001 && strike < xmax)
	{			
		index_s = (int) (log(strike / fwd) / log(meshx) + 0.5) + (*index);
		
		if (index_s == (*index))
		{
			index_s = 0;
		}
		else
		{		
			meshx = exp(log(strike / fwd) / (index_s - (*index)));
		}

		for (i=(*index)+1; i<=index_s; i++)
		{
			x[i] = x[i-1] * meshx;
		}

		meshx = exp(log(xmax / strike) / (nstp-1 - index_s));

		for (i=index_s+1; i<nstp-1; i++)
		{
			x[i] = x[i-1] * meshx;
		}
	}
	else
	{
		for (i=(*index)+1; i<nstp-1; i++)
		{
			x[i] = x[i-1] * meshx;
		}
	}

	x[nstp-1] = xmax;
}


void	disc_linleft_logright_bar2(	double	**x,
									long	*nstp,
									double	fwd,
									double	xmin,
									double	xmax,
									double	stdln,
									double	nstd,
									double	*barrier,
									long	nb_bar,
									double	*bar_pres,
									long	*index)
{
double	meshx, meshx2;
long	i, j;
double	*bar;
long	nb_bar2;
int		nb_dec;
long	temp;
int		exit;

	bar = malloc(nb_bar * sizeof(double));	

	if (fwd >= 10)
	{
		nb_dec = 10000;
	}
	else if (fwd >= 1)
	{
		nb_dec = 10000;
	}
	else
	{
		nb_dec = 10000;
	}
		
	for (i=0; i<nb_bar; i++)
	{		
		temp = (int) (barrier[i] * nb_dec);
		bar[i] = (temp * 1.0) / nb_dec;
	}

	num_f_add_number(&nb_bar, &bar, fwd);
	num_f_add_number(&nb_bar, &bar, xmax);
	

	nb_bar2 = nb_bar;
	num_f_sort_vector(nb_bar2, bar);
	num_f_unique_vector(&nb_bar2, bar);

	*x = realloc(*x, (*nstp + nb_bar2) * sizeof(double));

	meshx = log(xmax / xmin) / (*nstp - 1);
	(*index) = (int) (log(fwd / xmin) / meshx + 0.5);
	meshx = (fwd - xmin) / (*index);

	meshx2 = exp(log(xmax / fwd) / (*nstp-1 - (*index)));
				
	(*x)[0] = xmin;

	/* first Barrier */ 

	j = 0;
	while (j < nb_bar2 && bar[j] < xmin)
	{
		j++;
	}

	i = 1;

	exit = 1;
	while (exit)
	{
		/* discretisation of x */
		

		if (j < nb_bar2 && bar[j] < (*x)[i-1] + meshx + 1.0E-08)
		{
			/* we add a barrier point */
			
			(*x)[i] = bar[j];							
							
			if (fabs(bar[j] - fwd) < 1.0E-08)
			{		
				*index = i;
				exit = 0;
			}			

			i++;
			j++;
		}
		else
		{
			/* we add a regular point */
			(*x)[i] = (*x)[i-1] + meshx;
			i++;
		}
	}

		

	exit = 1;	
	while (exit)
	{
		/* discretisation of x */
		if (j <nb_bar2 && bar[j] < (*x)[i-1] * meshx2 + 1.0E-08)
		{
			/* we add a barrier point */
			
			(*x)[i] = bar[j];
						
			if (fabs(bar[j] - xmax) < 1.0E-08)
			{
				*nstp = i + 1;
				exit = 0;
			}			

			i++;
			j++;
		}
		else
		{
			/* we add a regular point */
			(*x)[i] = (*x)[i-1] * meshx2;
			i++;
		}
	}

	*bar_pres = 1.0 / nb_dec + 1.0E-08;

	free (bar);

}

double solve_for_next_coef(	double	**res_iter,
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

		if (vega < 0)
		{
			/* convergence problem */
			coef = 0.5 * (res_iter[i+1][0] - res_iter[i][0]);
		}
		else
		{
			coef = res_iter[i][0] + (premium_tgt - res_iter[i][1]) / vega;
		}
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
					

		a =		res_iter[i-1][1] / (res_iter[i-1][0] - res_iter[i][0]) / (res_iter[i-1][0] - res_iter[i+1][0])
			+	res_iter[i][1] / (res_iter[i][0] - res_iter[i-1][0]) / (res_iter[i][0] - res_iter[i+1][0])
			+	res_iter[i+1][1] / (res_iter[i+1][0] - res_iter[i-1][0]) / (res_iter[i+1][0] - res_iter[i][0]);

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
			coef = solve_for_next_coef(	res_iter,
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

#define BUTTERFLY	0.0005
#define ITER_MAX	50
#define	PRES_PERC	0.00001

double find_beta_lim(	double	forward,
						double	std_beta,
						double	beta,
						double	percent)
{
double	eps, pres, eps_lim;
double	level1, level2, perc1, perc2, perc_tgt, deriv;
int		i;

	eps = BUTTERFLY * forward;
	eps_lim = eps * 0.999;
	pres = PRES_PERC * (2.0 * eps);
	perc_tgt = percent * (2.0 * eps);

	level1 = forward;

	perc1 = - srt_f_optblkschbeta(forward, level1 + eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM)
		    + srt_f_optblkschbeta(forward, level1 - eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM);
	
	level2 = forward * 1.01;

	perc2 = - srt_f_optblkschbeta(forward, level2 + eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM)
		    + srt_f_optblkschbeta(forward, level2 - eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM);

	i = 0;

	while (i < ITER_MAX && fabs(perc2 - perc_tgt) > pres)
	{
		deriv = (perc2 - perc1) / (level2 - level1);
		level1 = level2;
		perc1 = perc2;
		level2 = level2  + (perc_tgt - perc2) / deriv;

		if (level2 < 0)
		{
			level2 = level1 * 0.9;
		}

		if (level2 > eps_lim)
		{				
			perc2 =	- srt_f_optblkschbeta(forward, level2 + eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM)
					+ srt_f_optblkschbeta(forward, level2 - eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM);
		}
		else
		{
			perc2 =	2.0 * (- srt_f_optblkschbeta(forward, level2 + eps, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM)
					+ srt_f_optblkschbeta(forward, level2, std_beta, 1.0, beta, 1.0, SRT_CALL, PREMIUM));
		}

		i++;
	}
	
	return level2;
}