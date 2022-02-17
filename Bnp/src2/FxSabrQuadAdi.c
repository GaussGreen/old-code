/* ==========================================================================
   FILE_NAME:	LGM2FSabr.c

   PURPOSE:		ADI implementation of the Sabr Fx model.
				Discretisation is ADI.

   DATE:		03/16/01
   
   AUTHOR:		L.C.
   ========================================================================== */

#include "FxSabrAdi.h"
#include "FxSabrQuadAdi.h"
#include "FxSabrQuadAdi.h"
#include "opfnctns.h"
#include "Fx3FUtils.h"
#include "FxSabrGrfn.h"
#include "opsabrcalib.h"
#include "math.h"

#define	NSTD_FX		5.0
#define	NSTD_VOL	5.0

#define MINNODE		25
#define	MAXNODE2	50

#define EPSALPHA	0.02
#define	EPSRHO		0.05
#define EPSCVXTY	0.001
#define	EPSBETA		0.05

#define EPSALPHAMIN	0.01
#define EPSALPHAMAX	0.10

#define	EPSRHOMIN	0.01
#define	EPSRHOMAX	0.05

#define	EPSCVXTYMIN	0.0001
#define	EPSCVXTYMAX	0.005

#define	EPSBETAMIN	0.01
#define	EPSBETAMAX	0.25

#define	NBITERMAX	50

void SabrQuadGetParams(	double	cash_fx,
						double	gamma,
						double	beta,
						double	*quad_param,
						double	*lin_param,
						double	*const_param)
{
	/*
	(*quad_param) = gamma;
	*/

	(*quad_param) = exp((beta - 2.0) * log(cash_fx)) * (gamma + beta * (beta - 1.0)) / 2.0;

	(*lin_param) = beta * exp((beta - 1.0) * log(cash_fx)) - 2.0 * (*quad_param) * cash_fx;
	(*const_param) = exp(beta * log(cash_fx)) - cash_fx * ((*quad_param) * cash_fx + (*lin_param));
}

Err FxSabrQuadPrecalculations(
							/*	Time data		*/
							int			nstp,
							double		*time,
							
							/*	Model data		*/
							double		sig0,
							double		*drift,
							double		alpha,
							double		a,
							double		b,
							double		c,
							double		rho,
							double		lambda,

							/*	Market data */
							double		spot_fx,

							/*	Output			*/							
							double		*expect_z,
							double		*drift_z,
							double		*std1)
{
double	mean, mean2, const1, const2, const3, tmid;
double	temp, exp_z, varxf, sig, dt;
int		i;

	/* Initialisation */	
	exp_z = rho * alpha * vol_quadra(spot_fx, a, b, c, 3);	
	mean = 2.0 * lambda - alpha * alpha;

	varxf = 0.0;
	sig = 1.0;
	temp = 0.0;

	if (fabs(mean) < 1.0E-10)
	{
		mean = 1.0E-05;
	}

	/*	for the next calculations we are going to use the Quad calcul
		with the new following new a */

	a = 2 * a * spot_fx + b;

	const3 = -0.5 * rho * alpha * a / mean * sig0 * sig0;
	const1 = const3 *  2.0 * lambda;
	const2 = -const3 * alpha * alpha;	

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

		/* Then update the expectation and drift at mid time */
		dt /= 2.0;
		tmid = time[i] + dt;

		if (fabs(drift[i]) > 1.0E-10)
		{
			exp_z += const1 * sig * sig / (2.0 * drift[i]) * (exp(2.0 * drift[i] * dt) - 1.0);			
		}
		else
		{
			exp_z += const1 * sig * sig * dt;			
		}

		mean2 = 2.0 * drift[i] - mean;

		if (fabs(mean2) > 1.0E-10)
		{
			exp_z += const2 * sig * sig / mean2 * (exp(2.0 * drift[i] * dt - mean * tmid) 
												- exp(-mean * time[i]));
		}
		else
		{
			exp_z += const2 * sig * sig * exp(-2.0 * drift[i] * time[i]) * dt;
		}

		temp += drift[i] * dt;
		sig = exp(temp);		

		expect_z[i] = exp_z - sig0 * sig;
		

		drift_z[i] = (const3 * sig * sig * (2.0 * lambda - alpha * alpha * exp(-mean * tmid))
					+ (lambda - drift[i]) * sig0 * sig) * 2.0 * dt;

		/* Eventually update expexctation at time [i+1] */

		if (fabs(drift[i]) > 1.0E-10)
		{
			exp_z += const1 * sig * sig / (2.0 * drift[i]) * (exp(2.0 * drift[i] * dt) - 1.0);			
		}
		else
		{
			exp_z += const1 * sig * sig * dt;			
		}

		mean2 = 2.0 * drift[i] - mean;

		if (fabs(mean2) > 1.0E-10)
		{
			exp_z += const2 * sig * sig / mean2 * (exp(2.0 * drift[i] * dt - mean * time[i+1]) 
												- exp(-mean * tmid));
		}
		else
		{
			exp_z += const2 * sig * sig * exp(-2.0 * drift[i] * tmid) * dt;
		}
		
		temp += drift[i] * dt;
		sig = exp(temp);
	}

	expect_z[nstp-1] = exp_z - sig0 * sig;

	*std1 = sig0 * sqrt(varxf);	

	return NULL;
}
				  
Err FxSabrQuad_adi(	
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
				double		a,
				double		b,
				double		c,
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
				double		*res,
				
				/* Additional informations */
				int			calc_greeks,
				double		**greeks,		/* array 6 * nprod containing delta, gamma, theta, vega, volga and vanna */
				
				/* For calibration purpose */
				int			calc_at_point,
				int			column,
				double		target,
				double		*vol,
				double		*res_at_point)
{
Err					err = NULL;

long				today;				
int					i, j, k, step;
int					index_x, index_z, nstepx, nstepz;
double				fux, flx, dt, dft;
double				std1, std3, stdln;
					
double				*expect_z		= NULL,
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
					*muzinit		= NULL,
					*sigBeta		= NULL,
					**r				= NULL;
					
double				*varxi,
					*varzi,
					*muzi,					
					**valuesi,
					**valuesip,
					**valuesim,
					*valuesij;

double				const_sigi, const_varz, const_muz, const_muz1, const_muz2, const_varxi, sigBetaij;
double				alpharho_a, alpharhoa5;

					
int					lx, ux, lz, uz;

double				bar_pres = 1.0E-08,
					bar_c1, bar_c2;
long				bar_index;

double				dxu, dxd, dzu, dzd, vega, coef, is_up, newz;
					
clock_t				t1, t2;

CNPDE_TEMP_2D_ADI	pdestr, *pde = NULL;

double				logx, logmax, logmin;

	/*	For stability reasons */

	if (alpha < 1.0e-05)
	{
		alpha = 1.0E-05;
	}
	
	t1 = clock();
	
	/* Constant	calculations */	

	alpharho_a = alpha * rho;
	alpharhoa5 = 0.5 * alpha * rho;
	const_varz = alpha * alpha * (1.0 - rho * rho);
	today = (long) (date[0] + 1.0e-06);	

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

	expect_z = dvector (0, nstp-1);	
	drift_z = dvector (0, nstp-1);
	
	if (!expect_z || !drift_z || !x || !z)
	{
		err = "Memory allocation error (1) in FxSabrQuad_Adi";
		goto FREE_RETURN;
	}

	/*	Calculate the standard deviation of X and the expectaion of Z	*/	
	err = FxSabrQuadPrecalculations(nstp,
								time,
								sig0,
								drift,
								alpha,
								a,
								b,
								c,
								rho,
								lambda,
								spot_fx,
								expect_z,
								drift_z,
								&std1);
	
	/*	Boundaries */
	
	std1 /= sqrt(time[nstp-1]);

	stdln = op_sabrgen(spot_fx,
						spot_fx,
						time[nstp-1],
						std1,
						alpha,
						a,
						b,
						c,
						rho,
						vol_quadra);

	stdln *= sqrt(time[nstp-1]);
	
	flx = spot_fx * exp(-0.5 * stdln * stdln - NSTD_FX * stdln);
	fux = spot_fx * exp(-0.5 * stdln * stdln + NSTD_FX * stdln);

	/*	Bound for the calculation of max and min drift and var */

	logmin = spot_fx * exp(-0.5 * stdln * stdln - floorstd * stdln);
	logmax = spot_fx * exp(-0.5 * stdln * stdln + floorstd * stdln);
	
	std3 = sig0 * sqrt(const_varz * time[nstp-1]);
													
	/*	Discretisation of space in the orthogonal system x / z */

	if (bar_lvl)
	{
		disc_SL_center_barrier(&x, &nstepx, spot_fx, flx, fux, 0.9999, 0.00001, stdln, NSTD_FX, bar_lvl, nstp, &index_x);
	}
	else
	{
		//disc_linleft_logright_center(x, nstepx, spot_fx, flx, fux, stdln, NSTD_FX, &index_x);
		//disc_linleft_linright_center(x, nstepx, spot_fx, flx, fux, stdln, NSTD_FX, &index_x);
		//disc_Quad_center(x, nstepx, spot_fx, flx, fux, 0.9999, 0.00001, stdln, NSTD_FX, &index_x);
		disc_SL_center(x, nstepx, spot_fx, flx, fux, 0.9999, 0.00001, stdln, NSTD_FX, &index_x);
	}

	/*
	disc_normal_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3, NSTD_VOL, &index_z);
	*/

	disc_linleft_linright_center(z, nstepz, 0, -NSTD_VOL * std3, NSTD_VOL * std3, std3, NSTD_VOL, &index_z);

	/* first allocate memory */

	values = f3tensor (0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor (0, nstepx-1, 0, nstepz-1, 0, nprod-1);

	mux = dmatrix (0, nstepx-1, 0, nstepz-1);
	muz = dmatrix (0, nstepx-1, 0, nstepz-1);
	varx = dmatrix (0, nstepx-1, 0, nstepz-1);
	varz = dmatrix (0, nstepx-1, 0, nstepz-1);
	r = dmatrix (0, nstepx-1, 0, nstepz-1);

	sigBeta = dvector (0, nstepx-1);	
	varxinit = dvector (0, nstepx-1);
	muzinit = dvector (0, nstepx-1);

	if (!x || !z || !values || !values_p1 || !mux || !muz || !sigBeta
		|| !varx || !varxinit || !muzinit || !varz  || !r)
	{
		err = "Memory allocation error (2) in FxSabrQuad_Adi";
		goto FREE_RETURN;
	}	

	/*	Precalculate variances and expectations				*/	
	for (i=0; i<nstepx; i++)
	{		
		logx = min(max(x[i], logmin), logmax);
		logx = logx * (a * logx + b) + c;
		varxinit[i] = logx * logx;

		logx = min(max(x[i], logmin), logmax);
		muzinit[i] = alpharhoa5 * (2.0 * a * logx + b);
		sigBeta[i] = alpharho_a * vol_quadra(logx, a, b, c, 3);
	}

	/*	Final payoff valuation					*/
	if (!eval_evt[nstp-1])
	{
		err = "No event at last step in FxSabrQuad_Adi";
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
		err = "Memory allocation error (3) in FxSabrQuad_Adi";
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

		const_muz1 = (drift[step] - lambda) * dt;						

		for (i=lx; i<=ux; i++)
		{			
			varxi = varx[i];
			varzi = varz[i];
			muzi = muz[i];
			const_muz = muzinit[i];
	
			const_sigi = sigBeta[i] - expect_z[step];
			const_varxi = varxinit[i];
						
			for (j=lz; j<=uz; j++)
			{
				/* reconstruction of the vol */
				sigBetaij = max(const_sigi - z[j], 1.0E-08);

				const_muz2 = const_muz1 * sigBetaij;
				
				sigBetaij *= sigBetaij * dt;
				
				varxi[j] = const_varxi * sigBetaij;
				muzi[j] = -const_muz * sigBetaij - const_muz2 - drift_z[step];
				varzi[j] = const_varz * sigBetaij;
			}
		}		

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

	/* compute greeks if needed */
	if (calc_greeks)
	{
		dxu = x[index_x+1] - x[index_x];
		dxd = x[index_x] - x[index_x-1];

		dzu = -(z[index_z+1] - z[index_z]);
		dzd = -(z[index_z] - z[index_z-1]);

		dt = time[1] - time[0];

		for (k=0; k<nprod; k++)
		{
			greeks[0][k] = (values_p1[index_x+1][index_z][k] - values_p1[index_x-1][index_z][k]) / (dxu + dxd);
			greeks[1][k] = 2.0 * (values_p1[index_x+1][index_z][k] * dxd + values_p1[index_x-1][index_z][k] * dxu
				- (dxd + dxu) * res[0]) / (dxu * dxd * (dxu + dxd));
			greeks[2][k] = (values[index_x][index_z][k] - values_p1[index_x][index_z][k]) / dt / DAYS_IN_YEAR;

			greeks[3][k] = (values_p1[index_x][index_z+1][k] - values_p1[index_x][index_z-1][k]) / (dzu + dzd);
			greeks[4][k] = 2.0 * (values_p1[index_x][index_z+1][k] * dzd + values_p1[index_x][index_z-1][k] * dzu
				- (dzu + dzd) * res[0]) / (dzu * dzd * (dzu + dzd));

			greeks[5][k] = ((values_p1[index_x+1][index_z+1][k] - values_p1[index_x-1][index_z+1][k]) / (dxu + dxd)
							- (values_p1[index_x+1][index_z-1][k] - values_p1[index_x-1][index_z-1][k]) / (dxu + dxd)) / (dzu + dzd);
		}
	}

	if (calc_at_point)
	{
		/* look on z for the closest value to target */

		if ((values_p1[index_x][index_z+1][column] - res[column]) / (z[index_z+1] - z[index_z]) > 0.0)
		{
			/* increasing function of index_z */
			is_up = 1;
		}
		else
		{
			is_up = -1;
		}

		if ((target - res[column]) * is_up > 0.0)
		{
			i = index_z;

			while (target < values_p1[index_x][i][column] && i < nstepz)
			{
				i++;
			}
			if (i == nstepz)
			{
				i = nstepz - 1;
			}

			vega = (values_p1[index_x][i][column] - values_p1[index_x][i-1][column]) / (z[i] - z[i-1]);

			newz = (target - values_p1[index_x][i-1][column]) / vega + z[i-1];

			*vol = sig0 - newz;
			
			coef = (newz - z[i-1]) / (z[i] - z[i-1]);

			for (k=0; k<nprod; k++)
			{
				res_at_point[k] = coef * values_p1[index_x][i][k] + (1.0 - coef) * values_p1[index_x][i-1][k];
			}
		}
		else
		{
			i = index_z;

			while (target > values_p1[index_x][i][column] && i >= 0)
			{
				i--;
			}
			if (i == -1)
			{
				i = 0;
			}

			vega = (values_p1[index_x][i+1][column] - values_p1[index_x][i][column]) / (z[i+1] - z[i]);

			newz = (target - values_p1[index_x][i][column]) / vega + z[i];

			*vol = sig0 - newz;
			
			coef = (newz - z[i]) / (z[i+1] - z[i]);

			for (k=0; k<nprod; k++)
			{
				res_at_point[k] = coef * values_p1[index_x][i+1][k] + (1.0 - coef) * values_p1[index_x][i][k];
			}
		}
	}

	t2 = clock();
	
	smessage ("Convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

	
	/* Allocation 1 */
	if (expect_z) free_dvector (expect_z, 0, nstp-1);
	if (drift_z) free_dvector (drift_z, 0, nstp-1);	

	if (x) free(x);
	if (z) free_dvector (z, 0, nstepz-1);

	/* Allocation 2 */
	if (values) free_f3tensor (values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor (values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	
	if (mux) free_dmatrix (mux, 0, nstepx-1, 0, nstepz-1);
	if (muz) free_dmatrix (muz, 0, nstepx-1, 0, nstepz-1);
	if (varx) free_dmatrix (varx, 0, nstepx-1, 0, nstepz-1);
	if (varz) free_dmatrix (varz, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix (r, 0, nstepx-1, 0, nstepz-1);

	if (varxinit) free_dvector (varxinit, 0, nstepx-1);
	if (sigBeta) free_dvector (sigBeta, 0, nstepx-1);
	if (muzinit) free_dvector (muzinit, 0, nstepx-1);

	/* Allocation 3 */
	if (pde) num_f_pde_free_2d_adi (pde, nstepx, nstepz, nprod);
				
	return err;
}

Err	FxSabrQuadCalibration(
						/*	Underlying	*/
						char	*dom_yc,
						char	*for_yc,
						long	today,
						double	spot_fx,

						/*	Model Parameters			*/						
						double	alpha,
						double	coefa,
						double	coefb,
						double	coefc,
						double	rho,
						double	lambda,

						double	floormu,

						/*	Options Parameters			*/
						double	*exercise,
						double	*maturity,
						double	*volatility,
						int		nbOpt,
						int		start_opt,
						int		is_straddle,

						/*	Discretisation Parameters	*/
						long	nstpt_tgt,
						int		nstpfx,
						int		nstpvol,
						int		nbIter,
						double	precision,

						double	*guess,

						/*	Result						*/
						double	*sigma						
						)
{

Err		err = NULL;

int		i, j, k, l, nb_iter;
double	dt, df, bsvol, coef, newCoef, var, vega; 
double	premium, premium2, premium_tgt;
double	newVar, newVar2;
double	error, errorv;
double	sig0, newSig0;
double	vol1, vol2, volQuad1, volQuad2, var_diff;
double	shift_spot;
double	coef_lin, coef_linab;
int		is_put;

long	spot_date;

double	*time		= NULL,
		*date		= NULL,
		*drift_res	= NULL,
		*drift_adi	= NULL,
		**res_iter	= NULL;

int		*eval_evt	= NULL;

int		maxstep, nstpt, old_index;

double	**func_parm	= NULL;
double	strike, shift_strike;

double	newCoef2, temp, alpha_eq;

time_t	t1, t2;

	t1 = clock();
	
	/*	For stability reasons */
	if (alpha < 1.0E-05)
	{
		alpha = 1.0E-05;
	}
	
	coef_lin = coefa * spot_fx + coefb;
	if (coef_lin > 0.0)
	{
		coef_linab = coef_lin;
		is_put = 0;
	}
	else
	{
		coef_linab = -coef_lin;
		is_put = 1;
	}

	/* look for the spot date */
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	spot_fx = spot_fx * swp_f_df (today, spot_date, dom_yc) / swp_f_df (today, spot_date, for_yc);

	shift_spot = spot_fx * coef_lin + coefc;
	
	maxstep = nstpt_tgt + nbOpt;
	
	date = dvector(0, maxstep-1);
	drift_res = dvector(0, nbOpt-1);
	drift_adi = dvector(0, maxstep-1);
	eval_evt = lvector(0, maxstep-1);
	func_parm = dmatrix(0, maxstep-1, 0, 0);
	res_iter = dmatrix(0, nbIter, 0, 1);	
		
	if (!date || !drift_res || !drift_adi || !eval_evt || !func_parm || !res_iter)
	{
		err = "Memory allocation error (1) in FxSabrCalibration";
		goto FREE_RETURN;
	}
				
	date[0] = today;
	var = 0;
	strike = spot_fx;
	shift_strike = shift_spot;

	smessage("Starting Calibration of Fx Quad Sabr");

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
	if (start_opt == 0)
	{
		alpha_eq = alpha;

		if (fabs(lambda) > 1.0E-08)
		{
			alpha_eq *= sqrt((1.0 - exp(-2.0 * lambda * exercise[0])) / (2.0 * lambda * exercise[0]));
		}

		time = calloc(2, sizeof(double));
		if (!time)
		{
			err = "Memory allocation error (2) in FxSabrCalibration";
			goto FREE_RETURN;
		}

		for (j=0; j<=nbIter; j++)
		{
			res_iter[j][0] = 0.0;
			res_iter[j][1] = 0.0;
		}

		nb_iter = 0;

		if (is_straddle)
		{
			strike = spot_fx * exp(-0.5 * volatility[0] * volatility[0] * exercise[0]);
			shift_strike = strike * coef_lin + coefc;
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
		sig0 = volatility[0] * spot_fx / (spot_fx * coef_lin + coefc);		

		premium_tgt = srt_f_optblksch(	spot_fx,
										strike,
										volatility[0],
										exercise[0],
										df,
										SRT_CALL,
										PREMIUM
										);

		err = srt_f_optimpvol(		
							premium_tgt * coef_linab,
							shift_spot,
							shift_strike,
							exercise[0],								
							df,
							is_put,
							SRT_LOGNORMAL,
							&sig0);

		sig0 /= coef_linab;
		vol2 = sig0;

		sig0 = op_sabrgen_calib(spot_fx,
								strike,
								exercise[0],
								volatility[0],
								alpha_eq,
								coefa,
								coefb,
								coefc,								
								rho,
								vol_quadra);

		volQuad2 = sig0;

		if (guess)
		{
			if (guess[0] != 0.0)
			{
				sig0 = guess[0];
			}
		}

		newVar = sig0 * sig0 * exercise[0];

		eval_evt[nstpt-1] = 1;
		func_parm[nstpt-1][0] = strike;

		err = FxSabrQuad_adi(nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							sig0,
							drift_adi,
							alpha,
							coefa,
							coefb,
							coefc,
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
							&premium,
							0,
							NULL,
							1,
							0,
							premium_tgt,
							&newSig0,
							&premium2
							);

		if (err)
		{
			goto FREE_RETURN;
		}

		error = fabs(premium - premium_tgt);
		err = srt_f_optimpvol(		
								premium,
								spot_fx,
								strike,
								exercise[0],						
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								&bsvol);

		errorv = fabs(bsvol - volatility[0]);

		if (errorv >= precision || err)
		{
			if (err || newSig0 < 0.0 || fabs(alpha) < 0.001)
			{
				if (err)
				{
					newSig0 = sig0 * 0.5;
				}
				else
				{
					newSig0 = sqrt(premium_tgt * premium_tgt / premium / premium * newVar / exercise[0]);
				}
			}
			
			if (!err)
			{
				res_iter[0][0] = sig0;
				res_iter[0][1] = premium;
				nb_iter++;
			}
		}
		else
		{
			newSig0 = sig0;
		}

		/* Now do a Newton algorithm	*/

		k = 0;
		while ((k < nbIter) && (errorv >= precision))
		{
			err = FxSabrQuad_adi(
							nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							newSig0,
							drift_adi,
							alpha,
							coefa,
							coefb,
							coefc,
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
							&premium2,
							0,
							NULL,
							0,
							0,
							0,
							NULL,
							NULL);
			if (err)
			{
				goto FREE_RETURN;
			}

			err = srt_f_optimpvol(		
								premium2,
								spot_fx,
								strike,
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
				/*
				newSig0 = sig0;
				*/

				vega = (premium2 - premium) / (newSig0 - sig0);
				if (vega < 0)
				{
					/* problem !!! */
				}
				else
				{
					l = 0;
					while (l < nb_iter && res_iter[l][0] < newSig0)
					{
						l++;
					}

					for (j=nb_iter-1; j>=l; j--)
					{
						res_iter[j+1][0] = res_iter[j][0];
						res_iter[j+1][1] = res_iter[j][1];
					}

					res_iter[l][0] = newSig0;
					res_iter[l][1] = premium2;
					nb_iter++;

					sig0 = newSig0;
					premium = premium2;

					newSig0 = solve_for_next_coef(	res_iter,
													nb_iter,
													premium_tgt,
													1);

					if (newSig0 < 0.0)
					{
						newSig0 = res_iter[0][0] * 0.25;
					}
				}
			}
			
			k++;
		}

		sig0 = newSig0;
		sigma[0] = sig0;

		if (errorv <= precision)
		{
			smessage("Option %d: Success at iteration %d", 1, k);
		}
		else
		{
			smessage("Option %d: May have failed", 1);
		}
	}
	else
	{
		sig0 = sigma[0];
	}	

	if (time)
	{
		free (time);
		time = NULL;
	}
	
	drift_res[0] = 0.0;

	var = sig0 * sig0 * exercise[0];

	/* First update var */
	for (i=1; i<start_opt; i++)
	{
		dt = exercise[i] - exercise[i-1];
		coef = log(sigma[i] / sigma[i-1]) / dt;
		var += sigma[i-1] * sigma[i-1] / (2.0 * coef) * (exp(2.0 * coef * dt) - 1.0);
		drift_res[i] = coef;
	}

	/* Update last calculation of vol */

	if (start_opt > 0)
	{
		alpha_eq = alpha;

		if (fabs(lambda) > 1.0E-08)
		{
			alpha_eq *= sqrt((1.0 - exp(-2.0 * lambda * exercise[start_opt-1])) / (2.0 * lambda * exercise[start_opt-1]));
		}

		if (is_straddle)
		{
			strike = spot_fx * exp(-0.5 * volatility[start_opt-1] * volatility[start_opt-1] * exercise[start_opt-1]);
			shift_strike = coefa * strike + coefb;
		}

		df = swp_f_df(today, (long) (today + exercise[start_opt-1] * 365.0 + 1.0E-08), dom_yc);

		premium_tgt = srt_f_optblksch(	spot_fx,
										strike,
										volatility[start_opt-1],
										exercise[start_opt-1],
										df,
										SRT_CALL,
										PREMIUM
										);

		err = srt_f_optimpvol(	premium_tgt * coef_linab,
								shift_spot,
								shift_strike,
								exercise[start_opt-1],								
								df,
								is_put,
								SRT_LOGNORMAL,
								&vol2);

		vol2 /= coef_linab;

		volQuad2 = op_sabrgen_calib(spot_fx,
									strike,
									exercise[start_opt-1],
									volatility[start_opt-1],
									alpha_eq,
									coefa,
									coefb,
									coefc,								
									rho,
									vol_quadra);		
	}

	/*	Next options to be calibrated */
	for (i=max(start_opt, 1); i<nbOpt; i++)
	{
		/* Initialisation	*/

		alpha_eq = alpha;

		if (fabs(lambda) > 1.0E-08)
		{
			alpha_eq *= sqrt((1.0 - exp(-2.0 * lambda * exercise[i])) / (2.0 * lambda * exercise[i]));
		}

		vol1 = vol2;
		volQuad1 = volQuad2;

		for (j=0; j<=nbIter; j++)
		{
			res_iter[j][0] = 0.0;
			res_iter[j][1] = 0.0;
		}

		nb_iter = 0;

		if (is_straddle)
		{
			strike = spot_fx * exp(-0.5 * volatility[i] * volatility[i] * exercise[i]);
			shift_strike = coefa * strike + coefb;
		}

		/* First construct the time steps */

		dt = exercise[i] - exercise[i-1];
		
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

		df = swp_f_df(today, (long) (today + exercise[i] * 365.0 + 1.0E-08), dom_yc);
		
		premium_tgt = srt_f_optblksch(	spot_fx,
										strike,
										volatility[i],
										exercise[i],
										df,
										SRT_CALL,
										PREMIUM
										);

		/* Calculation of the first guess */
		
		var_diff = volatility[i-1] * spot_fx / (spot_fx * (coefa * spot_fx + coefb) + coefc);
		vol2 = volatility[i] * spot_fx / (spot_fx * (coefa * spot_fx + coefb) + coefc);
		var_diff = (vol2 * vol2 * exercise[i] - var_diff * var_diff * exercise[i-1]);

		/* Other guess of var_diff */

		err = srt_f_optimpvol(	premium_tgt * coef_linab,
								shift_spot,
								shift_strike,
								exercise[i],								
								df,
								is_put,
								SRT_LOGNORMAL,
								&vol2);

		vol2 /= coef_linab;

		volQuad2 = op_sabrgen_calib(spot_fx,
									strike,
									exercise[i],
									volatility[i],
									alpha,
									coefa,
									coefb,
									coefc,
									rho,
									vol_quadra);

		var_diff = (volQuad2 * volQuad2 * exercise[i] - volQuad1 * volQuad1 * exercise[i-1]);
		
		if (var_diff < 0.0)
		{
			var_diff = (vol2 * vol2 * exercise[i] - vol1 * vol1 * exercise[i-1]);
		}

		if (var_diff < 0.0)
		{
			coef = 0.0;
		}
		else
		{						
			coef = solve_var_diff(dt, var_diff / (sigma[i-1] * sigma[i-1]));
		}

		/* Look if there is a guess provided */
		if (guess)
		{
			if (guess[i] != 0.0)
			{
				coef = log(guess[i] / guess[i-1]) / (exercise[i] - exercise[i-1]);
			}
		}
		
		/*	Fill term Structure	*/
		for (j=old_index; j<nstpt; j++)
		{
			drift_adi[j] = coef;
		}
		
		newVar = var + sigma[i-1] * sigma[i-1] / (2.0 * coef) * (exp(2.0 * coef * dt) - 1.0);
				
		
		/*	Adjustment in the PDE	*/		

		/*	Fill the param for the payoff funtion	*/		

		eval_evt[nstpt-1] = 1;
		func_parm[nstpt-1][0] = strike;
		
		j = 0;

		err = FxSabrQuad_adi(nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							sig0,
							drift_adi,
							alpha,
							coefa,
							coefb,
							coefc,
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
							&premium,
							0,
							NULL,
							1,
							0,
							premium_tgt,
							&newSig0,
							&temp
							);
				
		if (err)
		{
			goto FREE_RETURN;
		}

		error = fabs(premium - premium_tgt);
		err = srt_f_optimpvol(		
								premium,
								spot_fx,
								strike,
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
				
				if (newSig0 > 0 && fabs(alpha > 0.001))
				{
					newVar2 = newSig0 * newSig0 / sig0 / sig0 * newVar;

					if (newVar2 < var)
					{
						newVar2 = premium_tgt * premium_tgt / premium / premium * newVar;
					}
				}
				else
				{
					newVar2 = premium_tgt * premium_tgt / premium / premium * newVar;
				}
				
				dt = exercise[i] - exercise[i-1];
				
				if (newVar2 < var)
				{
					/* we have a problem */					
					newCoef = coef - 0.2;
				}
				else
				{
					newCoef = solve_var_diff(dt, (newVar2 - var) / (sigma[i-1] * sigma[i-1]));
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
			err = FxSabrQuad_adi(
							nstpt,
							time,
							date,
							nstpfx,
							nstpvol,
							sig0,
							drift_adi,
							alpha,
							coefa,
							coefb,
							coefc,
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
							&premium2,
							0,
							NULL,
							0,
							0,
							0,
							NULL,
							NULL
							);

			if (err)
			{
				goto FREE_RETURN;
			}

			err = srt_f_optimpvol(		
								premium2,
								spot_fx,
								strike,
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


				vega = (premium2 - premium) / (newCoef - coef);
				error = fabs(premium2 - premium_tgt);
				
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
					newCoef2 = solve_for_next_coef(	res_iter,
													nb_iter,
													premium_tgt,
													1);

					coef = newCoef;

					/*
					newCoef = coef - (premium2 - premium_tgt) / vega;
					*/

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
		

		sigma[i] = sigma[i-1] * exp(newCoef * dt);
		
		var += sigma[i-1] * sigma[i-1] / (2.0 * newCoef) * (exp(2.0 * newCoef * dt) - 1.0);
		
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


static Err GetStrikesVol(	double	fwd_fx,
							double	mat_smile,						
							double	df_dom,
							double	df_for,
							double	atm_vol,
							double	risk_reversal,
							double	butterfly,
							double	*atm_strike,
							double	*strike1,
							double	*vol1,
							double	*strike2,
							double	*vol2,
							double	*alpha,
							double	beta,
							double	*rho,
							double	lambda,
							int		is_straddle,
							int		use_input,
							double	*fitting_error)

{
double	temp;
double	*strike_sabr	= NULL,
		*vol_sabr		= NULL;

Err		err = NULL;
	
	if (is_straddle)
	{
		*atm_strike = fwd_fx * exp(-0.5 * atm_vol * atm_vol * mat_smile);
	}
	else
	{
		*atm_strike = fwd_fx;
	}

	if (*vol1 < 1.0E-08 || *vol2< 1.0E-08 || *strike1 < 1.0E-08 || *strike2 < 1.0E-08)
	{

		/* 25 PUT */
		*vol1 = butterfly + atm_vol - 0.5 * risk_reversal;
		//*strike1 = SolveForKPutPercent(-0.25, fwd_fx, *vol1 * sqrt(mat_smile), df_dom, df_for);

		/* 25 CALL */
		*vol2 = *vol1 + risk_reversal;
		//*strike2 = SolveForKCallOutPercent(0.25, fwd_fx, *vol2 * sqrt(mat_smile), df_dom, df_for);
	}
	else
	{
		if (*strike1 > *strike2)
		{
			temp = *strike1;
			*strike1 = *strike2;
			*strike2 = temp;

			temp = *vol1;
			*vol1 = *vol2;			
			*vol2 = temp;
		}
	}


	/* if alpha and rho are not provided we can do a first guess */
	if (!use_input)
	{
		/* Use SABR calibration */

		strike_sabr = dvector(1, 3);
		vol_sabr = dvector(1, 3);

		if (!strike_sabr || !vol_sabr)
		{
			err = "Memory allocation error in GetStrikesVol";
			goto FREE_RETURN;
		}		

		if (*atm_strike < *strike1)
		{
			strike_sabr[2] = *strike1;
			strike_sabr[3] = *strike2;
			strike_sabr[1] = *atm_strike;

			vol_sabr[2] = *vol1;
			vol_sabr[3] = *vol2;
			vol_sabr[1] = atm_vol;
		}
		else if (*atm_strike > *strike2)
		{
			strike_sabr[1] = *strike1;
			strike_sabr[2] = *strike2;
			strike_sabr[3] = *atm_strike;

			vol_sabr[1] = *vol1;
			vol_sabr[2] = *vol2;
			vol_sabr[3] = atm_vol;
		}
		else
		{
			strike_sabr[1] = *strike1;
			strike_sabr[3] = *strike2;
			strike_sabr[2] = *atm_strike;

			vol_sabr[1] = *vol1;
			vol_sabr[3] = *vol2;
			vol_sabr[2] = atm_vol;
		}

		*alpha = 0.05;
		*rho = 0.0;

		err = opsabrcalib(	fwd_fx,
							mat_smile,
							3,
							strike_sabr,
							vol_sabr,
							&atm_vol,
							alpha,
							1,
							&beta,
							0,
							rho,
							1,
							fitting_error);

		if (err)
		{
			*alpha = 0.05;
			*rho = 0.0;
		}
		else
		{
			if (fabs(lambda) > 1.0E-10)
			{
				/* adjust for mean- reversion */
				*alpha /= sqrt((1.0 - exp(-2.0 * lambda * mat_smile)) / (2.0 * lambda * mat_smile));
			}
		}
	}

FREE_RETURN:

	if (strike_sabr) free_dvector(strike_sabr, 1, 3);
	if (vol_sabr) free_dvector(vol_sabr, 1, 3);

	return err;
}

static Err	FxSabrQuadCalibSmilePricer(
							 /* Model params */
							double	today,
							double	spot_fx,
							double	fwd_fx,
							double	alpha,						
							double	a,
							double	b,
							double	c,
							double	rho,						
							double	lambda,
							double	floormu,

							char	*dom_yc,
							char	*for_yc,

							double	*sigma_time_fx,
							double	*sigma_fx,
							int		sigma_n_fx,

							/* Discretisation params */
							double	*time,
							double	*date,
							double	*drift,
							int		nstp,

							double	**func_parm,
							int		*is_event,

							long	nstpt_tgt,
							int		nstpfx,
							int		nstpvol,

							/* Product infos */
							double	df,
							double	mat_smile,
							double	strike1,
							double	strike2,

							double	*res,
							double	*impvol1,
							double	*impvol2)

{
int		i, index;
double	coef;
Err		err = NULL;

	/* fill the drift vector */	
	i = 0;
	index = 0;
	while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
	{
		drift[i] = 0.0;
		i++;
	}
	
	while (i < nstp)
	{
		if (index < sigma_n_fx-1)
		{			
			coef = log(sigma_fx[index+1] / sigma_fx[index]) / (sigma_time_fx[index+1] - sigma_time_fx[index]);
			index++;
			while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
			{
				drift[i] = coef;
				i++;
			}
		}
		else
		{
			coef = 0.0;
			while (i < nstp)
			{
				drift[i] = coef;
				i++;
			}
		}
	}	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;			
	}

	/* pricing */

	err = FxSabrQuad_adi(nstp,
						time,
						date,
						nstpfx,
						nstpvol,
						sigma_fx[0],
						drift,
						alpha,
						a,
						b,
						c,
						rho,
						lambda,
						floormu,
						func_parm,
						is_event,
						NULL,
						NULL,
						NULL,
						spot_fx,
						dom_yc,
						for_yc,
						payoff_fx_sabr_adi_opt,
						2,
						res,
						0,
						NULL,
						0,
						0,
						0,
						NULL,
						NULL
						);

	if (err)
	{
		return err;
	}

	/* get the corresponding imp vol */

	err = srt_f_optimpvol(		
							res[0],
							fwd_fx,
							strike1,
							mat_smile,								
							df,
							SRT_CALL,
							SRT_LOGNORMAL,
							impvol1);

	err = srt_f_optimpvol(		
							res[1],
							fwd_fx,
							strike2,
							mat_smile,								
							df,
							SRT_CALL,
							SRT_LOGNORMAL,
							impvol2);

	return err;
}

Err	FxSabrQuadCalibSmile(
						/*	Underlying						*/
						char	*dom_yc,
						char	*for_yc,
						long	today,
						double	spot_fx,

						/*	Model Parameters				*/						
						double	*alpha_out,
						double	gamma,
						double	beta,						
						double	*rho_out,
						double	lambda,
						double	floormu,

						/*	User input of starting points	*/
						int		calib_smile,
						int		use_input,
						int		use_total_ts,

						/*	Options Parameters				*/
						double	*exercise,
						double	*maturity,
						double	*volatility,
						int		nbOpt,
						int		is_straddle,

						long	mat_date,
						double	risk_reversal,
						double	butterfly,						

						/*	Discretisation Parameters		*/
						long	nstpt_tgt,
						int		nstpfx,
						int		nstpvol,
						int		nbIterATM,
						double	precisionATM,
						int		nbIterSmile,
						double	precisionSmile,

						/*	Result							*/
						double	*strike1_out,
						double	*vol1_out,
						double	*strike2_out,
						double	*vol2_out,
						
						double	*sigma		
						)
{
double	coefa, coefb;
double	alpha1, alpha2, alpha3, rho1, rho2, rho3, 
		impvol11, impvol21, impvol12, impvol22, impvol13, impvol23;
double	sensalpha1, sensalpha2, sensrho1, sensrho2, delta, error1, error2, delta_alpha, delta_rho, fitting_error;
double	shift_spot;

double	mat_smile, setlt_smile, fwd_fx, df, df_for, spot_fx_cash;
long	exe_date, spot_date;

double	*time			= NULL,
		*date			= NULL,
		*drift			= NULL,
		**func_parm		= NULL,
		*res			= NULL,
		*sigma_alpha	= NULL,
		*sigma_rho		= NULL,
		*guess			= NULL;

double	a, b, c;

double	*exercise_cal, *maturity_cal, *volatility_cal;

double	bt_theo, rr_theo, bt, rr;

int		*is_event	= NULL;
int		start_opt;

double	strike1, vol1, strike2, vol2, atm_strike, atm_vol, coef1, coef2, coef3;
int		i, nstp, k;

int		nb_calib;

double	error;

Err		err = NULL;
time_t	t1, t2;

	t1 = clock();

		
	smessage("Starting Calibration of Fx Sabr Smile");	

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	spot_fx_cash = spot_fx * swp_f_df (today, spot_date, dom_yc) / swp_f_df (today, spot_date, for_yc);

	/* Transformation from beta and convexity to b and c */
	SabrQuadGetParams(spot_fx_cash, gamma, beta, &a, &b, &c);

	coefa = a * spot_fx + b;
	coefb = c;

	shift_spot = coefa * spot_fx_cash + coefb;

	exe_date = add_unit (mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	mat_smile = (exe_date - today) * YEARS_IN_DAY;
	setlt_smile = (mat_date - today) * YEARS_IN_DAY;

	fwd_fx = spot_fx_cash * swp_f_df (today, mat_date, for_yc)
								/ swp_f_df (today, mat_date, dom_yc);
	
	df = swp_f_df(today, exe_date, dom_yc);
	df_for = swp_f_df(today, exe_date, for_yc);

	alpha1 = *alpha_out;
	rho1 = *rho_out;
	strike1 = *strike1_out;
	strike2 = *strike2_out;
	vol1 = *vol1_out;
	vol2 = *vol2_out;

	start_opt = 0;

	if (calib_smile)
	{
		/* First find the corresponding strikes */
		i = 0;
		while (fabs(exercise[i] - mat_smile) > 1.0E-08 && i<nbOpt)
		{
			i++;
		}

		if (i == nbOpt)
		{
			err = serror("Maturity %d of calibration is not present in the market ATM term struct", mat_date);
			goto FREE_RETURN;
		}

		start_opt = min(i + 1, nbOpt);

		atm_vol = volatility[i];	

		err = GetStrikesVol(fwd_fx,
							mat_smile,
							df,
							df_for,
							atm_vol,
							risk_reversal,
							butterfly,
							&atm_strike,
							&strike1,
							&vol1,
							&strike2,
							&vol2,
							&alpha1,
							0,
							&rho1,
							lambda,
							is_straddle,
							use_input,
							&fitting_error);
		

		if (err)
		{
			goto FREE_RETURN;
		}

		*strike1_out = strike1;
		*vol1_out = vol1;
		*strike2_out = strike2;
		*vol2_out = vol2;

		if (use_total_ts)
		{
			nb_calib = start_opt;
			exercise_cal = exercise;
			maturity_cal = maturity;
			volatility_cal = volatility;
		}
		else
		{
			nb_calib = 1;
			exercise_cal = &(exercise[start_opt-1]);
			maturity_cal = &(maturity[start_opt-1]);
			volatility_cal = &(volatility[start_opt-1]);
		}
		

		/* discretise in time			*/		
		nstp = 2;
		time = (double*) calloc (nstp, sizeof (double));
		if (!time)
		{
			err = "Memory allocation error (1) in FxSabrQuadCalibSmile";
			goto FREE_RETURN;
		}

		time[0] = 0.0;
		time[1] = mat_smile;
		
		/*	Add revelant vol times	*/		
		for (i=0; i<nb_calib; i++)
		{
			num_f_add_number (&nstp, &time, exercise_cal[i]);
		}		


		num_f_sort_vector(nstp, time);
		num_f_unique_vector(&nstp, time);

		/*	Fill the time vector */
		err = fill_time_vector (	&time, 
									&nstp,
									0,
									NULL,
									0,
									NULL, 
									nstpt_tgt);				
		if (err)
		{
			goto FREE_RETURN;
		}	

		date = dvector (0, nstp - 1);
		drift = dvector (0, nstp - 1);
		func_parm = dmatrix(0, nstp-1, 0, 1);
		res = dvector(0, 1);
		is_event = calloc(nstp, sizeof(int));

		sigma_alpha = dvector(0, nbOpt-1);
		sigma_rho = dvector(0, nbOpt-1);
		guess = dvector(0, nbOpt-1);
			
		if (!date || !drift || !func_parm || !res || !sigma_alpha || !sigma_rho || !guess)
		{
			err = "Memory allocation failure (2) in FxSabrQuadCalibSmile";
			goto FREE_RETURN;
		}	
		
		func_parm[nstp-1][0] = strike1;
		func_parm[nstp-1][1] = strike2;	
		is_event[nstp-1] = 1;	

		/* first run */

		/* calibration */
		err = FxSabrQuadCalibration(dom_yc, for_yc, today, spot_fx, alpha1, a, b, c, rho1, lambda, floormu,
									exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		/* pricer */
		err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a, b, c, rho1, lambda, floormu,
										dom_yc, for_yc, exercise_cal, sigma, nb_calib,
										time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
										df, mat_smile, strike1, strike2, res, &impvol11, &impvol21);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		error1 = impvol11 - vol1;
		error2 = impvol21 - vol2;

		error = max(fabs(error1), fabs(error2));

		if (error > precisionSmile)
		{
			/*	First guess failed but to have sensitivity to rho
				alpha must be gretear than 0.0						*/
			if ((use_input || (!use_input && fabs(fitting_error) > 0.025)) && alpha1 < 0.05)
			{
				alpha1 = 0.05;
			

				/* calibration */
				err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a, b, c, rho1, lambda, floormu,
											exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
											nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
			
				if (err)
				{
					goto FREE_RETURN;
				}

				/* pricer */
				err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a, b, c, rho1, lambda, floormu,
												dom_yc, for_yc, exercise_cal, sigma, nb_calib,
												time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
												df, mat_smile, strike1, strike2, res, &impvol11, &impvol21);
				
				if (err)
				{
					goto FREE_RETURN;
				}

				error1 = impvol11 - vol1;
				error2 = impvol21 - vol2;

				error = max(fabs(error1), fabs(error2));
			}
		}

		/* we need to adjust the alpha and rho */
		/* it is a very stupid algorithm at the moment */

		k = 0;

		delta_alpha = EPSALPHA;
		delta_rho = EPSRHO;

		bt_theo = 0.5 * (vol1 + vol2) - atm_vol;
		rr_theo = vol2 - vol1;
		
		bt = 0.5 * (impvol11 + impvol21) - atm_vol;						
		rr = impvol21 - impvol11;

		if (bt_theo < bt)
		{
			delta_alpha *= -1.0;
		}

		if (rr_theo < rr)
		{
			delta_rho *= -1.0;
		}
		
		while ((k < nbIterSmile) && (error >= precisionSmile))
		{

			/* Compute sensi */

			if (alpha1 + delta_alpha > 0.0)
			{
				alpha2  = alpha1 + delta_alpha;
			}
			else
			{
				alpha2  = alpha1 / 2.0;
			}

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha2, a, b, c, rho1, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle, 
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, sigma, sigma_alpha);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha2, a, b, c, rho1, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma_alpha, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol12, &impvol22);

			sensalpha1 = (impvol12 - impvol11) / (alpha2 - alpha1);
			sensalpha2 = (impvol22 - impvol21) / (alpha2 - alpha1);

			if (fabs(rho1 + delta_rho) < 0.99)
			{
				rho2  = rho1 + delta_rho;
			}
			else
			{
				if (rho1 + delta_rho > 0.99)
				{
					rho2 = 0.5 * (rho1 + 1.0);
				}
				else
				{
					rho2 = 0.5 * (rho1 - 1.0);
				}
			}

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a, b, c, rho2, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, sigma, sigma_rho);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a, b, c, rho2, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma_rho, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol12, &impvol22);

			sensrho1 = (impvol12 - impvol11) / (rho2 - rho1);
			sensrho2 = (impvol22 - impvol21) / (rho2 - rho1);

			if (fabs(sensrho1) < 10E-06 || fabs(sensrho2) < 10E-06)
			{
				/*	No more sensitivity to rho, it means that alpha is too small */

				if (error < 5.0 * precisionSmile)
				{
					/* we stop here before the algorithm crashs */
					smessage("Nb Iter %d: Smile Calib error %.0f bp", k, error * 10000);
					*alpha_out = alpha1;
					*rho_out = rho1;

					goto FREE_RETURN;
				}
				else
				{
					/* try to put some alpha */
					alpha3 += 0.05;
				}
			}
			else
			{
				/* solve the matrix */

				delta = sensalpha1 * sensrho2 - sensalpha2 * sensrho1;

				if (fabs(delta) < 1.0E-10)
				{
					err = "cannot calibrate";
					goto FREE_RETURN;
				}

				alpha3 = alpha1 - (sensrho2 * error1 - sensrho1 * error2) / delta; 
				rho3 = rho1 - (-sensalpha2 * error1 + sensalpha1 * error2) / delta;

				while (alpha3 < 0 || fabs(rho3) > 0.9999)
				{
					if (alpha3 < 0)
					{
						coef1 = (alpha1 - alpha1 * 0.5) / ((sensrho2 * error1 - sensrho1 * error2) / delta);
					}

					if (rho3 < -0.9999)
					{
						coef2 = (rho1 + 0.99) / ((-sensalpha2 * error1 + sensalpha1 * error2) / delta);
					}
					else if (rho3 > 0.9999)
					{
						coef2 = (rho1 - 0.99) / ((-sensalpha2 * error1 + sensalpha1 * error2) / delta);
					}

					coef3 = min(coef1, coef2);

					alpha3 = alpha1 - coef3 * (sensrho2 * error1 - sensrho1 * error2) / delta; 
					rho3 = rho1 - coef3 * (-sensalpha2 * error1 + sensalpha1 * error2) / delta;
				}				
			}

			/* Fill the guess */
			for (i=0; i<nb_calib; i++)
			{
				guess[i] = sigma[i] 
						+ (alpha3 - alpha1) * (sigma_alpha[i] - sigma[i]) / (alpha2 - alpha1)
						+ (rho3 - rho1) * (sigma_rho[i] - sigma[i]) / (rho2 - rho1);
			}

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha3, a, b, c, rho3, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, guess, sigma);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha3, a, b, c, rho3, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol13, &impvol23);

			error1 = impvol13 - vol1;
			error2 = impvol23 - vol2;

			error = max(fabs(error1), fabs(error2));

			/* Update alpha, rho, implied vols and the shift */

			delta_alpha = max(min(fabs(alpha3 - alpha1), EPSALPHAMAX), EPSALPHAMIN);
			delta_rho = max(min(fabs(rho3 - rho1), EPSRHOMAX), EPSRHOMIN);

			bt = 0.5 * (impvol13 + impvol23) - atm_vol;						
			rr = impvol23 - impvol13;

			if (bt_theo < bt)
			{
				delta_alpha *= -1.0;
			}

			if (rr_theo < rr)
			{
				delta_rho *= -1.0;
			}

			alpha1 = alpha3;
			rho1 = rho3;

			impvol11 = impvol13;
			impvol21 = impvol23;

			k++;
		}

		smessage("Nb Iter %d: Smile Calib error %.0f bp", k, error * 10000);
	}

	/* calibration of the last sigma */

	if (use_total_ts || nbOpt == 1)
	{
		err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a, b, c, rho1, lambda, floormu,
									exercise, maturity, volatility, nbOpt, start_opt, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
	}
	else
	{
		err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a, b, c, rho1, lambda, floormu,
									exercise, maturity, volatility, nbOpt, 0, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
	}


	*alpha_out = alpha1;
	*rho_out = rho1;


FREE_RETURN:

	t2 = clock();

	if (calib_smile)
	{
		smessage ("Total calibration time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	}

	if (time) free (time);
	if (date) free_dvector (date, 0, nstp - 1);
	if (drift)free_dvector (drift, 0, nstp - 1);
	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, 1);
	if (res) free_dvector(res, 0, 1);
	if (is_event) free (is_event);

	if (sigma_alpha) free_dvector(sigma_alpha, 0, nbOpt-1);
	if (sigma_rho) free_dvector(sigma_rho, 0, nbOpt-1);
	if (guess) free_dvector(guess, 0, nbOpt-1);	

	return err;
}

Err	FxSabrQuadCalibSmile2(
						/*	Underlying						*/
						char	*dom_yc,
						char	*for_yc,
						long	today,
						double	spot_fx,

						/*	Model Parameters				*/						
						double	*alpha_out,
						double	*gamma_out,
						double	*beta_out,						
						double	*rho_out,
						double	lambda,
						double	floormu,

						/*	User input of starting points	*/
						int		calib_smile,
						int		use_input,
						int		use_total_ts,
						int		use_alpharho,

						/*	Options Parameters				*/
						double	*exercise,
						double	*maturity,
						double	*volatility,
						int		nbOpt,
						int		is_straddle,

						long	mat_date,
						double	risk_reversal,
						double	butterfly,						

						/*	Discretisation Parameters		*/
						long	nstpt_tgt,
						int		nstpfx,
						int		nstpvol,
						int		nbIterATM,
						double	precisionATM,
						int		nbIterSmile,
						double	precisionSmile,

						/*	Result							*/
						double	*strike1_out,
						double	*vol1_out,
						double	*strike2_out,
						double	*vol2_out,
						
						double	*sigma		
						)
{
double	alpha1, alpha2, alpha3, rho1, rho2, rho3, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
		a1, a2, a3, b1, b2, b3, c1, c2, c3, 
		impvol11, impvol21, impvol12, impvol22, impvol13, impvol23;
double	spotbeta1, spotbeta2, spotbeta3;
double	sensbt1, sensbt2, sensrr1, sensrr2, delta, error1, error2, delta_bt, delta_rr, fitting_error;

double	mat_smile, setlt_smile, fwd_fx, df, df_for, spot_fx_cash;
long	exe_date, spot_date;

double	*time			= NULL,
		*date			= NULL,
		*drift			= NULL,
		**func_parm		= NULL,
		*res			= NULL,
		*sigma_alpha	= NULL,
		*sigma_rho		= NULL,
		*guess			= NULL;

double	*param_RR1, *param_BT1;
double	*param_RR2, *param_BT2;
double	*param_RR3, *param_BT3;

double	*exercise_cal, *maturity_cal, *volatility_cal;

double	bt_theo, rr_theo, bt, rr;

int		*is_event	= NULL;
int		start_opt;

double	strike1, vol1, strike2, vol2, atm_strike, atm_vol, coef1, coef2, coef3;
int		i, nstp, k;

int		nb_calib;

double	error;

Err		err = NULL;
time_t	t1, t2;

	t1 = clock();

		
	smessage("Starting Calibration of Fx Sabr Smile");	

	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	spot_fx_cash = spot_fx * swp_f_df (today, spot_date, dom_yc) / swp_f_df (today, spot_date, for_yc);	

	exe_date = add_unit (mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	mat_smile = (exe_date - today) * YEARS_IN_DAY;
	setlt_smile = (mat_date - today) * YEARS_IN_DAY;

	fwd_fx = spot_fx_cash * swp_f_df (today, mat_date, for_yc)
								/ swp_f_df (today, mat_date, dom_yc);
	
	df = swp_f_df(today, exe_date, dom_yc);
	df_for = swp_f_df(today, exe_date, for_yc);

	alpha1 = *alpha_out;
	rho1 = *rho_out;
	gamma1 = *gamma_out;
	beta1 = *beta_out;

	SabrQuadGetParams(spot_fx_cash, gamma1, beta1, &a1, &b1, &c1);
	
	strike1 = *strike1_out;
	strike2 = *strike2_out;
	vol1 = *vol1_out;
	vol2 = *vol2_out;

	start_opt = 0;

	if (calib_smile)
	{
		/* First find the corresponding strikes */
		i = 0;
		while (fabs(exercise[i] - mat_smile) > 1.0E-08 && i<nbOpt)
		{
			i++;
		}

		if (i == nbOpt)
		{
			err = serror("Maturity %d of calibration is not present in the market ATM term struct", mat_date);
			goto FREE_RETURN;
		}

		start_opt = min(i + 1, nbOpt);

		atm_vol = volatility[i];	

		err = GetStrikesVol(fwd_fx,
							mat_smile,
							df,
							df_for,
							atm_vol,
							risk_reversal,
							butterfly,
							&atm_strike,
							&strike1,
							&vol1,
							&strike2,
							&vol2,
							&alpha1,
							0,
							&rho1,
							lambda,
							is_straddle,
							use_input,
							&fitting_error);
		

		if (err)
		{
			goto FREE_RETURN;
		}

		*strike1_out = strike1;
		*vol1_out = vol1;
		*strike2_out = strike2;
		*vol2_out = vol2;

		if (use_total_ts)
		{
			nb_calib = start_opt;
			exercise_cal = exercise;
			maturity_cal = maturity;
			volatility_cal = volatility;
		}
		else
		{
			nb_calib = 1;
			exercise_cal = &(exercise[start_opt-1]);
			maturity_cal = &(maturity[start_opt-1]);
			volatility_cal = &(volatility[start_opt-1]);
		}
		

		/* discretise in time			*/		
		nstp = 2;
		time = (double*) calloc (nstp, sizeof (double));
		if (!time)
		{
			err = "Memory allocation error (1) in FxSabrQuadCalibSmile";
			goto FREE_RETURN;
		}

		time[0] = 0.0;
		time[1] = mat_smile;
		
		/*	Add revelant vol times	*/		
		for (i=0; i<nb_calib; i++)
		{
			num_f_add_number (&nstp, &time, exercise_cal[i]);
		}		


		num_f_sort_vector(nstp, time);
		num_f_unique_vector(&nstp, time);

		/*	Fill the time vector */
		err = fill_time_vector (	&time, 
									&nstp,
									0,
									NULL,
									0,
									NULL, 
									nstpt_tgt);				
		if (err)
		{
			goto FREE_RETURN;
		}	

		date = dvector (0, nstp - 1);
		drift = dvector (0, nstp - 1);
		func_parm = dmatrix(0, nstp-1, 0, 1);
		res = dvector(0, 1);
		is_event = calloc(nstp, sizeof(int));

		sigma_alpha = dvector(0, nbOpt-1);
		sigma_rho = dvector(0, nbOpt-1);
		guess = dvector(0, nbOpt-1);
			
		if (!date || !drift || !func_parm || !res || !sigma_alpha || !sigma_rho || !guess)
		{
			err = "Memory allocation failure (2) in FxSabrQuadCalibSmile";
			goto FREE_RETURN;
		}	
		
		func_parm[nstp-1][0] = strike1;
		func_parm[nstp-1][1] = strike2;	
		is_event[nstp-1] = 1;	

		/* first run */

		/* calibration */
		err = FxSabrQuadCalibration(dom_yc, for_yc, today, spot_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
									exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		/* pricer */
		err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
										dom_yc, for_yc, exercise_cal, sigma, nb_calib,
										time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
										df, mat_smile, strike1, strike2, res, &impvol11, &impvol21);
		
		if (err)
		{
			goto FREE_RETURN;
		}

		error1 = impvol11 - vol1;
		error2 = impvol21 - vol2;

		error = max(fabs(error1), fabs(error2));

		if (error > precisionSmile)
		{
			/*	First guess failed but to have sensitivity to rho
				alpha must be gretear than 0.0						*/
			if ((use_input || (!use_input && fabs(fitting_error) > 0.025)) && ((use_alpharho == 1 && alpha1 < 0.05) || (use_alpharho == 0 && gamma1 < 0.0005)))
			{
				if (use_alpharho)
				{
					alpha1 = 0.05;
				}
				else
				{
					gamma1 = 0.005;
					SabrQuadGetParams(spot_fx_cash, gamma1, beta1, &a1, &b1, &c1);
				}

				/* calibration */
				err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
											exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
											nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
			
				if (err)
				{
					goto FREE_RETURN;
				}

				/* pricer */
				err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
												dom_yc, for_yc, exercise_cal, sigma, nb_calib,
												time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
												df, mat_smile, strike1, strike2, res, &impvol11, &impvol21);
				
				if (err)
				{
					goto FREE_RETURN;
				}

				error1 = impvol11 - vol1;
				error2 = impvol21 - vol2;

				error = max(fabs(error1), fabs(error2));
			}
		}
		
		/* set the calibration variables */

		if (use_alpharho)
		{
			param_RR1 = &rho1;
			param_BT1 = &alpha1;
			param_RR2 = &rho2;
			param_BT2 = &alpha2;
			param_RR3 = &rho3;
			param_BT3 = &alpha3;
			delta_bt = EPSALPHA;
			delta_rr = EPSRHO;
			
			beta2 = beta1;
			beta3 = beta1;
			gamma2 = gamma1;
			gamma3 = gamma1;
		}
		else
		{
			param_RR1 = &beta1;
			param_BT1 = &gamma1;
			param_RR2 = &beta2;
			param_BT2 = &gamma2;
			param_RR3 = &beta3;
			param_BT3 = &gamma3;
			delta_bt = EPSCVXTY;
			delta_rr = EPSBETA;
			
			alpha2 = alpha1;
			alpha3 = alpha1;
			rho2 = rho1;
			rho3 = rho1;
		}

		/* we need to adjust the alpha and rho */
		/* it is a very stupid algorithm at the moment */

		k = 0;

		bt_theo = 0.5 * (vol1 + vol2) - atm_vol;
		rr_theo = vol2 - vol1;
		
		bt = 0.5 * (impvol11 + impvol21) - atm_vol;						
		rr = impvol21 - impvol11;

		if (bt_theo < bt)
		{
			delta_bt *= -1.0;
		}

		if (rr_theo < rr)
		{
			delta_rr *= -1.0;
		}
		
		while ((k < nbIterSmile) && (error >= precisionSmile))
		{

			/* Compute sensi */

			if ((*param_BT1) + delta_bt > 0.0)
			{
				(*param_BT2)  = (*param_BT1) + delta_bt;
			}
			else
			{
				(*param_BT2)  = (*param_BT1) / 2.0;
			}

			SabrQuadGetParams(spot_fx_cash, gamma2, beta1, &a2, &b2, &c2);			

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha2, a2, b2, c2, rho1, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle, 
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, sigma, sigma_alpha);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha2, a2, b2, c2, rho1, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma_alpha, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol12, &impvol22);

			sensbt1 = (impvol12 - impvol11) / ((*param_BT2) - (*param_BT1));
			sensbt2 = (impvol22 - impvol21) / ((*param_BT2) - (*param_BT1));

			if (use_alpharho)
			{
				if (fabs(rho1 + delta_rr) < 0.99)
				{
					rho2  = rho1 + delta_rr;
				}
				else
				{
					if (rho1 + delta_rr > 0.99)
					{
						rho2 = 0.5 * (rho1 + 1.0);
					}
					else
					{
						rho2 = 0.5 * (rho1 - 1.0);
					}
				}
			}
			else
			{
				beta2 = beta1 + delta_rr;
			}

			SabrQuadGetParams(spot_fx_cash, gamma1, beta2, &a2, &b2, &c2);

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a2, b2, c2, rho2, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, sigma, sigma_rho);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha1, a2, b2, c2, rho2, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma_rho, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol12, &impvol22);

			sensrr1 = (impvol12 - impvol11) / ((*param_RR2) - (*param_RR1));
			sensrr2 = (impvol22 - impvol21) / ((*param_RR2) - (*param_RR1));

			if (fabs(sensrr1) < 10E-06 || fabs(sensrr2) < 10E-06)
			{
				/*	No more sensitivity to rho, it means that alpha is too small */

				if (error < 5.0 * precisionSmile)
				{
					/* we stop here before the algorithm crashs */
					smessage("Nb Iter %d: Smile Calib error %.0f bp", k, error * 10000);

					if (use_alpharho)
					{
						*alpha_out = alpha1;
						*rho_out = rho1;
					}
					else
					{
						*gamma_out = gamma1;
						*beta_out = beta1;
					}

					goto FREE_RETURN;
				}
				else
				{
					/* try to put some alpha */
					if (use_alpharho)
					{
						alpha3 += 0.05;
					}
					else
					{
						gamma3 += 0.005;
					}
				}
			}
			else
			{
				/* solve the matrix */

				delta = sensbt1 * sensrr2 - sensbt2 * sensrr1;

				if (fabs(delta) < 1.0E-10)
				{
					err = "cannot calibrate";
					goto FREE_RETURN;
				}

				(*param_BT3) = (*param_BT1) - (sensrr2 * error1 - sensrr1 * error2) / delta; 
				(*param_RR3) = (*param_RR1) - (-sensbt2 * error1 + sensbt1 * error2) / delta;

				while ((*param_BT3) < 0 || (use_alpharho && fabs(rho3) > 0.9999))
				{
					if ((*param_BT3) < 0)
					{
						coef1 = ((*param_BT1) - (*param_BT1) * 0.5) / ((sensrr2 * error1 - sensrr1 * error2) / delta);
					}

					if (use_alpharho && rho3 < -0.9999)
					{
						coef2 = (rho1 + 0.99) / ((-sensbt2 * error1 + sensbt1 * error2) / delta);
					}
					else if (use_alpharho && rho3 > 0.9999)
					{
						coef2 = (rho1 - 0.99) / ((-sensbt2 * error1 + sensbt1 * error2) / delta);
					}
					else
					{
						coef2 = (beta1 - delta_rr) / ((-sensbt2 * error1 + sensbt1 * error2) / delta);
					}

					coef3 = min(coef1, coef2);

					(*param_BT3) = (*param_BT1) - coef3 * (sensrr2 * error1 - sensrr1 * error2) / delta; 
					(*param_RR3) = (*param_RR1) - coef3 * (-sensbt2 * error1 + sensbt1 * error2) / delta;
				}				
			}

			/* Fill the guess */
			for (i=0; i<nb_calib; i++)
			{
				if (use_alpharho)
				{
					guess[i] = sigma[i] 
							+ ((*param_BT3) - (*param_BT1)) * (sigma_alpha[i] - sigma[i]) / ((*param_BT2) - (*param_BT1))
							+ ((*param_RR3) - (*param_RR1)) * (sigma_rho[i] - sigma[i]) / ((*param_RR2) - (*param_RR1));
				}
				else
				{
					spotbeta1 = exp(beta1 * log(spot_fx_cash));
					spotbeta2 = exp(beta2 * log(spot_fx_cash));
					spotbeta3 = exp(beta3 * log(spot_fx_cash));					

					guess[i] = (sigma[i] * spotbeta1 							
							+ ((*param_RR3) - (*param_RR1)) * (sigma_rho[i] * spotbeta2 - sigma[i] * spotbeta1) / ((*param_RR2) - (*param_RR1))) / spotbeta3;

					guess[i] += ((*param_BT3) - (*param_BT1)) * (sigma_alpha[i] - sigma[i]) / ((*param_BT2) - (*param_BT1));
				}
			}

			SabrQuadGetParams(spot_fx_cash, gamma3, beta3, &a3, &b3, &c3);

			/* calibration */
			err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha3, a3, b3, c3, rho3, lambda, floormu,
										exercise_cal, maturity_cal, volatility_cal, nb_calib, 0, is_straddle,
										nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, guess, sigma);
		
			if (err)
			{
				goto FREE_RETURN;
			}

			/* pricer */
			err	= FxSabrQuadCalibSmilePricer(	today, spot_fx_cash, fwd_fx, alpha3, a3, b3, c3, rho3, lambda, floormu,
											dom_yc, for_yc, exercise_cal, sigma, nb_calib,
											time, date, drift, nstp, func_parm, is_event, nstpt_tgt, nstpfx, nstpvol,
											df, mat_smile, strike1, strike2, res, &impvol13, &impvol23);

			error1 = impvol13 - vol1;
			error2 = impvol23 - vol2;

			error = max(fabs(error1), fabs(error2));

			/* Update alpha, rho, implied vols and the shift */
			if (use_alpharho)
			{
				delta_bt = max(min(fabs((*param_BT3) - (*param_BT1)), EPSALPHAMAX), EPSALPHAMIN);
				delta_rr = max(min(fabs((*param_RR3) - (*param_RR1)), EPSRHOMAX), EPSRHOMIN);
			}
			else
			{
				delta_bt = max(min(fabs((*param_BT3) - (*param_BT1)), EPSCVXTYMAX), EPSCVXTYMIN);
				delta_rr = max(min(fabs((*param_RR3) - (*param_RR1)), EPSBETAMAX), EPSBETAMIN);
			}

			bt = 0.5 * (impvol13 + impvol23) - atm_vol;						
			rr = impvol23 - impvol13;

			if (bt_theo < bt)
			{
				delta_bt *= -1.0;
			}

			if (rr_theo < rr)
			{
				delta_rr *= -1.0;
			}

			(*param_BT1) = (*param_BT3);
			(*param_RR1) = (*param_RR3);

			impvol11 = impvol13;
			impvol21 = impvol23;

			k++;
		}

		smessage("Nb Iter %d: Smile Calib error %.0f bp", k, error * 10000);
	}

	/* calibration of the last sigma */

	SabrQuadGetParams(spot_fx_cash, gamma1, beta1, &a1, &b1, &c1);
	
	if (use_total_ts || nbOpt == 1)
	{
		err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
									exercise, maturity, volatility, nbOpt, start_opt, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
	}
	else
	{
		err = FxSabrQuadCalibration(	dom_yc, for_yc, today, spot_fx, alpha1, a1, b1, c1, rho1, lambda, floormu,
									exercise, maturity, volatility, nbOpt, 0, is_straddle,
									nstpt_tgt, nstpfx, nstpvol, nbIterATM, precisionATM, NULL, sigma);
	}


	*alpha_out = alpha1;
	*rho_out = rho1;
	*gamma_out = gamma1;
	*beta_out = beta1;


FREE_RETURN:

	t2 = clock();

	if (calib_smile)
	{
		smessage ("Total calibration time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	}

	if (time) free (time);
	if (date) free_dvector (date, 0, nstp - 1);
	if (drift)free_dvector (drift, 0, nstp - 1);
	if (func_parm) free_dmatrix(func_parm, 0, nstp-1, 0, 1);
	if (res) free_dvector(res, 0, 1);
	if (is_event) free (is_event);

	if (sigma_alpha) free_dvector(sigma_alpha, 0, nbOpt-1);
	if (sigma_rho) free_dvector(sigma_rho, 0, nbOpt-1);
	if (guess) free_dvector(guess, 0, nbOpt-1);	

	return err;
}