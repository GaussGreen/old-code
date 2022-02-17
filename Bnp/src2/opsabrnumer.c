/*	New SABR pricing functions	*/
/*	Dec99	*/

#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"
#include		"BGMUtils.h"

#include "num_h_cn_pde.h"
#define NSTD 5.0

static void disc_normal_center(double *x, long nstp, double fwd, double xmin, double xmax, double stdn, double nstd, long *index);
static void disc_linleft_logright(double *x, long nstp, double fwd, double xmin, double xmax, double stdln, double nstd, long *index);
static void disc_linleft_logright_strike(double *x, long nstp, double fwd, double xmin, double xmax, double stdln, double nstd, double strike, long *index);

Err op_sabr_adi(
				double							forward,
				double							strike,
				double							maturity,
				double							disc,
				int								call_put_var,	/*	0:	(F - K)+
																	1:	(K - F)+
																	2:	F^2	*/
				double							sigma_beta,
				double							alpha,
				double							beta,
				double							rho,
				int								nt,
				int								nx,
				int								nz,
				double							*result)
{
int							cx, cz;
double						fl, fu, z0, zl, zu;
double						mx, mz;
double						dt, std;

int							i, j, stp;

double						*x		= NULL,
							*z		= NULL,
							**mux	= NULL,
							**muz	= NULL,
							**r		= NULL,
							**varx	= NULL,
							**varz	= NULL,
							***v1	= NULL,
							***v2	= NULL,
							***vtmp = NULL;

double						*varx_init	= NULL,
							*muz_init	= NULL,
							**sig_init	= NULL;

double						*exp_z		= NULL,
							*drift_z	= NULL;

double						*muxi,
							*muzi,
							*varxi,
							*varzi,
							*ri,
							*sig_initi,
							**v1i;

double						sig, pay;

double						F, temp, temp2, temp3, tempi, tempi2, tempi3;

double						const_exp1, const_exp2, exp_temp, t;
CNPDE_TEMP_2D_ADI			adistr, *adi = &adistr;


Err							err = NULL;

	/*	For stability reasons */
	if (alpha < 1.0e-8) alpha = 1.0e-8;

	/*	Force even number of nodes */
	nx = 2 * (nx / 2) + 1;
	nz = 2 * (nz / 2) + 1;

	/*	Precalculate dt */
	dt = maturity / (nt - 1);

	if (fabs(beta - 1.0) < 1.0e-08)
	{
		beta = 1.0 - 1.0E-08;
	}
		
	/*	Calculate grids */
	
	/*	Case 1: beta = 0 Normal grids on F-F0 and log[sigma/sigma0] */
	if (fabs(beta) < 1.0e-08)
	{
		/*	Limits */
		std = sigma_beta * exp (alpha * sqrt (maturity)) * sqrt (maturity);
		fl = - NSTD * std;
		fu = + NSTD * std;
		mx = (fu - fl) / (nx - 1);
		cx = nx / 2;

		z0 = - sigma_beta;
		std = alpha * sigma_beta * sqrt (1.0 - rho * rho) * sqrt (maturity);
		zl = z0 - NSTD * std;
		zu = z0 + NSTD * std;
		mz = (zu - zl) / (nz - 1);
		cz = nz / 2;

		/*	Allocate */
		x = (double*) calloc (nx, sizeof (double));
		z = (double*) calloc (nz, sizeof (double));
		mux = dmatrix (0, nx-1, 0, nz-1);
		muz = dmatrix (0, nx-1, 0, nz-1);
		r = dmatrix (0, nx-1, 0, nz-1);
		varx = dmatrix (0, nx-1, 0, nz-1);
		varz = dmatrix (0, nx-1, 0, nz-1);
		v1 = f3tensor (0, nx-1, 0, nz-1, 0, 0);
		v2 = f3tensor (0, nx-1, 0, nz-1, 0, 0);

		if (!x || !z || !mux || !muz || !r || !varx || !varz || !v1 || !v2)
		{
			err = "Memory allocation error in op_sabr_adi";
			goto FREE_RETURN;
		}


		/*	x */
		for (i=0; i<nx; i++)
		{
			x[i] = fl + i * mx;
		}
		
		/*	z */
		for (j=0; j<nz; j++)
		{
			z[j] = zl + j * mz;
		}

		/*	Final payoff */
		for (i=0; i<nx; i++)
		{
			if (call_put_var == 0)
			{
				pay = max (0, forward + x[i] - strike);
			}
			else
			if (call_put_var == 1)
			{
				pay = max (0, strike - forward - x[i]);
			}
			else
			if (call_put_var == 2)
			{
				pay = (forward + x[i]) * (forward + x[i]);
			}

			for (j=0; j<nz; j++)
			{
				v1[i][j][0] = pay;
			}
		}

		/*	Drifts and variances */
		temp = (1.0 - rho * rho) * alpha * alpha;
		for (i=0; i<nx; i++)
		{
			for (j=0; j<nz; j++)
			{
				mux[i][j] = muz[i][j] = r[i][j] = 0.0;
				sig = max (0.0, rho * alpha * x[i] - z[j]);
				sig *= sig * dt;
				varx[i][j] = sig ;
				varz[i][j] = temp * sig;
			}
		}

		smessage("Convolution: nstept = %d nstepx = %d nstepy = %d", nt, nx, nz);
		/*	Prepare the adi */
		num_f_pde_init_2d_adi (adi, nx, nz, 1);

		/*	Precalculate the probas	*/
		calculate_proba_backward_2f_adi(
										adi,
										nx,
										x,
										nz,
										z,
										1,
										v1,
										mux,
										muz,
										varx,
										varz,
										r,
										v2,
										0,
										nx-1,
										0,
										nz-1);
		
		/*	Do the adi */
		for (stp = nt-2; stp >=0; stp--)
		{
			num_f_pde_one_step_backward_2f_adi_precalc_proba(	
															adi,
															nx,
															x,
															nz,
															z,
															1,
															v1,
															mux,
															muz,
															varx,
															varz,
															r,
															v2,
															0,
															nx-1,
															0,
															nz-1);
		
			vtmp = v2;
			v2 = v1;
			v1 = vtmp;
		}
	}	
	else
	/*	Case 2: Beta != 0 Beta grid on F and normal grid on log[sigma/sigma0]*/
	{
		/*	Limits */
								
		/*	Allocate */
		x = (double*) calloc (nx, sizeof (double));
		z = (double*) calloc (nz, sizeof (double));
		mux = dmatrix (0, nx-1, 0, nz-1);
		muz = dmatrix (0, nx-1, 0, nz-1);
		r = dmatrix (0, nx-1, 0, nz-1);
		varx = dmatrix (0, nx-1, 0, nz-1);
		varz = dmatrix (0, nx-1, 0, nz-1);
		v1 = f3tensor (0, nx-1, 0, nz-1, 0, 0);
		v2 = f3tensor (0, nx-1, 0, nz-1, 0, 0);

		sig_init = dmatrix (0, nx-1, 0, nz-1);
		varx_init = dvector(0, nx-1);
		muz_init = dvector(0, nx-1);
		
		exp_z = dvector(0, nt-1);
		drift_z = dvector(0, nt-1);

		if (!x || !z || !mux || !muz || !r || !varx || !varz || !v1 || !v2 
			|| !sig_init || !varx_init || !muz_init || !exp_z || !drift_z)
		{
			err = "Memory allocation error in op_sabr_adi";
			goto FREE_RETURN;
		}

		/*	Limits of F */

		/* first approximation for variance of the Fx */
		std = sigma_beta * pow (forward, (beta - 1.0)) * sqrt (maturity);
		
		
		fl = log(forward) -std * std / 2.0 - NSTD * std;
		fu = fl  + 2.0 * NSTD * std;
		disc_linleft_logright_strike(x, nx, forward, exp(fl), exp(fu), std, NSTD, strike, &cx);		

		/*	Limits of Z */
		z0 = rho * alpha / (1.0 - beta) * (pow (forward, 1.0 - beta) - 1.0) - sigma_beta;
		std = alpha * sigma_beta * sqrt (1.0 - rho * rho) * sqrt (maturity);

		zl = -NSTD * std;
		zu = NSTD * std;
		disc_normal_center(z, nz, 0, zl, zu, std, NSTD, &cz);

		/*	Estimation of the expectation of Z */
				
		const_exp1 = -0.5 * rho * alpha * beta * pow (forward, beta - 1.0) * sigma_beta * sigma_beta;
		const_exp2 = const_exp1 / alpha / alpha;
		const_exp1 *= dt;

		t = dt / 2;

		for (i=0; i<nt; i++)
		{
			exp_temp = exp(alpha * alpha * t);

			drift_z[i] = const_exp1 * exp_temp;
			exp_z[i] = z0 + const_exp2 * (exp_temp - 1.0);

			t += dt;
		}

		/*	Final payoff */
		for (i=0; i<nx; i++)
		{
			if (call_put_var == 0)
			{
				pay = max (0, x[i] - strike);
			}
			else
			if (call_put_var == 1)
			{
				pay = max (0, strike - x[i]);
			}
			else
			if (call_put_var == 2)
			{
				pay = x[i] * x[i];
			}			

			v1i = v1[i];
			for (j=0; j<nz; j++)
			{
				v1i[j][0] = pay;
			}
		}

		/*	Precalculation of drifts and variances */
		temp = (1.0 - rho * rho) * alpha * alpha * dt;
		temp2 = -0.5 * rho * alpha * beta * dt;
		temp3 = rho * alpha / (1.0 - beta);

		for (i=0; i<nx; i++)
		{
			F = x[i];
			tempi = pow(F, beta);							/*	F ^ beta								*/
			tempi2 = tempi / F;								/*	F ^ (beta - 1)							*/
			tempi3 = temp3 / tempi2 - temp3;	
			varx_init[i] = tempi * tempi * dt;				/*	F ^ (2 * beta)							*/
			muz_init[i] = temp2 * tempi2;					/*	temp2 * F ^ (beta - 1)					*/
			
			sig_initi = sig_init[i];
			ri = r[i];
			muxi = mux[i];

			for (j=0; j<nz; j++)
			{
				sig_initi[j] = tempi3 - z[j];
				ri[j] = 0.0;
				muxi[j] = 0.0;
			}
		}
		
		smessage("Convolution: nstept = %d nstepx = %d nstepy = %d", nt, nx, nz);

		/*	Initialize the adi */
		num_f_pde_init_2d_adi (adi, nx, nz, 1);
			
		/*	Do the adi */
		for (stp = nt-2; stp >=0; stp--)
		{
			for (i=0; i<nx; i++)
			{			
				varxi = varx[i];
				varzi = varz[i];
				muzi = muz[i];

				sig_initi = sig_init[i];
				
				for (j=0; j<nz; j++)
				{
					sig = max(sig_initi[j] - exp_z[stp], 1.0E-08);
					sig *= sig;

					varxi[j] = varx_init[i] * sig;
					muzi[j] = muz_init[i] * sig - drift_z[stp];
					varzi[j] = temp * sig;
				}
			}


			num_f_pde_one_step_backward_2f_adi(	
												adi,
												nx,
												x,
												nz,
												z,
												0,
												0,
												v1,
												mux,
												muz,
												varx,
												varz,
												r,
												v2,
												0,
												nx-1,
												0,
												nz-1);
		
			vtmp = v2;
			v2 = v1;
			v1 = vtmp;
		}
	}

	
	*result = disc * v1[cx][cz][0];

FREE_RETURN:
		
	/*	Free */
	
	if (x) free (x);
	if (z) free (z);
	if (mux) free_dmatrix (mux, 0, nx-1, 0, nz-1);
	if (muz) free_dmatrix (muz, 0, nx-1, 0, nz-1);
	if (r) free_dmatrix (r, 0, nx-1, 0, nz-1);
	if (varx) free_dmatrix (varx, 0, nx-1, 0, nz-1);
	if (varz) free_dmatrix (varz, 0, nx-1, 0, nz-1);
	if (v1) free_f3tensor (v1, 0, nx-1, 0, nz-1, 0, 0);
	if (v2) free_f3tensor (v2, 0, nx-1, 0, nz-1, 0, 0);

	if (sig_init) free_dmatrix (sig_init, 0, nx-1, 0, nz-1);
	if (varx_init) free_dvector(varx_init, 0, nx-1);
	if (muz_init) free_dvector(muz_init, 0, nx-1);

	if (exp_z) free_dvector(exp_z, 0, nt-1);
	if (drift_z) free_dvector(drift_z, 0, nt-1);


	num_f_pde_free_2d_adi (adi, nx, nz, 1);

	return err;
}

double op_sabr_adi_simple(
double							forward,
double							strike,
double							maturity,
double							disc,
int								call_put_var,	/*	0:	(F - K)+
													1:	(K - F)+
													2:	F^2	*/
double							sigma_beta,
double							alpha,
double							beta,
double							rho,
int								nt,
int								nx,
int								nz)
{
	int							cx, cz;
	double						fl, fu, z0, zl, zu;
	double						mx, mz;
	double						dt, std;

	int							i, j, stp;
	double						*x, *z, **mux, **muz, **r, **varx, **varz, ***v1, ***v2, ***vtmp;
	double						*muxi, *muzi, *varxi, *varzi, *ri, **v1i;
	double						sig, pay;

	double						F, temp, temp2, temp3, tempi, tempi2, tempi3, tempi4, tempi5;

	CNPDE_TEMP_2D_ADI			adistr, *adi = &adistr;
	double						result, bs_vol;
	
	Err							err = NULL;

	/*	For stability reasons */
	if (alpha < 1.0e-8) alpha = 1.0e-8;

	/*	Force even number of nodes */
	nx = 2 * (nx / 2) + 1;
	nz = 2 * (nz / 2) + 1;

	/*	Precalculate dt */
	dt = maturity / (nt - 1);

	if (fabs(beta - 1.0) < 1.0e-08)
	{
		beta = 1.0 - 1.0E-08;
	}
		

	/*	Calculate grids */
	
	/*	Case 1: beta = 0 Normal grids on F-F0 and log[sigma/sigma0] */
	if (fabs(beta) < 1.0e-08)
	{
		/*	Limits */
		std = sigma_beta * exp (alpha * sqrt (maturity)) * sqrt (maturity);
		fl = - NSTD * std;
		fu = + NSTD * std;
		mx = (fu - fl) / (nx - 1);
		cx = nx / 2;

		z0 = - sigma_beta;
		std = alpha * sigma_beta * sqrt (1.0 - rho * rho) * sqrt (maturity);
		zl = z0 - NSTD * std;
		zu = z0 + NSTD * std;
		mz = (zu - zl) / (nz - 1);
		cz = nz / 2;

		/*	Allocate */
		x = (double*) calloc (nx, sizeof (double));
		z = (double*) calloc (nz, sizeof (double));
		mux = dmatrix (0, nx-1, 0, nz-1);
		muz = dmatrix (0, nx-1, 0, nz-1);
		r = dmatrix (0, nx-1, 0, nz-1);
		varx = dmatrix (0, nx-1, 0, nz-1);
		varz = dmatrix (0, nx-1, 0, nz-1);
		v1 = f3tensor (0, nx-1, 0, nz-1, 0, 0);
		v2 = f3tensor (0, nx-1, 0, nz-1, 0, 0);

		/*	x */
		for (i=0; i<nx; i++)
		{
			x[i] = fl + i * mx;
		}
		
		/*	z */
		for (j=0; j<nz; j++)
		{
			z[j] = zl + j * mz;
		}

		/*	Final payoff */
		for (i=0; i<nx; i++)
		{
			if (call_put_var == 0)
			{
				pay = max (0, forward + x[i] - strike);
			}
			else
			if (call_put_var == 1)
			{
				pay = max (0, strike - forward - x[i]);
			}
			else
			if (call_put_var == 2)
			{
				pay = (forward + x[i]) * (forward + x[i]);
			}

			for (j=0; j<nz; j++)
			{
				v1[i][j][0] = pay;
			}
		}

		/*	Drifts and variances */
		temp = (1.0 - rho * rho) * alpha * alpha;
		for (i=0; i<nx; i++)
		{
			for (j=0; j<nz; j++)
			{
				mux[i][j] = muz[i][j] = r[i][j] = 0.0;
				sig = max (0.0, rho * alpha * x[i] - z[j]);
				sig *= sig * dt;
				varx[i][j] = sig ;
				varz[i][j] = temp * sig;
			}
		}
	}	
	else
	/*	Case 2: Beta != 0 Beta grid on F and normal grid on log[sigma/sigma0]*/
	{
		/*	Limits */
								
		/*	Allocate */
		x = (double*) calloc (nx, sizeof (double));
		z = (double*) calloc (nz, sizeof (double));
		mux = dmatrix (0, nx-1, 0, nz-1);
		muz = dmatrix (0, nx-1, 0, nz-1);
		r = dmatrix (0, nx-1, 0, nz-1);
		varx = dmatrix (0, nx-1, 0, nz-1);
		varz = dmatrix (0, nx-1, 0, nz-1);
		v1 = f3tensor (0, nx-1, 0, nz-1, 0, 0);
		v2 = f3tensor (0, nx-1, 0, nz-1, 0, 0);

		/*	Limits of log(F) */

		/* first approximation */
		std = sigma_beta * pow (forward, (beta - 1.0)) * sqrt (maturity);

		/* we use SABR Closed Form */
		if (beta <= 1.0 && beta >= 0.0)
		{
			err = srt_f_optsarbvol(
									forward,							
									forward, 
									maturity, 
									sigma_beta,  
									alpha,
									beta,
									rho,
									SRT_BETAVOL,
									SRT_LOGNORMAL,
									&bs_vol);
			
			std = bs_vol * sqrt(maturity);
		}
								
		fl = log(forward) -std * std / 2.0 - NSTD * std;
		fu = fl  + 2.0 * NSTD * std;

		mx = (fu - fl) / (nx - 1);

		/*	We want to have a point at the forward */
		cx = (int) ((log(forward) - fl) / mx + 0.5);
		fl = log(forward) - cx * mx;

		/*	x */
		temp = fl;
		x[0] = exp(temp);
		for (i=1; i<nx; i++)
		{
			temp += mx;
			x[i] = exp(temp);
		}

		/*	We also want a point at the strike of the option */
		cz = (int) ((log(strike) - fl) / mx + 0.5);
		if (cz >= 0 && cz < nx &&  cz != cx)
		{
			x[cz] = strike;
		}
		
		z0 = rho * alpha / (1.0 - beta) * (pow (forward, 1.0 - beta) - 1.0) - sigma_beta;
		/*
		std = alpha * sigma_beta * sqrt (1.0 - rho * rho) * sqrt (maturity);
		*/

		std = sigma_beta * sqrt (1.0 - rho * rho) * sqrt (exp(alpha * alpha * maturity) - 1.0);


		zl = z0 - NSTD * std;
		zu = z0 + NSTD * std;
		mz = (zu - zl) / (nz - 1);
		cz = nz / 2;
		
		/*	z */
		z[0] = zl;
		for (j=1; j<nz; j++)
		{
			z[j] =z[j-1] + mz;
		}

		/*	Final payoff */
		for (i=0; i<nx; i++)
		{
			if (call_put_var == 0)
			{
				pay = max (0, x[i] - strike);
			}
			else
			if (call_put_var == 1)
			{
				pay = max (0, strike - x[i]);
			}
			else
			if (call_put_var == 2)
			{
				pay = x[i] * x[i];
			}

			v1i = v1[i];
			for (j=0; j<nz; j++)
			{
				v1i[j][0] = pay;
			}
		}

		/*	Drifts and variances */
		temp = (1.0 - rho * rho) * alpha * alpha;
		temp2 = -0.5 * rho * alpha * beta;
		temp3 = rho * alpha / (1.0 - beta);

		for (i=0; i<nx; i++)
		{
			F = fabs (x[i]) + 1.0e-06;						/*	F has to be strictly greater than 0.0	*/
			tempi = pow(F, beta);							/*	F ^ beta								*/
			tempi2 = tempi / F;								/*	F ^ (beta - 1)							*/
			tempi3 = temp3 / tempi2 - temp3;	
			tempi4 = tempi * tempi;							/*	F ^ (2 * beta)							*/
			tempi5 = temp2 * tempi2;						/*	temp2 * F ^ (beta - 1)					*/

			ri = r[i];
			muxi = mux[i];
			muzi = muz[i];
			varxi = varx[i];
			varzi = varz[i];

			for (j=0; j<nz; j++)
			{
				ri[j] = 0.0;
				sig = max (0.0, tempi3 - z[j]);
				sig *= sig * dt;
				muxi[j] = 0.0;
				muzi[j] = tempi5 * sig;
				varxi[j] = sig * tempi4;
				varzi[j] = temp * sig;
			}
		}
	}

	smessage("Convolution: nstept = %d nstepx = %d nstepy = %d", nt, nx, nz);
	/*	Prepare the adi */
	num_f_pde_init_2d_adi (adi, nx, nz, 1);

	/*	Precalculate the probas	*/
	calculate_proba_backward_2f_adi(
									adi,
									nx,
									x,
									nz,
									z,
									1,
									v1,
									mux,
									muz,
									varx,
									varz,
									r,
									v2,
									0,
									nx-1,
									0,
									nz-1);
	
	/*	Do the adi */
	for (stp = nt-2; stp >=0; stp--)
	{
		num_f_pde_one_step_backward_2f_adi_precalc_proba(	
														adi,
														nx,
														x,
														nz,
														z,
														1,
														v1,
														mux,
														muz,
														varx,
														varz,
														r,
														v2,
														0,
														nx-1,
														0,
														nz-1);
	
		vtmp = v2;
		v2 = v1;
		v1 = vtmp;
	}
	
	result = disc * v1[cx][cz][0];
		
	/*	Free */
	
	free (x);
	free (z);
	free_dmatrix (mux, 0, nx-1, 0, nz-1);
	free_dmatrix (muz, 0, nx-1, 0, nz-1);
	free_dmatrix (r, 0, nx-1, 0, nz-1);
	free_dmatrix (varx, 0, nx-1, 0, nz-1);
	free_dmatrix (varz, 0, nx-1, 0, nz-1);
	free_f3tensor (v1, 0, nx-1, 0, nz-1, 0, 0);
	free_f3tensor (v2, 0, nx-1, 0, nz-1, 0, 0);

	num_f_pde_free_2d_adi (adi, nx, nz, 1);

	return result;
}


Err op_sabr_mc_ln(
					double							forward,
					int								nb_strike,
					double							*strike,
					double							maturity,
					double							disc,
					double							sigma,
					double							alpha,
					double							rho,
					int								nt,
					long							npaths,
					double							*res,
					double							*volres,
					double							*fwdres)
{
	double	sqdt, dt;
	double	var, var2, sqvar, vol;
	double	fwd, price;
	double	rho2;
	long	path, i;
	long	seed = -123456789;

	double	d1, d2;
	double	lnvol, lnvol0;
	double	drift, rhoratio, rhoratio2;
	double	*lnct = NULL,
			*res1 = NULL,
			*res2 = NULL;

	double	res3;

	double	muSig, sigSig;
	Err		err = NULL;
		
	/*	Time discretisation	*/
	dt = maturity / nt;
	sqdt = sqrt(dt);

	/*	Sigma dynamic		*/
	sigSig = 2.0 * alpha * sqdt;
	muSig = -alpha * alpha * dt;
	
	/*	Initial values		*/
	lnct = dvector(0, nb_strike - 1);
	res1 = dvector(0, nb_strike - 1);
	res2 = dvector(0, nb_strike - 1);

	if (!lnct || !res1 || !res2)
	{
		err = "memory allocation faillure in op_sabr_mc_ln";
		return err;
	}

	for (i=0; i<nb_strike; i++)
	{
		lnct[i] = log(forward / strike[i]);
	}

	rho2 = (1.0 - rho * rho) * dt;
	rhoratio = 0.5 * (1.0 - 2.0 * rho * rho) / (1.0 - rho * rho);
	rhoratio2 = 0.5 * rho * rho / (1.0 - rho * rho);	

	lnvol0 = 2.0 * log(sigma);
		
	res3 = 0;

	/* We discretise the log of the vol			*/
	for (path=0; path<npaths; path++)
	{
		lnvol = lnvol0;
		var = 0.0;

		for (i=0; i<nt; i++)
		{		
			var += exp(lnvol);
			lnvol += muSig + sigSig * gauss_sample (&seed);
		}

		var = var * rho2;
		var2 = var * rhoratio;
		sqvar = sqrt(var);

		vol = exp(0.5 * lnvol);
		drift = rho * (vol - sigma) / alpha;
		
		fwd = forward * exp(-rhoratio2 * var + drift);

		for (i=0; i<nb_strike; i++)
		{
			d1 = (lnct[i] + drift + var2) / sqvar;
			d2 = d1 - sqvar;
			price = fwd * norm(d1) - strike[i] * norm(d2);							
			res1[i] += price / npaths;
			res2[i] += price * price / npaths;				
		}
	}

	for (i=0; i<nb_strike; i++)
	{
		res[i] = disc * res1[i];
		volres[i] = disc * sqrt((res2[i] - res1[i] * res1[i]) / npaths);		
	}

	*fwdres = disc * res3;

	if (lnct) free_dvector (lnct, 0, nb_strike - 1);
	if (res1) free_dvector (res1, 0, nb_strike - 1);
	if (res2) free_dvector (res2, 0, nb_strike - 1);

	return NULL;
}

Err op_sabr_mc_nor(
					double							forward,
					int								nb_strike,
					double							*strike,
					double							maturity,
					double							disc,					
					double							sigma,
					double							alpha,
					double							rho,
					int								nt,
					long							npaths,
					double							*res,
					double							*volres)
{
	double	dt, lnsig0, lnsig, var, vol, price;
	double	std, sqvar;
	double	rho2dt, rhosqdt;
	double	brow, drift, driftvol;	
	long	path, i;
	long	seed = -123456789;

	double	*diff = NULL,
			*res1 = NULL,
			*res2 = NULL;

	double	d1;
	double	pirac;

	Err		err = NULL;
	
	dt = maturity / nt;
	std = 2.0 * alpha * sqrt(dt);	
	rho2dt = (1.0 - rho * rho) * dt;
	rhosqdt = rho * sqrt(dt);
	driftvol = -alpha * alpha * dt;

	lnsig0 = 2.0 * log(sigma);

	diff = dvector(0, nb_strike - 1);
	res1 = dvector(0, nb_strike - 1);
	res2 = dvector(0, nb_strike - 1);

	if (!diff || !res1 || !res2)
	{
		err = "memory allocation faillure in op_sabr_mc_ln";
		return err;
	}

	for (i=0; i<nb_strike; i++)
	{
		diff[i] = forward - strike[i];
	}	
	
	pirac = 0.398942280401433;

	for (path=0; path<npaths; path++)
	{
		lnsig = lnsig0;
		var = 0.0;
		drift = 0.0;

		for (i=0; i<nt; i++)
		{			
			brow = gauss_sample (&seed);			
			var += exp(lnsig);						
			lnsig += driftvol + std * brow;
		}

		var *= rho2dt;
		sqvar = sqrt(var);
		vol = exp(0.5 * lnsig);
		drift = rho * (vol - sigma) / alpha;

		for (i=0; i<nb_strike; i++)
		{
			d1 = (diff[i] + drift) / sqvar;
			price = (diff[i] + drift) * norm(d1) + sqvar * exp(-d1 * d1 / 2.0) * pirac;
			res1[i] += price / npaths;
			res2[i] += price * price / npaths;
		}
	}

	for (i=0; i<nb_strike; i++)
	{
		res[i] = disc * res1[i];
		volres[i] = disc * sqrt((res2[i] - res1[i] * res1[i]) / npaths);
	}

	if (diff) free_dvector (diff, 0, nb_strike - 1);
	if (res1) free_dvector (res1, 0, nb_strike - 1);
	if (res2) free_dvector (res2, 0, nb_strike - 1);

	return NULL;
}
	
/* This routine currently uses a very simple Euler discretization of the SABR equations */
Err op_sabr_mc_beta_Euler(
    double  forward,
    int     nb_strike,
    double  *strike,    // [in]
    double  maturity,   // Option tenor in years
    int     optionType, // 1 for call, -1 for put
    double  df,         // discount factor
    double  sigma,      // "sigma beta"
    double  alpha,      // annualized lognormal vol of sigma
    double  beta,       // in [0,1]
    double  rho,        // in (-1, 1)
    int     nt,         // number of time points
    long    npaths,     // number of paths
    long    seed,
    double  *value,
    double  *sterr)
{
    Err err = NULL;
    int    i, k;
    long   path;
    double dpaths = (double) npaths;
    double dt, rdt;
    double z1, z2;
    double Lf, Ly, h;
    double c1, c2;
    double payoff;

    // Issue a warning
    // err = "op_sabr_mc_beta_Euler: This function is only for testing purposes";
    // return err;

    // Time discretisation
    dt = maturity / (double) nt;
    rdt = sqrt(dt);
    
    // Constants
    c1 = sqrt(1.0 - rho*rho);
    if( seed == 0 ) seed = -123456789;

    // Main Monte Carlo loop
    for( path = 0; path < npaths; ++path )
    {
        // Initialize logarithmic variables
        Lf = log(forward);
        Ly = log(sigma);

        // Path of forward and vol
        for( k = 0; k < nt; ++k )
        {
            // Get two normal deviates
            z1 = gauss_sample (&seed);
            z2 = gauss_sample (&seed);
            
            // Update the processes
            h = exp((beta-1)*Lf + Ly);

            // Put a cap on the local volatility of f
            if( h > 5.0 ) h = 5.0;
            
            Lf += h * (z1 * rdt - 0.5 * h * dt);
            Ly += alpha * ((rho*z1 + c1*z2) * rdt - 0.5 * alpha * dt);

            // Debug output
            if( (path == 8528) && (nt == 10) )
            {
                printf("f=%le, y=%le\n", exp(Lf), exp(Ly));
            }
        }
        
        // Compute option payoff
        for( i = 0; i < nb_strike; ++i )
        {
            payoff = optionType * (exp(Lf) - strike[i]);
            if( payoff < 0.0 ) payoff = 0.0;
            value[i] += payoff;
            sterr[i] += payoff * payoff;
        }
    }

    // Renormalize results
    for( i = 0; i < nb_strike; ++i )
    {
        c2 = value[i] / dpaths;
        value[i] = df * c2;
        sterr[i] = df * sqrt((sterr[i]/dpaths - c2*c2) / dpaths);
    }

    return err;
}


Err	op_sabr_calib_adi(
								double							forward,
								double							strike,
								double							maturity,							
								double							tgt_vol,
								SrtDiffusionType				input_vol_type,
								double							alpha,
								double							beta,
								double							rho,
								int								nt,
								int								nx,
								int								nz,
								int								nbiter,
								double							precision,
								double							*res)
{
double	tgt_price;
double	vol1, vol2, price1, price2, impvol1, impvol2;
double	Sabrp1, Sabrp2;
double	errorp, errorv, shift, vega;
int		k;
Err		err = NULL;

	tgt_price = srt_f_optblksch(forward,
								strike,
								tgt_vol,
								maturity,
								1,
								SRT_CALL,
								PREMIUM
								);

	/*	First Guess: classical formula: sigBeta * K^Beta = sigLN * K = sigNormal */

	if (input_vol_type == SRT_LOGNORMAL)
	{
		vol1 = tgt_vol * exp(log(forward) * (1.0 - beta));
	}
	else
	{
		vol1 = tgt_vol * exp(-log(forward) * beta);
	}

	/*	Calculate the price and impvol1 */

	k = -1;

	err = "Init err";

	while (err && k < nbiter)
	{

		err = op_sabr_adi(	forward,
							strike,
							maturity,
							1.0,
							0,
							vol1,
							alpha,
							beta,
							rho,
							nt,
							nx,
							nz,
							&price1);

		if (err)
		{
			return err;
		}
	
		err = srt_f_optimpvol(		
								price1,
								forward,
								strike,
								maturity,
								1.0,
								SRT_CALL,
								input_vol_type,
								&impvol1);

		if (err)
		{
			vol1 *= 0.8;
		}

		k += 1;
	}

	errorp = (price1 - tgt_price);
	errorv = fabs(impvol1 - tgt_vol);

	if (errorv < precision)
	{
		*res = vol1;
		return err;
	}

	/*	If we have an error, calculate a second point using SABR closed form vega */
	
	if (beta >=0.0 && beta <=1)
	{
		err = op_sabr_pricing(	forward,
								strike,
								maturity,
								1.0,
								SRT_CALL,
								alpha,
								beta,
								rho,
								vol1,
								SABR_ATM_BETA,
								&Sabrp1
								);
		if (err)
		{
			return err;
		}

		shift = vol1 * 0.05;
		if (errorp < 0)
		{
			shift *= -1.0;
		}

		err = op_sabr_pricing(	forward,
								strike,
								maturity,
								1.0,
								SRT_CALL,
								alpha,
								beta,
								rho,
								(vol1 + shift),
								SABR_ATM_BETA,
								&Sabrp2
								);
		if (err)
		{
			return err;
		}

		vega = (Sabrp2 - Sabrp1) / shift;
	}
	else
	{
		shift = vol1 * 0.05;
		if (errorp < 0)
		{
			shift *= -1.0;
		}

		err = op_sabr_adi(	forward,
								strike,
								maturity,
								1.0,
								0,
								(vol1 + shift),
								alpha,
								beta,
								rho,
								nt,
								nx,
								nz,
								&price2);
		if (err)
		{
			return err;
		}

		vega = (price2 - price1) / shift;
	}


	/*	Second guess: vol1 adjusted by the vega */

	vol2 = vol1 - errorp / vega;

	/*	Now do a Newton to calibrate */
	
	k = 0;

	while ((k < nbiter) && (errorv > precision))
	{
		err = op_sabr_adi(	forward,
								strike,
								maturity,
								1.0,
								0,
								vol2,
								alpha,
								beta,
								rho,
								nt,
								nx,
								nz,
								&price2);

		if (err)
		{
			return err;
		}

		err = srt_f_optimpvol(		
							price2,
							forward,
							forward,
							maturity,
							1.0,
							SRT_CALL,
							input_vol_type,
							&impvol2);

		if (err)
		{
			return err;
		}

		errorp = (price2 - tgt_price);
		errorv = fabs(impvol2 - tgt_vol);

		if (errorv < precision)
		{
			*res = vol2;
			return err;
		}

		vega = (price2 - price1) / (vol2 - vol1);

		if (vega < 0)
		{
			/* convergence problem... */
			nx = (int) (nx * 1.5);
			nz = (int) (nz * 1.5);
		}
		else
		{
			/* adjust */
			vol1 = vol2;
			price1 = price2;

			vol2 = vol1 - errorp / vega;
		}

		k++;
	}

	*res = vol2;

	return err;
}

static void	disc_normal_center(		double	*x,
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

static void	disc_linleft_logright(	double	*x,
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

static void	disc_linleft_logright_strike(	double	*x,
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

	if (strike > fwd || strike < xmin || index_s == (*index))
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

	if (strike <fwd || strike > xmax || index_s == (*index))
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

/* --------------------------------------------------------------------------------------------------------------------------------- */
/* Recalibrates a SABR smile by changing the beta or the rho using ADI */
/* --------------------------------------------------------------------------------------------------------------------------------- */
#define R 0.61803399
#define C (1.0-R)
#define SHIFT2(a,b,c) (a)=(b);(b)=(c);
#define SHIFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

Err srt_ADINewRhoSABR(
					 double forward,
					 double maturity,
					 double ATMVol,
					 double alpha,
					 double beta,
					 double rho,
					 double newbeta,
					 int nt, int nx, int nz,
					 int iNumStrikes,
					 double *pdNewRho
					 )
{

	Err err=NULL;
	double max_std=2.5; /* for the smile */
	int    num_strikes; /* has to be odd */
	double dPrice;
	double *strikes = NULL;
	double *vols = NULL;
	double dTempVol, dSigmaBeta;
	double ax, bx, cx, f1, f2, x0, x1, x2, x3;
	double tol = 1E-05;
	double dJump = 0.005;
	int i;

/* Set local variables */
	num_strikes = iNumStrikes;


/* Allocates memory */
	strikes=dvector(1,num_strikes);
	vols=dvector(1,num_strikes);

/* Calculates the strikes */
	for (i=1;i<=num_strikes;i++)
		strikes[i]=forward * exp( (i-(iNumStrikes+1)/2)*ATMVol*sqrt(maturity) );

/* Calculate the sigmabeta for the given ATM vol (using the old value of beta and rho */ 
	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,beta,rho,nt,nx,nz,128,1e-04,&dSigmaBeta);


/* Calculates the vols at strikes */
	for (i=1;i<num_strikes+1 && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, beta, rho, nt, nx, nz, &dPrice );
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &vols[i] );
	}

/* Use the SABR vol approximation to get a good first guess */
	err=srt_BGMNewRhoSABR(forward,maturity,ATMVol,alpha,beta,rho,newbeta,&bx);



/* ----------------------------------------------------------------------------------------------------------------------------------*/
/* 1-D Calibration using 'golden'.  Assumes that the error function is quadratic in rho.  We also need a point that has an error     */
/* less than the error at +/- 1  Try to bracket the minimum by expanding about the analytic approximation */

	ax = bx;
	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,newbeta,ax,nt,nx,nz,128,1.0e-04,&dSigmaBeta);
	f1 = 0.0;
	for (i=1;i<num_strikes+1 && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, ax, nt, nx, nz, &dPrice);
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f1 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
	}
	f2 = f1 - 1.0;
	while ( f2 < f1 && !err)
	{
		if ( (ax -= dJump) < -1.0 )
			err = "ax got too small";
		else
		{
			f2 = 0.0;
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, ax, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );

			}
		}
	}

	cx = bx;
	f2 = f1 - 1.0;
	while ( f2 < f1 && !err)
	{
		if ( (cx += dJump) > 1.0 )
			err = "cx got too big";
		else
		{
			f2 = 0.0;
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, cx, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );

			}
		}
	}



/* Start of the NR code */
	x0=ax; 
	x3=cx;
	if(fabs(cx-bx)>fabs(bx-ax))
	{
		x1=bx;
		x2=bx+C*(cx-bx);
    }
    else
    {
		x2=bx; 
		x1=bx-C*(bx-ax);
	}

/* Calculate the function values at rho = x1, x2 */
	f1 = 0.0;
	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,newbeta,x1,nt,nx,nz,128,1.0e-04,&dSigmaBeta);
	for (i=1;i<num_strikes+1 && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, x1, nt, nx, nz, &dPrice );
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f1 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
	}

	f2 = 0.0;
	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,newbeta,x2,nt,nx,nz,128,1.0e-04,&dSigmaBeta);
	for (i=1;i<num_strikes+1 && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, x2, nt, nx, nz, &dPrice );
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
	}
	
/* Loop until we have a sufficiently small bracket on the minimum */
	while ( fabs(x3-x0) > tol && !err )
    {
		if(f2<f1) 
		{ 
			SHIFT3(x0,x1,x2,R*x1+C*x3); 
	  		f1=f2;
			f2 = 0.0;
			err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,newbeta,x2,nt,nx,nz,128,1.0e-04,&dSigmaBeta);
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, x2, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
			}
		}
		else
		{ 
			SHIFT3(x3,x2,x1,R*x2+C*x0); 
			f2=f1;
			f1 = 0.0;
			err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,newbeta,x1,nt,nx,nz,128,1.0e-04,&dSigmaBeta);
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, alpha, newbeta, x1, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f1 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
			}
		}
    }

/* Set the value */
    if(f1<f2)
		*pdNewRho=x1; 
    else
		*pdNewRho=x2; 

/* free memory */
	free_dvector(strikes,1,num_strikes);
	strikes = NULL;
	free_dvector(vols,1,num_strikes);
	vols= NULL;

	return err;
}

#undef R 
#undef C 
#undef SHIFT2
#undef SHIFT3



/* --------------------------------------------------------------------------------------------------------------------------------- */
/* Given a set of SABR parameters, calculates the price, and then calculates the alpha necessary in the ADI model to match the price */
/* --------------------------------------------------------------------------------------------------------------------------------- */
#define R 0.61803399
#define C (1.0-R)
#define SHIFT2(a,b,c) (a)=(b);(b)=(c);
#define SHIFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

Err srt_AlphaFromADI(
					 double forward,
					 double maturity,
					 double ATMVol,
					 double alpha,
					 double beta,
					 double rho,
					 double Disc,

					 int nt, 
					 int nx, 
					 int nz,
					 int iNumStrikes,
					
					 double strikeMin,
					 double strikeMax,

					 double *pdNewAlpha
					 )
{

	Err err=NULL;
	double max_std=2.5; /* for the smile */
	int    num_strikes; /* has to be odd */
	double dPrice;
	double *strikes = NULL;
	double *vols = NULL;
	double dTempVol;
	double ax, bx, cx, f1, f2, x0, x1, x2, x3;
	double tol = 1E-07;
	double dJump = 0.005,gap,atm_SABR_vol,atm_ADI_vol;
	int i;

/* Set local variables */
	num_strikes = iNumStrikes;


/* Allocates memory */
	strikes=dvector(1,num_strikes);
	vols=dvector(1,num_strikes);

/* Calculates the strikes */
/*
		dApproxStdev = ATMVol / pow(forward, 1.0 - beta) * sqrt(maturity);		
		n_StdDev=3;

		strikeMin = forward * exp( - 0.5 * dApproxStdev * dApproxStdev - n_StdDev * dApproxStdev); 
		strikeMax = forward * exp( - 0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);
*/
	for (i=1;i<=num_strikes;i++) 
		strikes[i]=strikeMin+(i-1)*(strikeMax-strikeMin)/(num_strikes-1);

/* Calculate the sigmabeta for the given ATM vol (using the old value of beta and rho */ 
//	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,alpha,beta,rho,nt,nx,nz,128,&dSigmaBeta);

/* Calculates the SABR Log smile */
	for (i=1;i<num_strikes+1 && !err;i++)
		err = srt_f_optsarbvol(	forward,strikes[i],maturity,ATMVol,alpha,beta,rho,SABR_BETAVOL,SABR_LOGVOL,&vols[i]);


/* ----------------------------------------------------------------------------------------------------------*/
/* Step1: We try to "bracket the minimum of the objective function. We find three points ax<bx<cx such that
	      the value of the objective function is minimum at f(bx). These points will be the inputs of the next 
	      function        
*/

/*
	bx=alpha+0.07; 

	ax = bx;
	err=op_sabr_calib_adi( forward,forward,maturity,ATMVol,SRT_LOGNORMAL,ax,beta,rho,nt,nx,nz,128,&dSigmaBeta);
	f1 = 0.0;
	for (i=1;i<num_strikes+1 && !err;i++)
	{
		dPrice = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, ax, beta, rho, nt, nx, nz );
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f1 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );
	}
	f2 = f1 - 1.0;
	while ( f2 < f1 && !err)
	{
		if ( (ax -= dJump) < -1.0 )
			err = "ax got too small";
		else
		{
			f2 = 0.0;
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				dPrice = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, ax, beta, rho, nt, nx, nz );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );

			}
		}
	}

	cx = bx;
	f2 = f1 - 1.0;
	while ( f2 < f1 && !err)
	{
		if ( (cx += dJump) > 1.0 )
			err = "cx got too big";
		else
		{
			f2 = 0.0;
			for (i=1;i<num_strikes+1 && !err;i++)
			{
				dPrice = op_sabr_adi( forward, strikes[i], maturity, 1.0, 0, dSigmaBeta, cx, beta, rho, nt, nx, nz );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol-vols[i] ) * ( dTempVol-vols[i] );

			}
		}
	}

*/

	ax=alpha;
	bx=alpha+0.1;
	cx=alpha+0.2;

/* Start of the NR code */
	x0=ax; 
	x3=cx;
	if(fabs(cx-bx)>fabs(bx-ax))
	{
		x1=bx;
		x2=bx+C*(cx-bx);
    }
    else
    {
		x2=bx; 
		x1=bx-C*(bx-ax);
	}

/* Calculate the function values at alpha = x1, x2 */

	err = srt_f_optsarbvol(	forward,forward,maturity,ATMVol,alpha,beta,rho,SABR_BETAVOL,SABR_LOGVOL,&atm_SABR_vol);

	err = op_sabr_adi(forward,forward,maturity,Disc,SRT_CALL,ATMVol,x1,beta,rho,nt,nx,nz,&dPrice);
	err = srt_f_optimpvol(dPrice,forward,forward,maturity,Disc,SRT_CALL,SABR_LOGVOL,&atm_ADI_vol);
	gap=atm_SABR_vol-atm_ADI_vol;

	f1=0.0;
	for (i=1;i<=num_strikes && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, Disc, SRT_CALL, ATMVol,x1,beta,rho,nt,nx,nz,&dPrice);
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, Disc, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f1 += ( dTempVol+gap-vols[i] ) * ( dTempVol+gap-vols[i] );
	}

	err = op_sabr_adi(forward,forward,maturity,Disc,SRT_CALL,ATMVol,x2,beta,rho,nt,nx,nz,&dPrice);
	err = srt_f_optimpvol(dPrice,forward,forward,maturity,Disc,SRT_CALL,SABR_LOGVOL,&atm_ADI_vol);
	gap=atm_SABR_vol-atm_ADI_vol;

	f2=0.0;
	for (i=1;i<=num_strikes && !err;i++)
	{
		err = op_sabr_adi( forward, strikes[i], maturity, Disc, SRT_CALL, ATMVol, x2, beta, rho, nt, nx, nz, &dPrice );
		err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, Disc, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
		f2 += ( dTempVol+gap-vols[i] ) * ( dTempVol+gap-vols[i] );
	}
	
/* Loop until we have a sufficiently small bracket on the minimum */
	while ( fabs(x3-x0) > tol && !err )
    {
		if(f2<f1) 
		{ 
			SHIFT3(x0,x1,x2,R*x1+C*x3); 
	  		f1=f2;
			f2 = 0.0;

			err = op_sabr_adi(forward,forward,maturity,Disc,SRT_CALL,ATMVol,x2,beta,rho,nt,nx,nz,&dPrice);
			err = srt_f_optimpvol(dPrice,forward,forward,maturity,Disc,SRT_CALL,SABR_LOGVOL,&atm_ADI_vol);
			gap=atm_SABR_vol-atm_ADI_vol;

			for (i=1;i<=num_strikes && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, Disc, SRT_CALL, ATMVol, x2, beta, rho, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, Disc, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f2 += ( dTempVol+gap-vols[i] ) * ( dTempVol+gap-vols[i] );
			}
		}
		else
		{ 
			SHIFT3(x3,x2,x1,R*x2+C*x0); 
			f2=f1;
			f1 = 0.0;

			err = op_sabr_adi(forward,forward,maturity,Disc,SRT_CALL,ATMVol,x1,beta,rho,nt,nx,nz,&dPrice);
			err = srt_f_optimpvol(dPrice,forward,forward,maturity,Disc,SRT_CALL,SABR_LOGVOL,&atm_ADI_vol);
			gap=atm_SABR_vol-atm_ADI_vol;

			for (i=1;i<num_strikes+1 && !err;i++)
			{
				err = op_sabr_adi( forward, strikes[i], maturity, Disc, SRT_CALL, ATMVol, x1, beta, rho, nt, nx, nz, &dPrice );
				err = srt_f_optimpvol( dPrice, forward, strikes[i], maturity, Disc, SRT_CALL, SRT_LOGNORMAL, &dTempVol );
				f1 += ( dTempVol+gap-vols[i] ) * ( dTempVol+gap-vols[i] );
			}
		}
    }

/* Set the value */
    if(f1<f2)
		*pdNewAlpha=x1; 
    else
		*pdNewAlpha=x2; 

/* free memory */
	free_dvector(strikes,1,num_strikes);
	strikes = NULL;
	free_dvector(vols,1,num_strikes);
	vols= NULL;

	return err;
}

#undef R 
#undef C 
#undef SHIFT2
#undef SHIFT3


