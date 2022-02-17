#include "srt_h_all.h"
#include "CheyBetaGrfn.h"
#include "Fx3FUtils.h"
#include "math.h"


Err	 cheybeta_pricing_pde(	
					/*	Time data		*/
					int			nstept,					
					double		*time,
					double		*date,

					/*	Discretisation	*/
					int			nx,
					int         nphi,
										
					/*	Model data		*/
					SrtUndPtr   und,					

					/*	Product data */
					void		**func_parm_tab, 
					int			*eval_evt,
					
					/*	Market data */
					double		*ifr,					
					char		*yc,
					
					/*	Payoff function */
					Err (*payoff_func)(
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										
										/* Market data	*/										
										void	*yc,
										double  ifr,
										/* Model data	*/
										TermStruct	*ts,					
										/* Gride data	*/
										int		lx,
										int		ux,
										int		lphi,
										int		uphi,
										double	*x,
										double	*phi,
																
										/* Vector of results to be updated */
										int		nprod,
										double	***prod_val
										),
					/*	Result */
					int			nprod, 
					double		*res)
{
Err						err = NULL;
CHEYBETA_MDL			mdl;
CHEYBETA_GRID			grid;
CHEYBETA_BCK_PDE_TEMP   tmp;
CNPDE_TEMP			    tmppde;

int                     i,t;
int						indexx, indexphi;
double                  meshx;
double                  stubb;

	/*	nstep has to be a odd nuber			*/
	nx = ((int) (nx / 2)) * 2 + 1;
	nphi = ((int) (nphi / 2)) * 2 + 1;


	/*	Init chey beta bck grids and structures	*/
	chey_beta_mdl_init (&mdl);
	chey_beta_mdl_build_from_und (&mdl, und);
	chey_beta_grid_init (&grid);
	chey_beta_grid_time (&mdl, nstept, time, nstept, &grid);
	chey_beta_grid_phi_x (&mdl, &grid, nx, nphi);

	/* force on x to be 0 */
	meshx = ( - grid.x[0] + grid.x[grid.nx-1]) / (grid.nx - 1);
	indexx = (int) (( - 1.0 *  grid.x[0]) / meshx + 0.5);
	stubb = grid.x[indexx];
	if (fabs(stubb) > 1e-9)
	{
		for (i=0; i<grid.nx; i++)
			grid.x[i] -= stubb;
	}
	indexphi = 0;

	chey_beta_grid_prod_val(&grid,nprod);
	srt_f_chey_beta_bck_pde_alloc_temp(&tmp,nx,nphi,nprod);

	/*	Init PDE	*/
	num_f_pde_init (&tmppde, grid.nx, grid.num_prod);
	
	/*	Final payoff valuation					*/
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in cheybeta_pricing";
		goto FREE_RETURN;
	}

	/*	Eval payoff */
	err = payoff_func(
					date[nstept-1],
					time[nstept-1],
					func_parm_tab[nstept-1], yc, ifr[0],
					/* Model data	*/
					((SrtIrDesc*) (und->spec_desc))->ts,					
					/* Gride data	*/
					0,grid.nx-1,0,grid.nphi[nstept-1]-1,grid.x,grid.phi[nstept-1],				
					/* Vector of results to be updated */
					nprod,grid.vt2);
	if (err)
	{
		goto FREE_RETURN;
	}

	/* initialise the loop */
	for (t=nstept-2; t>=0; t--)
	{

		srt_f_chey_beta_bck_pde_one_step(
						&tmppde,&tmp,&mdl,&grid,t,t+1);

				/*	Eval payoff */
		if (eval_evt[t])
		{

			err = payoff_func(
					date[t],
					time[t],
					func_parm_tab[t], yc, ifr[t],
					/* Model data	*/
					((SrtIrDesc*) (und->spec_desc))->ts,					
					/* Gride data	*/
					0,grid.nx-1,0,grid.nphi[t]-1,grid.x,grid.phi[t],				
					/* Vector of results to be updated */
					nprod,grid.vt2);

			if (err)
			{
				goto FREE_RETURN;
			}
		}
	}

	/* copy the result */
	for (i=0; i<nprod; i++)
		res[i] = grid.vt2[indexphi][indexx][i];


	/* Free struct PDE */
	num_f_pde_free (&tmppde, grid.nx, grid.num_prod);
	srt_f_chey_beta_bck_pde_free_temp (&tmp);
	chey_beta_mdl_free (&mdl);
	chey_beta_grid_free (&grid);

FREE_RETURN:

	return err;
}

Err	 cheybeta_pricing_mc(	
					/*	Time data		*/
					int			nstept,					
					double		*time,
					double		*date,

					/*	Discretisation	*/
					int			numpaths,
					int			method, /* 0 balantisam, 1 balantisam adjusted */
					SrtMCSamType	gen_method,
					
					/*	Model data		*/
					SrtUndPtr   und,					

					/*	Product data */
					void		**func_parm_tab, 
					int			*eval_evt,
					
					/*	Market data */
					double		*ifr,					
					char		*yc,
					int		    *stop_vol,
					int         *stop_lambda,
					
					/*	Payoff function */
					Err (*payoff_func)(
										double	evt_date,
										double	evt_time,
										void	*func_parm, 
										
										/* Market data	*/										
										void	*yc,
										/* var of model	*/
										double  r,
										double	x,
										double	phi,				
										/* Vector of results to be updated */
										int		ncols,
										double	*cols_val
										),
					/*	Result */
					int			nprod, 
					double		**res)
{
Err						err = NULL;
int                     n,i,t;

double					***cube = NULL;
double                  **matrix = NULL;
double					*currentvalues = NULL;
double					*temppayoffs = NULL;
double                  *sum_payoff = NULL;
double                  *sum_payoff2 = NULL;

double					t1, t2, dt, num, x, phi, stdv;

GRFNPARMCHEYBETA		total;
CHEYBETA_MDL			mdl;
double					beta, tau, betavol;

long                    seed;

double                  strike = 1.0;
double                  payout = 0.0;

	/* initialise the curvalues */
	currentvalues = dvector(0,nprod-1);
	temppayoffs = dvector(0,nprod-1);
	sum_payoff = dvector(0,nprod-1);
	sum_payoff2 = dvector(0,nprod-1);
	
	/* in MC case we pre-compute the coefficient of the reconstruction formula */
	for (t=0; t<=nstept-1; t++)
	{	
		if (eval_evt[t])
		{
			total = (GRFNPARMCHEYBETA) func_parm_tab[t];
	
			if (total->num_df > 0)
			{		
				for (i=0; i<total->num_df; i++)
				{
					if (total->df_tms[i] >= YEARS_IN_DAY)
					{
						total -> dff[i] = swp_f_df (date[t], total->df_dts[i], yc);
						total -> gam[i] = Lambda_func(time[t],time[t]+total->df_tms[i],((SrtIrDesc*) (und->spec_desc))->ts);
						total -> gam_sqr[i] = 0.5 *total -> gam[i]*total -> gam[i];
					}
					else
					{
						total -> dff[i] = 1.0;
						total -> gam[i] = 0.0;
						total -> gam_sqr[i] = 0.0;
					}
			
				}
			}
	
		}
	}
	
	chey_beta_mdl_init(&mdl);
	chey_beta_mdl_build_from_und(&mdl,und);
	
	/* fill the Brownian cube */

	/* npaths has to be odd */
	numpaths = 2 * ((long) (numpaths / 2)) + 1;
	
	/* by default use Balantisam */
	cube = (double ***) malloc(numpaths*sizeof(double **));
	for (i=0; i<numpaths; i++)
	{
		cube[i] = (double **) malloc(1*sizeof(double *));
		cube[i][0] = (double *) malloc((nstept-1)*sizeof(double));
	}

	seed = -1234567;

	err = srt_f_New_gen_random(
				cube, 
				0, 
				numpaths-1, 
				0, 
				nstept-2, 
				1, 
				gen_method, 
				&seed, 
				time);
	if (err)
		goto FREE_RETURN;

	if (method && ((gen_method == BAL_ANTI) || (gen_method == RANDOM_GAUSS) || (gen_method == ANTITHETIC)) )
	{
		/* indeed, it's ugly ! */
		matrix = dmatrix(0,numpaths-1,0,nstept-2);
		for (n=0; n< numpaths; n++)
			for (t=0; t<nstept-1; t++)
				matrix[n][t] = cube[n][0][t];

		adjust_emp_covar (matrix, numpaths,  nstept - 1);

		for (n=0; n< numpaths; n++)
			for (t=0; t<nstept-1; t++)
				cube[n][0][t] = matrix[n][t];
	}

	/* loop on paths */
	for (n=0; n< numpaths; n++)
	{
		t1 = time[0];
		x = 0.0;
		phi = 0.0;
		num = 1.0 / numpaths;

		betavol = mdl.sigma[0];
		tau = mdl.lambda[0];
		beta = mdl.beta[0];

		if (eval_evt[0])
		{
			/*	Eval payoff at the first evt payoff */
			err = payoff_func(
							date[0],
							time[0],
							func_parm_tab[0], yc,
							/* Model data	*/
							ifr[0],x,phi,				
							/* Vector of results to be updated */
							nprod,currentvalues);
			if (err)
				goto FREE_RETURN;
		}

		/* loop on time  */
		for (t=1; t<=nstept-1; t++)
		{
			t2 = time[t];
			dt = t2 - t1;

			if (stop_vol[t-1])
				betavol = mdl.sigma[stop_vol[t-1]];

			if (stop_lambda[t-1])
				tau = mdl.lambda[stop_lambda[t-1]];
	
			/* implementation of the cheyette beta dynamic */
			stdv = betavol * sqrt(dt) * pow(x + ifr[t-1], beta);
			
			x += ( phi   - x * tau ) * dt 
				+ stdv * cube[n][0][t-1];
			
			phi +=   stdv * stdv - 2.0 * tau * phi * dt;

			num *= exp ( - (x + ifr[t-1]) * dt);

			/*	Eval payoff */
			if (eval_evt[t])
			{
				err = payoff_func(
							date[t],
							time[t],
							func_parm_tab[t], yc,
							/* Model data	*/
							x+ifr[t],x,phi,				
							/* Vector of results to be updated */
							nprod,temppayoffs);

				if (err)
					goto FREE_RETURN;

				for (i=0; i<nprod; i++)
					currentvalues[i] += temppayoffs[i] * num;
			}
			t1 = t2;
		}

		for (i=0; i<nprod; i++)
		{
			sum_payoff[i] += currentvalues[i];
			sum_payoff2[i] += currentvalues[i] * currentvalues[i];
			currentvalues[i] = 0.0;
		}
	}

	/* copy the result */
	for (i=0; i<nprod; i++)
	{
		/* the mean */
		res[0][i] = sum_payoff[i];
		/* the confidence interval */
		res[1][i] = sqrt((sum_payoff2[i] * numpaths - sum_payoff[i] * sum_payoff[i]) / numpaths);
	}

FREE_RETURN:

	if (sum_payoff) free_dvector(sum_payoff,0,nprod);
	if (sum_payoff2) free_dvector(sum_payoff2,0,nprod);
	if (currentvalues) free_dvector(currentvalues,0,nprod);
	if (temppayoffs) free_dvector(temppayoffs,0,nprod);
	if (cube)
	{
		for (i=0; i<numpaths; i++)
		{
			free(cube[i][0]);
			free(cube[i]);
		}
		free(cube);
	}

	sum_payoff = NULL; sum_payoff2 = NULL;
	currentvalues = NULL; temppayoffs = NULL; cube = NULL;

	if (&mdl) chey_beta_mdl_free(&mdl);

	return err;
}