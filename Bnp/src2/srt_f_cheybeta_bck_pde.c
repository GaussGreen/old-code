/* ==============================================================================
   
   FILE NAME:      srt_f_cheybeta_bck_pde.c 
   
   OBJECT:         functions implementing the local vol Cheyette backward PDE

  =============================================================================== */

#include "srt_h_all.h"

/*	Calculate the coefficients of a linear interpolation,
	i.e. find c1..cn such that f(x0) = sum ci * f(xi)	*/
static void lininter_quick_coef(
double							*x,
int								n,
double							x0,
double							*cup,
double							*cdown,
long                            *indexup,
long							*indexdown)
{
	int i;

	/*	Left extrapolation flat	*/
	if (x0 <= x[0])
	{
		*cdown = 1.0;
		*cup = 0.0;
		*indexup = 0;
		*indexdown = 0;
		return;
	}
	 
	/*	Right extrapolation flat	*/
	if (x0 >= x[n-1])
	{
		*cup = 1.0;
		*cdown = 0.0;
		*indexup = n-1;
		*indexdown = n-1;
		return;
	}

	/*	Find i such that xi <= x0 < xi+1	*/
	i = 0;
	while (x0 >= x[i+1])
	{
		i++;
	}

	/* fill the index */
	*indexup = i+1;
	*indexdown = i;

	/*	Calculate linear interpolation coefficients of x0 in terms of xi and xi+1	*/
	*cdown = (x[i+1] - x0) / (x[i+1] - x[i]);
	*cup = 1.0 - (*cdown);
}

static void lininter_coef(
double							*x,
int								n,
double							x0,
double							*c)
{
	int i;

	/*	Set all ci to 0	*/
	memset (c, 0, n * sizeof (double));

	/*	Left extrapolation flat	*/
	if (x0 <= x[0])
	{
		c[0] = 1.0;
		return;
	}
	 
	/*	Right extrapolation flat	*/
	if (x0 >= x[n-1])
	{
		c[n-1] = 1.0;
		return;
	}

	/*	Find i such that xi <= x0 < xi+1	*/
	i = 0;
	while (x0 >= x[i+1])
	{
		i++;
	}


	/*	Calculate linear interpolation coefficients of x0 in terms of xi and xi+1	*/
	c[i] = (x[i+1] - x0) / (x[i+1] - x[i]);
	c[i+1] = 1.0 - c[i];
}

/*	Backward PDE Functions	*/

/*	Allocate a few temp vectors	*/
void srt_f_chey_beta_bck_pde_alloc_temp(
CHEYBETA_BCK_PDE_TEMP			*tmp,
int								nx,
int								nphi,
int								num_prod)
{
	tmp->coef = (double*) calloc (nphi, sizeof (double));
	tmp->val = dmatrix (0, nx-1 , 0, num_prod-1);
	tmp->mu = (double*) calloc (nx, sizeof (double));
	tmp->var = (double*) calloc (nx, sizeof (double));
	tmp->r = (double*) calloc (nx, sizeof (double));
}

/*	Free temp vectors	*/
void srt_f_chey_beta_bck_pde_free_temp(
CHEYBETA_BCK_PDE_TEMP			*tmp)
{
	free (tmp->coef);
	free_dmatrix (tmp->val, 0, 1, 0, 1);
	free (tmp->mu);
	free (tmp->var);
	free (tmp->r);
}

/*	Get backward phi	*/
static double chey_beta_mdl_phi_back(
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Param	*/
CHEYBETA_PARAM					*param,
/*	Statevars	*/
double							x,
double							phi,
/*	Output from chey_beta_mdl_var	*/
double							var)
{
	return (phi - var) / ( 1.0 - 2 * param->lambda * param->dt) ;
}

/*	Go backward one step	*/
void srt_f_chey_beta_bck_pde_one_step(
CNPDE_TEMP						*tmppde,
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								t1,
/*	Next time step idx	*/
int								t2)
{
	int k;
	int phi_slice;
	double nvar;
	CHEYBETA_PARAM p;
	int xi;
	double x;
	double ***temp;

	double  cup,cdown;
	long    indexup, indexdown;

	/*	Get local diffusion parameters	*/
	chey_beta_mdl_param (mdl, grid->t[t1], grid->t[t2], &p);

	/* before to enter into the loop on phi pre-compute a few things */
	for (xi=0; xi<grid->nx; xi++)
	{
		/* Store the x */
		x = grid->x[xi];

		/*	Calculate r	*/
		tmpcb->r[xi] = (grid->fwdrt[t1] + x) * p.dt;
		
		/*	Calculate x coefficients var : independent of phi */
		nvar = chey_beta_mdl_norm_var_at_t (	
				&p, 
				x, 
				0.0,
				grid->fwdrt[t1],
				grid->maxvar[t1],
				grid->minvar[t1]);
	
		tmpcb->var[xi] = chey_beta_mdl_var (mdl, 
											&p, 
											x, 
											0.0, 
											nvar);
	}

	/*	Go phi slice by phi slice	*/
	for (phi_slice=0; phi_slice<grid->nphi[t1]; phi_slice++)
	{
		/*	Calculate drifts, variances and discount factors of x and y	*/
		for (xi=0; xi<grid->nx; xi++)
		{											
			x = grid->x[xi];
			
			/*	Calculate mu of x	*/
			tmpcb->mu[xi] = chey_beta_mdl_drift (	mdl, 
													&p, 
													x, 
													grid->phi[t1][phi_slice]);

			/*	Do explicit convolution in phi	*/
			lininter_quick_coef (
							grid->phi[t2], 
							grid->nphi[t2], 
							chey_beta_mdl_phi (	mdl, 
												&p, 
												x, 
												grid->phi[t1][phi_slice],
												tmpcb->var[xi] ),
							&cup,
							&cdown,
							&indexup,
							&indexdown);

			for (k=0; k<grid->num_prod; k++)
			{
				tmpcb->val[xi][k] = 0.0;
			
				tmpcb->val[xi][k] += cup * grid->vt2[indexup][xi][k]
								+ cdown * grid->vt2[indexdown][xi][k];
			}
		}

		/*	Convolve slice into a temp vector using Crank-Nicholson	*/
		num_f_pde_one_step_backward	(	tmppde,
										grid->nx,
										grid->x,
										grid->num_prod,
										tmpcb->val,
										tmpcb->mu,
										tmpcb->var,
										tmpcb->r,
										0.55,
										grid->vt1[phi_slice]	);


	}/* end of the loop phi_slice */

	/*	Switch grids	*/
	temp = grid->vt1;
	grid->vt1 = grid->vt2;
	grid->vt2 = temp;
}

/*	Go backward one step	*/
void srt_f_chey_beta_bck_pde_one_step2(
CNPDE_TEMP						*tmppde,
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								t1,
/*	Next time step idx	*/
int								t2)
{
	int j, k;
	int phi_slice;
	double nvar;
	CHEYBETA_PARAM p;
	int xi;
	double x;
	double ***temp;

	double  cup,cdown;
	long    indexup, indexdown;

	/*	Get local diffusion parameters	*/
	chey_beta_mdl_param (mdl, grid->t[t1], grid->t[t2], &p);

	/* before to enter into the loop on phi pre-compute a few things */
	for (xi=0; xi<grid->nx; xi++)
	{
		/* Store the x */
		x = grid->x[xi];

		/*	Calculate r	*/
		tmpcb->r[xi] = (grid->fwdrt[t1] + x) * p.dt;
		
		/*	Calculate x coefficients var : independent of phi */
		nvar = chey_beta_mdl_norm_var_at_t (	
				&p, 
				x, 
				0.0,
				grid->fwdrt[t1],
				grid->maxvar[t1],
				grid->minvar[t1]);
	
		tmpcb->var[xi] = chey_beta_mdl_var (mdl, 
											&p, 
											x, 
											0.0, 
											nvar);

		/* reset the value to zero */
		for (j=0; j<grid->nphi[t1]; j++)
		{
			for (k=0; k<grid->num_prod; k++)
				grid->vt1[j][xi][k] = 0.0;
		}
	}

	/*	Go phi slice by phi slice	*/
	for (phi_slice=0; phi_slice<grid->nphi[t1]; phi_slice++)
	{
		/*	Calculate drifts, variances and discount factors of x and y	*/
		for (xi=0; xi<grid->nx; xi++)
		{											
			x = grid->x[xi];
			
			/*	Calculate mu of x	*/
			tmpcb->mu[xi] = chey_beta_mdl_drift (	mdl, 
													&p, 
													x, 
													grid->phi[t1][phi_slice]);

			
		}

		/* convole in x for a phi_slice */
		num_f_pde_one_step_backward	(	tmppde,
										grid->nx,
										grid->x,
										grid->num_prod,
										grid->vt2[phi_slice],
										tmpcb->mu,
										tmpcb->var,
										tmpcb->r,
										0.55,
										tmpcb->val	);


		for (xi=0; xi<grid->nx; xi++)
		{
			x = grid->x[xi];

			lininter_quick_coef (	grid->phi[t1], 
									grid->nphi[t1], 
									chey_beta_mdl_phi_back (	mdl, 
												&p, 
												x, 
												grid->phi[t2][phi_slice],
												tmpcb->var[xi] ),
									&cup,
									&cdown,
									&indexup,
									&indexdown);

			for (k=0; k<grid->num_prod; k++)
			{			
					grid->vt1[indexup][xi][k] += cup * tmpcb->val[xi][k];
					grid->vt1[indexdown][xi][k] += cdown * tmpcb->val[xi][k];
			}
		}

	}/* end of the loop phi_slice */

	/*	Switch grids	*/
	temp = grid->vt1;
	grid->vt1 = grid->vt2;
	grid->vt2 = temp;
}

/*	Go backward one step	*/
void srt_f_chey_beta_bck_pde_one_step1(
CNPDE_TEMP						*tmppde,
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								t1,
/*	Next time step idx	*/
int								t2)
{
	int j, k;
	int phi_slice;
	double nvar;
	CHEYBETA_PARAM p;
	int xi;
	double x;
	double ***temp;

	/*	Get local diffusion parameters	*/
	chey_beta_mdl_param (mdl, grid->t[t1], grid->t[t2], &p);

	/*	Go phi slice by phi slice	*/
	for (phi_slice=0; phi_slice<grid->nphi[t1]; phi_slice++)
	{
		/*	Calculate drifts, variances and discount factors of x and y	*/
		for (xi=0; xi<grid->nx; xi++)
		{											
			x = grid->x[xi];
			
			/*	Calculate x coefficients	*/
			nvar = chey_beta_mdl_norm_var_at_t (	
				&p, 
				x, 
				grid->phi[t1][phi_slice],
				grid->fwdrt[t1],
				grid->maxvar[t1],
				grid->minvar[t1]	);

			/*	Calculate mu and var of x	*/

			tmpcb->var[xi] = chey_beta_mdl_var (	mdl, 
													&p, 
													x, 
													grid->phi[t1][phi_slice], 
													nvar);

			tmpcb->mu[xi] = chey_beta_mdl_drift (	mdl, 
													&p, 
													x, 
													grid->phi[t1][phi_slice]);

			/*	Calculate r	*/
			tmpcb->r[xi] = (grid->fwdrt[t1] + x) * p.dt;

			/*	Do explicit convolution in phi	*/
			lininter_coef (	grid->phi[t2], 
							grid->nphi[t2], 
							chey_beta_mdl_phi (	mdl, 
												&p, 
												x, 
												grid->phi[t1][phi_slice],
												tmpcb->var[xi] ),
							tmpcb->coef);

			for (k=0; k<grid->num_prod; k++)
			{
				tmpcb->val[xi][k] = 0.0;
			
				for (j=0; j<grid->nphi[t2]; j++)
				{
					tmpcb->val[xi][k] += tmpcb->coef[j] * grid->vt2[j][xi][k];
				}
			}
		}

		/*	Convolve slice into a temp vector using Crank-Nicholson	*/
		num_f_pde_one_step_backward	(	tmppde,
										grid->nx,
										grid->x,
										grid->num_prod,
										tmpcb->val,
										tmpcb->mu,
										tmpcb->var,
										tmpcb->r,
										0.55,
										grid->vt1[phi_slice]	);

	}

	/*	Switch grids	*/
	temp = grid->vt1;
	grid->vt1 = grid->vt2;
	grid->vt2 = temp;
}

/*	Go backward n steps	*/
void srt_f_chey_beta_bck_pde_n_steps(
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								*ti,
/*	Number of steps backward	*/
int								n)
{
	int i;
	CNPDE_TEMP tmppde;

	/*	Init PDE	*/
	num_f_pde_init (&tmppde, grid->nx, grid->num_prod);

	for (i=0; i<n; i++)
	{
		srt_f_chey_beta_bck_pde_one_step (&tmppde, tmpcb, mdl, grid, *ti-i-1, *ti-i);
	}

	*ti -= n;

	/*	Free PDE	*/
	num_f_pde_free (&tmppde, grid->nx, grid->num_prod);
}