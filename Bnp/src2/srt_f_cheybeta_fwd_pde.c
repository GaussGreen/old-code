/* ==============================================================================
   
   FILE NAME:      srt_f_cheybeta_fwd_pde.c 
   
   OBJECT:         functions implementing the local vol Cheyette forward PDE

  =============================================================================== */

#include "srt_h_all.h"
#include "math.h"

/*	Calculate the coefficients of a linear interpolation,
	i.e. find c1..cn such that f(x0) = sum ci * f(xi)	*/
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

/*	Forward PDE Functions	*/

/*	Allocate a few temp vectors	*/
void srt_f_chey_beta_fwd_pde_alloc_temp(
CHEYBETA_FWD_PDE_TEMP			*tmp,
int								nx,
int								nphi)
{
	tmp->coef = (double*) calloc (nphi+1, sizeof (double));
	tmp->ad = (double*) calloc (nx+1, sizeof (double));
	tmp->mu = (double*) calloc (nx+1, sizeof (double));
	tmp->var = (double*) calloc (nx+1, sizeof (double));
	tmp->r = (double*) calloc (nx+1, sizeof (double));
}

/*	Free temp vectors	*/
void srt_f_chey_beta_fwd_pde_free_temp(
CHEYBETA_FWD_PDE_TEMP			*tmp)
{
	free (tmp->coef);
	free (tmp->ad);
	free (tmp->mu);
	free (tmp->var);
	free (tmp->r);
}

/*	Get Arrow-Debreu prices @ first step	*/
void srt_f_chey_beta_fwd_pde_first_step(
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid)
{
	double std;
	CHEYBETA_PARAM p;
	double df;
	int i;
	double sum;
	double nvar;

	/*	Get params at time 0	*/
	chey_beta_mdl_param (mdl, grid->t[0], grid->t[1], &p);
	
	/*	Calculate std of X to time 1	*/
	nvar = chey_beta_mdl_norm_var_at_t (	&p,
											0.0,
											0.0,
											grid->fwdrt[0],
											DBL_MAX,
											1.0e-16	);

	std = sqrt (chey_beta_mdl_var			(	mdl, 
												&p, 
												0.0, 
												0.0, 
												nvar));
	/*	Calculate df to time 1	*/
	df = exp (-grid->fwdrt[0] * (grid->t[1] - grid->t[0]));

	/*	Calculate densitites at time 1	*/
	sum = 0.0;
	for (i=0; i<grid->nx; i++)
	{
		grid->adt1[0][i] = gauss (grid->x[i] / std);
		sum += grid->adt1[0][i];
	}

	/*	Renormalise	*/
	for (i=0; i<grid->nx; i++)
	{
		grid->adt1[0][i] *= df / sum;
	}
}

/*	Go forward one step	*/
void srt_f_chey_beta_fwd_pde_one_step(
CNPDE_TEMP						*tmppde,
CHEYBETA_FWD_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								t1,
/*	Next time step idx	*/
int								t2)
{
	int i, j;
	int phi_slice;
	double nvar;
	CHEYBETA_PARAM p;
	int xi;
	double x;
	double **temp;

	/*	Get local diffusion parameters	*/
	chey_beta_mdl_param (mdl, grid->t[t1], grid->t[t2], &p);

	/*	Set new AD grid to 0	*/
	for (i=0; i<grid->nphi[t2]; i++)
	{
		memset (grid->adt2[i], 0, grid->nx * sizeof (double));
	}

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
		}

		/*	Convolve slice into a temp vector using Crank-Nicholson	*/
		num_f_pde_one_step_forward	(	tmppde,
										grid->nx,
										grid->x,
										grid->adt1[phi_slice],
										tmpcb->mu,
										tmpcb->var,
										tmpcb->r,
										0.55,
										tmpcb->ad	);

		/*	Share densities between next phis	*/
		for (xi=0; xi<grid->nx; xi++)
		{
			/*	Take x in the middle	*/
			x = grid->x[xi];

			lininter_coef (	grid->phi[t2], 
							grid->nphi[t2], 
							chey_beta_mdl_phi (	mdl, 
												&p, 
												x, 
												grid->phi[t1][phi_slice],
												tmpcb->var[xi] ),
							tmpcb->coef);

			for (j=0; j<grid->nphi[t2]; j++)
			{
				grid->adt2[j][xi] += tmpcb->ad[xi] * tmpcb->coef[j];				
			}
		}
	}

	/*	Switch grids	*/
	temp = grid->adt1;
	grid->adt1 = grid->adt2;
	grid->adt2 = temp;
}

/*	Go forward n steps	*/
void srt_f_chey_beta_fwd_pde_n_steps(
CHEYBETA_FWD_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								*ti,
/*	Number of steps forward	*/
int								n)
{
	int i;
	CNPDE_TEMP tmppde;

	/*	Init PDE	*/
	num_f_pde_init (&tmppde, grid->nx, 1);

	for (i=0; i<n; i++)
	{
		srt_f_chey_beta_fwd_pde_one_step (&tmppde, tmpcb, mdl, grid, *ti+i, *ti+i+1);
	}

	*ti += n;

	/*	Free PDE	*/
	num_f_pde_free (&tmppde, grid->nx, 1);
}