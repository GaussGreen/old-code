
/* ==========================================================================
   FILE_NAME:	Fx3FUtils.c

   DATE:		12/20/00
   
   AUTHOR:		L.C.
   ========================================================================== */

#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "num_h_diagonalise.h"
#include "math.h"

long Get_Index(double T, double *Maturity, long nbrMat)
{
static int i;

	T -= 1.0E-10;

	for (i=0 ; (i < nbrMat) && (Maturity[i] < T) ; i++);

	if (i >= nbrMat)
		i = i-1;

	return i;
}

long Get_Index_Rec(double T, double *Maturity, long nbrMat, long StartIndex)
{
static int i;

	for (i=StartIndex ; (i < nbrMat) && (Maturity[i] < T) ; i++);

	if (i >= nbrMat)
		i = i-1;

	return i;
}

Err get_transf_matrix(
							double *sigma, 
							double *correl, 
							double **transf_matrix, 
							double **transf_matrix_inv, 
							double **covar_matrix,
							double **eigen_vector, 
							double **eigen_vector_transp, 
							double *eigen_val)
{
static	unsigned short i, j, gi, gj, gk;
static	double gsum;
static	Err err = NULL;

	/*	Construct Covariance Matrix */
	for (i=0; i<3; i++)
	{
		covar_matrix[i][i] = sigma[i] * sigma[i];
		for (j=i+1; j<3; j++)
		{
			covar_matrix[i][j] = sigma[i] * sigma[j] * correl[i+j-1];
			covar_matrix[j][i] = covar_matrix[i][j];
		}
	}
	
	/*	Diagonalise Covariance Matrix */
	err = diagonalise_symmetric_matrix(
										covar_matrix,
										3,
										eigen_val,
										eigen_vector);
				
	if (err)
	{
		return err;
	}

	/*	Calculation of the transformation matrix: from Dual to Original Base */
	for (i=0; i<3 ; i++)
	{
		if (eigen_val[i] < 1.0e-16)
		{
			err = "Covariance matrix is NOT positive definite - check your correlations";
			return err;
		}

		eigen_val[i] = sqrt (eigen_val[i]);
	}
		
	prod_matrix_diago(	eigen_vector,
						eigen_val,
						transf_matrix,
						gi,
						gj,
						gsum);
		
	/*	Calculation of the transformation matrix inverse : from Original to Dual Base */
	for (i=0; i<3 ; i++)
	{
		eigen_val[i] = 1.0 / eigen_val[i];
	}
	
	transpo_matrix(	eigen_vector,
					eigen_vector_transp,
					gi,
					gj);	

	prod_diago_matrix(	eigen_vector_transp,
						eigen_val,
						transf_matrix_inv,
						gi,
						gj,
						gsum);
			
	return err;
}

Err get_lambda_from_ir_ts(
						TermStruct *ts,
						double *lambda)
{				
	SrtLst			*l;
	long			n;
	double			tau;	

	l = ts -> head;
	n = 0;
	tau = 0.0;

	while (l) 
	{	
		if ((((IrmTermStructVal*)l->element->val.pval)->val_origin == TAU_DATE) ||
		(((IrmTermStructVal*)l->element->val.pval)->val_origin == BOTH_DATE))
		{	
			if (fabs (((IrmTermStructVal*)l->element->val.pval) -> tau - tau) > EPS)
			{
				tau = ((IrmTermStructVal*)l->element->val.pval) -> tau;
				n++;
			}
		}
		
		l = l -> next;
	}

	if (n != 1)
	{
		return "Term structure contains maore than one tau in get_lambda_from_ir_ts";
	}

	*lambda = 1.0 / tau;
	return NULL;
}

Err compute_vol_times (
						char	*und3dfx, 
						int		*num_vol_times, 
						double	**vol_times,
						double	last_time)
{

	int sig_n_dom, tau_n_dom, sig_n_for, tau_n_for, sig_n_fx, nb_corr;
	
	double *sig_date_dom = NULL,
		   *sig_dom = NULL,
		   *tau_date_dom = NULL,
		   *tau_dom = NULL,
		   *sig_date_for = NULL,
		   *sig_for = NULL,
		   *tau_date_for = NULL,
		   *tau_for = NULL,
		   *sig_date_fx = NULL,
		   *sig_fx = NULL,
		   *mat_c = NULL,
		   *c1 = NULL,
		   *c2 = NULL,
		   *c3 = NULL;	

	int i;
	Err err = NULL;

	err = Get_FX_StochRate_TermStructures_corr (
		und3dfx,
		&sig_date_dom,
		&sig_dom,
		&sig_n_dom,
		&tau_date_dom,
		&tau_dom,
		&tau_n_dom,
		&sig_date_for,
		&sig_for,
		&sig_n_for,
		&tau_date_for,
		&tau_for,
		&tau_n_for,
		&sig_date_fx,
		&sig_fx,
		&sig_n_fx,
		&mat_c,
		&c1,
		&c2,
		&c3,
		&nb_corr);

	if (err)
	{
		goto FREE_RETURN;
	}

	*vol_times = (double*) calloc (sig_n_dom + sig_n_for + sig_n_fx, sizeof (double));
	if (!vol_times)
	{
		err = "Memory allocation error in compute_vol_times";
		goto FREE_RETURN;
	}

	*num_vol_times = 0;

	i = 1;

	while ((i<sig_n_dom) && (sig_date_dom[i-1] < last_time))
	{
		if (fabs (sig_dom[i] - sig_dom[i-1]) > EPS)
		{
			(*vol_times)[*num_vol_times] = sig_date_dom[i-1];
			(*num_vol_times)++;
		}
		i++;
	}

	i = 1;
	while ((i<sig_n_for) && (sig_date_for[i-1] < last_time))
	{
		if (fabs (sig_for[i] - sig_for[i-1]) > EPS)
		{
			(*vol_times)[*num_vol_times] = sig_date_for[i-1];
			(*num_vol_times)++;
		}
		i++;
	}

	i = 1;
	while ((i<sig_n_fx) && (sig_date_fx[i-1] < last_time))
	{
		if (fabs (sig_fx[i] - sig_fx[i-1]) > EPS)
		{
			(*vol_times)[*num_vol_times] = sig_date_fx[i-1];
			(*num_vol_times)++;
		}
		i++;
	}

	if (!(*num_vol_times))
	{
		free (*vol_times);
		*vol_times = NULL;
	}

FREE_RETURN:
	
	if (sig_date_dom) free (sig_date_dom);
	if (sig_dom) free (sig_dom);
	if (tau_date_dom) free (tau_date_dom);
	if (tau_dom) free (tau_dom);
	if (sig_date_for) free (sig_date_for);
	if (sig_for) free (sig_for);
	if (tau_date_for) free (tau_date_for);
	if (tau_for) free (tau_for);
	if (sig_date_fx) free (sig_date_fx);
	if (sig_fx) free (sig_fx);
	if (mat_c) free (mat_c);
	if (c1) free (c1);
	if (c2) free (c2);
	if (c3) free (c3);

	return err;
}

/*	Fill the time vector */
Err fill_time_vector(
					double				**time, 
					int					*nstp, 
					int					num_bar_times,
					double				*bar_times, 
					int					num_vol_times, 
					double				*vol_times, 
					int					target_nstp)
{
	Err				err = NULL;

	/*	Add today if required */
	if ((*time)[0] < -EPS)
	{
		err = "Past event date in SrtGrfn3DFXTree";
		goto FREE_RETURN;
	}
	if ((*time)[0] > EPS)
	{
		num_f_add_number (nstp, time, 0.0);
		num_f_sort_vector (*nstp, *time);
		num_f_unique_vector (nstp, *time);
	}

	/*	If only one event today, add empty event */
	if (*nstp == 1)
	{
		num_f_add_number (nstp, time, 1.0);		
	}

	/*	Fill the vector */		
	/*	New algorithm */
	num_f_fill_vector_newalgo (nstp, time, target_nstp);

	/*	Old algorithm */
/*	num_f_fill_vector (*nstp, time, target_nstp, 0);
	*nstp = (*nstp > target_nstp? *nstp: target_nstp);
	num_f_fill_vector_maxtime (nstp, time, 1.1 * (*time)[*nstp-1]/(*nstp-1)); */

FREE_RETURN:

	return err;
}

Err find_phi(
			double						*date, 
			double						*phi, 
			int							nstp,
			double						theDate, 
			double						*phi_at_t)
{
	int i;

	for (i=0; i<nstp; i++)
	{
		if (fabs (date[i] - theDate) < 1.0e-08)
		{
			*phi_at_t = phi[i];
			return NULL;
		}
	}

	return serror ("Cannot find phi at date %d", theDate);
}

/*	Empirical covariance adjustment */
void adjust_emp_covar(
							double		**g, 
							int			npth, 
							int			n)
{
	int						i, j, k;
	double					**cov		= NULL, 
							**chol		= NULL,
							*temp		= NULL,
							**invchol	= NULL;
	double					sum;

	/*	Allocate */
	cov = dmatrix (0, n-1, 0, n-1);
	chol = dmatrix (0, n-1, 0, n-1);
	invchol = dmatrix (0, n-1, 0, n-1);
	temp = dvector (0, n-1);

	if (!cov || !chol || !invchol || !temp)
	{
		goto FREE_RETURN;
	}

	/*	Calculate empirical covar */
	for (i=0; i<n; i++)
	{
		for (j=i; j<n; j++)
		{
			sum = 0.0;
			for (k=0; k<npth; k++)
			{
				sum += g[k][i] * g[k][j];
			}
			cov[i][j] = cov[j][i] = sum / npth;
		}
	}

	/*	Calculate its Cholesky */
	nr_choldc (n, cov, chol);

	/*	Calculate the inverse */
	inverse_lower_triangular_matrix (chol, n, invchol);

	/*	Multiply on the right by chol-1' (lower triangular) */
	for (i=0; i<npth; i++)
	{
		for (j=0; j<n; j++)
		{
			sum = 0.0;
			for (k=0; k<=j; k++)
			{
				sum += g[i][k] * invchol[j][k];
			}
			temp[j] = sum;
		}
		memcpy (g[i], temp, n * sizeof (double));
	}

FREE_RETURN:

	if (cov) free_dmatrix (cov, 0, n-1, 0, n-1);
	if (chol) free_dmatrix (chol, 0, n-1, 0, n-1);
	if (temp) free_dvector (temp, 0, n-1);
	if (invchol) free_dmatrix (invchol, 0, n-1, 0, n-1);
}

/*	Generate BalSam numbers */
Err balsam_generation(long nbPaths, /* Must be odd */
					  long nbSteps, 
					  double **matrix)
{
double	*init_gauss = NULL;
double	step, prob;
long	i, j, seed;
int		rand;

	init_gauss = dvector(0, nbPaths-1);
	nbPaths -= 1;
	nbPaths /= 2;

	if (!init_gauss)
		return "Memory allocation failure in balsam_generation";

	step = 0.5 / (nbPaths + 1);
	prob = step;

	seed = -123456789;

	/* Generation of the fractiles of the gaussian */

	for (i=0; i<nbPaths; i++)
	{
		init_gauss[i] = inv_cumnorm_fast(prob);
		init_gauss[nbPaths+i+1] = -init_gauss[i];
		prob += step;
	}
	init_gauss[nbPaths] = 0.0;

	/* shuffle it */
	nbPaths *= 2;
	nbPaths += 1;
	for (j=0; j<nbSteps; j++)
	{
		for (i=0; i<nbPaths-1; i++)
		{
			/* rand = random_int(nbPaths-1-i, &seed) + i; */
			rand = i + (int) ((nbPaths-i) * uniform (&seed));
			matrix[i][j] = init_gauss[rand];
			init_gauss[rand] = init_gauss[i];
			init_gauss[i] = matrix[i][j];
		}
		matrix[nbPaths-1][j] = init_gauss[nbPaths-1];
	}	
	
	if (init_gauss)
	{
		free_dvector(init_gauss, 0, nbPaths-1);
	}

	return NULL;
}

Err correlate_variable_corr(
							double	*dom_std,
							double	*for_std,
							double	*fx_std,
							double	*corr12,
							double	*corr13,
							double	*corr23,
							long	nbPaths,
							long	nbSteps,
							double	**matrix)
{
long	i, j;
double	*a = NULL,
		*b = NULL,
		*c = NULL,
		*d = NULL,
		*e = NULL,
		*f = NULL;
double  x, y, z;

Err		err = NULL;

	a = (double*) calloc (nbSteps, sizeof (double));
	b = (double*) calloc (nbSteps, sizeof (double));
	c = (double*) calloc (nbSteps, sizeof (double));
	d = (double*) calloc (nbSteps, sizeof (double));
	e = (double*) calloc (nbSteps, sizeof (double));
	f = (double*) calloc (nbSteps, sizeof (double));

	if (!a || !b || !c || !d || !e || !f)
	{
		err = "Memory allocation error in Covarate variables";
		goto FREE_RETURN;
	}

	for (j=0; j<nbSteps; j++)
	{
		a[j] = dom_std[j];
		b[j] = corr12[j] * dom_std[j] * for_std[j] / a[j];
		c[j] = for_std[j] * for_std[j] - b[j] * b[j];
		if (c[j] < 0)
		{
			return "Check your correlations";
		}
		c[j] = sqrt(c[j]);

		d[j] = corr13[j] * dom_std[j] * fx_std[j] / a[j];
		
		if (c[j] != 0)
		{
			e[j] = (corr23[j] * for_std[j] * fx_std[j] - b[j] * d[j]) / c[j];
			f[j] = fx_std[j] * fx_std[j] - d[j] * d[j] - e[j] * e[j];
			if (f[j] < 0)
			{
				return "Check your correlations";
			}
			f[j] = sqrt(f[j]);
		}
		else
		{
			if (corr13[j] * dom_std[j] * fx_std[j] == corr23[j] * for_std[j] * fx_std[j])
			{
				e[j] = fx_std[j] * fx_std[j] - d[j] * d[j];
				f[j] = 0;
				if (e[j] < 0)
				{
					return "Check your correlations";
				}
				e[j] = sqrt(e[j]);
			}
			else
			{
				return "Check your correlations";
			}
		}
	}
	
	for (i=0; i<nbPaths; i++)
	{
		for (j=0; j<nbSteps; j++)
		{
			x = matrix[i][3*j];
			y = matrix[i][3*j+1];
			z = matrix[i][3*j+2];

			matrix[i][3*j]= a[j] * x;
			matrix[i][3*j+1] = b[j] * x + c[j] * y;
			matrix[i][3*j+2] = d[j] * x + e[j] * y + f[j] * z;
		}
	}

FREE_RETURN:

	if (a) free(a);
	if (b) free(b);
	if (c) free(c);
	if (d) free(d);
	if (e) free(e);
	if (f) free(f);

	return err;
}

Err correlate_variable_covar(	double	*dom_std,
								double	*for_std,
								double	*fx_std,
								double	*covar12,
								double	*covar13,
								double	*covar23,
								long	nbPaths,
								long	nbSteps,
								double	**matrix)
{
long	i, j;
double	*a = NULL,
		*b = NULL,
		*c = NULL,
		*d = NULL,
		*e = NULL,
		*f = NULL;
double  x, y, z;

Err		err = NULL;

	a = (double*) calloc (nbSteps, sizeof (double));
	b = (double*) calloc (nbSteps, sizeof (double));
	c = (double*) calloc (nbSteps, sizeof (double));
	d = (double*) calloc (nbSteps, sizeof (double));
	e = (double*) calloc (nbSteps, sizeof (double));
	f = (double*) calloc (nbSteps, sizeof (double));

	if (!a || !b || !c || !d || !e || !f)
	{
		err = "Memory allocation error in Covarate variables";
		goto FREE_RETURN;
	}

	for (j=0; j<nbSteps; j++)
	{
		a[j] = dom_std[j];
		b[j] = covar12[j] / a[j];
		c[j] = for_std[j] * for_std[j] - b[j] * b[j];
		if (c[j] < 0)
		{
			return "Check your correlations";
		}
		c[j] = sqrt(c[j]);

		d[j] = covar13[j] / a[j];
		
		if (c[j] != 0)
		{
			e[j] = (covar23[j] - b[j] * d[j]) / c[j];
			f[j] = fx_std[j] * fx_std[j] - d[j] * d[j] - e[j] * e[j];
			if (f[j] < 0)
			{
				return "Check your correlations";
			}
			f[j] = sqrt(f[j]);
		}
		else
		{
			if (covar13[j] == covar23[j])
			{
				e[j] = fx_std[j] * fx_std[j] - d[j] * d[j];
				f[j] = 0;
				if (e[j] < 0)
				{
					return "Check your correlations";
				}
				e[j] = sqrt(e[j]);
			}
			else
			{
				return "Check your correlations";
			}
		}
	}
	
	for (i=0; i<nbPaths; i++)
	{
		for (j=0; j<nbSteps; j++)
		{
			x = matrix[i][3*j];
			y = matrix[i][3*j+1];
			z = matrix[i][3*j+2];

			matrix[i][3*j]= a[j] * x;
			matrix[i][3*j+1] = b[j] * x + c[j] * y;
			matrix[i][3*j+2] = d[j] * x + e[j] * y + f[j] * z;
		}
	}

FREE_RETURN:

	if (a) free(a);
	if (b) free(b);
	if (c) free(c);
	if (d) free(d);
	if (e) free(e);
	if (f) free(f);

	return err;
}