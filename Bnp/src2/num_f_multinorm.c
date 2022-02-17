/* ========================================================
   FILENAME:  num_f_multinorm.c
   
   PURPOSE:   Calculate multinomial Gaussian distributions
   ======================================================== */

#include "num_h_allhdr.h"
#include "time.h"

/*	Generation of Sobol Matrix */
double ***gen_sobol_mat(
int							n,				/*	Dimension */
double						**cov,			/*	Covariance matrix */
int							*dir,			/*	-1: below, +1: above */
int							npth)			/*	Number of Sobol paths */
{
	int						i, j, k;
	double					***g = NULL, **chol = NULL, *temp = NULL;
	double					sum;
	
	/*	Allocate memory */
	g = (double***) f3tensor (0, npth-1, 0, 0, 0, n-1);
	/*	Calculate cholesky */
	chol = dmatrix (0, n-1, 0, n-1);
	temp = dvector (0, n-1);
	nr_choldc (n, cov, chol);

	if (!g || !chol)
	{
		goto FREE_RETURN;
	}

	/*	Initialise Sobol sequences */
	sobol_init (0, npth-1, 0, 0, 0, n-1);
	sobol_cube (g, 0, npth-1, 0, 0, 0, n-1);

	/*	Adjust chol for direction */
	for (j=0; j<n; j++)
	{
		for (k=0; k<=j; k++)
		{
			chol[j][k] *= dir[j];
		}
	}	

	/*	Incorporate covariance */
	for (i=0; i<npth; i++)
	{
		for (j=0; j<n; j++)
		{
			sum = 0.0;
			for (k=0; k<=j; k++)
			{
				sum += chol[j][k] * g[i][0][k];
			}
			temp[j] = sum;
		}
		memcpy (g[i][0], temp, n * sizeof (double));
	}

FREE_RETURN:
	
	if (chol) free_dmatrix (chol, 0, n-1, 0, n-1);
	if (temp) free_dvector (temp, 0, n-1);

	return g;
}

/*	Multinormal */
double num_f_multinorm(
int							n,				/*	Dimension */
double						*e,				/*	Expectations */
double						**cov,			/*	Covariance matrix */
double						*x,				/*	Strikes */
int							*dir,			/*	-1: below, +1: above */
int							npth)			/*	Number of Sobol paths */
{
	int						i, j;
	double					*ux = NULL, ***g = NULL;
	double					in, ave = 0.0;
	clock_t					t1, t2;

	/*	Allocate memory */
	ux = (double*) calloc (n, sizeof (double));
	
	t1 = clock();
	g = gen_sobol_mat (n, cov, dir, npth);
	t2 = clock();
	smessage ("Generation time: %.2f sec", ((double) (t2-t1)) / CLOCKS_PER_SEC);

	if (!ux || !g)
	{
		goto FREE_RETURN;
	}

	/*	Initialise strikes */
	for (j=0; j<n; j++)
	{
		ux[j] = dir[j] * (x[j] - e[j]);
	}	

	/*	Go through the paths */
	ave = 0.0;
	for (i=0; i<npth-1; i++)
	{		
		/*	Calculate proba */
		in = 1.0;
		for (j=0; j<n; j++)
		{
			if (g[i][0][j] < ux[j])
			{
				in = 0.0;
				break;
			}
		}

		ave += in;
	}

	ave /= npth;

FREE_RETURN:

	if (ux) free (ux);
	if (g) free_f3tensor (g, 0, npth-1, 0, 0, 0, n-1);
	sobol_free ();

 	return ave;
}