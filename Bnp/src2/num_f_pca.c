/* ====================================================================================

   FILE NAME : 	   NUM_F_PCA.C 

   PURPOSE:        Perform PCA on a covariance matrix
                   Project the factors on a two dimensional space

   ==================================================================================== */


#include "num_h_allhdr.h"             
#include "num_h_pca.h"
#include "math.h"


#define 	TINY 			1.0e-20;
#define 	MAX_MATRIX_SIZE 	50
#define 	TINY_EPSILON		1e-50

/* ==========================================================================
			A useful function for a square matrix 
   ========================================================================== */

/*
Matrixes are described as follows:

	                     col=i
	                   _____________
	                  |		|
     	                  | data[i][j]  |
	                  |		|
	                  |		|
	                  |		|
	                  |-------------|
			      vector[i]
*/

static void copy_matrix( double **input,int n, double **output)
/* Returns the copy of  matrix input into the matrix ouput */
{
	int i,j;	

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
			output[i][j] = input[i][j];
	}
	return;
}



/* See paper Principal Component Analysis - SORT Fixed Income */
	
Err principal_component_analysis(double **covar, int n,
			double eigen_val[], double **vol_func)
{
	int i,j, nrot;
	Err err=NULL;
	double sqrt_eg_val,**eigen_vec,**copy;
	
	copy = dmatrix(0,n-1,0,n-1);   
	copy_matrix(covar, n, copy);
	eigen_vec = dmatrix(0,n-1,0,n-1);   

/* 
   We need the eigen_values and associated eigen vectors from the 
   covariance matrix
*/
	if ( err=jacobi_diagonalisation(copy, n, eigen_val, eigen_vec, &nrot))
		return err;
/* 
   From the sorted eigen vectors, we build a possible set of 
   volatlity functions (defined up to a rotation) by multiplying 
   the eigen vectors by the square_root of the corresponding eigen value
*/
	for (j=0;j<n;j++)
	{
		sqrt_eg_val = sqrt(eigen_val[j]);
		for(i=0;i<n;i++)
			vol_func[i][j] = sqrt_eg_val * eigen_vec[i][j];
	}
		
	free_dmatrix(copy,0,n-1,0,n-1);
	free_dmatrix(eigen_vec,0,n-1,0,n-1);
	return err;
 
}

/* ========================================================================= */

/*
   This function takes as input a REAL COVARIANCE matrix (and its size)
   It returns several results when projecting the covariance information
   on the main two axes (given by the main two eigen values):
	- the explained percentage of variance expl_per 
	- the first two volatilities of the factors (sqrt of the eigen values)
	- the function that will be used to fit a two factor model,
	  that defines the correlation (it is linked to correlation by the
	  following relationship:
			         
		               1 + func(x)*func(y)
		cor(x,y) = ------------------------------------
			   sqrt(1+func(x)^2)*sqrt(1+func(y)^2)

	- the "new" two factor covariance matrix after projection on the 
	  two dimensionnal plane ( two_fac_covar)                
	  given by the previous equation for each value of x and y
*/

Err two_factor_projection(double **real_covar,int n,
		double *expl_perc, double *vol_fac1, double *vol_fac2,
		double mkt_fit_func[], double **two_fac_correl)
{
	Err err = NULL;
	int i,j;
 	double sgn;	
	double *eigen_val;
	double **vol_func;
	double sum = 0.0;

	eigen_val = dvector(0,n-1);
	vol_func = dmatrix(0,n-1,0,n-1);

	principal_component_analysis(real_covar, n, eigen_val, vol_func);

	if (eigen_val[0] >= 0.0)
		*vol_fac1 = sqrt(eigen_val[0]);
	else                                      
		return serror("First eigen value is negative...");
	if (eigen_val[1] >= 0.0)
		*vol_fac2 = sqrt(eigen_val[1]);
	else
		return serror("Second eigen value is negative...");

	for (i=0;i<n;i++)
		sum += eigen_val[i];
	*expl_perc = (eigen_val[0] + eigen_val[1])  / sum;
	
        if (vol_func[1][0] < 0.00)
		sgn = -1.0;
	else
		sgn = 1.0;
	for(i=0;i<n;i++)
	{	
	    if ( fabs(vol_func[0][i]) > TINY_EPSILON )
		mkt_fit_func[i] = sgn * vol_func[1][i]/vol_func[0][i];
	    else
		return serror("The first factor has zero vol on %d th date",i+1);

	}
	for (i=0;i<n;i++)
	{
	    for (j=0;j<n;j++)
		two_fac_correl[i][j] = (1+ mkt_fit_func[i]*mkt_fit_func[j]) /
			sqrt( (1+mkt_fit_func[i]*mkt_fit_func[i]) *
			       (1+mkt_fit_func[j]*mkt_fit_func[j]) )	;
        }

	free_dmatrix(vol_func,0,n-1,0,n-1);
	free_dvector(eigen_val,0,n-1);

	return err;

}

