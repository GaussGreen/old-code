/* ====================================================================================

   FILE NAME : 	   NUM_H_PCA.H

   PURPOSE:        Perform PCA on a covariance matrix
                   Project the factors on a two dimensional space

   ====================================================================================
 */

#ifndef NUM_H_PCA_H
#define NUM_H_PCA_H

Err principal_component_analysis(double **covar, int n, double eigen_val[],
                                 double **vol_func);

/*
   This function takes as input a REAL COVARIANCE matrix (and its size)
   It returns several results when projecting the covariance information
   on the main two axes (given by the main two eigen values):
        - the explained percentage of variance expl_per
        - the first two volatilities of the factors (sqrt of the eigen values)
        - the function that will be used to fit a two factor model      ,
          that defines the correlation (it is linked to correlation by the
          following relationship:

                               1 + func(x)*func(y)
                cor(x      ,y) = ------------------------------------
                           sqrt(1+func(x)^2)*sqrt(1+func(y)^2)

        - the "new" two factor correlation matrix after projection on the
          two dimensionnal plane (two_fac_corr)
          given by the previous equation for each value of x and y
*/

Err two_factor_projection(double **real_covar, int n, double *expl_perc,
                          double *vol_fac1, double *vol_fac2, double fit_func[],
                          double **two_fac_corr);

#endif
