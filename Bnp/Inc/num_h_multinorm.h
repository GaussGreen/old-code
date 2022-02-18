/* ========================================================
   FILENAME:  num_h_multinorm.h

   PURPOSE:   Calculate multinomial Gaussian distributions
   ======================================================== */

#ifndef NUM_H_MULTINORM_H
#define NUM_H_MULTINORM_H

/*	Generation of Sobol Matrix */
double*** gen_sobol_mat(
    int      n,   /*	Dimension */
    double** cov, /*	Covariance matrix */
    int*     dir, /*	-1: below, +1: above */
    int      npth);    /*	Number of Sobol paths */

/*	Multinormal */
double num_f_multinorm(
    int      n,   /*	Dimension */
    double*  e,   /*	Expectations */
    double** cov, /*	Covariance matrix */
    double*  x,   /*	Strikes */
    int*     dir, /*	-1: below, +1: above */
    int      npth);    /*	Number of Sobol paths */

#endif