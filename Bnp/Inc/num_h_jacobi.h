/* ========================================================================

  MODULE:	    num_h_jacobi.h

  LIBRARY:	    NUM_LIB

  FUNCTION:	    jacobi_diagonalisation

  DESCRIPTION:	From Numerical Recipes 
                to diagonalise a symmetrix matrix

  ========================================================================== */

#ifndef NUM_H_JACOBI_H
#define NUM_H_JACOBI_H

Err jacobi_diagonalisation(double **sym_mat, int n, 
			double *eigen_val, double **eigen_vec, int *nrot);

Err jacobi_diagonalisation2(
				double  **sym_mat, 
				int     n, 
				double  *eigen_val, 
				double  **eigen_vec, 
				int     *nrot);

static Err sort_eigen_values2(double *eigen_val, double **eigen_vec, int n);

#endif
