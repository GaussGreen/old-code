/* ======================================================
   FILENAME:  num_h_matrix.h
   
   PURPOSE:   Matrix functions
   ====================================================== */

#ifndef NUM_H_MATRIX_H
#define NUM_H_MATRIX_H

/* returns tranpose (x) */
double **transpose_matrix (double **old_matrix, long old_rl, long old_rh, long old_cl, long old_ch);

/* returns xy */
double **product_matrix (double **x, long rxl, long rxh, long cxl, long cxh, 
						 double **y, long ryl, long ryh, long cyl, long cyh);

/* returns inv(x) */
double **inverse_matrix (double **x, long rl, long rh);

/*	Cholesky decomposition */
void nr_choldc (
int							n, /* 0.. n-1*/
double						**cov,
double						**chol);

/*	Cholesky decomposition with error checking */
Err choldc (int n, double **a, double **chol);

/*	Inverse Lower Trinagular Matrix */
void inverse_lower_triangular_matrix (
							double **mat, 
							int n,	/*	0..n-1 */ 
							double **res);

Err PositiveMatrix(double **Matrix,
				   int n
				   );

#endif
