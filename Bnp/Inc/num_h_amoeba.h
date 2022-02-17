#ifndef NUM_H_AMOEBA_H
#define NUM_H_AMOEBA_H

Err amoeba(
		double   **p,                /* [1..ndim+1][1..ndim]: starting points*/ 
		double   y[],                /* f(x_i) with x_i = row of p */
		int      ndim,               /* Dimension of x */
		double   ftol,               /* Fractionnal convergence tolerance */
		int      maxiter,            /* Maximum number of function eval */ 
		double   (*funk)(double []), /* Function to minimise */
		int      *nfunk);            /* Number of function evaluations */

#endif
