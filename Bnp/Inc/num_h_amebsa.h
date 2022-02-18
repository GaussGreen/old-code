/* ------------------------------------------------------------------------
                            SIMULATED ANNEALING
   ------------------------------------------------------------------------ */

#ifndef NUM_H_AMEBA_H
#define NUM_H_AMEBA_H

void amebsa(
    double** p,    /* Starting points of the simplex */
    double   y[],  /* Values at the starting points */
    int      ndim, /* Number of dimension */
    double   pb[], /* Best point ever encountered */
    double*  yb,   /* Best function value encountered */
    double   ftol, /* Fractional convergence tolerance */
    double (*funk)(double[]),
    int*   iter,    /* Number of function evaluations */
    double temptr); /* Temperature */

#endif
