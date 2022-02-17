
/* ------------------------------------------------------------------------
        FILENAME: 	num_h_annealing.c

        FUNCTION:   simulated_annealing
   ------------------------------------------------------------------------ */

#ifndef NUM_H_ANNEALING_H
#define NUM_H_ANNEALING_H

typedef enum AnnealMethod { PROPORTIONAL_EPSILON, POWER_ALPHA } AnnealMethod;

Err interp_annealing_method(String sMethod, AnnealMethod *pMethod);

Err simulated_annealing(double **p, /* From [1] to [ndim+1]  * [1] to [ndim] */
                        double *y,  /* From [1] to [ndim] */
                        long ndim, double ftol, long niter,
                        double (*funcs)(double *), double max_temp,
                        double min_temp, double decr_factor);

#endif
