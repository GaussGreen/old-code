/* ======================================================
   FILENAME:  num_h_sobol.h
   
   PURPOSE:   Geberate Sobol sequences (quasi random)
   ====================================================== */

#ifndef NUM_H_SOBOL_H
#define NUM_H_SOBOL_H


Err sobol_init (long pl, long ph, long rl, long rh, long sl, long sh);

Err sobol_vector(double *v, long sl, long sh);

Err sobol_matrix(double **v, long rl, long rh, long sl, long sh);

Err sobol_cube(double ***v, long pl, long ph, long rl, long rh, long sl, long sh);

Err sobol_free (void);

double **GetUnifDev(long p,long d);

void GetSobolMatrix(long p,long d,double **mtx);


#endif