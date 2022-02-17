/* ======================================================
   FILENAME:  num_h_random.h

   PURPOSE:   Random number generations
   ====================================================== */

#ifndef NUM_H_RANDOM_H
#define NUM_H_RANDOM_H

#define RANDINIT -123456789

double uniform(long *seed);
double uniform_fast(long *seed);

/* To generate random independent Gaussian variables */
double gauss_sample(long *seed);
Err gauss_vector(double *v, long vl, long vh, long *seed);
Err gauss_matrix(double **v, long rl, long rh, long cl, long ch, long *seed);
Err gauss_cube(double ***v, long pl, long ph, long bl, long bh, long sl,
               long sh, long *seed);

Err gauss_anti_cube(double ***v, long pl, long ph, long bl, long bh, long sl,
                    long sh, long *seed);

Err gauss_box_muller_init(void);
double gauss_box_muller(long *seed);

int random_int(int n, long *seed);

/* ======================================================================== */

/* Correlation of independent Gaussian increments */

Err correl_random(double **r, long sl, long sh, long nbr, double **coeff);
Err compute_coeff_from_correl(double **correl, long nbr, double **coeff);

Err correl_random_eigen(double **r, long sl, long sh, long nbr, double **coeff);
Err compute_eigen_from_correl(double **correl, long nbr, double **coeff);

/* ======================================================================== */

/* ------------------------------------------------------------------------
   BALANCED SAMPLING
   ------------------------------------------------------------------------ */

Err BalSampMatrix(double **m, int pathindexstart, int pathindexend,
                  int stepindexstart, int stepindexend, long *seed);

Err BalSampCube(double ***rand, long pathindexstart, long pathindexend,
                long browindexstart, long browindexend, long stepindexstart,
                long stepindexend, long *seed);

/* ======================================================================== */

#endif
