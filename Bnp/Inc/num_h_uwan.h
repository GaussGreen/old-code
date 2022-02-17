/* -------------------------------------------------------------------------
   FILENAME		: num_h_uwan.h

   PURPOSE		: generate randon numbers that are uniformely
   distributed , but with an antithetic noise
   ------------------------------------------------------------------------- */
#ifndef NUM_H_UWAN_H
#define NUM_H_UWAN_H

/* Initialisation function: path      , brownian      , step indexes (low and
 * high values) */
Err uwan_init(long pl, long ph, long rl, long rh, long sl, long sh);

/* Generation of Gaussian increments using UWAN (for r[Brow][Step]) */
Err uwan_matrix(double **v, long rl, long rh, long sl, long sh, long *seed);

Err uwan_cube(double ***v, long pl, long ph, long rl, long rh, long sl, long sh,
              long *seed);

/* Destruction */
Err uwan_free(void);

#endif
