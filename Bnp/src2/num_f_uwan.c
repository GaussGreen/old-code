/* -------------------------------------------------------------------------
   FILENAME		: num_f_uwan.c

   PURPOSE		: generate randon numbers that are uniformely distributed
   , but with an antithetic noise
   ------------------------------------------------------------------------- */
#include "math.h"
#include "num_h_proba.h"
#include "num_h_random.h"
#include "num_h_uwan.h"
#include "utallhdr.h"

static long ipow(long i, long j) { return DTOL(pow(i, j)); }

/* -------------------------------------------------------------------------- */

static long ipow_inv(long i, long j) {

  long r = DTOL(pow((double)i, 1.00 / j));

  /* Since 2^32 is the biggest long  , set num segments to 1 if dimension >= 32
   */
  if (j >= 32 && i > 1)
    return 1;
  while (ipow(r, j) <= i)
    r++;

  return r - 1;
}

/* -------------------------------------------------------------------------- */

static long UniformVector(double *x, long dimension, long *seed) {

  long i;

  for (i = 0; i < dimension; i++)
    x[i] = uniform_fast(seed);

  return 1;
}

/* -------------------------------------------------------------------------- */

static long Symetrise(double *x, long dimension) {

  long i;

  for (i = 0; i < dimension; i++)
    x[i] = 1 - x[i];

  return 1;
}

/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
   Private functions to perform all the operations
   -------------------------------------------------------------------------- */

static long new_path, dim, num_segment, num_square, square_num,
    *axis_index = NULL;
static double width, *random_vector = NULL;

static long UwanInit(long num_path, long dimension) {

  new_path = 0;
  dim = dimension;
  square_num = 0;

  axis_index = (long *)calloc(dim, sizeof(long));
  if (!axis_index)
    return 0;

  axis_index[0] = -1;

  random_vector = (double *)calloc(dim, sizeof(double));
  if (!random_vector) {
    free(axis_index);
    return 0;
  }

  num_square = num_path / 2;
  num_segment = ipow_inv(num_square, dim);
  num_square = ipow(num_segment, dim);

  width = 1.00 / num_segment;

  return 1;

} /* END long UwanInit(...) */

/* -------------------------------------------------------------------------- */

static long UwanNext(double *x, long *seed) {

  long i, err;

  new_path = (new_path + 1) % 2;
  square_num += new_path;

  if (square_num > num_square) {

    err = new_path ? UniformVector(random_vector, dim, seed)
                   : Symetrise(random_vector, dim);

    if (!err)
      return err;

    memcpy(x, random_vector, dim * sizeof(double));

    return 1;
  }

  if (new_path) {

    for (i = 0;; i++) {

      axis_index[i]++;

      if (axis_index[i] < num_segment)
        break;
      else
        axis_index[i] = 0;
    }

    err = UniformVector(random_vector, dim, seed);
    if (!err)
      return err;

  }

  else {

    err = Symetrise(random_vector, dim);
    if (!err)
      return err;
  }

  for (i = 0; i < dim; i++)
    x[i] = ((double)axis_index[i] + random_vector[i]) * width;

  return err;

} /* END UwanNext(...) */

/* -------------------------------------------------------------------------- */

static Err UwanFree(void) {

  if (axis_index)
    free(axis_index);
  axis_index = NULL;
  if (random_vector)
    free(random_vector);
  random_vector = NULL;
  return NULL;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Reference to UWAN sequence generator UTUWAN.C
   (FUNCTIONS FOR EXTERNAL USE)
   _________________________________________________________________________ */

static double *UWAN_RAN_NUM = NULL;
static long UWAN_BRW_NUM, UWAN_PTH_NUM, UWAN_STP_NUM, UWAN_DIM;

/* Initialisation function: path  , brownian  , step indexes (low and high
 * values) */
Err uwan_init(long pl, long ph, long rl, long rh, long sl, long sh) {

  UWAN_PTH_NUM = ph - pl + 1;
  UWAN_BRW_NUM = rh - rl + 1;
  UWAN_STP_NUM = sh - sl + 1;
  UWAN_DIM = UWAN_BRW_NUM * UWAN_STP_NUM;

  if (!UwanInit(UWAN_PTH_NUM, UWAN_DIM))
    return serror("Uwan initialisation failed");

  UWAN_RAN_NUM = dvector(0, UWAN_DIM);

  if (UWAN_RAN_NUM == NULL)
    return serror("Memory allocation error in Uwan vector");

  return NULL;

} /* END Err uwan_init() */

/* Generation of Gaussian increments using UWAN (for r[Brow][Step]) */
Err uwan_matrix(double **v, long rl, long rh, long sl, long sh, long *seed) {

  long i, j;

  /* Collect the next full path (all Brownians and steps: one point in dim DIM
   * )*/
  if (!UwanNext(UWAN_RAN_NUM, seed))
    return serror("Uwan sequence failed");

  /* Transfer the full path */
  for (i = 0; i < UWAN_BRW_NUM; i++)
    for (j = 0; j < UWAN_STP_NUM; j++)
      v[rl + i][sl + j] = inv_cumnorm_fast(UWAN_RAN_NUM[i + UWAN_BRW_NUM * j]);

  return NULL;

} /* END Err uwan_matrix() */

/* Generation of Gaussian increments using UWAN (for r[Path][Brow][Step]) */
Err uwan_cube(double ***v, long pl, long ph, long rl, long rh, long sl, long sh,
              long *seed) {

  long i, j, k;

  for (k = 0; k < UWAN_PTH_NUM; k++) {
    /* Collect the next full path (all Brownians and steps: one point in dim DIM
     * )*/
    if (!UwanNext(UWAN_RAN_NUM, seed))
      return serror("Uwan sequence failed");

    /* Transfer the full path */
    for (i = 0; i < UWAN_BRW_NUM; i++)
      for (j = 0; j < UWAN_STP_NUM; j++)
        v[rl + k][rl + i][sl + j] =
            inv_cumnorm_fast(UWAN_RAN_NUM[i + UWAN_BRW_NUM * j]);
  }

  return NULL;

} /* END Err uwan_cube() */

/* To properly free everything */

Err uwan_free(void) {

  if (!UwanFree())
    return serror("Uwan free failed");

  if (UWAN_RAN_NUM)
    free_dvector(UWAN_RAN_NUM, 0, UWAN_DIM);

  UWAN_RAN_NUM = NULL;

  return NULL;
}
