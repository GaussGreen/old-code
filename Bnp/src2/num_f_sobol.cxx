/* ===============================================================
   FILE_NAME:	num_f_sobol.cxx

   PURPOSE:     Generation of Sobol sequences for quasi-random
   =============================================================== */

// S.Galluccio: 21 November 2000
// Added functions GetUnifDev and GetUnifCplDev

#include "num_h_allhdr.h"
#include "utallhdr.h"
#ifdef _DEBUG
#include "string.h"
#endif

#define MAXBIT 30

/*  Static variables		*/

static long *pol, *otpol, *deg, *npol, *cnpol, *rinit = NULL;
static long in, *ix = NULL, **iu = NULL;
static double fac;
extern long *PrimitivePolynomials[];
extern int number_of_prim_pols_of_degree[];
/* ----------------------------------------------------------------------------
 */
/* Private Functions (not to be used externally ) for the Polynomial operations
 */

/*
 *	The previous version of Sobol was using irreducible
 *	polynomials and not the primitive polynomials. We
 *  correct this here.
 *  This change however does not impact the QRS effiency
 *  and convergence        , which was observed by Antoine Savine
 *  when he implemented the first version and
 *   which underlines the difficulty of
 *  using QRS in high dimensions. Recall that all results
 *  for QRS are asympotically proved and a rough approximation
 *  of the number of points needed to get the N log(N) convergence
 *  is of order d^d where d is the dimension. Their convergence
 *  can be poorer the Pseudo-Random Sequences when truncated too early.
 *	Author: C. Godart
 *  Date: 02/01/2001
 */
static long PolInit(long dimension) {
  long polcount, poldeg;
  long nb_prim_pol_deg_smaller_than_poldeg;

  deg = (long *)calloc(dimension + 1, sizeof(long));
  otpol = (long *)calloc(dimension + 1, sizeof(long));

  if ((deg == NULL) || (otpol == NULL))

  {

    srt_free(deg);
    srt_free(otpol);

    return 0;
  }

  polcount = 1;
  poldeg = 1;
  nb_prim_pol_deg_smaller_than_poldeg = 0;

  while (polcount <= dimension)

  {
    if (polcount > nb_prim_pol_deg_smaller_than_poldeg +
                       number_of_prim_pols_of_degree[poldeg - 1]) {
      nb_prim_pol_deg_smaller_than_poldeg +=
          number_of_prim_pols_of_degree[poldeg - 1];
      poldeg++;
    }
    otpol[polcount] =
        PrimitivePolynomials[poldeg - 1][polcount - 1 -
                                         nb_prim_pol_deg_smaller_than_poldeg];
    deg[polcount] = poldeg;
    polcount++;
  }

  deg[0] = 0;

  return poldeg;
}

void select_all_init_sequences(long *init_sequence, long dimension,
                               long maxdeg);
static long RecInit(long dimension) {

  long maxdeg, j, k;
  long val, lim;
  /*
  #ifdef _DEBUG
          char init_seq_string[5096];
          char tmp[4];
  #endif
  */
  maxdeg = PolInit(dimension);

  if (!maxdeg)
    return 0;

  rinit = (long *)calloc(dimension * MAXBIT + 1, sizeof(long));
  if (rinit == NULL)

  {

    srt_free(deg);
    srt_free(otpol);
    srt_free(rinit);

    return 0;
  }

  /*select_all_init_sequences(rinit        , dimension        , maxdeg);*/

  lim = 2;

  for (j = 0; j < maxdeg; j++)

  {

    val = 1;

    for (k = 0; k < dimension; k++)

    {

      val += 2;
      rinit[1 + dimension * j + k] = val % lim;
    }

    lim <<= 1;
  }
  /*
  #ifdef _DEBUG
          for (k=0;k<maxdeg;k++) {
                  for (j=0;j<dimension;j++) {
                          sprintf(tmp        ,"%d "        ,
  rinit[j+k*dimension]); strcat(init_seq_string        , tmp);
                  }
                  sprintf(init_seq_string        ,"\0");
          }
  #endif
  */
  return maxdeg;
}
/* ----------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
   Private functions for Initialisation        , Extraction and Memory  Free
   ----------------------------------------------------------------------------
 */

static long SobInit(long dimension) {

  long j, k, l, poldeg, memerror;
  long i, ipp;

  memerror = RecInit(dimension);

  if (!memerror)
    return 0;

  ix = (long *)calloc(dimension + 1, sizeof(long));
  iu = (long **)calloc(MAXBIT + 1, sizeof(long *));

  if ((ix == NULL) || (iu == NULL))

  {

    srt_free(deg);
    srt_free(otpol);
    srt_free(rinit);

    return 0;
  }

  for (j = 1, k = 0; j <= MAXBIT; j++, k += dimension)
    iu[j] = &rinit[k];

  for (k = 1; k <= dimension; k++)

  {

    poldeg = deg[k];

    for (j = 1; j <= poldeg; j++)
      iu[j][k] <<= (MAXBIT - j);

    for (j = poldeg + 1; j <= MAXBIT; j++)

    {

      ipp = otpol[k];
      i = iu[j - poldeg][k];
      i ^= (i >> poldeg);

      for (l = poldeg - 1; l >= 1; l--)

      {

        if (ipp & 1)
          i ^= iu[j - l][k];

        ipp >>= 1;
      }

      iu[j][k] = i;
    }
  }

  fac = 1.00 / (1L << MAXBIT);
  in = dimension * dimension * dimension;

  srt_free(otpol);
  srt_free(deg);

  return 1;
}

/* ----------------------------------------------------------------------------
 */

static long SobNext(long dimension, double *otp) {

  long j, k;
  long im;

  im = in;

  for (j = 1; j <= MAXBIT; j++)

  {

    if (!(im & 1))
      break;

    im >>= 1;
  }

  im = (j - 1) * dimension;

  for (k = 1; k <= dimension; k++)

  {

    ix[k] ^= rinit[im + k];
    otp[k - 1] = (double)ix[k] * fac;
  }

  in++;

  return 1;
}

/*
 *
 *	sobol_next: we need this function as an externally callable
 *				function.
 *
 *
 */
void sobol_next(long dimension, double *otp) { SobNext(dimension, otp); }

/* ----------------------------------------------------------------------------
 */

static long SobFree(void) {

  if (rinit)
    srt_free(rinit);
  rinit = NULL;
  if (ix)
    srt_free(ix);
  ix = NULL;
  if (iu)
    srt_free(iu);
  iu = NULL;

  return 1;
}

/* ----------------------------------------------------------------------------
 */

/* Public functions (to be used from the outside world...) */

/* -------------------------------------------------------------------------
   Reference to Sobol sequence generator
   _________________________________________________________________________ */

static double *SOBOL_RAN_NUM = NULL;
static long SOBOL_BRW_NUM, SOBOL_STP_NUM, SOBOL_PTH_NUM, SOBOL_DIM;

Err sobol_init(long pl, long ph, long rl, long rh, long sl, long sh) {

  SOBOL_PTH_NUM = ph - pl + 1;
  SOBOL_BRW_NUM = rh - rl + 1;
  SOBOL_STP_NUM = sh - sl + 1;
  SOBOL_DIM = SOBOL_BRW_NUM * SOBOL_STP_NUM;

  if (!SobInit(SOBOL_DIM))
    return serror("Sobol initialisation failed");

  SOBOL_RAN_NUM = dvector(0, SOBOL_DIM);

  if (SOBOL_RAN_NUM == NULL)
    return serror("Memory allocation error in Sobol vector");

  return NULL;
}

/* --------------------------------------------------------------
   Fills in a vector of random numbers generated using Sobol
   The vector correspond to one rndom draw in the n dimensional
   space (nsteps)
   This will fill in v[sl..sh]
   -------------------------------------------------------------- */
Err sobol_vector(double *v, long sl, long sh) {

  long j;

  /* Checks the indexes are in line with the initialisation */
  if (sh - sl + 1 != SOBOL_STP_NUM)
    return ("Cannot get next Sobol: inconsitency in dimensions");

  /* Generates the next full path (all Brownians and steps: one point in dim DIM
   * )*/
  if (!SobNext(SOBOL_DIM, SOBOL_RAN_NUM))
    return serror("Sobol sequence failed");

  /* Collect and transfer the full path */
  for (j = 0; j < SOBOL_STP_NUM; j++)
    v[sl + j] = inv_cumnorm_fast(SOBOL_RAN_NUM[j]);

  return NULL;
}

/* --------------------------------------------------------------
   Fills in a matrix of random numbers generated using Sobol
   Each point correspond to one "path": all steps / all Brownians
   The order in terms of v[r][s] for the fill in is:
                take the first r        ,
                        fill in all the s        ,
                move to the next r ...
   -------------------------------------------------------------- */
Err sobol_matrix(double **v, long rl, long rh, long sl, long sh) {

  long i, j;

  /* Checks the indexes are in line with the initialisation */
  if ((sh - sl + 1 != SOBOL_STP_NUM) || (rh - rl + 1 != SOBOL_BRW_NUM))
    return ("Cannot get next Sobol: inconsitency in dimensions");

  /* Generates the next full path (all Brownians and steps: one point in dim DIM
   * )*/
  if (!SobNext(SOBOL_DIM, SOBOL_RAN_NUM))
    return serror("Sobol sequence failed");

  /* Collect and transfer the full path */
  for (i = 0; i < SOBOL_BRW_NUM; i++)
    for (j = 0; j < SOBOL_STP_NUM; j++)
      v[rl + i][sl + j] =
          inv_cumnorm_fast(SOBOL_RAN_NUM[i + SOBOL_BRW_NUM * j]);

  return NULL;
}

/* --------------------------------------------------------------
   Fills in a cube of random numbers generated using Sobol
   The cube is represented as follows:
                        rand[Path][Brow][Step]
        with:
                        pl  <= Path <= ph
                        rl  <= Brow <= rh
                        sl  <= Step <= sh

   More importance is givent to the Brownian compare to the Steps
   (for better balance between the various Brownians...)
   -------------------------------------------------------------- */

Err sobol_cube(double ***v, long pl, long ph, long rl, long rh, long sl,
               long sh) {

  long i, j, k;

  /* Checks the indexes are in line with the initialisation */
  if ((sh - sl + 1 != SOBOL_STP_NUM) || (rh - rl + 1 != SOBOL_BRW_NUM) ||
      (sh - sl + 1 != SOBOL_STP_NUM))
    return ("Cannot get next Sobol: inconsitency in dimensions");

  /* Loops on all paths to generate */
  for (k = 0; k < SOBOL_PTH_NUM; k++) {
    /* Generates the next full path (all Brownians and steps: one point in dim
     * DIM )*/
    if (!SobNext(SOBOL_DIM, SOBOL_RAN_NUM))
      return serror("Sobol sequence failed");

    /* Collect and transfer the full path */
    for (i = 0; i < SOBOL_BRW_NUM; i++)
      for (j = 0; j < SOBOL_STP_NUM; j++)
        v[pl + k][rl + i][sl + j] =
            inv_cumnorm_fast(SOBOL_RAN_NUM[i + SOBOL_BRW_NUM * j]);
  }

  return NULL;
}

/* -------------------------------------------------------------------------------
 */
Err sobol_free(void) {

  if (!SobFree())
    return serror("Sobol free failed");

  if (SOBOL_RAN_NUM)
    free_dvector(SOBOL_RAN_NUM, 0, SOBOL_DIM);

  SOBOL_RAN_NUM = NULL;

  return NULL;
}

/* ------------------------------------------------------------------------------------------
 */

double **GetUnifDev(long p, long d)

{
  long i;
  double **res = dmatrix(0, p - 1, 0, d - 1);

  SobInit(d);

  for (i = 0; i < p; i++)
    SobNext(d, &(res[i][0]));

  SobFree();

  return (res);
}

void GetSobolMatrix(long p, long d, double **mtx)

{
  long i;

  SobInit(d);
  for (i = 0; i < p; i++)
    SobNext(d, &(mtx[i][0]));
  SobFree();
}
