/* ---------------------------------------------------------------------------------
   FILENAME:         srt_f_mc_genrand.cxx

   MAIN FUNCTIONS:   Err SrtMCRanSrc_start(...)
                     Err SrtMCRanSrc_correl(...)
                     Err SrtMCRanSrc_next(SrtMCRanSrc *mcrs)
                     Err SrtMCRanSrc_end(...)
   PRIVATE FUNCTION: Err srt_f_gen_random(...)


   PURPOSE:          The main core for the random number generation when using
                     Monte Carlo. The numbers finally output correspond to
   normal variables        , that should be translated into Brownian increments
   for each time step.
   ---------------------------------------------------------------------------------
 */

#include <math.h"
#include <stdio.h"

#include "srt_h_all.h"

#include "srt_h_mc_genrand.h"
#include "utallhdr.h"

static Err srt_f_extract_time_at_steps(SrtStpPtr top, double **time_at_steps,
                                       long *size);

/* -----------------------------------------------------------------------
   PRIVATE FUNCTION : NOT TO BE CALLED FROM THE OUTSIDE
   Function which calls relevent routines        , depending on generation
   method to fill a cube
                        rand[first_path..last_path][0..nbr-1][first_step..last_step]
   This function should be used only when the entire cube of random needs
   to be stored (non dynamic methods)
   Memory has already been allocated
   ----------------------------------------------------------------------- */

static Err srt_f_gen_random(double ***rand, int first_path, int last_path,
                            int first_step, int last_step, int nbr,
                            SrtMCSamType gen_method, long *seed,
                            SrtStpPtr top) {
  int i, w, n, segment, extra;
  int num_path;
  int num_step;
  double **matrix = NULL;
  double *time_at_step = NULL;
  long size;
  Err err = NULL;

  num_path = last_path - first_path + 1;
  num_step = last_step - first_step + 1;

  if (gen_method == BAL_SAMPLING) {
    /* BalSam returns a cube [path][brownian][step] */
    err = BalSampCube(rand, first_path, last_path, 0, nbr - 1, first_step,
                      last_step, seed);
    if (err)
      return err;

  } /* END if (method == BAL_SAMPLING) */

  else if (gen_method == ABS) {
    /* ABS returns a cube [path][brownian][step] */
    err = ABSCube(rand, first_path, last_path, 0, nbr - 1, first_step,
                  last_step, seed);
    if (err)
      return err;

  } /* END if (method == ABS) */

  else if (gen_method == RANDOM_GAUSS) {
    err = gauss_cube(rand, first_path, last_path, 0, nbr - 1, first_step,
                     last_step, seed);
    if (err)
      return err;

  } /* END if (method == RANDOM_GAUSS) */

  else if (gen_method == ANTITHETIC) {
    /* See if there is an even or odd number of paths */
    if (num_path % 2 != 0)
      extra = 1;
    else
      extra = 0; /** extra= 1 indicates an odd no. of paths **/

    segment = (num_path + extra) / 2;

    /* Generates the first half ( + 1 point if odd) of rand[i] */
    err = gauss_cube(rand, first_path, first_path + segment, 0, nbr - 1,
                     first_step, last_step, seed);
    if (err)
      return err;

    /* Makes the symetry (path wise) */
    for (w = first_path; w < first_path + segment - extra; w++) {
      for (i = 0; i < nbr; i++) {
        for (n = first_step; n <= last_step; n++) {
          rand[w + segment][i][n] = -rand[w][i][n];

        } /* END of loop on steps */
      }   /* END of loop on Brownian */
    }     /* END of loop on paths */

  } /* END if (method == ANTITHETIC) */

  else

      if (gen_method == BAL_ANTI) {
    /* See if there is an even or odd number of paths */
    if (num_path % 2 != 0)
      extra = 1;
    else
      extra = 0; /** extra= 1 indicates an odd no. of paths **/
    segment = (num_path + extra) / 2;

    /* Generates the first half ( + 1 point if odd) of rand[i] with Balance
     * Sampling */
    err = BalSampCube(rand, first_path, first_path + segment, 0, nbr - 1,
                      first_step, last_step, seed);
    if (err)
      return err;

    /* Makes the symetry (path wise) */
    for (w = first_path; w < first_path + segment - extra; w++) {
      for (i = 0; i < nbr; i++) {
        for (n = first_step; n <= last_step; n++) {
          rand[w + segment][i][n] = -rand[w][i][n];

        } /* END of loop on steps */
      }   /* END of loop on Brownian */
    }     /* END of loop on paths */

  } /* END of if (method = BAL_ANTI)*/

  else

      if (gen_method == SPECTRUNC) {
    /* Extracts from the SrtStpPtr the list of all time steps  */
    err = srt_f_extract_time_at_steps(top, &time_at_step, &size);
    if (err)
      return err;

    /* Generate the non correlated Brownian paths according to the covar matrix
     */
    err = SpecTruncCube(rand, time_at_step, 1.0e-04, first_path, last_path, 0,
                        nbr - 1, first_step, last_step, seed);
    if (err)
      return err;

    /* Free Memory */
    free_dvector(time_at_step, 0, size - 1);

  } /* END of if (gen_method == SPECTRUNC) */

  else
    return serror("Unknown method requested in gen_random");

  return NULL;

} /* END of static srt_f_gen_random(...)*/

/*---------------------------------------------------------------------------------------------------------*/
/*New function that does the same as the previous one but which does not take a
SrtStpPtr as an input but only the vector of times_at_steps  which is the only
thing really needed*/

Err srt_f_New_gen_random(double ***rand, int first_path, int last_path,
                         int first_step, int last_step, int nbr,
                         SrtMCSamType gen_method, long *seed,
                         double *times_at_steps) {
  int i, w, n, segment, extra;
  int num_path;
  int num_step;
  double **matrix = NULL;
  double *time_at_step = NULL;
  Err err = NULL;

  num_path = last_path - first_path + 1;
  num_step = last_step - first_step + 1;

  if (gen_method == BAL_SAMPLING) {
    /* BalSam returns a cube [path][brownian][step] */
    err = BalSampCube(rand, first_path, last_path, 0, nbr - 1, first_step,
                      last_step, seed);
    if (err)
      return err;

  } /* END if (method == BAL_SAMPLING) */

  else if (gen_method == ABS) {
    /* ABS returns a cube [path][brownian][step] */
    err = ABSCube(rand, first_path, last_path, 0, nbr - 1, first_step,
                  last_step, seed);
    if (err)
      return err;

  } /* END if (method == ABS) */

  else if (gen_method == RANDOM_GAUSS) {
    err = gauss_cube(rand, first_path, last_path, 0, nbr - 1, first_step,
                     last_step, seed);
    if (err)
      return err;

  } /* END if (method == RANDOM_GAUSS) */

  else if (gen_method == ANTITHETIC) {
    /* See if there is an even or odd number of paths */
    if (num_path % 2 != 0)
      extra = 1;
    else
      extra = 0; /** extra= 1 indicates an odd no. of paths **/

    segment = (num_path + extra) / 2;

    /* Generates the first half ( + 1 point if odd) of rand[i] */
    err = gauss_cube(rand, first_path, first_path + segment, 0, nbr - 1,
                     first_step, last_step, seed);
    if (err)
      return err;

    /* Makes the symetry (path wise) */
    for (w = first_path; w < first_path + segment - extra; w++) {
      for (i = 0; i < nbr; i++) {
        for (n = first_step; n <= last_step; n++) {
          rand[w + segment][i][n] = -rand[w][i][n];

        } /* END of loop on steps */
      }   /* END of loop on Brownian */
    }     /* END of loop on paths */

  } /* END if (method == ANTITHETIC) */

  else

      if (gen_method == BAL_ANTI) {
    /* See if there is an even or odd number of paths */
    if (num_path % 2 != 0)
      extra = 1;
    else
      extra = 0; /** extra= 1 indicates an odd no. of paths **/
    segment = (num_path + extra) / 2;

    /* Generates the first half ( + 1 point if odd) of rand[i] with Balance
     * Sampling */
    err = BalSampCube(rand, first_path, first_path + segment, 0, nbr - 1,
                      first_step, last_step, seed);
    if (err)
      return err;

    /* Makes the symetry (path wise) */
    for (w = first_path; w < first_path + segment - extra; w++) {
      for (i = 0; i < nbr; i++) {
        for (n = first_step; n <= last_step; n++) {
          rand[w + segment][i][n] = -rand[w][i][n];

        } /* END of loop on steps */
      }   /* END of loop on Brownian */
    }     /* END of loop on paths */

  } /* END of if (method = BAL_ANTI)*/

  else

      if (gen_method == SPECTRUNC) {

    /* Generate the non correlated Brownian paths according to the covar matrix
     */
    err = SpecTruncCube(rand, times_at_steps, 1.0e-04, first_path, last_path, 0,
                        nbr - 1, first_step, last_step, seed);
    if (err)
      return err;

  } /* END of if (gen_method == SPECTRUNC) */

  else
    return serror("Unknown method requested in gen_random");

  return NULL;

} /* END of srt_f_New_gen_random(...)*/

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
        SrtMCRanSrc

        USE THIS STRUCTURE AND IT'S FUNCTIONS TO GET YOUR RANDOM SAMPLES , FOR
   MONTE CARLO SIMULATION. RIGHT NOW IT GENERATES ALL SAMPLES AT ONCE        ,
   BUT IT COULD CONCEIVABLY BE READING NUMBERS FROM A FILE        , OR
   GENERATING SAMPLES ON THE FLY.
   ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Initialises mcrs according to inputs:
                - for DYNAMIC methods        , sets mcrs->pl = mcrs->ph = 0
                - makes memory allocation for mcrs->r[pl..ph][0..nbr-1][sl..sh]
                - initialises the seed (or equivalent)
                - for NON DYNAMIC generates the non correlated random numbers
   Please note that the strucuture of random numbers is the following:
                                                mcrs->r[path][brownian][step]
   ------------------------------------------------------------------------ */

Err SrtMCRanSrc_start(SrtMCRanSrc *mcrs, long pl, long ph, long sl, long sh,
                      long nbr, SrtMCSamType gen_method, long seed,
                      SrtStpPtr top) {
  Err err = NULL;

  /* Initialisation of the fields in the structure */
  mcrs->seed = seed;
  mcrs->sample_type = gen_method;
  mcrs->pl = pl;
  mcrs->ph = ph;
  mcrs->lastindex = pl;
  mcrs->sl = sl;
  mcrs->sh = sh;
  mcrs->nbr = nbr;

  /* For the moment        , the current path is one before the first one
   * (increment after in MC...next) */
  mcrs->cur_path = mcrs->pl - 1;

  /* Check the inputs: something ha to be generated (for Evolve functions to
   * work)  */
  if (mcrs->sl > mcrs->sh)
    mcrs->sl = mcrs->sh;

  /* Dynamic methods: just initialisation and memory allocation FOR ONE PATH is
   * done */
  if (gen_method == RANDOM_GAUSS_DYNAMIC || gen_method == ANTITHETIC_DYNAMIC ||
      gen_method == SOBOL || gen_method == UWAN) {
    if (gen_method == SOBOL) {
      /* Initialises Sobol for nbr Paths        , Brownian and all the steps */

      if (err = sobol_init(mcrs->pl, mcrs->ph, 0, nbr - 1, mcrs->sl, mcrs->sh))
        return err;
    }

    else if (gen_method == UWAN) {
      /* Initialises Uwan for nbr Paths        , Brownian and all the steps */

      if (err = uwan_init(mcrs->pl, mcrs->ph, 0, nbr - 1, mcrs->sl, mcrs->sh))
        return err;
    }

    /* In this case        , only one path will be used in the mcrs->r cube */
    mcrs->pl = 0;
    mcrs->ph = 0;
    mcrs->lastindex = 0;

    /* Allocate memory for one path only ( path number mcrs->pl = 0 ) in mcrs->r
     */
    mcrs->r = dcube(mcrs->pl, mcrs->ph, 0, mcrs->nbr - 1, mcrs->sl, mcrs->sh);
    if (!mcrs->r)
      return serror("Couldn't allocate space for random numbers");

    /* Sets all the values to 0 for safety reasons */
    memset(&(mcrs->r[0][0][0]), 0, nbr * (sh - sl + 1) * sizeof(double));

  } /* END if method = RANDOM_GAUSS_DYNAMIC        , SOBOL        , UWAN ,
       ANTI_DYNAMIC */

  else

  /* Non dynamic methods: full memory allocation and generation is done here*/
  {
    /* Allocate some initial memory for all the random numbers
     * [path][brow][step] */
    mcrs->r = dcube(mcrs->pl, mcrs->ph, 0, mcrs->nbr - 1, mcrs->sl, mcrs->sh);
    if (!mcrs->r)
      return serror("Couldn't allocate space for random numbers");

    /* Sets all the values to 0.0 for safety reasons */
    memset(&(mcrs->r[0][0][0]), 0,
           nbr * (sh - sl + 1) * (ph - ph + 1) * sizeof(double));

    /* Generates and stores the numbers in the large
     * mcrs->r[path][brownian][step] array */
    err = srt_f_gen_random(mcrs->r, mcrs->pl, mcrs->ph, mcrs->sl, mcrs->sh,
                           mcrs->nbr, gen_method, &mcrs->seed, top);
    if (err)
      return err;

  } /* END if method != DYNAMIC or SOBOL */

  return NULL;

} /* END SrtMCRanSrc_start(...) */

/*----------------------------------------------------------------------------------------------*/
/*New Function which is the same as the previous one only taking times_at_steps
instead as a SrtStpPtr as an input*/

Err SrtNewMCRanSrc_start(SrtMCRanSrc *mcrs, long pl, long ph, long sl, long sh,
                         long nbr, SrtMCSamType gen_method, long seed,
                         double *times_at_steps) {
  Err err = NULL;

  /* Initialisation of the fields in the structure */
  mcrs->seed = seed;
  mcrs->sample_type = gen_method;
  mcrs->pl = pl;
  mcrs->ph = ph;
  mcrs->lastindex = pl;
  mcrs->sl = sl;
  mcrs->sh = sh;
  mcrs->nbr = nbr;

  /* For the moment        , the current path is one before the first one
   * (increment after in MC...next) */
  mcrs->cur_path = mcrs->pl - 1;

  /* Check the inputs: something ha to be generated (for Evolve functions to
   * work)  */
  if (mcrs->sl > mcrs->sh)
    mcrs->sl = mcrs->sh;

  /* Dynamic methods: just initialisation and memory allocation FOR ONE PATH is
   * done */
  if (gen_method == RANDOM_GAUSS_DYNAMIC || gen_method == ANTITHETIC_DYNAMIC ||
      gen_method == SOBOL || gen_method == UWAN) {
    if (gen_method == SOBOL) {
      /* Initialises Sobol for nbr Paths        , Brownian and all the steps */

      if (err = sobol_init(mcrs->pl, mcrs->ph, 0, nbr - 1, mcrs->sl, mcrs->sh))
        return err;
    }

    else if (gen_method == UWAN) {
      /* Initialises Uwan for nbr Paths        , Brownian and all the steps */

      if (err = uwan_init(mcrs->pl, mcrs->ph, 0, nbr - 1, mcrs->sl, mcrs->sh))
        return err;
    }

    /* In this case        , only one path will be used in the mcrs->r cube */
    mcrs->pl = 0;
    mcrs->ph = 0;
    mcrs->lastindex = 0;

    /* Allocate memory for one path only ( path number mcrs->pl = 0 ) in mcrs->r
     */
    mcrs->r = dcube(mcrs->pl, mcrs->ph, 0, mcrs->nbr - 1, mcrs->sl, mcrs->sh);
    if (!mcrs->r)
      return serror("Couldn't allocate space for random numbers");

    /* Sets all the values to 0 for safety reasons */
    memset(&(mcrs->r[0][0][0]), 0, nbr * (sh - sl + 1) * sizeof(double));

  } /* END if method = RANDOM_GAUSS_DYNAMIC        , SOBOL        , UWAN ,
       ANTI_DYNAMIC */

  else

  /* Non dynamic methods: full memory allocation and generation is done here*/
  {
    /* Allocate some initial memory for all the random numbers
     * [path][brow][step] */
    mcrs->r = dcube(mcrs->pl, mcrs->ph, 0, mcrs->nbr - 1, mcrs->sl, mcrs->sh);
    if (!mcrs->r)
      return serror("Couldn't allocate space for random numbers");

    /* Sets all the values to 0.0 for safety reasons */
    memset(&(mcrs->r[0][0][0]), 0,
           nbr * (sh - sl + 1) * (ph - ph + 1) * sizeof(double));

    /* Generates and stores the numbers in the large
     * mcrs->r[path][brownian][step] array */
    err = srt_f_New_gen_random(mcrs->r, mcrs->pl, mcrs->ph, mcrs->sl, mcrs->sh,
                               mcrs->nbr, gen_method, &mcrs->seed,
                               times_at_steps);
    if (err)
      return err;

  } /* END if method != DYNAMIC or SOBOL */

  return NULL;

} /* END SrtNewMCRanSrc_start(...) */

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
   This function fills in the mcrs->r[mcrs->cur_path][0...nbr-1][sl...sh] matrix
         , according to the method required        , and the PATH NUMBER If
   numbers have been generated before (non DYNAMIC methods)        , this
   function makes sure that mcrs->r[i] points to the right path. Memory
   allocation has to be made before for the mcrs->r[pathindex][][] matrix ,and
   for the mcrs Once this function has been called        , mcrs->r[i] could be
   used for i = pl ... ph [paths].
   ------------------------------------------------------------------------ */

Err SrtMCRanSrc_init(SrtMCRanSrc *mcrs) {
  Err err = NULL;
  long i, j;
  long index;

  /* Check that there is at least one time step ; if not: nothing */
  if (mcrs->sl > mcrs->sh)
    return NULL;

  /* set the index to the first path */
  for (index = mcrs->pl; index < mcrs->ph; index++) {

    /* DYNAMIC methods first: no initial storage of the full cube */
    if (mcrs->sample_type == RANDOM_GAUSS_DYNAMIC) {

      /* Fills in the random numbers for this path */
      if (err = gauss_matrix(mcrs->r[index], 0, mcrs->nbr - 1, mcrs->sl,
                             mcrs->sh, &mcrs->seed))
        return err;

      /* Sets the index for the last generated path to 0 (not really used) */
      mcrs->lastindex = 0;

    } /* END if (sample_type == RANDOM_GAUSS_DYNAMIC) */

    else if (mcrs->sample_type == SOBOL) {

      /* Fills in the random numbers for this path */
      if (err = sobol_matrix(mcrs->r[index], 0, mcrs->nbr - 1, mcrs->sl,
                             mcrs->sh))
        return err;

      /* Sets the index for the last generated path to 0 (not really used) */
      mcrs->lastindex = 0;

    } /* END if (sample_type == SOBOL) */

    else if (mcrs->sample_type == UWAN) {

      /* Fills in the random numbers for this path */
      if (err = uwan_matrix(mcrs->r[index], 0, mcrs->nbr - 1, mcrs->sl,
                            mcrs->sh, &mcrs->seed))
        return err;

      /* Sets the index for the last generated path to 0 (not really used) */
      mcrs->lastindex = 0;

    } /* END if (sample_type == UWAN) */

    else if (mcrs->sample_type == ANTITHETIC_DYNAMIC) {
      if ((index - mcrs->pl) % 2 == 0)
      /* If last generated path number is even        , generate a new path
         (Brow * steps) */
      {
        if (err = gauss_matrix(mcrs->r[index], 0, mcrs->nbr - 1, mcrs->sl,
                               mcrs->sh, &mcrs->seed))
          return err;
      } else /* (mcrs->lastindex - mcrs->pl) % 2 == 1 */
      /* If path number is odd        , use the previous one and make it
         symetrical*/
      {
        for (i = 0; i < mcrs->nbr; i++) {
          for (j = mcrs->sl; j <= mcrs->sh; j++)
            mcrs->r[index][i][j] *= -1;
        }
      }

    } /* END if ANTITHETIC_DYNAMIC */

    else
    /* Non DYNAMIC methods: all the numbers were stored by SrtMCRanSrc_init
       have been already generate by SrtMCRanSrc_start*/
    {
      /* Increment the current path index (for storage of mcrs->r[...])  */
    }
  }
  return err;

} /* END of SrtMCRanSrc_init(...) */

/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Once the non correlated numbers have been generated in
   mcrs->r[mcrs->cur_path]        ,
   this is used to correlate them.
   The correlated numbers overwrite the ones in
   mcrs->r[mcrs->cur_path][...][...] This function uses SrtStrPtr to associate
   the right correlation matrix to the right step        , assuming that one
   step correspond to one move in the list and one indexation in the random
   numbers matrix The correlation is made for:
                - all the Brownian (mcrs->nbr)        ,
        - all the steps (from mcrs->sl to mcrs->sh)
                - one path only (the one that corresponds to the rndm matrix)
   ------------------------------------------------------------------------ */

Err SrtMCRanSrc_correl(SrtMCRanSrc *mcrs, SrtStpPtr stp) {
  Err err = NULL;
  int rndm_size = (mcrs->sh - mcrs->sl + 1) * mcrs->nbr;
  int step_index = mcrs->sl;

  /* Check that there is at least one time step ; if not: nothing */
  if (rndm_size <= 0)
    return NULL;

  /* If only one Brownian motion        , no work has to be done */
  if (mcrs->nbr == 1)
    return NULL;
  else
      /* If the need_to_correl flag is set to SRT_NO        , no work has to be
       * done
       */
      if (mcrs->need_to_correl == SRT_NO)
    return NULL;
  else {
    /* Test if stp is well defined */
    if (!stp)
      return serror("SrtStpPtr not defined in SrtMCRanSrc_correl");

    while (stp->next) {
      if (step_index > mcrs->sh)
        return serror(
            "SrtStpPtr and mcrs do not correspond in SrtMCRanSrc_correl");
      if (!stp->coeff)
        return serror("Null coefficients in SrtMCRanSrc_correl");
      if (1) {
        /* The original function call was to "correl_random" with
                the same argumrnts at the eigenvalue equivalent call below.*/
        err = correl_random_eigen(mcrs->r[mcrs->cur_path], step_index,
                                  step_index, mcrs->nbr, stp->coeff);
      } else {
        err = correl_random(mcrs->r[mcrs->cur_path], step_index, step_index,
                            mcrs->nbr, stp->coeff);
      }
      if (err)
        return err;

      /* Increment the step index by one        , and move to next one in the
       * list*/
      step_index++;
      stp = stp->next;

    } /* END of loop on the SrtStpPtr */

  } /* END if ( nbr != 1 ) */

  return err;

} /* END of SrtMCRanSrc_correl(...) */

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
   This function fills in the mcrs->r[mcrs->cur_path][0...nbr-1][sl...sh] matrix
         , according to the method required        , and the PATH NUMBER If
   numbers have been generated before (non DYNAMIC methods)        , this
   function makes sure that mcrs->r[mcrs->cur_path] points to the right path.
   Memory allocation has to be made before for the mcrs->r[pathindex][][] matrix
   ,and for the mcrs Once this function has been called        , mcrs->lastindex
   refers to the index of the path that has just been generated.
   ------------------------------------------------------------------------ */

Err SrtMCRanSrc_next(SrtMCRanSrc *mcrs) {
  Err err = NULL;
  long i, j;

  /* Check that there is at least one time step ; if not: nothing */
  if (mcrs->sl > mcrs->sh)
    return NULL;

  /* DYNAMIC methods first: no initial storage of the full cube */
  if (mcrs->sample_type == RANDOM_GAUSS_DYNAMIC) {
    /* The current path index (for storage) is always the same : the first one
     */
    mcrs->cur_path = mcrs->pl;

    /* Fills in the random numbers for this path */
    if (err = gauss_matrix(mcrs->r[mcrs->pl], 0, mcrs->nbr - 1, mcrs->sl,
                           mcrs->sh, &mcrs->seed))
      return err;

    /* Sets the index for the last generated path to 0 (not really used) */
    mcrs->lastindex = 0;

  } /* END if (sample_type == RANDOM_GAUSS_DYNAMIC) */

  else if (mcrs->sample_type == SOBOL) {
    /* The current path index (for storage) is always the same : the first one
     */
    mcrs->cur_path = mcrs->pl;

    /* Fills in the random numbers for this path */
    if (err = sobol_matrix(mcrs->r[mcrs->pl], 0, mcrs->nbr - 1, mcrs->sl,
                           mcrs->sh))
      return err;

    /* Sets the index for the last generated path to 0 (not really used) */
    mcrs->lastindex = 0;

  } /* END if (sample_type == SOBOL) */

  else if (mcrs->sample_type == UWAN) {
    /* The current path index (for storage) is always the same : the first one
     */
    mcrs->cur_path = mcrs->pl;

    /* Fills in the random numbers for this path */
    if (err = uwan_matrix(mcrs->r[mcrs->pl], 0, mcrs->nbr - 1, mcrs->sl,
                          mcrs->sh, &mcrs->seed))
      return err;

    /* Sets the index for the last generated path to 0 (not really used) */
    mcrs->lastindex = 0;

  } /* END if (sample_type == UWAN) */

  else if (mcrs->sample_type == ANTITHETIC_DYNAMIC) {
    /* The current path index (for storage) is always the same : the first one
     */
    mcrs->cur_path = mcrs->pl;

    if ((mcrs->lastindex - mcrs->pl) % 2 == 0)
    /* If last generated path number is even        , generate a new path (Brow
       * steps) */
    {
      if (err = gauss_matrix(mcrs->r[mcrs->pl], 0, mcrs->nbr - 1, mcrs->sl,
                             mcrs->sh, &mcrs->seed))
        return err;
    } else /* (mcrs->lastindex - mcrs->pl) % 2 == 1 */
    /* If path number is odd        , use the previous one and make it
       symetrical*/
    {
      for (i = 0; i < mcrs->nbr; i++) {
        for (j = mcrs->sl; j <= mcrs->sh; j++)
          mcrs->r[mcrs->pl][i][j] *= -1;
      }
    }

    /* Increment the index of the last generated path (for tracking purposes) */
    mcrs->lastindex++;

  } /* END if ANTITHETIC_DYNAMIC */

  else

  /* Non DYNAMIC methods: all the numbers were stored by SrtMCRanSrc_init */
  {
    /* Increment the current path index (for storage of mcrs->r[...])  */
    mcrs->cur_path++;
  }

  return err;

} /* END of SrtMCRanSrc_next(...) */

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Undo whatever SrtMCRanSrc_init() did:
                - free memory
                - ....
  -------------------------------------------------------------------------- */

Err SrtMCRanSrc_end(SrtMCRanSrc *mcrs) {
  Err err = NULL;

  /* Check that there was at least one time step ; if not: nothing */
  if (mcrs->sl > mcrs->sh)
    return NULL;

  /* Free memory for SOBOL or UWAN */
  if (mcrs->sample_type == SOBOL) {
    err = sobol_free();
    if (err)
      return err;
  } else if (mcrs->sample_type == UWAN) {
    err = uwan_free();
    if (err)
      return err;
  }

  /* Free the full cube of random points (for dynamic: only one path...)  */
  if (mcrs->r)
    free_dcube(mcrs->r, mcrs->pl, mcrs->ph, 0, mcrs->nbr - 1, mcrs->sl,
               mcrs->sh);

  /* Sets the pointer to NULL for safety reasons (to prevent double free...) */
  mcrs->r = NULL;

  return err;

} /* END of SrtMCRanSrc_end(...) */

/* --------------------------------------------------------------------------------
 */

/* --------------------------------------------------------------------------------
               UTILITY FUNCTIONS FOR SPECTRAL TRUNCATION
   --------------------------------------------------------------------------------
 */

/* Build from a SrtStpPtr a vector of time at steps (allocation is done inside)
 */

static Err srt_f_extract_time_at_steps(SrtStpPtr top, double **time_at_steps,
                                       long *size) {
  Err err = NULL;
  long i;
  SrtStpPtr stp;

  /* Compute the number of time steps (hopefully        , time is 0 at first
   * step) */
  i = 0;
  stp = top;
  while (stp) {
    i++;
    stp = stp->next;
  }
  *size = i;

  /* Memory allocation */
  (*time_at_steps) = dvector(0, *size);

  /* Sets all the times (from today to step date) into a single vector for all
   * the step*/
  stp = top;
  for (i = 0; i < *size; i++) {
    (*time_at_steps)[i] = stp->time;
    stp = stp->next;
  }

  /* Return a success string */
  return NULL;

} /* END srt_f_extract_time_at_steps(...) */

/* --------------------------------------------------------------------------------
 */
