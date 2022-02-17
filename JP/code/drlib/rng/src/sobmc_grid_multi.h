/***************************************************************************
 
  File       : sobmc_grid_multi.h
  Written by : Viral Acharya
  Date       : 6th August 1997
 
  This is the header file accompaying the library implementation of random
  number generator for Monte Carlo
  simulations using Sobol sequences and one-dimensional grids, along with
  Brownian bridges. The full
  implementation is contained in sobmc_grid_multi.c
 
  Accompanying files : Normal.c InvertCumNormal.c sobol.c randomgen.c
                       Normal.h InvertCumNormal.h sobol.h randomgen.h
 
  M.Huq : Took out references to A-lib
  M.Huq : Thu Aug 16 14:20:40 EDT 2001
          Added pointer to Sequence to pass in low-discrepancy sequence as desired
 
****************************************************************************/

#ifndef _SC_SOBMC_GRID_MULTI_H
#define _SC_SOBMC_GRID_MULTI_H

#include "edginc/Sequence.h"
#include "edginc/IRNG.h"
#include <stdlib.h>

CORE_BEGIN_NAMESPACE

#define FREE(ptr)         free(ptr)

#define FREE_ARRAY(ptr)         free(ptr)

#define NEW_ARRAY(t,n)          (t *) malloc(sizeof(t)*(n))


/* Standard definitions   */

#ifndef TRUE
#define TRUE true
#endif

#ifndef FALSE
#define FALSE false
#endif

#ifndef BBSUCCESS
#define BBSUCCESS 0
#endif

#ifndef BBFAILURE
#define BBFAILURE -1
#endif

#ifndef IS
#define IS ==
#endif

#ifndef IS_NOT
#define IS_NOT !=
#endif


/**************************************************************************
 
  The following are variables defined for the implementation.
 
  *************************************************************************/

#define SGM_SOB_MAX       391   /* Maximum number of sobol dimensions is restricted to 50 */
#define SGM_POW_MAX       15   /* Used internally for computing powers of 2 */

/**************************************************************************
 
  The following are the possible types of dimension (dimensionType argument) :-
 
  Sobol dimension can take on any value from 1..SOB_MAX (defined above).
  A Sobol dimension cannot be repeated across dimensions.
 
  **************************************************************************/

#define SGM_PSEUDO_RANDOM -1L /* A regular pseudo-random deviate generated using ran2/gasdev2 */

/* A one-dimensional grid. Can be used on exactly
    one dimension among all
    dimensions (across all factors). Strongly
    recommended
    that this be used
    on the most important dimension
    - usually last time point in single factor
    implementation.
    - usually last time point in most important
    dimension in a multi-factor
    implementation.
    */

#define SGM_GRID           0L
/********************************************************************
 
  The following are the types of random number fetch (fetchType) :-
 
  ***********************************************************************/

#define SGM_FETCH_ONE_FACTOR          0L    /* To fetch all dimensions of one
factor */
#define SGM_FETCH_ONE_DIM_ONE_FACTOR  1L    /* To fetch deviates along a
particular
dimension of a
factor */
#define SGM_FETCH_ONE_DIM_ALL_FACTORS 2L    /* To fetch a cross-section of path
along one dimension
for all factors. */

/************************************************************************
 
  The following is used to start Sobol sequences at their 32nd point,
  as recommended by Arnon Levy.
 
  ***********************************************************************/

#define SGM_SOBOL_START  32
#define SGM_DEFAULT_SEED -3L    /* In case negative seed is not passed in
initialization call */

/* Initialization routine for monte carlo simulation.
   This does all memory allocation and must be associated with a closing
   free call to get rid of
   all memory it has allocated. It also updates data structures for
   building bridges.
   */

int
SGM_Initialize_Sobol_Bridge(long numFactorsPerPath,    /* Number of factors
                                                                 in a simulation */
                            long *numDimsPerFactor,    /* Number of dimensions
                                                          in each factor */
                            long numPathsPerRun,       /* Number of paths
                                                          (simulations) per
                                                          run*/
                            long numRuns,              /* Number of runs */
                            long numPathsPerFetch,     /* Number of paths
                                                          accessed
                                                          in one call
                                                          to fetch_Sobol_Bridge.
                                                          e.g. In our IR Monte
                                                          Carlo simulation,
                                                          we might want to
                                                          access 1000 points
                                                          per
                                                          fetch due to memory
                                                          reasons when the
                                                          total numer of paths
                                                          per
                                                          run is actually 5000.
                                                          */
                            long *dimensionType[],     /* Type of each dimension
                                                          in each factor -
                                                          could be
                                                          SGM_PSEUDO_RANDOM,
                                                          SGM_GRID or
                                                          1..SGM_GRID */
                            INormalSuperCubeRNGSP    rng,                /* Seed for random
                                                          number
                                                          generator */
                            long sampleTailsFlag,      /* TRUE/FALSE to sample
                                                          tails
                                                          or not in each
                                                          dimension */
                            /* The following are the three arguments to
                               Julia/Arnon's routine
                               for tail sampling. */
                            long numSampDistPts,       /* Number of points for
                                                          sampling distribution
                                                          */
                            double sampIntvlBound,     /* Determines sampling
                                                          distr interval
                                                          (-sampIntvlBound,
                                                          +sampIntvlBound)
                                                          in terms of standard
                                                          deviation */
                            double maxProbRatio,       /* Max ratio between
                                                          path
                                                          probabilities. If 1.,
                                                          the samples
                                                          are normal */
                            double *ptr2GaussianGrid,   /* Grid specification.
                                                         If NULL internal grid used.
                                                         User must specify grid in
                                                         non-uniform format. */
                            SequenceSP ldSeqPointer)
    ;   /* Pointer to Sequence class.
  This is a pointer to a low-discrpenacy
  sequence. This is used in the supports
  for the brownian bridge */

/* The following function should be called at the beginning of each fetch in
   Monte Carlo simulation.
   It prepares the context of monte carlo for the fetch, generates bridge
   end-points and performs other initializations.
   */

int
SGM_Prepare_Sobol_Bridge(void) ;

/* The following is used to fetch random gaussian deviates from this generator.
   It can be used to access numPathsPerFetch points along a given dimension
   in a given factor OR
   It can be used to access a cross-section of a path along a given dimension
   (timepoint) for all factors.
   */

int
SGM_Fetch_Sobol_Bridge(long typeFetch,         /* SGM_FETCH_ONE_FACTOR,
                                                      SGM_FETCH_ONE_DIM_ONE_FACTOR
                                                      or
                                                      SGM_FETCH_ONE_DIM_ALL_FACTORS
                                                      */
                       long factor,            /* factor to be fetched */
                       long dimension,         /* dimension to be fetched -
                                              ignored
                                              if SGM_FETCH_ONE_FACTOR */
                       long pathIndex,         /* path being fetched.
                                              This is used only for
                                              SGM_FETCH_ONE_DIM_ALL_FACTORS
                                              and must be between 0 and
                                              numPathsPerFetch-1. */
                       // unused         INormalSuperCubeRNGSP    rng,             /* seed for random number  generator */
                       double *gaussians,      /* Set of gaussian deviates --
                                              must be allocated in memory
                                              in calling routine */
                       double *probabilities); /* The probabilities
corresponding
to each point returned in
"gaussians" */

/* Final call to Sobol bridge.
   Cleans up all the static memory allocated in initialize_Sobol_Bridge(.)
   Must be called at the end of multi-factor Monte Carlo simulation to avoid
   any memory leakages.
   */

int
SGM_Free_Sobol_Bridge(void) ;

CORE_END_NAMESPACE

#endif // _SC_SOBMC_GRID_MULTI_H
