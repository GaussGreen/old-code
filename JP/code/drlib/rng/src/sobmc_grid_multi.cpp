
/*****************************************************************************
 
  File       : sobmc_grid_multi.c
  Written by : Viral Acharya
  Date       : 29th July 1997
  Modified by: Julia Chislenko, September 1997, January 1998
 
  This file contains the library implementation of random number generator
  for Monte Carlo
  simulations using Sobol sequences and one-dimensional grids, along with
  Brownian bridges.
 
  Accompanying files : Normal.c InvertCumNormal.c sobol.c randomgen.c
                       Normal.h InvertCumNormal.h sobol.h randomgen.h
 
 
  M.Huq : Added feature that allows a grid to be passed into the Brownian Bridge.
          Had to change interface to enable this.
  M.Huq : Aug 16th, 2001
          Update of access to low-discrepancy sequence for the bridge.
          Pass in a pointer to a Sequence that has been appropriately
   initialized. The Sequence can be a SobolSequence or can be a
   SuperCube sequence for instance that scrambles Sobol.
 
 
***************************************************************************/


#define  xSGM_DEBUG
#define  xSGM_DEBUG_SOB_NUM

#include "edginc/Sequence.h"
#include "SCErrorHandlingInterface.h"
#include "edginc/gasdev2.h"
#include "InvertCumNormal.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

CORE_BEGIN_NAMESPACE

// Needed for GtoErrMsg replacement.
// note that this file must be compiled as a c++ file
char logBuffer[8192];



/*****************************************************************************
 
  The following header files are used for function prototypes from the
  Accompanying files noted above.
 
  *****************************************************************************/

#include "Normal.h"          /* Cumulative normal density - numerical
approximation */
#include "InvertCumNormal.h" /* Invert Cumulative normal density -
numerical approximation */
#include "sobol.h"

#include "edginc/gasdev2.h" // INormalSuperCubeRNGSP

/****************************************************************************
 
  The following header files contains constants and prototypes specific to the
  multi-factor Sobol and Brownian Bridge implementation.
 
  *************************************************************************/

#include "sobmc_grid_multi.h" /* contains ALIB cerror.h and macros.h*/


/****************************************************************************
 
  The following structures are used to implement the multi-factor
  random-number generator, the
  one-dimensional grid and the brownian bridge.
 
 ****************************************************************************/

/* Specifications of the grid used for sampling uniform distribution */

typedef struct Grid
{
    double base ;               /* The grid is constructed from base .. top
           in [0,1] */
    double top  ;
    double grid_size ;          /* This is the width of grid in each run */
    double step_size ;          /* This is the size of step within a grid.
           A grid is spanned using steps of size
           "step_size" as we progress through runs. */
    double * gaussianGridPtr;   /* Pointer to grid. If NULL uses above parameters */

}
Grid ;

/* All specifications of a given monte carlo simulation */

typedef struct sobmcGridContext
{
    long numFactorsPerPath ;    /* Number of factor (e.g currencies) in a
           simulation */
    long *numDimsPerFactor ;    /* Number of dimensions in each factor */
    long *numBridgesPerFactor ; /* Number of sobol dimensions per factor --
           To be calculated */
    long maxDimsPerFactor  ;    /* Maximum number of dimensions per factor
           across all factors */
    long maxSobDim ;            /* Maximum sobol dimension used across all
           factors.
           Is checked to be <= SGM_SOB_MAX */
    long numPathsPerRun    ;    /* Number of paths(simulations) per run */
    long numRuns           ;    /* Number of runs */
    long numPathsPerFetch  ;    /* Number of paths accessed in one call to
           fetch_Sobol_Bridge.
           e.g. In our IR Monte Carlo simulation, we
           might want to
           access 1000 points per fetch due to memory
           reasons when the total numer of paths per
           run is actually 5000. */
    long **dimensionType   ;    /* Type of each dimension in each factor -
           could be SGM_PSEUDO_RANDOM, SGM_GRID or
           1..SGM_SOB_MAX.
           This is not used after initialization and
           hence is set only as
           a pointer to passed dimensionType array. */
    long default_dimTypes  ;    /* Flag used internally in getting default
           dimension Types */
    long *curBridge        ;    /* Used internally during generation of
           deviates */
    long curRunNum         ;    /* Current run number */
    long curPathNum        ;    /* Path number in current run number */
    long current_dimension ;      /* Used only in
             SGM_FETCH_ONE_DIM_ALL_FACTORS */


    Grid uniformGrid       ;    /* Grid specifications used for generating
           perfectly uniform points */
    double **nextsobs      ;    /* Next set of sobol sequences.
           nextsobs[i][j] gives sobol point for ith
           simulation in current run
           for jth sobol dimension. */
    INormalSuperCubeRNGSP rng; // it seems that we use only normal generators
    /*  long seed              ; */   /* Seed for pseudo-random number generation.
                             Must be < 0 in initialization call.
                             Or it is set internally. */
    SequenceSP LDSequencePtr;    /* Pointer to LDS for BB supports */
}
sobmcGridContext ;

typedef struct bridgeEnd
{
    double *deviates ;          /* The deviates generated for the bridgeEnd for
           a given run */
    long   dimensionType  ;     /* The Sobol dimension used for this particular
           deviate */
    long   t       ;            /* Timepoint t corresponding to this bridgeEnd */
    long   curt    ;            /* A temporary used while computing order table
           and while generating deviates. */
}
bridgeEnd ;

/* For computing the order of building brownian bridges */

typedef struct Order
{
    bridgeEnd *left  ;          /* The  left point of the bridge */
    bridgeEnd *right ;          /* The right point of the bridge */
    bridgeEnd *mid   ;          /* The   mid point of the bridge - To be
           computed. */
}
Order ;

/****************************************************************************
 
  The following static structures are used all through the program and
  across calls.
  This is the main source of any memory overhead from using this library
  in simulation software.
  The overhead however is limited by
  SGM_SOB_MAX * number of paths per run * sizeof(double).
 
  ***********************************************************************/

static sobmcGridContext context={0}; /* All about current simulation */
static bridgeEnd **bridges=NULL;  /* A two-dimensional array of current bridge
         end-points.
         bridges[i][j] corresponds to jth end
         point for ith factor. */
static Order  **orderTable=NULL ; /* Stores the order for computing the
         brownian bridges.
         orderTable[i][j] contains the order for
         building jth brownian
         bridge for factor i.
         Each factor has number of bridges
         (including terminal point)
         equal to numSobDims. However the order is
         computed only for
         numSobDims-1 bridges (excluding the first
         bridge from 0 to
         last sobol dimension in that factor). */
static int start[SGM_SOB_MAX] ;       /* Needed to initialize Arnon's sobol
      sequence generator */
// Variables needed for Sobol sequence generator
// Initialization done in SGM_Initialize_Sobol_Bridge
// static int **Sobol_gen;
// static int  *Sobol_xn;
// static int  *Sobol_n;
// static int  Sobol_nstream;



/******************************************************************************
 
  The following are constants declared static to avoid redundant computations
  on repeated function calls.
 
  ****************************************************************************/

static double ln2inv          ; /* 1.0 / log(2.0) */
static long   pow2[SGM_POW_MAX+1] ; /* Powers of 2 needed for building
           brownian bridges */

/****************************************************************************
 
  The main functions start here.
 
 ***************************************************************************/

/* Computes powers of an integer upto a given limit.
   Limited by max_powers.
   */

int
compute_powers(long base,
               long *powers,
               long max_powers,
               long upper_limit)
{
    long i ;

    powers[0] = 1;

    for (i=1; (i < max_powers) && (powers[i-1] < upper_limit); i++)
        powers[i] = powers[i-1] * base ;

    return BBSUCCESS ;
}

/* The following is used for initialization of is_Computed array to FALSE. */

int
init_int_array(int  *array,
               long num_points,
               int  value)
{
    long i ;

    for (i=0; i < num_points; i++)
        array[i] = value ;

    return BBSUCCESS ;

}

/* Gets next set of "paths_per_fetch" sequences for sobol dimensions used
   for monte carlo simulation */

int
get_next_sobs(double **next_sobs, long paths_per_fetch, long max_sob)
{
    long i ;
    long j;

    if(max_sob > 0 ) {
        /* Assumes that Arnon Levy's Sobol sequence generator has been initialized */

        for (i=0; i<paths_per_fetch; i++)
        {
            /* Get next sequence of sobol points for all dimensions */

            //SobolSequenceAL(0,max_sob-1,start,next_sobs[i],0) ;
            //     sobvect(max_sob,
            //      next_sobs[i],
            //      Sobol_gen,
            //      Sobol_xn,
            //      Sobol_n,
            //      &Sobol_nstream );
            // Copy over vector
            context.LDSequencePtr->populateVector();
            double *RNvect = context.LDSequencePtr->getVector();
            for(j=0; j<max_sob; j++) {
                next_sobs[i][j]=RNvect[j];
            }//j


#ifdef SGM_DEBUG_SOB_NUM
            for(j=0; j<max_sob; j++)
            {
                printf("%14.12f\t", next_sobs[i][j]);
            }
            printf("\n");
#endif

        }
    }
    return BBSUCCESS ;

}

/* Returns whether a dimensiontype is Sobol dimension or not */

int
is_Sobol(long dimtype)
{
    if ((dimtype >= 1) && (dimtype <= SGM_SOB_MAX))
        return TRUE ;
    else
        return FALSE ;
}

/* Initialize the bridges.
   Allocate memory for bridges. Bridges are updated to reflect right bridge
   end-points (along time
   dimension) and to also reflect the Sobol dimension used to compute them.
   */

int
initialize_Bridges(void)
{
    int  status ;    /* All temporaries */
    long factor ;
    long sobdim ;
    /*  long run ;
     */
    long timepoint=0 ; /* Used to obtain the timepoint for each sobol dimension */

    status = BBSUCCESS ;

    if ((bridges = NEW_ARRAY(bridgeEnd*, context.numFactorsPerPath))
            /*(bridgeEnd **) malloc (context.numFactorsPerPath * sizeof(bridgeEnd *)))*/
            IS (bridgeEnd **) NULL) {

        SCLogMessage("Initialization of Bridges : Memory allocation failed for bridges",FAILURE) ;
        status = BBFAILURE ;
        goto done ;
    }

    for (factor=0; factor < context.numFactorsPerPath; factor++) {

#ifdef SGM_DEBUG
        printf("Allocating memory to bridges[%ld]\n",factor) ;
#endif

        if ((bridges[factor] = NEW_ARRAY(bridgeEnd,
                                         context.numBridgesPerFactor[factor]+1)
                               /* (bridgeEnd *)
                                malloc ((context.numBridgesPerFactor[factor]+1) * sizeof(bridgeEnd))*/
            ) IS  NULL) {

            sprintf(logBuffer, "Initialization of Bridges : Memory allocation failed "
                    "for bridges[%ld]\n",factor) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        for (sobdim=0; sobdim <= context.numBridgesPerFactor[factor]; sobdim++) {

#ifdef SGM_DEBUG
            printf("Allocating memory to bridges[%ld][%ld].deviates\n",factor,
                   sobdim) ;
#endif

            if ((bridges[factor][sobdim].deviates =
                        NEW_ARRAY(double,context.numPathsPerFetch)
                        /* (double *)
                           malloc (context.numPathsPerFetch * sizeof(double))*/
                ) IS (double *) NULL) {

                sprintf(logBuffer,"Initialization of Bridges : Memory allocation failed "
                        "for bridges[%ld][%ld]\n",
                        factor, sobdim) ;
                SCLogMessage(logBuffer,FAILURE);

                status = BBFAILURE ;
                goto done ;
            }

            if (sobdim IS 0) {

                /* bridges[factor][0] is the time point 0 -- always to be
                   initialized to 0.0  */

                bridges[factor][sobdim].t    = 0 ;
                bridges[factor][sobdim].curt = 0 ;
                bridges[factor][sobdim].dimensionType = SGM_PSEUDO_RANDOM ;
                /* Doesnt matter - Not used */

                /* Initialize for bridges now */

                timepoint = 1 ;
            }
            else /* Look for a bridge along the dimensions of the factor --
                                skip pseudorandom points */ 
            {

                while (context.dimensionType[factor][timepoint-1] IS SGM_PSEUDO_RANDOM)
                    timepoint++ ;

                /* Identified the bridgeEnd. Update bridge information in bridges
                   data structure */

                bridges[factor][sobdim].t    = timepoint ;
                bridges[factor][sobdim].curt = 0 ;
                bridges[factor][sobdim].dimensionType =
                    context.dimensionType[factor][timepoint-1] ;

                timepoint++ ;
            }

#ifdef SGM_DEBUG
            printf("Initialized bridges[%ld][%ld] t = %ld DimensionType = %ld\n",
                   factor,sobdim,bridges[factor][sobdim].t,
                   bridges[factor][sobdim].dimensionType) ;
#endif

        }
    }

done :

    return status ;

}

/* Frees the memory allocated to bridges */

void
free_Bridges(void)
{
    long factor, sobdim/*, run */;

    if (bridges IS_NOT (bridgeEnd **) NULL) {

        for (factor=0; factor < context.numFactorsPerPath; factor++) {

            if (bridges[factor] IS_NOT (bridgeEnd *) NULL) {

                for (sobdim=0;sobdim<= context.numBridgesPerFactor[factor];sobdim++) {

                    if (bridges[factor][sobdim].deviates IS_NOT (double *) NULL)
                        FREE(bridges[factor][sobdim].deviates) ;

                }
                FREE(bridges[factor]) ;
            }
        }
        FREE(bridges) ;
        bridges = NULL;
    }
}

/* Prints Bridges. Used only for SGM_DEBUGging purposes */

void
print_Bridges(void)
{
    long i, j ;

    for (i=0; i<context.numFactorsPerPath; i++) {
        printf("Factor %ld : \n",i) ;
        for (j=0; j<=context.numBridgesPerFactor[i]; j++)
            printf("Bridge %ld : t = %ld DimensionType = %ld\n",
                   j, bridges[i][j].t, bridges[i][j].dimensionType) ;
    }


}

/* Allocate memory for order_table */

int
initialize_Order_Table(void)
{
    int  status ;
    long factor ;

    status = BBSUCCESS ;

    if ((orderTable = NEW_ARRAY(Order*, context.numFactorsPerPath)
                      /* (Order **) malloc (context.numFactorsPerPath * sizeof(Order *)) */
        )
            IS (Order **) NULL) {

        sprintf(logBuffer,"Initialization of Order Table : Memory allocation failed "
                "for orderTable\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    for (factor=0; factor < context.numFactorsPerPath; factor++) {

        if (context.numBridgesPerFactor[factor] > 1) {

            if ((orderTable[factor] =
                        NEW_ARRAY(Order, context.numBridgesPerFactor[factor]-1)
                        /* (Order *) malloc ((context.numBridgesPerFactor[factor]-1) *
                                sizeof(Order)) */
                ) IS (Order *) NULL) {

                sprintf(logBuffer,"Initialization of Order Table : Memory allocation "
                        "failed for orderTable[%ld]\n",
                        factor) ;
                SCLogMessage(logBuffer,FAILURE);
                status = BBFAILURE ;
                goto done ;
            }
        } else
            orderTable[factor] = (Order *) NULL ;
    }

done :

    return status ;
}

/* Frees the memory allocated to order_table */

void
free_Order_Table(void)
{
    long factor ;

    if (orderTable IS_NOT (Order **) NULL) {

        for (factor=0; factor < context.numFactorsPerPath; factor++)
            if (orderTable[factor] IS_NOT (Order *) NULL)
                FREE(orderTable[factor]) ;

        FREE(orderTable);
        orderTable = NULL;
    }
}

/* The following sets the order in orderTable for one point of a brownian
   bridge. We cheat here and use "curt" field in bridges to check whether
   we are visiting a bridgeEnd again.
   */

int
aux_set_Order(long      left,
              long      right,
              long      index,
              Order     *order_table, /* Should be order_table of current
                            factor */
              bridgeEnd *bridges)     /* Should be bridges of current factor */
{
    long mid ;

    mid = (long) ceil((right+left)*0.5) ;

    /* Check if mid has already been computed.
       If not computed, set the order for computing mid from left and right.
       Else do nothing.
       */

    if (bridges[mid].curt IS 0) /* mid not yet computed */
    {

        order_table[index].left  = &bridges[left]  ;
        order_table[index].right = &bridges[right] ;
        order_table[index].mid   = &bridges[mid]   ;
        bridges[mid].curt        = 1               ; /* Note that mid has been
                      computed */

        return BBSUCCESS ;

    }
    else

        return BBFAILURE ;

}

/* The following sets the order correct for a given bridge specification.
   It generates the order for building brownian bridges specified by bridges.
   Note that the bridge is formed by halving the bridge-length at each step.
   In an ordinary call to this routine,
 
   start_point = 0
   end_point   = numBridgesPerFactor
 
   The generated order of building the bridge is built into order_table,
   starting at start_index entry.
   In an ordinary call,
 
   start_index = 0
 
   */ /* This algorithm is wrong - it creates a bridge order that is
       * inconsistent with ascent of Sobol dimentions. Use instead another
       * function which is ordering bridges according to min Sobol dims
       * - JC, October 1997 */

int
set_Order_old(long             start_point,
              long             end_point,
              long             start_index,
              Order            *order_table, /* Should be order_table of current
                                  factor */
              bridgeEnd        *bridges)     /* Should be bridges of current
      factor */
{
    long    num_points ;                          /* Number of bridges */
    long    log2_num_points ;                     /* Temporary -- used to avoid
               extra computation */
    long    num_iter, count, index, left, right ; /* Auxiliaries for building
               bridges */
    double  width ;                               /* Width of current window in
               bridge formation */
    long    i ;                                   /* Temporaries */

    /* Start by initializing timepoint 0 and ending bridge as "computed" */

    bridges[start_point].curt = 1 ;
    bridges[end_point].curt   = 1 ;

    /* This is an unravelled recursion of building the bridge between
       start_point and end_point, halving */

    index           = start_index ;
    num_points      = end_point - start_point ;
    log2_num_points = (long) ceil(ln2inv*log((double)num_points)) ;
    width           = num_points * 1.0 ;

    for (num_iter=0; num_iter < log2_num_points; num_iter++)
    {

        /* Initialize for this iteration */

        count = pow2[num_iter] ;
        left  = start_point ;
        right = start_point ;

        /* Build "count" no. of bridges between left and right */

        for (i=0; (i < count) && (index < num_points-1); i++) {

            right = (long) ceil((i+1)*width) ;
            while (bridges[right].curt IS 0)
                right++ ;

            /* Try building a bridge between current left and right */

            if (aux_set_Order(left,
                              right,
                              index,
                              order_table,
                              bridges)
                    IS BBSUCCESS)
                index++ ;
            left = right ;
        }

        /* Prepare for next iteration */

        width = width * 0.5 ;
    }

    return BBSUCCESS ;

}

/*  Julia Chislenko, October 1997 - use this function instead of the old
 *  Set order according to the order of sobol dims - */
int
set_Order(long             start_point,
          long             end_point,
          long             start_index,
          Order            *order_table, /* Should be order_table of current
                          factor */
          bridgeEnd        *bridges)     /* Should be bridges of current
      factor */
{
    long    num_points ;                          /* Number of bridges */
    long    count, index, left, right ;           /* Auxiliaries for building
               bridges */
    long    min, mid=0;

    /* Start by initializing timepoint 0 and ending bridge as "computed" */

    num_points      = end_point - start_point ;
    index           = start_index ;

    bridges[start_point].curt = 1 ; /* to have always the first one
                built between the first and
                the last points */
    bridges[end_point].curt   = 1 ;

    if(num_points == 1)
        return BBSUCCESS;

    left  = start_point ;
    right = end_point ;

    for (index = start_index; index<num_points-1; )
    {
        /* Find the smallest dimType between left and right
         */
        min = SGM_SOB_MAX+1;
        for(count=left+1; count < right; count++) {
            if(bridges[count].dimensionType<min) {
                min = bridges[count].dimensionType;
                mid = count;
            }
        }

        if(min <= SGM_SOB_MAX) /* if left+1<right */
        {
            /* instead of aux_set_Order()
             */
            if (bridges[mid].curt IS 0) /* mid not yet computed */
            {

                order_table[index].left  = &bridges[left]  ;
                order_table[index].right = &bridges[right] ;
                order_table[index].mid   = &bridges[mid]   ;
                bridges[mid].curt        = 1               ; /* mid has been
                                      computed */
                index++;
            } else
                return BBFAILURE ;
        }

        left = (right==end_point ? start_point :right);
        right = left+1;
        while(right<end_point &&
                bridges[right].curt<1)
            right++;
    }

    return BBSUCCESS ;

}

/* Prints the order table - Used for debugging only */

void
print_Order_Table(void)
{
    long i, j ;

    for (i=0; i<context.numFactorsPerPath; i++) {
        printf("Factor %ld : \n",i) ;
        for (j=0; j<context.numBridgesPerFactor[i]-1; j++)
            printf("Order %ld : Left = %ld (Type = %ld) Right = %ld "
                   "(Type = %ld) Mid = %ld (Type = %ld)\n",
                   j,
                   orderTable[i][j].left->t,  orderTable[i][j].left->dimensionType,
                   orderTable[i][j].right->t, orderTable[i][j].right->dimensionType,
                   orderTable[i][j].mid->t,   orderTable[i][j].mid->dimensionType) ;
    }

}

/* This routine sets the dimension types across all factors when no
   dimensionType is specified in input.
   It goes across factors, starting with GRID at last dimension of first
   factor, 1st Sobol dimension at last
   dimension of second factor and so on .... and then halves the interval
   across all factors .......
   This routine assumes that all factors have same number of dimensions.
   If not, the dimensionTypes cannot be obtained internally.
   */

int
set_Default_Dimension_Types(void)
{
    long    num_points ;                          /* Number of dimensions along
               each factor */
    long    log2_num_points ;                     /* Temporary -- used to avoid
               extra computation */
    long    num_iter, count, index ;              /* Auxiliaries for building
               bridges */
    long    left, right, mid ;
    double  width ;                               /* Width of current window in
               bridge formation */
    long    i, j ;                                /* Temporaries */
    long    cur_sobdim ;                          /* Current sob dim to be
               assigned to next dimension*/

    int     *is_Computed = NULL;                  /* Used to avoid visiting a
               dimension type */
    int     status = BBSUCCESS ;

    /* Initialization */

    num_points = context.maxDimsPerFactor ;

    /* Allocate memory to is_Computed */

    if ((is_Computed = NEW_ARRAY(int, num_points+1)
                       /* (int *) malloc ((num_points+1) * sizeof(int)) */
        ) IS (int *) NULL) {
        sprintf(logBuffer,"set_Default_Dimension_Types : Memory allocation failed for "
                "is_Computed\n") ;
        SCLogMessage(logBuffer, FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* Initialize is_Computed */

    is_Computed[0]          = TRUE ;
    is_Computed[num_points] = TRUE ;

    for (i=1; i < num_points; i++)
        is_Computed[i] = FALSE ;

    /* Allocate last dimension first to set off the recursion */

    cur_sobdim = 1 ;

    /* Set dimension type to PSEUDO_RANDOM by default */

    for (i=0; i < context.numFactorsPerPath; i++)
        for (j=0; j < num_points; j++)
            context.dimensionType[i][j] = SGM_PSEUDO_RANDOM ;

    context.dimensionType[0][num_points-1] = SGM_GRID ;

    for (j=1; j < context.numFactorsPerPath; j++) {
        if (cur_sobdim <= SGM_SOB_MAX) {
            context.dimensionType[j][num_points-1] = cur_sobdim ;
            cur_sobdim++ ;
        }
    }

    /* This is an unravelled recursion of building the bridge between
       start_point and end_point,
       halving */

    index           = 0 ;
    log2_num_points = (long) ceil(ln2inv*log((double)num_points)) ;
    width           = num_points * 1.0 ;

#ifdef SGM_DEBUG

    printf("Initializing default dim type matrix: ln2inv=%f log()=%f"
           " ln*log()=%f ceil()=%f log2_num_points=%ld\n",
           ln2inv,log(num_points), ln2inv*log(num_points),
           ceil(ln2inv*log(num_points)), log2_num_points) ;
#endif

    for (num_iter=0; num_iter <= log2_num_points; num_iter++) {

#ifdef SGM_DEBUG
        printf("Num Iter = %ld\n",num_iter) ;
#endif

        /* Initialize for this iteration */

        count = pow2[num_iter] ;
        left  = 0 ;
        right = 0 ;

        /* Build "count" no. of bridges between left and right */

        for (i=0; (i < count) && (index < num_points-1); i++) {

            right = (long) ceil((i+1)*width) ;
            while (is_Computed[right] IS FALSE)
                right++ ;

            /* Try building a bridge between current left and right */

            mid = (long) (0.5 * (left + right)) ;

#ifdef SGM_DEBUG

            printf("left = %ld right = %ld mid = %ld width = %f\n",
                   left, right, mid, width) ;
#endif

            if (is_Computed[mid] IS FALSE) {

                /* Set dimensions across all factors */

                for (j=0; j < context.numFactorsPerPath; j++) {
                    if (cur_sobdim <= SGM_SOB_MAX) {
                        context.dimensionType[j][mid-1] = cur_sobdim ;
                        cur_sobdim++ ;
                    }
                }

                is_Computed[mid] = TRUE ;
                index++ ;

            }

            left = right ;
        }


        /* Prepare for next iteration */

        width = width * 0.5 ;
    }

    /* Swap half of the types - JC */
    /*  mid = num_points/2;
    for (j=0; j<mid-2; j+=3)
    {
        i = context.dimensionType[0][j];
        context.dimensionType[0][j] = context.dimensionType[0][j+mid];
        context.dimensionType[0][j+mid] = i;
    }
    */
done :

#ifdef SGM_DEBUG

    for (i=0; i < context.numFactorsPerPath; i++)
        for (j=0; j < num_points; j++)
            printf("Factor = %ld Dimension = %ld Type = %ld\n",
                   i, j, context.dimensionType[i][j]) ;
#endif

    if (is_Computed IS_NOT (int *) NULL)
        FREE(is_Computed) ;

    return status ;
}


/* *******************************************************************
   Returns a gaussian deviate based on dimension of the timepoint and
   current runNumber and pathNumber
   MFH : Modification. Allowed use of prespecifed grid in uniformGrid.
         Note the grid is passed in as a gaussian grid already.
   ********************************************************************/
inline double get_Gaussian_Deviate(long             curFactor, /* Current factor being
                                                                                processed */
                                   long             dimType)   /* Dimension Type of the
   factor   */
{
    long   runNum  ;  /* Used in grid and sobol points */
    long   pathNum ;

    double uniform ;  /* Deviate on uniform [0,1] support */
    double z       ;  /* Gaussian deviate */

    if (dimType IS SGM_PSEUDO_RANDOM)
    {
        //z = SC_gasdev2/*2*/(&context.seed) ;
        z = context.rng->fetch();
    }
    else if (dimType IS SGM_GRID)
    {

        pathNum = context.curPathNum ;
        runNum  = context.curRunNum  ;
        // Check to see if we have prespecified grid
        if(!context.uniformGrid.gaussianGridPtr) {
            uniform = context.uniformGrid.base +
                      context.uniformGrid.grid_size * (pathNum);

            // To avoid uniform deviate being zero, set it to 1/10th of step_size.

            if (uniform IS 0.0)
                uniform = 0.1 * context.uniformGrid.step_size ;
            z       = SC_InvertCumNormal(uniform) ;
        } else { // If a grid is specified. Use it!
            z = context.uniformGrid.gaussianGridPtr[pathNum];
        }


    }
    else /* dimType is SOBOL */
    {

        pathNum = context.curPathNum % context.numPathsPerFetch ;
        uniform = context.nextsobs[pathNum][dimType-1] ; /* In array, dimension
                   1 is in 0th entry */
        z       = SC_InvertCumNormal(uniform) ;
    }

#ifdef SGM_DEBUG
    if (z IS 0.0)
    {
        printf("Zero gaussian \n") ;
    }
#endif

    return z ;
}

/* Builds a brownian bridge between left and right at mid.
   Look at documentation on brownian bridges for the formula employed here.
   */
inline int
get_Conditional_Deviate(long   left,          /* Left  timepoint */
                        long   right,         /* Right timepoint */
                        long   mid,           /* Mid   timepoint */
                        double left_deviate,  /* Left  deviate   */
                        double right_deviate, /* Right deviate   */
                        double *mid_deviate,  /* Mid   deviate - to be
                                                    computed */
                        double z)             /* Random shock at mid */
{
    int    status ;
    double t, T ; /* The bridge is built with 0 at left, T at right and t at
           mid */
    double invT;

    status = BBSUCCESS ;

    /* The following check should be removed once the code is tested */

    //   if ((mid <= left) || (mid >= right) || (left >= right)) {

    //     sprintf(logBuffer,"get_Conditional_Deviate : Invalid left = %ld mid = %ld "
    //        "right = %ld\n",
    //      left, mid, right) ;
    //     SCLogMessage(logBuffer,FAILURE);
    //     status = BBFAILURE ;
    //     goto done ;
    //   }

    t = (mid - left)  ;
    T = (right-left)  ;

    /* Brownian bridge formula */

    /*
    (*mid_deviate) = (1.0 - t/T) * left_deviate + t/T * right_deviate +
        sqrt(t - t*t/T) * z ;
    */
    // Reduce operation count by predividing.

    invT = t/T;
    (*mid_deviate) = (1.0 - invT) * left_deviate + invT * right_deviate +
                     sqrt(t - t*invT) * z ;


#ifdef SGM_DEBUG

    if ((*mid_deviate) IS 0.0)
    {
        printf("Zero conditional gaussian : left = %f right = %f z = %f\n",
               left_deviate, right_deviate, z) ;
    }
#endif

done :

    return status ;
}

/* The following routine is one of the most important ones.
   It generates the deviates at bridges by using the conditional
   distribution of lognormal given
   fixed initial and ending point.
   Once the bridge deviates are computed, the intermediate points
   can be generated sequentially
   again by building bridges.
   */

int
build_Bridges(bridgeEnd *bridges,        /* Bridges in the current factor */
              Order     *order_table,    /* Order Table for computing bridges
                              in the current factor */
              long      factor)          /* Current factor */
{
    long order_index ; /* Current index into order_table */
    /*  long bridgeNum   ; *//* Number of bridge being currently processed */
    long dimType     ; /* Dimension type of current bridge being processed */
    long lastBridge  ; /* Number of last bridge */
    long pathNum     ; /* Number of path for which deviates are being generated*/
    long startPathNum; /* Original number of current path in context -
         context.curPathNum to be reset to this value at the
         end of the routine */
    double z         ; /* Gaussian shock for current bridgePoint - could be from
         a uniform grid or from Sobol dimension */

    bridgeEnd *leftptr  ; /* Pointers to bridgeEnd structures for left, right
            and mid points */
    bridgeEnd *rightptr ;
    bridgeEnd *midptr   ;

    int status ;

    status = BBSUCCESS ;

    startPathNum = context.curPathNum ;

#ifdef SGM_DEBUG

    printf("Generating bridge points for factor %ld\n",factor) ;
#endif

    for (pathNum=0; pathNum < context.numPathsPerFetch; pathNum++)
    {

        context.curPathNum = startPathNum + pathNum ;

        /* Initialize timepoint 0 */

        bridges[0].deviates[pathNum] = 0.0 ;

        /* First set the final bridge point using standard lognormal formula */

        lastBridge = context.numBridgesPerFactor[factor] ;
        dimType    = bridges[lastBridge].dimensionType ;
        z          = get_Gaussian_Deviate(factor,
                                          dimType) ;

        /* left = 0 for first bridge, right = T */

        bridges[lastBridge].deviates[pathNum] =
            sqrt(bridges[lastBridge].t*1.0) * z ;

        /* Now set the rest of the bridges in place by Brownian bridge
           construction */

        for (order_index=0; order_index < lastBridge-1; order_index++) {

            leftptr  = order_table[order_index].left ;
            rightptr = order_table[order_index].right ;
            midptr   = order_table[order_index].mid ;

            //       cout << leftptr->t << " " << midptr->t << " " << rightptr->t << " ";
            //       cout << leftptr->dimensionType << " "
            //     << midptr->dimensionType << " "
            //     << rightptr->dimensionType << " : ";
            //       cout << leftptr->deviates[pathNum] << " "
            //     << rightptr->deviates[pathNum] << endl;

            dimType  = midptr->dimensionType ;
            z        = get_Gaussian_Deviate(factor,
                                            dimType) ;

            if (get_Conditional_Deviate(leftptr->t, /* left  time point of bridge  */
                                        rightptr->t,/* right time point of bridge  */
                                        midptr->t,  /* mid   time point of bridge  */
                                        leftptr->deviates[pathNum],/* left deviate */
                                        rightptr->deviates[pathNum],/* right deviate*/
                                        &midptr->deviates[pathNum],/*mid  deviate
                                                                                        to be computed */
                                        z)                      /* gaussian shock */
                    IS BBFAILURE) {

                sprintf(logBuffer,"build_Bridges : get_Conditional_Deviate failed\n") ;
                SCLogMessage(logBuffer,FAILURE);
                status = BBFAILURE ;
                goto done ;
            }
        }
    }

done :

    /* Set pathNum back to startPathNum in Context */

    context.curPathNum = startPathNum ;

    return status ;
}

/* Initializes the grid specification inside the context.
   Kept separate from initialize_Context for greater visibility.
 
   MFH : Changed so that grid can be passed in.
   */

void
initialize_Grid(double *pointer2Grid)
{
    /* JFC modified 5/24/99 - There is only 1 run */
    /* context.uniformGrid.grid_size = 10.0 / context.numPathsPerRun ;
     context.uniformGrid.step_size = context.uniformGrid.grid_size ;
     context.uniformGrid.base = -5.0 ;
     context.uniformGrid.top  = 5.0 ;
    */
    // Set these parameters anyway
    context.uniformGrid.grid_size = 1.0 / (context.numPathsPerRun) ;
    context.uniformGrid.step_size = context.uniformGrid.grid_size / context.numRuns ;
    context.uniformGrid.base = context.uniformGrid.step_size*.5;
    context.uniformGrid.top  = 1.0-context.uniformGrid.base;

    // Set up the grid pointer. Later if defined it is used over the internal generation
    // of the grid based on the above parameters.
    context.uniformGrid.gaussianGridPtr = pointer2Grid;

    //  printf("%lf %lf %i\n",context.uniformGrid.base,context.uniformGrid.top,
    //    context.numPathsPerRun);
    //  context.uniformGrid.base =0.;
    //  context.uniformGrid.base = 1.;


}

/* Allocate memory for nextsobs. For each sobol dimension, we need
   upto numPathsPerRun sequences stored at any point of time.
   */

int
initialize_Nextsobs(void)
{
    int  status = BBSUCCESS ;
    long i ;

    if (context.maxSobDim <= 0) {
        context.nextsobs = NULL;
        return BBSUCCESS;
    }

    if ((context.nextsobs = NEW_ARRAY(double*, context.numPathsPerFetch)
                            /*(double **) malloc (context.numPathsPerFetch * sizeof(double))*/
        )
            IS (double **) NULL) {

        sprintf(logBuffer,"Initialization : Memory allocation failed for nextsobs\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;

    }

    for (i=0; i < context.numPathsPerFetch; i++)

        if ((context.nextsobs[i] = NEW_ARRAY(double, context.maxSobDim)
                                   /*(double *) malloc (context.maxSobDim * sizeof(double))*/
            )
                IS (double *) NULL) {

            sprintf(logBuffer,"Initialization : Memory allocation failed for "
                    "nextsobs[%ld]\n",i) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

done :

    return status ;

}

/* Frees memory allocated to nextsobs */

void
free_Nextsobs(void)
{
    long i ;

    if (context.nextsobs IS_NOT (double **) NULL) {
        for (i=0; i < context.numPathsPerFetch; i++)
            if (context.nextsobs[i] IS_NOT (double *) NULL)
                FREE(context.nextsobs[i]) ;
        FREE(context.nextsobs) ;
        context.nextsobs = NULL;
    }

}

/* The following checks if dimensions passed for multi-factor
   monte carlo are valid.
   If yes, it updates the other entries in context for further processing.
   Else returns FAILURE.
   */

int
initialize_Context(long numFactorsPerPath,    /* Number of factors in a
                                             simulation */
                   long *numDimsPerFactor,    /* Number of dimensions in
                                        each factor */
                   long numPathsPerRun,       /* Number of paths(simulations)
                                        per run */
                   long numRuns,              /* Number of runs */
                   long numPathsPerFetch,     /* Number of paths accessed in
                                        one call to fetch_Sobol_Bridge.
                                        e.g. In our IR Monte Carlo
                                        simulation, we might want to
                                        access 1000 points per fetch
                                        due to memory reasons when the
                                        total numer of paths per run
                                        is actually 5000.
                                        */
                   long *dimensionType[])     /* Type of each dimension in each
  factor -
  could be SGM_PSEUDO_RANDOM,
  SGM_GRID or 1..SGM_SOB_MAX */
{
    int  used_Sobol_Dimension[SGM_SOB_MAX+1] ;  /* Used to test that sobol
             dimension is not repeated */
    int  used_Grid ;                          /* Used to test that grid has been
                  used is only once */

    long dim, factor, dimtype ;               /* Temporaries */
    int  status =0;
    long i ;

    context.numDimsPerFactor = NULL;
    context.numBridgesPerFactor = NULL;
    context.curBridge = NULL;
    context.dimensionType = NULL;

    /* Check for valid arguments before initializing context */

#ifdef SGM_DEBUG

    printf("Checking valid argument for initializing Context\n") ;

#endif

    if (numFactorsPerPath <= 0)
    {
        sprintf(logBuffer,"Initialization : Invalid numFactorsPerPath = %ld "
                "(Must be > 0)\n",
                numFactorsPerPath) ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    if (numPathsPerRun <= 0)
    {
        sprintf(logBuffer,"Initialization : Invalid numPathsPerRun = %ld (Must be > 0)\n",
                numPathsPerRun) ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    if (numRuns <= 0)
    {
        sprintf(logBuffer,"Initialization : Invalid numRuns = %ld (Must be > 0)\n",
                numRuns) ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    if (numPathsPerFetch <= 0)
    {
        sprintf(logBuffer,"Initialization : Invalid numPathsPerFetch = %ld (Must be > 0)\n",
                numPathsPerFetch) ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    if (numDimsPerFactor IS (long *) NULL)
    {
        sprintf(logBuffer,"Initialization : Invalid numDimsPerFactor Array "
                "(Null pointer)\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* If no dimensionType is passed, set up context.dimensionType
       to get default dimension types */

    if (dimensionType IS (long **) NULL)
    {

        context.default_dimTypes = TRUE ;
        if ((context.dimensionType = NEW_ARRAY(long*, numFactorsPerPath)
                                     /*(long **) malloc (numFactorsPerPath * sizeof (long *))*/
            ) IS (long **) NULL) {
            sprintf(logBuffer,"Initialization : Memory allocation failed for default "
                    "dimension types\n") ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        for (i=0; i < numFactorsPerPath; i++) {
            if ((context.dimensionType[i] = NEW_ARRAY(long, numDimsPerFactor[i])
                                            /*(long *) malloc (numDimsPerFactor[i] * sizeof(long))*/
                ) IS (long *) NULL) {
                sprintf(logBuffer,"Initialization : Memory allocation failed for default "
                        "dimension types\n") ;
                SCLogMessage(logBuffer,FAILURE);
                status = BBFAILURE ;
                goto done ;
            }

        }

#ifdef SGM_DEBUG
        printf("Allocated memory to dimension Type since null passed\n") ;

#endif

    } else
        context.default_dimTypes = FALSE ;

    /* Initialize */

    context.numFactorsPerPath   = numFactorsPerPath ;
    context.numPathsPerRun      = numPathsPerRun ;
    context.numRuns             = numRuns ;
    context.numPathsPerFetch    = numPathsPerFetch ;
    context.maxDimsPerFactor    = 0 ;
    context.maxSobDim           = 0 ;
    context.curRunNum           = 0 ;
    context.curPathNum          = 0 ;
    context.current_dimension   = 1;

    context.numDimsPerFactor    = (long *) NULL ;
    context.numBridgesPerFactor = (long *) NULL ;
    context.curBridge           = (long *) NULL ;
    context.nextsobs            = (double **) NULL ;

    /* Allocate memory to numDimsPerFactor */

#ifdef SGM_DEBUG

    printf("Allocating memory to context.numDimsPerFactor\n") ;

#endif

    if ((context.numDimsPerFactor = NEW_ARRAY(long, numFactorsPerPath)
                                    /*(long *) malloc (numFactorsPerPath * sizeof(long))*/
        )
            IS (long *) NULL)
    {
        sprintf(logBuffer,"Initialization : Memory allocation failed for "
                "context.numDimsPerFactor\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* Allocate memory to numBridgesPerFactor array */

#ifdef SGM_DEBUG
    printf("Allocating memory to context.numBridgesPerFactor\n") ;

#endif

    if ((context.numBridgesPerFactor = NEW_ARRAY(long, numFactorsPerPath)
                                       /*(long *) malloc (numFactorsPerPath * sizeof(long))*/
        )
            IS (long *) NULL)
    {
        sprintf(logBuffer,"Initialization : Memory allocation failed for "
                "context.numBridgesPerFactor\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    for (factor=0; factor < context.numFactorsPerPath; factor++)
    {

        context.numDimsPerFactor[factor] = numDimsPerFactor[factor] ;
        context.maxDimsPerFactor = std::max(context.maxDimsPerFactor,
                                            context.numDimsPerFactor[factor]) ;
    }

    if(context.maxDimsPerFactor == 0)
    {
        sprintf(logBuffer,"Initialization : Dimensions for all factors are 0."
                "Attempt to initialize Sobol bridges for 0 simulations.\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* Initialize static constants */

    ln2inv = 1.0 / log(2.0) ;
    compute_powers(2,
                   pow2,
                   SGM_POW_MAX,
                   context.maxDimsPerFactor) ;

    /* Check and initialize context */

    /* Allocate memory to curBridge array */

#ifdef SGM_DEBUG

    printf("Allocating memory to context.curBridge\n") ;

#endif

    if ((context.curBridge = NEW_ARRAY(long, numFactorsPerPath)
                             /*(long *) malloc (numFactorsPerPath * sizeof(long))*/
        )
            IS (long *) NULL)
    {
        sprintf(logBuffer,"Initialization : Memory allocation failed for "
                "context.curBridge\n") ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* Set dimensionType in context to passed dimensionType array.
       This is safe since dimensionType is never used after initialization. */

    if (context.default_dimTypes IS FALSE)
        context.dimensionType = dimensionType ;
    else
    {

#ifdef SGM_DEBUG
        printf("Assigning default dimension types\n") ;

#endif

        /* Use internally generated default dimension Types */

        if (set_Default_Dimension_Types() IS BBFAILURE)
            goto done ;
    }

#ifdef SGM_DEBUG
    printf("Initializing used_Sobol_Dimension array\n") ;

#endif

    init_int_array(used_Sobol_Dimension,
                   SGM_SOB_MAX+1,
                   FALSE) ;
    used_Grid = FALSE ;
    status    = BBSUCCESS ;

    /* Go through dimensions specified for each factor and check for
       repetitions */

    for (factor=0; factor < context.numFactorsPerPath; factor++)
    {

        /* Initialize uninitialized parts of the context */

        context.curBridge[factor]        = 0 ;

        if (context.dimensionType[factor] IS (long *) NULL) {
            sprintf(logBuffer,"Initialization : Invalid dimensionType[factor=%ld] "
                    "Array (Null pointer)\n",
                    factor) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        /* Initialize numBridgesPerFactor before processing dimension types */

        context.numBridgesPerFactor[factor] = 0 ;

        /* Check for valid dimensions and dimension types across all factors */

        for (dim=0; dim < numDimsPerFactor[factor]; dim++) {

            dimtype = context.dimensionType[factor][dim] ;

            /* Check if dimtype is out of bounds */

            if ((dimtype < SGM_PSEUDO_RANDOM) || (dimtype > SGM_SOB_MAX)) {

                sprintf(logBuffer,"Initialization : Invalid dimension type = %ld "
                        "Factor = %ld Dimension = %ld\n",
                        dimtype,factor,dim) ;
                SCLogMessage(logBuffer,FAILURE);
                status = BBFAILURE ;
                goto done ;
            }

            /* If Sobol dimension - Check repetition.
               If valid Sobol dimension, update maxSobDim and numBridgesPerFactor
            entries. */

            else if (is_Sobol(dimtype) IS TRUE) {

                if (used_Sobol_Dimension[dimtype] IS FALSE) {

                    used_Sobol_Dimension[dimtype] = TRUE ;
                    context.numBridgesPerFactor[factor]++ ;
                    context.maxSobDim = std::max(context.maxSobDim,
                                                 dimtype) ;
                } else {

                    sprintf(logBuffer,"Initialization : Sobol dimension %ld used for two points :"
                            " Second Usage - Factor = %ld Dimension = %ld\n",
                            dimtype,factor,dim) ;
                    SCLogMessage(logBuffer,FAILURE);
                    status = BBFAILURE ;
                    goto done ;
                }
            }

            /* If Grid - Check repetition */

            else if (dimtype IS SGM_GRID) {

                if (used_Grid IS FALSE) {
                    used_Grid = TRUE ;
                    context.numBridgesPerFactor[factor]++ ;
                } else {
                    sprintf(logBuffer,"Initialization : Grid used for two points : "
                            "Second Usage - Factor = %ld Dimension = %ld\n",
                            factor,dim) ;
                    SCLogMessage(logBuffer,FAILURE);
                    status = BBFAILURE ;
                    goto done ;
                }
            }
        }

    }

    /* If all above successful, allocate memory for storage of sobol sequences */

#ifdef SGM_DEBUG
    printf("Context : maxDimsPerFactor = %ld maxSobDim = %ld\n",
           context.maxDimsPerFactor, context.maxSobDim) ;
    for (i=0; i<context.numFactorsPerPath; i++)
        printf("numBridgesPerFactor[%ld] = %ld\n",
               i, context.numBridgesPerFactor[i]) ;
    printf("Initializing internal storage of Sobol points\n") ;

#endif

    if ((status = initialize_Nextsobs()) IS BBFAILURE)
        goto done ;

    /* Finished all checking or found errors -- Quit now */

done :

    return status ;
}


/* Frees the memory allocated to context.
   Should be called on free call as well as if an error occurs during
   initialization to ensure
   that there is no memory leakage.
   */

void
free_Context(void)
{
    long factor ;

    if (context.numDimsPerFactor IS_NOT (long *) NULL)
        FREE(context.numDimsPerFactor) ;
    context.numDimsPerFactor    = NULL ;

    if (context.numBridgesPerFactor IS_NOT (long *) NULL)
        FREE(context.numBridgesPerFactor) ;
    context.numBridgesPerFactor = NULL ;

    if (context.curBridge IS_NOT (long *) NULL)
        FREE(context.curBridge) ;
    context.curBridge           = NULL ;

    if ((context.default_dimTypes IS TRUE) && (context.dimensionType IS_NOT
            (long **) NULL)) {
        for (factor=0; factor < context.numFactorsPerPath; factor++) {
            if (context.dimensionType[factor] IS_NOT (long *) NULL)
                FREE(context.dimensionType[factor]) ;
        }
        FREE(context.dimensionType) ;
        context.dimensionType = NULL ;
    }

    free_Nextsobs() ;
}

/* Final call to Sobol bridge.
   Cleans up all the static memory allocated in initialize_Sobol_Bridge(.)
   Must be called at the end of multi-factor Monte Carlo simulation
   to avoid any memory leakages.
   */

int
SGM_Free_Sobol_Bridge(void)
{
    free_Bridges() ;
    free_Order_Table() ;
    free_Context() ;
    return BBSUCCESS;
}

/* Initialization routine for monte carlo simulation.
   This does all memory allocation and must be associated with a closing free
   call to get rid of
   all memory it has allocated. It also updates data structures for building
   bridges.
   */

int
SGM_Initialize_Sobol_Bridge(long numFactorsPerPath,    /* Number of factors in
                                                                 a simulation */
                            long *numDimsPerFactor,    /* Number of dimensions
                                                          in each factor */
                            long numPathsPerRun,       /* Number of paths
                                                          (simulations)
                                                          per run */
                            long numRuns,              /* Number of runs */
                            long numPathsPerFetch,     /* Number of paths
                                                          accessed in one call
                                                          to fetch_Sobol_Bridge.
                                                          e.g. In our IR Monte
                                                          Carlo simulation, we
                                                          might want to
                                                          access 1000 points
                                                          per fetch due to
                                                          memory reasons when
                                                          the total numer of
                                                          paths per run is
                                                          actually 5000.
                                                          */
                            long *dimensionType[], /* Type of each dimension in
                                                             each factor -
                                                             could be
                                                             SGM_PSEUDO_RANDOM,
                                                             SGM_GRID or
                                                             1..SGM_SOB_MAX */
                            INormalSuperCubeRNGSP rng,
                            /*       long *seed,         */   /* Seed for random number
                                                                        generator of normally dist
                                                      generator */
                            long sampleTailsFlag,  /* TRUE/FALSE to sample
                                                             tails or not in each
                                                             dimension */

                            /* The following are the three arguments to Julia/Arnon's routine
                               for tail sampling. */

                            long   numSampDistPts,/* Number of points for
                                                              sampling distribution */
                            double sampIntvlBound,/* Determines sampling
                                                            distr interval
                                                            (-sampIntvlBound,
                                                            +sampIntvlBound) in
                                                            terms of standard
                                                            deviation */
                            double maxProbRatio,  /* Max ratio between path
                                                            probabilities. If 1.,
                                                            the samples  are normal */
                            double *ptr2GaussianGrid, /* Grid specification.
                                                         If NULL internal grid used.
                                                         User must specify grid in
                                                         non-uniform format. */
                            SequenceSP ldSeqPointer) /* pointer to a Sequence class.
      To be used for supports in BB */
{
    int status = BBSUCCESS ;

    long factor/*, dim */;                               /* Temporaries */

#ifdef SGM_DEBUG

    long i,j;

    printf("Initializing Monte Carlo Context\n") ;
    printf("NumFactors = %ld\n",numFactorsPerPath) ;
    printf("NumPathsPerRun = %ld\n",numPathsPerRun) ;
    printf("NumRuns = %ld\n",numRuns) ;
    printf("NumPathsPerFetch = %ld\n",numPathsPerFetch) ;
    for (i=0; i<numFactorsPerPath; i++)
    {
        printf("numDimsPerFactor[%ld] = %ld\n",i,numDimsPerFactor[i]) ;
        if (dimensionType IS_NOT (long **) NULL) {
            for (j=0; j<numDimsPerFactor[i]; j++) {
                if (dimensionType[i] IS_NOT (long *) NULL) {
                    printf("Type of dimension[%ld][%ld] = %ld\n",i,j,
                           dimensionType[i][j]) ;
                }
            }
        }
    }
#endif

    /* Ensure that low discrepancy sequence for brownian bridge is properly initialized */
    if(!ldSeqPointer)
    {
        sprintf(logBuffer,"SGM_initialize_Sobol_Bridge : ldSeqPointer not defined. Must pass a valid sequence pointer\n");
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    } else
    {
        context.LDSequencePtr = ldSeqPointer;
    }

    /* Initialize static memory pointers for careful exits */

    bridges    = (bridgeEnd **) NULL ;
    orderTable = (Order     **) NULL ;

    if ((status = initialize_Context(numFactorsPerPath,
                                     numDimsPerFactor,
                                     numPathsPerRun,
                                     numRuns,
                                     numPathsPerFetch,
                                     dimensionType)) IS BBFAILURE)
        goto done ;

    /* Initialize bridges, order_tables and nextsobs arrays */

#ifdef SGM_DEBUG

    printf("Initializing Bridges\n") ;
#endif

    if ((status = initialize_Bridges()) IS BBFAILURE)
        goto done ;

#ifdef SGM_DEBUG

    print_Bridges() ;
    printf("Initializing Order Table\n") ;
#endif

    if ((status = initialize_Order_Table()) IS BBFAILURE)
        goto done ;

    /* Create order for building bridges along each factor to avoid this
       calculation on a repeated basis */

#ifdef SGM_DEBUG

    printf("Setting Order to build bridges once and for all\n") ;

#endif

    for (factor=0; factor < context.numFactorsPerPath; factor++)
        if (context.numBridgesPerFactor[factor] > 1)
            if((status = set_Order(0,
                                   /* start_point for
                                   bridges */
                                   context.numBridgesPerFactor[factor],  /* end_point   for
                                                                            bridges */
                                   0,                                    /* start_index for
                                                                            order_table */
                                   orderTable[factor],                   /* order_table of
                                                                            current factor */
                                   bridges[factor])) == BBFAILURE)          /* bridges     of
                                      current factor */
                goto done;

#ifdef SGM_DEBUG

    print_Order_Table() ;
    printf("Grid parameters\n") ;
#endif

    //=======================================================================
    // SOBOL generator should have been initialized coming into this routine.
    //=======================================================================

    //   /* Initialize arrays for sobol sequences and initialize grid
    //      specification */
    //   Sobol_gen = (int **)malloc(sizeof(int*)*MAXBIT);
    //   for(dim=0; dim < MAXBIT; dim++){
    //     Sobol_gen[dim] = (int*)malloc(sizeof(int)*DIMENSION);
    //   }
    //   Sobol_xn = (int *)malloc(sizeof(int)*DIMENSION);
    //   Sobol_n  = (int *)malloc(sizeof(int)*DIMENSION);


    //   /* Dummy call to Arnon's sobol sequence generator for initialization */
    //   if(context.maxSobDim > 0){
    //     sobinit(context.maxSobDim,
    //      Sobol_gen,
    //      Sobol_xn,
    //      Sobol_n,
    //      &Sobol_nstream);
    //     //SobolSequenceAL(0,context.maxSobDim-1,start,context.nextsobs[0],1) ;

    //     // Generate the first 50 sobol
    //     for(dim=0; dim<50; dim++){
    //       sobvect(context.maxSobDim,
    //        context.nextsobs[0],
    //        Sobol_gen,
    //        Sobol_xn,
    //        Sobol_n,
    //        &Sobol_nstream );
    //     }
    //   }
    // Initialize the grid
    initialize_Grid(ptr2GaussianGrid) ;

#ifdef SGM_DEBUG

    printf("Uniform grid : Base = %f Top = %f grid_size = %f step_size = %f\n",
           context.uniformGrid.base, context.uniformGrid.top,
           context.uniformGrid.grid_size, context.uniformGrid.step_size) ;

#endif

    /* Set seed */
    //rng->init(); moved init() into superCubeAdjustment
    rng->superCubeAdjustment();
    // This is buggy, as we really want to use the same generator;
    // See bug description in SuperCubeBBridge.cpp:788
    context.rng = DYNAMIC_POINTER_CAST<INormalSuperCubeRNG> (rng->clone()); // FIXME FIXME

    //    context.rng = INormalSuperCubeRNGSP(rng);
#if 0 // original code

    SC_init_gasdev2();
    if (*seed < 0)
        context.seed = (*seed) ;
    else
    {

        // CHANGE this to now set the seed to be -1 * input seed
        /* If initialization seed is not negative, get one using system clock */
        //context.seed = - ( (long) (time( (time_t *) NULL)) % 1000000 ) ;
        //SC_gasdev2/*2*/(&context.seed) ;
        //(*seed) = context.seed ;
        context.seed = - (*seed);
    }
    (*seed) = context.seed;
#endif

#ifdef SGM_DEBUG

    printf("Seed set to %ld\n",(*seed)) ;

#endif

    /* All done or found errors -- Quit now */

done :

    if (status IS BBFAILURE)
        SGM_Free_Sobol_Bridge() ;

    return status ;
}


/* The following function should be called at the beginning of each fetch
   of Monte Carlo simulation.
   One fetch -> numPathsPerFetch have been accessed for all dimensions for
   all factors.
   It prepares the context of monte carlo for the current run for the current
   fetch, generates bridge end-points and performs other initializations.
   */

int
SGM_Prepare_Sobol_Bridge(void)
{
    int  status ;
    long factor ;
    long bridgeNum ;

    status = BBSUCCESS ;

    /* Get next set of Sobol points in context for this run */

#ifdef SGM_DEBUG

    printf("Run %ld : Getting next set of Sobol sequences\n",
           context.curRunNum) ;

#endif

    get_next_sobs(context.nextsobs,
                  context.numPathsPerFetch,
                  context.maxSobDim) ;

#ifdef SGM_DEBUG

    printf("Run %ld : Generating bridge points for all factors\n",
           context.curRunNum) ;

#endif

    for (factor=0; factor < context.numFactorsPerPath; factor++) {

        /* Build bridges along each factor */

        if (build_Bridges(bridges[factor],
                          orderTable[factor],
                          factor) IS BBFAILURE) {
            status = BBFAILURE ;
            goto done ;
        }

        /* Initialize curBridge in each factor to 0 */

        context.curBridge[factor] = 0 ;

        /* Set all curt's in bridge to t's */

        for (bridgeNum=0; bridgeNum <= context.numBridgesPerFactor[factor];
                bridgeNum++)
            bridges[factor][bridgeNum].curt = bridges[factor][bridgeNum].t ;
    }

done :

    return status ;
}

/* The following is used to fetch random gaussian deviates from this generator.
   It can be used to access one entire factor Or
   It can be used to access points for an entire run along a given dimension in
   a given factor.
   */

int
SGM_Fetch_Sobol_Bridge(long typeFetch,       /* SGM_FETCH_ONE_FACTOR or
                                                    SGM_FETCH_ONE_DIM_ONE_FACTOR
                                                    or
                                                    SGM_FETCH_ONE_DIM_ALL_FACTORS
                                                    */
                       long factor,          /* factor to be fetched */
                       long dimension,       /* dimension to be fetched -
                                              ignored if
                                              SGM_FETCH_ONE_FACTOR */
                       long pathIndex,       /* path being fetched. Used only
                                              for
                                              SGM_FETCH_ONE_DIM_ALL_FACTORS
                                              and must be between 0 and
                                              numPathsPerFetch-1. */
                       //         long *seed,           /* seed for random number
                       //            generator */
                       // unused               INormalSuperCubeRNGSP    rng,
                       double *gaussians,    /* Set of gaussian deviates --
                                              must be allocated in memory
                                              in calling routine */
                       double *probabilities)/* The probabilities corresponding
to each point returned in
"gaussians". */
{
    int status = BBSUCCESS ;

    long lastBridge ;                        /* lastBridge in given factor */
    long curBridge  ;                        /* left endpoint of current bridge */
    long nextBridge ;                        /* right endpoint of current bridge*/

    long pathNum    ;                        /* to run through all points for a
                 dimension in
                 SGM_FETCH_ONE_DIM_ALL_FACTORS
                 or the current path in case of
                 SGM_FETCH_ONE_FACTOR. */
    long bridgeNum  ;

    long left ;                              /* All for construction of Brownian
                 bridge */
    long right ;
    double leftdeviate ;
    double rightdeviate ;
    double thisdeviate ;                     /* deviate generated for this
                 dimension */
    double z        ;                        /* gaussian shock generated for
                 this dimension */

    /******************** Fetch an entire factor *****************************/

    if (typeFetch IS SGM_FETCH_ONE_FACTOR)
    {


        /* Check if factor is within bounds */

        if ((factor < 0) || (factor >= context.numFactorsPerPath)) {

            sprintf(logBuffer,"fetch_Sobol_Bridge : Invalid factor %ld : Must be "
                    "between %d and %ld\n",
                    factor, 0, context.numFactorsPerPath-1) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        /* All arguments ok -- Process now */


        /* The following line modified by Viral Acharya - 6th Feb 1998 -- % .... is the stuff added */

        pathNum = context.curPathNum % context.numPathsPerFetch;

        lastBridge = context.numBridgesPerFactor[factor] ;
        curBridge  = context.curBridge[factor] ;
        nextBridge = curBridge + 1 ;

        /* Process each dimension of this factor */

        for (dimension=1; dimension <= context.numDimsPerFactor[factor];
                dimension++) {

            if (dimension > bridges[factor][lastBridge].t) /* All bridges built */
            {

                gaussians[dimension-1] = get_Gaussian_Deviate(factor,
                                         SGM_PSEUDO_RANDOM) ;

#ifdef SGM_DEBUG

                printf("Factor =%ld Dimension = %ld CurPathNum = %ld ... "
                       "Past all bridges\n",
                       factor, dimension, context.curPathNum) ;

#endif

            }
            /* Dimension is a bridge endpoint */
            else if (bridges[factor][nextBridge].t IS dimension)  {

                gaussians[dimension-1] = bridges[factor][nextBridge].deviates[pathNum]
                                         - bridges[factor][curBridge].deviates[pathNum] ;

                curBridge  = nextBridge ;
                nextBridge = curBridge + 1 ;

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Hit bridge %ld\n",
                       factor, dimension, context.curPathNum, nextBridge) ;

#endif

            }
            else  { /* The timepoint is within a bridge. So need to build the
                               bridge and shift the left-point of
                               bridge to this timepoint. */

                left  = bridges[factor][curBridge].curt ;
                right = bridges[factor][nextBridge].curt ;

                leftdeviate  = bridges[factor][curBridge].deviates[pathNum] ;
                rightdeviate = bridges[factor][nextBridge].deviates[pathNum] ;
                z = get_Gaussian_Deviate(factor,
                                         SGM_PSEUDO_RANDOM) ;
                if (get_Conditional_Deviate(left,     /* left  time point of bridge  */
                                            right,    /* right time point of bridge  */
                                            dimension,/* mid   time point of bridge  */
                                            leftdeviate,  /* left  deviate */
                                            rightdeviate, /* right deviate */
                                            &thisdeviate, /* mid deviate - to be
                                                                                           computed */
                                            z)            /* gaussian shock */
                        IS BBFAILURE) {

                    sprintf(logBuffer,"fetch_Sobol_Bridge : get_Conditional_Deviate failed : "
                            "factor = %ld left = %ld right = %ld mid = %ld\n",
                            factor, left, right, dimension) ;
                    SCLogMessage(logBuffer,FAILURE);
                    status = BBFAILURE ;
                    goto done ;
                }

                gaussians[dimension-1] = thisdeviate - leftdeviate ;

                /* Update deviates at left endpoint of bridge and move current bridge
                   endpoint to the right */

                bridges[factor][curBridge].deviates[pathNum] = thisdeviate ;
                bridges[factor][curBridge].curt = dimension ;

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Within bridges %ld and %ld\n",
                       factor, dimension, context.curPathNum, curBridge, nextBridge) ;

#endif

            }

        }

        /* Check if pathNum needs to be incremented by one */

        if (factor IS (context.numFactorsPerPath-1))
            context.curPathNum++ ;

        /* Check if run is completed */

        if (context.curPathNum IS context.numPathsPerRun) {
            context.curRunNum++ ;
            context.curPathNum = 0 ;
        }

        /* Set all curt's in bridges back to t's */

        for (bridgeNum=0; bridgeNum <= context.numBridgesPerFactor[factor];
                bridgeNum++)
            bridges[factor][bridgeNum].curt = bridges[factor][bridgeNum].t ;
    }

    /*************** Fetch an entire dimension ********************************/

    else if (typeFetch IS SGM_FETCH_ONE_DIM_ONE_FACTOR)
    {

        /* Check if factor is within bounds */

        if ((factor < 0) || (factor >= context.numFactorsPerPath)) {

            sprintf(logBuffer,"fetch_Sobol_Bridge : Invalid factor %ld : Must be between "
                    "%d and %ld\n",
                    factor, 0, context.numFactorsPerPath-1) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        /* Check if dimension is within bounds for the given factor */

        if ((dimension <= 0) || (dimension > context.numDimsPerFactor[factor])) {

            sprintf(logBuffer,"fetch_Sobol_Bridge : Invalid dimension %ld for factor %ld : "
                    "Must be between %d and %ld\n",
                    dimension, factor, 1, context.numDimsPerFactor[factor]) ;
            SCLogMessage(logBuffer,FAILURE);
            status = BBFAILURE ;
            goto done ;
        }

        /* All arguments ok -- Process now */

        lastBridge = context.numBridgesPerFactor[factor] ;
        curBridge  = context.curBridge[factor] ;
        nextBridge = curBridge + 1 ;

        if (dimension > bridges[factor][lastBridge].t) /* All bridges built */
        {

            for (pathNum=0; pathNum < context.numPathsPerFetch; pathNum++)
            {

                gaussians[pathNum] = get_Gaussian_Deviate(factor,
                                     SGM_PSEUDO_RANDOM) ;
                context.curPathNum++ ;

#ifdef SGM_DEBUG

                printf("Factor =%ld Dimension = %ld CurPathNum = %ld ... Past all "
                       "bridges\n",
                       factor, dimension, context.curPathNum) ;

#endif

            }

        }
        /* Dimension is a bridge endpoint */
        else if (bridges[factor][nextBridge].t IS dimension)  {

            for (pathNum=0; pathNum < context.numPathsPerFetch; pathNum++) {

                gaussians[pathNum] = bridges[factor][nextBridge].deviates[pathNum] -
                                     bridges[factor][curBridge].deviates[pathNum] ;
                context.curPathNum++ ;

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Hit bridge %ld\n",
                       factor, dimension, context.curPathNum, nextBridge) ;
                printf("Left bridge = %f Right bridge = %f\n",
                       bridges[factor][curBridge].deviates[pathNum],
                       bridges[factor][nextBridge].deviates[pathNum]) ;

#endif

            }

            context.curBridge[factor] = nextBridge ;

        }
        else  { /* The timepoint is within a bridge. So need to build the bridge
                           and shift the left-point of
                           bridge to this timepoint. */

            left  = bridges[factor][curBridge].curt ;
            right = bridges[factor][nextBridge].curt ;

            for (pathNum=0; pathNum < context.numPathsPerFetch; pathNum++) {

                leftdeviate  = bridges[factor][curBridge].deviates[pathNum] ;
                rightdeviate = bridges[factor][nextBridge].deviates[pathNum] ;
                z = get_Gaussian_Deviate(factor,
                                         SGM_PSEUDO_RANDOM) ;
                if (get_Conditional_Deviate(left,       /* left  time point of bridge */
                                            right,      /* right time point of bridge */
                                            dimension,  /* mid   time point of bridge */
                                            leftdeviate, /* left  deviate */
                                            rightdeviate,/* right deviate */
                                            &thisdeviate,/* mid   deviate - to be
                                                                                          computed */
                                            z)           /* gaussian shock */
                        IS BBFAILURE) {

                    sprintf(logBuffer,"fetch_Sobol_Bridge : get_Conditional_Deviate failed :"
                            " factor = %ld left = %ld right = %ld mid = %ld\n",
                            factor, left, right, dimension) ;
                    SCLogMessage(logBuffer,FAILURE);
                    status = BBFAILURE ;
                    goto done ;
                }

                gaussians[pathNum] = thisdeviate - leftdeviate ;

                /* Update deviates at left endpoint of bridge */

                bridges[factor][curBridge].deviates[pathNum] = thisdeviate ;

                context.curPathNum++ ;

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Within bridges %ld and %ld\n",
                       factor, dimension, context.curPathNum, curBridge, nextBridge) ;

#endif

            }

            /* Move current bridge endpoint to the right */

            bridges[factor][curBridge].curt = dimension ;

        }

        /* Check if this is the end of a Fetch and end of a run */

        if ((factor IS (context.numFactorsPerPath-1)) &&
                (dimension IS context.numDimsPerFactor[factor])) {
            if (context.curPathNum IS context.numPathsPerRun) {
                context.curRunNum++ ;
                context.curPathNum = 0 ;
            }
        } else
            context.curPathNum -= context.numPathsPerFetch ;
    }

    /*********** Fetch one dimension for all factors ***********************/

    else if (typeFetch IS SGM_FETCH_ONE_DIM_ALL_FACTORS)
    {

        pathNum = context.curPathNum % context.numPathsPerFetch ;

        for (factor=0; factor < context.numFactorsPerPath; factor++) {

            /* All arguments ok -- Process now */

            lastBridge = context.numBridgesPerFactor[factor] ;
            curBridge  = context.curBridge[factor] ;
            nextBridge = curBridge + 1 ;

            if (context.current_dimension > bridges[factor][lastBridge].t) {           /* All bridges built */

                gaussians[factor] = get_Gaussian_Deviate(factor,
                                    SGM_PSEUDO_RANDOM) ;

#ifdef SGM_DEBUG

                printf("Factor =%ld Dimension = %ld CurPathNum = %ld ... "
                       "Past all bridges\n",
                       factor, context.current_dimension, context.curPathNum) ;

#endif

            }
            /* Dimension is a bridge endpoint */
            else if (bridges[factor][nextBridge].t IS context.current_dimension)  {

                gaussians[factor] = bridges[factor][nextBridge].deviates[pathNum] -
                                    bridges[factor][curBridge].deviates[pathNum] ;
                //  printf("next : %lf prev : %lf gaussian: %lf\n",
                //         bridges[factor][nextBridge].deviates[pathNum],
                //         bridges[factor][curBridge].deviates[pathNum],
                //         gaussians[factor]);

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Hit bridge %ld\n",
                       factor, context.current_dimension,
                       context.curPathNum, nextBridge) ;
                printf("Left bridge = %f Right bridge = %f Gaussian = %f\n",
                       bridges[factor][curBridge].deviates[pathNum],
                       bridges[factor][nextBridge].deviates[pathNum],
                       gaussians[factor] ) ;

#endif

            }
            else  { /* The timepoint is within a bridge. So need to build the bridge
                               and shift the left-point of
                               bridge to this timepoint. */

                left  = bridges[factor][curBridge].curt ;
                right = bridges[factor][nextBridge].curt ;

                leftdeviate  = bridges[factor][curBridge].deviates[pathNum] ;
                rightdeviate = bridges[factor][nextBridge].deviates[pathNum] ;

#ifdef SGM_DEBUG

                printf("Factor = %ld Dimension = %ld CurPathNum = %ld ... "
                       "Within bridges %ld and %ld\n",
                       factor, context.current_dimension, context.curPathNum,
                       curBridge,
                       nextBridge) ;

#endif

                z = get_Gaussian_Deviate(factor,
                                         SGM_PSEUDO_RANDOM) ;
                if (get_Conditional_Deviate(left,       /* left  time point of bridge */
                                            right,      /* right time point of bridge */
                                            context.current_dimension, /* mid time point of
                                                                                         bridge  */
                                            leftdeviate,  /* left  deviate */
                                            rightdeviate, /* right deviate */
                                            &thisdeviate, /* mid   deviate - to be
                                                                                           computed */
                                            z)            /* gaussian shock */
                        IS BBFAILURE) {

                    sprintf(logBuffer,"fetch_Sobol_Bridge : get_Conditional_Deviate failed : "
                            "factor = %ld left = %ld right = %ld mid = %ld\n",
                            factor, left, right, context.current_dimension) ;
                    SCLogMessage(logBuffer,FAILURE);
                    status = BBFAILURE ;
                    goto done ;
                }

                gaussians[factor] = thisdeviate - leftdeviate ;

#ifdef SGM_DEBUG

                printf("Gaussian = %f\n",
                       gaussians[factor]) ;

#endif
                /* Update deviates at left endpoint of bridge */

                bridges[factor][curBridge].deviates[pathNum] = thisdeviate ;

            }

        }

        /* Check if all pathsPerFetch accessed for current_dimension.
           If yes, increment current_dimension and update current bridge
           end-point */

        context.curPathNum++ ;

        if (pathNum IS (context.numPathsPerFetch-1)) {

            /* Update all current Bridges since dimension has been accessed
            completely */

            for (factor=0; factor < context.numFactorsPerPath; factor++) {

                curBridge  = context.curBridge[factor] ;
                nextBridge = curBridge + 1 ;

                if (bridges[factor][nextBridge].t IS context.current_dimension)
                    /* Dimension is a bridge endpoint */
                    context.curBridge[factor] = nextBridge ;
                else /* Move the leftpoint of the current bridge */
                    bridges[factor][curBridge].curt = context.current_dimension ;

            }

            /* Move to next dimension */

            context.current_dimension++  ;

            /* Check if this is the end of a run */

            if (context.current_dimension > context.maxDimsPerFactor) {

                context.current_dimension = 1 ;
                if (context.curPathNum IS context.numPathsPerRun) {
                    context.curRunNum++ ;
                    context.curPathNum = 0 ;
                }

            }
            else

                context.curPathNum -= context.numPathsPerFetch ;

        }

    }
    else
    { /* fetchType is invalid */

        sprintf(logBuffer,"fetch_Sobol_Bridge : Invalid fetchType %ld : Must be %ld "
                "(entire factor) or %ld (one dimension, one factor) or %ld "
                "(one dimension, all factors)\n",
                typeFetch, SGM_FETCH_ONE_FACTOR, SGM_FETCH_ONE_DIM_ONE_FACTOR,
                SGM_FETCH_ONE_DIM_ALL_FACTORS) ;
        SCLogMessage(logBuffer,FAILURE);
        status = BBFAILURE ;
        goto done ;
    }

    /* Return the seed for random number generation */

    //  (*seed) = context.seed ; Not needed now

done :

    return status ;
}

CORE_END_NAMESPACE
