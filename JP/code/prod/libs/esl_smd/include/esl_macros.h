#ifndef ESL_MACROS_DOT_H
#define ESL_MACROS_DOT_H

/** NOTE: This file should be only included through 'esl.h'
 */


#ifdef  __cplusplus
extern "C" {
#endif

#ifndef PI
extern const double PI; /* valued in esl_macros.h */
#endif

#define  TRUE          1
#define  FALSE         0
#define  SUCCESS       0
#define  FAILURE       -1
#define  MAXBUFF       250

/* global constant variables are valued in esl_util.h */
extern const double TINY;
extern const double QCUTOFF;   /* Normal model for |q|<QCUTOFF  */

#define  BIG           exp (15.)
#define  ERROR         1E-7           /* Error tolerance                */
#define  BARRIER_TOL   1E-4           /* Error tolerance for barriers   */
#define  JUMPCOEFF     3.0            /* Jump size coefficient          */
#define  NBCRITDATE    50             /* Number of critical date arrays */
#define  MAX_INST      NBCRITDATE
#define  Nb_Daily_Pts  0              /* Number of daily points         */

/** the first zerobank event type  */
#define  ZbkEVENT      (NBCRITDATE-3) 

#define  MAXNBDATE     500    /* Max nb of elements in input date array */
#define  MAXNBSTATES   200    /* Max nb of state variables              */
#define  MAXNBEVCURVES 5      /* Max nb of event curves in EVENT_LIST   */
#define  NBSWAPTION    25     /* Maximum number of volatilities in swaption matrix  */
#define  MAXNBDATE     500    /* Maximum number of elements in input date array     */

#define   MAX_ITERATIONS    50
#define   MAXINDEX        8   /* Maximum number of characters in index name         */

/* Types for memory allocation */
#define     CRITDATE        0
#define     INT             1
#define     INT_PTR         2
#define     INT_D_PTR       3
#define     LONG            4
#define     LONG_PTR        5
#define     LONG_D_PTR      6
#define     DOUBLE          7
#define     DOUBLE_PTR      8
#define     DOUBLE_D_PTR    9
#define     CHAR            10
#define     CHAR_PTR        11
#define     MIN_SMD_CORR    -0.95
#define     MAX_SMD_CORR    0.95





/** different mode choices for GetDLOffset function */
#define  CbkEXACT   0
#define  CbkLOWER   -1
#define  CbkHIGHER  1

/* total nb bytes in char block is this +1 */
#ifndef GTO_MAX_STR_LEN
#define GTO_MAX_STR_LEN    127 
#endif

 /* nb bytes for array count: 1 => max 255 elements in char block array */
#ifndef GTO_NUM_BYTES_FOR_STR_CNT
#define GTO_NUM_BYTES_FOR_STR_CNT    1
#endif

/* Macros */
#ifndef NEAR_INT
#define NEAR_INT(a) (int) ((a) + ((a) >= 0. ? 0.5 : -0.5));
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a) > 0 ? (a) : -(a))
#endif

#ifndef SIGN
#define SIGN(a) ((a) >= 0. ? 1 : -1)
#endif

#ifndef SQUARE
#define SQUARE(a) ((a)*(a))
#endif

#ifndef CUBE
#define CUBE(a) ((a)*(a)*(a))
#endif

#ifndef COLLAR
#define COLLAR(amt,cap,flr) \
                MAX(MIN((amt),(cap)),(flr))
#endif

#ifndef ACC_FN
#define ACC_FN(rate,dcf,isSIMPLE)            \
              ( (isSIMPLE) ?                 \
                ((rate)*(dcf)) :             \
                (pow((1+(rate)),(dcf))-1.0)  \
              )
#endif

#ifndef ACC_FN_DIFF  /* equals ACC_FN(r1,dcf,type) - ACC_FN(r2,dcf,type) */
#define ACC_FN_DIFF(r1,r2,dcf,isSIMPLE)                    \
              ( (isSIMPLE) ?                               \
                (((r1)-(r2))*(dcf)) :                      \
                (pow((1+(r1)),(dcf))-pow((1+(r2)),(dcf)))  \
              )
#endif

#ifndef IS_EQUAL
#define IS_EQUAL(m,n) (fabs((m)-(n)) < TINY)
#endif

#ifndef DR_REALLOC
#define DR_REALLOC(dl,s) (((dl)==NULL) ? malloc((s)) : realloc((dl),(s)))
#endif

#ifndef IS_Q
#define IS_Q(q) (fabs((q)) > QCUTOFF)
#endif

#ifndef SIZEOFARRAY
#define SIZEOFARRAY(a)	(sizeof(a)/sizeof(a[0]))
#endif

/* positioning of first element of each char block in array */
#ifndef GTO_STR_POS
#define GTO_STR_POS(idx)                          \
            (                                     \
                (idx-1)*(GTO_MAX_STR_LEN+1)       \
                      + GTO_NUM_BYTES_FOR_STR_CNT \
            )
#endif

/* total nb bytes in char block array  */
#ifndef DR_CHAR_BLOCK_ARRAY_SIZE
#define DR_CHAR_BLOCK_ARRAY_SIZE(NbBlocks)    \
            (                                 \
                GTO_STR_POS(NbBlocks +1) -1   \
            )
#endif


#ifdef  __cplusplus
}
#endif


#endif

