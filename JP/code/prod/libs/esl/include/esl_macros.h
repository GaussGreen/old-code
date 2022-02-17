#ifndef ESL_MACROS_DOT_H
#define ESL_MACROS_DOT_H

/** NOTE: This file should be only included through 'esl.h'
 */

#ifdef ESL_NEW_DATE
#include "irx/macros.h"
#endif




#if defined (_WIN32) || defined (WIN32) || defined (__WINDOWS__)
    #include <float.h>
#else 
    #if defined(SOLARIS) 
        #include <ieeefp.h>
    #endif
#endif


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
#ifdef  ERROR
#undef  ERROR
#endif
#define  ERROR         1E-7           /* Error tolerance                */
#define  BARRIER_TOL   1E-4           /* Error tolerance for barriers   */
#define  STRIKE_TOL    0.0            /* Error tolerance for strikes    */
#define  JUMPCOEFF     3.0            /* Jump size coefficient          */
#define  NBCRITDATE    50             /* Number of critical date arrays */
#define  MAX_INST      NBCRITDATE
#define  Nb_Daily_Pts  0              /* Number of daily points         */

/** the first zerobank event type  */
#define  ZbkEVENT      (NBCRITDATE-3) 
/* TMX: numeraire event */
#define  NMREVENT         (NBCRITDATE-4) /* Numeraire date event type   */
#define  MAXPRODEVENT     (NBCRITDATE-5) /* Max product critical type   */

#define  MAXNBDATE     1000    /* Max nb of elements in input date array */
#define  MAXNBSTATES   200    /* Max nb of state variables              */
#define  MAXNBEVCURVES 5      /* Max nb of event curves in EVENT_LIST   */
#define  NBSWAPTION    25     /* Maximum number of volatilities in swaption matrix  */
#define  MAXNBTD       50     /* Maximum number of time dependent quantities */
/* #define  MAX_INPUT_ZRATES MAXNBDATE */

#define   MAX_ITERATIONS    50
#define   MAXINDEX        8   /* Maximum number of characters in index name         */

/* TMX constants */   
#define  MAXNBNMR          256  
#define  NMR_CUTOFF        1E-16      /* Low numeraire cutoff (> 0)     */
#define  MAXNBEXPIRY       30
#define  MAXNBFWDMAT       20
#define  NBVOLPARS         9          /* sig,q,vvol,bbV,bbR,dL,tl,dR,tR */
#define  NBSTRIKE          5          /* number of strikes for calib    */
#define  TMX_MAX_VOL_VOL   1.5
#define  TMX_MAX_SKEW      4

/*****************************************************/
/* Types for memory allocation.                      */
/*                                                   */
/* Warning: change the order at your own risk,       */
/* since these macros are duplicated in other header */
/* files (Q3, Q3TMX) that are included by products   */
/* that call esl functions such as DR_Array.         */
/*                                                   */
/* Warning: change to an "enum" at your own risk,    */
/* since you would then find conflicts with types    */
/* defined in MSVS/VC98.                             */
/*****************************************************/

/* Types for memory allocation */
#define     CRITDATE             0
#define     INT                  1
#define     INT_PTR              2
#define     INT_D_PTR            3
#define     LONG                 4
#define     LONG_PTR             5
#define     LONG_D_PTR           6
#define     DOUBLE               7
#define     DOUBLE_PTR           8
#define     DOUBLE_D_PTR         9
#define     CHAR                10
#define     CHAR_PTR            11
#define     INT_T_PTR           12
#define     TYPE_TPROB_0        13
#define     TYPE_TPROB_1        14
#define     TYPE_TPROB_2        15
#define     TYPE_SWAPRATE_DATA  16

/* IRDate is currently represented as a long, but allow us flexibility to change. */
#define     IDATE       LONG
#define     IDATE_PTR   LONG_PTR
#define     IDATE_D_PTR LONG_D_PTR

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

#ifndef MAXMIN
#define MAXMIN(amt,cap,flr) (                                                    \
                                (flr) < (cap) ? (                                \
                                    (amt) < (flr) ? (flr) : (                    \
                                        (amt) > (cap) ? (cap) : (amt)            \
                                                            )                    \
                                                ) : (flr)                        \
                            )
#endif

#ifndef MINMAX
#define MINMAX(amt,flr,cap) (                                                    \
                                (flr) < (cap) ? (                                \
                                    (amt) < (flr) ? (flr) : (                    \
                                        (amt) > (cap) ? (cap) : (amt)            \
                                                            )                    \
                                                ) : (cap)                        \
                            )
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

#ifndef IS_ZERO
#define IS_ZERO(m) ( IS_EQUAL((m),0.0) )
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


#ifndef DR_IS_FINITE
#if defined (_WIN32) || defined (WIN32) || defined (__WINDOWS__)
#define DR_IS_FINITE(a) (_finite((double)(a)))
#else
#define DR_IS_FINITE(a) (finite((double)(a)))
#endif
#endif


/* returns MAX(a,b) if a & b both finite, else +/- infinity, a if both are infinite */
#ifndef MAX_INF
#define MAX_INF(a,b)  ((DR_IS_FINITE(a)) ? ( DR_IS_FINITE(b) ? MAX(a,b) : b ) : a )
#endif

#ifndef MIN_INF
#define MIN_INF(a,b)  ((DR_IS_FINITE(a)) ? ( DR_IS_FINITE(b) ? MIN(a,b) : b ) : a )
#endif

#ifndef DR_FREE
#define DR_FREE(a,deleteFunc) if(a){deleteFunc(a);a=NULL;}
#endif


#define DR_CHECK(COND, FAILVALUE, ERRMSG, ROUTINE) \
    if ((COND) == (FAILVALUE)) {\
        DR_Error("%s : %s! ", (ROUTINE), (ERRMSG) ); \
        goto RETURN; \
    } \

#define DR_ASSERT(CONDITON,  ERRMSG, ROUTINE) \
    if (!(CONDITON)) {\
        DR_Error("%s : %s! ", (ROUTINE), (ERRMSG) ); \
        goto RETURN; \
    } \

#ifdef  __cplusplus
}
#endif


#endif

