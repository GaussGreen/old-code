#ifndef _drmacros_h
#define _drmacros_h

#ifdef  __cplusplus
extern "C" {
#endif


#define     TRUE          1
#define     FALSE         0

#define     SUCCESS       0
#define     FAILURE      -1
#define     MAXBUFF       250

#define    MAXNBCHAR    256    /* Maximum string size           */
#define    MAXNBSV      200    /* Maximum nb of state variables */

#define MAXPOINTS 500
#define TINY 1.0e-14

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


#ifndef DR_REALLOC
#define DR_REALLOC(dl,s) (((dl)==NULL) ? malloc((s)) : realloc((dl),(s)))
#endif

#ifndef IS_EQUAL
#define IS_EQUAL(m,n) (fabs((m)-(n)) < TINY)
#endif

#ifndef DR_FREE
#define DR_FREE(a,deleteFunc) if(a){deleteFunc(a);a=NULL;}
#endif

#ifndef NUM_BUSIDAYS_PER_YEAR
#define NUM_BUSIDAYS_PER_YEAR 252
#endif


#ifdef  __cplusplus
}
#endif


#endif
