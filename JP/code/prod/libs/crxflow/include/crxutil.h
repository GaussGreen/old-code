/*********************************************************************************
 * CRXUTIL.H 
 * crx misc utils
 *
 ********************************************************************************/

#ifndef __CRXUTIL_H__
#define __CRXUTIL_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <alib/gtobf.h>
    
#include "t_curve.h"
#include "l_date.h"
#include <common/include/drmacros.h>
#include "crxerror.h"
    
#define DR_FREE(a,deleteFunc) if(a){deleteFunc(a);a=NULL;}



    
/********************************************************************************
 * Date routines : make TDate.
 *
 * Creates and returns a date corresponding to
 * month "m", day "d" and year "y".
 */
TDate CrxTDateMake(int m, int d, int y);



/********************************************************************************
 * Date routines : convert from YYYYMMDD long format.
 *                              
 * Converts a DR date type (encoded as YYYYMMDD in a long)
 * to a TDate. Returns -1 if failed.
 ********************************************************************************/
TDate CrxDrDate2TDate(long drDate);

/********************************************************************************
 * Date routines : convert to YYYYMMDD long format.
 *                              
 * Converts a date in TDate type to DR date (encoded as YYYYMMDD in a long)
 * Returns -1 if failed.
 ********************************************************************************/
long CrxTDate2DrDate(TDate tDate);

/********************************************************************************
 * Convert Freq (char) to Alib Basis (double)
 ********************************************************************************/
char   CrxBasis2Freq(double Basis);

/********************************************************************************
 * Convert Dr T_CURVE to Alib TCurve
 ********************************************************************************/
int CrxDrCurve2TCurve(TCurve  **AlibCurve,
                      T_CURVE * DrCurve);




/**
 * Linear interp inside the range, and flat outside
 */
int CrxTDateLinearInterp1d(
    TDate *xa,  /* (I) array of X(i) (i=0,..,m-1) */
    double *ya, /* (I) array of Y(i) (i=0,..,m-1) */
    int m,      /* (I) arrays size */
    TDate x,    /* (I) point to intepolate */
    double *y); /* (O) interpolated value */
        


/**----------------------------------------------------
 * Merge and sort two date lists
 */
int CrxSortMergeDateList(
    int   *NbMerge,        /* (O) Num of points in merged list   */
    long  *DLMerge,        /* (O) Merged & sorted date list      */
    int   Nb1,             /* (I) Num of points in 1st list */
    long  *DL1,            /* (I) 1st date list */
    int   Nb2,             /* (I) Num of points in 2nd list */
    long  *DL2);            /* (I) 2nd date list */



/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
