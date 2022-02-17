/************************************************************************
 * File:        dllutil.h
 * Function:    Utility functions for wrapping DLL
 ************************************************************************/
/* Header Files */ 
#ifndef _dllutil_H
#define _dllutil_H

#include <stdlib.h> 
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "esl_types.h"
#include "esl_date.h"

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * Define the platform compatibility for DLL export function
 */
# if defined(_WINDLL) && (!defined(_WIN32) && !defined(WIN32))
#  define DLL_EXPORT(type)      type __export
# elif defined(_WIN32) || defined(WIN32)
#  define DLL_EXPORT(type)      __declspec(dllexport) type
# else
#  define DLL_EXPORT(type)      type
# endif


/* Task: Macros for counted array 
 */
#define ASSERT_OR_FAIL(cond)   {\
    if (!(cond)) {DR_Error("ASSERT_OR_FAIL: \
    `%s' failed.",\
    #cond); return(FAILURE);}}

#define IS_VALID_ARRAY_L(arrayL)  \
    ((arrayL) != NULL && (int)(arrayL)[0] != 0)

#define CHECK_ARRAY_LEN_L(arrayL, size)    {\
    if ((int)arrayL[0] != (int)(size)) {DR_Error("CHECK_ARRAY_LEN_L: \
    argument `%s' size (%d) is not %d.", #arrayL, (int)arrayL[0], (int)(size)); return(FAILURE);}}

/* Task: Convert counted array to pointer 
 */
#define COPY_FROM_ARRAY_L(arrayL, MAX_SIZE, array)    {\
    int i;\
    if ((int)arrayL[0] >= (MAX_SIZE)) {DR_Error("COPY_FROM_ARRAY_L: \
    %s too large for target.", #arrayL); return(FAILURE);}\
    for (i=1;i<=(int)arrayL[0];i++) {array[i-1]=arrayL[i];}}

/* Task: Convert pointer to counted array
 */
#define COPY_TO_ARRAY_L(arrayL, n, array)    {\
    int i;\
    if (n >= (MAXNBDATE)) {DR_Error("COPY_TO_ARRAY_L: \
    %s too large for target.", #array); return(FAILURE);}\
    arrayL[0]=n;\
    for (i=1;i<=n;i++) {arrayL[i]=array[i-1];}}

/* Set date conversion function -- IRDateFromYMDDate  by default */
int SetDateConv (long (*f) (long));
long DRDate (long); /* output is in internal format */

#define COPY_FROM_DATES_L(arrayL, MAX_SIZE, array)    {\
    int i;\
    if ((int)arrayL[0] >= (MAX_SIZE)) {DR_Error("COPY_FROM_ARRAY_L: \
    %s too large for target.", #arrayL); return(FAILURE);}\
    for (i=1;i<=(int)arrayL[0];i++) {array[i-1]=DRDate(arrayL[i]);}}

/* Task: Convert to counted array of dates (no conversion)
 */
#define COPY_TO_DATES_L(arrayL, n, array)    {\
    int i;\
    if (n >= (MAXNBDATE)) {DR_Error("COPY_TO_ARRAY_L: \
    %s too large for target.", #array); return(FAILURE);}\
    arrayL[0]=n;\
    for (i=1;i<=n;i++) {arrayL[i]=YMDDateFromIRDate(array[i-1]);}}

/* Task: Convert array of string to array of chars
 */
#define COPY_FROM_CHARS_L(arrayL, MAX_SIZE, array)    {\
    int i; \
    if ((int)arrayL[0] >= (MAX_SIZE)) {DR_Error("COPY_FROM_CHARS_L: \
    %s too large for target.", #arrayL); return(FAILURE);}\
    for (i=1;i<=(int)arrayL[0];i++) {array[i-1]=(char)toupper(arrayL[i][0]);}}

/* Task: Convert array of chars to counted array of strings
 */
#define COPY_TO_CHARS_L(arrayL, n, array)    {\
    int i; \
    if (n >= (MAXNBDATE)) {DR_Error("COPY_TO_CHARS_L: \
    %s too large for target.", #array); return(FAILURE);}\
    arrayL[0]=(char *)n;\
    for (i=1;i<=n;i++) {sprintf(arrayL[i],"%c",array[i-1]);}}

/* Task: Convert counted array of strings to pointer 
 */
#define COPY_FROM_STRNG_L(arrayL, MAX_SIZE, array)    {\
    int i; \
    if ((int)arrayL[0] >= (MAX_SIZE)) {DR_Error("COPY_FROM_STRNG_L: \
    %s too large for target.", #arrayL); return(FAILURE);}\
    for (i=1;i<=(int)arrayL[0];i++) {strcpy(array[i-1],arrayL[i]);}}

/* Task: Convert pointer to counted array of strings
 */
#define COPY_TO_STRNG_L(arrayL, n, array)    {\
    int i; \
    if (n >= (MAXNBDATE)) {DR_Error("COPY_TO_STRNG_L: \
    %s too large for target.", #array); return(FAILURE);}\
    arrayL[0]=(char *)n;\
    for (i=1;i<=n;i++) {strcpy(arrayL[i],array[i-1]);}}

/* Task: allocate and free counted array of strings (for regression test)
 */
char** Alloc_Strng_L (long, long);
int    Free_Strng_L (char **, long, long);


typedef enum  
{  
   L_ArrayCount_BASEDATE,    
   L_Today,  
   L_ValueDate0,  
   L_ValueDate1,  
   L_ValueDate2,  
   L_VolBaseDate,  
  
   L_BASEDATE_SIZE  
}  
   BASEDATE;  

  
typedef enum  
{  
   L_ArrayCount_SMILEPARAM,
   L_QLeft,  
   L_QRight,  
   L_FwdShift,  
   L_CetNbIter,  
  
   L_SMILEPARAM_SIZE  
}  
   SMILEPARAM;  


/* Original FIX3 parameters */
typedef enum  
{  
   L_ArrayCount_TREEPARAM,
   L_Ppy,
   L_NbSigmaMax,  
   L_CvDiff,  
   L_CvDisc,  
  
   L_TREEPARAM_SIZE  
}  
   TREEPARAM_ORIGINAL;  


typedef enum
{
   L_ArrayCount_CALIBFLAG,
   L_VolSkipFlag,
   L_VolCalibFlag,
   L_VolFilterSpotVolFlag,
   L_VolSmoothingFlag,

   L_CALIBFLAG_SIZE
}
   CALIBFLAG_ORIGINAL;


typedef enum
{
   L_ArrayCount_ZCURVE,
   L_Zero0_curve,
   L_Zero1_curve,
   L_Zero2_curve,

   L_ZCURVE_SIZE
}
   ZCURVE;


/* Unified -- include model choice */
typedef enum
{
    NUM_ArrayCount,
    NUM_Ppy,
    NUM_NbSigmaMax,
    NUM_CetNbIter,
    NUM_CvDiff,
    NUM_CvDisc,
    NUM_InterpType,
    NUM_NbSigmaMQ,
    NUM_NckMQ,

    NUMERICS_SIZE
}
    NUMERICS;


typedef enum
{
    VOL_ArrayCount,
    VOL_ModelChoice,
    VOL_SkipFlag,
    VOL_CalibFlag,
    VOL_SmoothingFlag,
    VOL_Freq,
    VOL_DCC,
    VOL_Type,
    VOL_Backbone,

    VOLCALIB_SIZE
}
    VOLCALIB;


typedef enum  
{  
   CLS_ArrayCount,
   CLS_QLeft,  
   CLS_QRight,  
   CLS_FwdShift,

   CLASSIC_SMILEPARAM_SIZE  
}  
   SMILEPARAM_CLASSIC; 


typedef enum  
{  
   SMD_ArrayCount,
   SMD_QLeft,  
   SMD_QRight,  
   SMD_FwdShift,
   SMD_Afac,
   SMD_Bfac,
   SMD_Cfac,
   SMD_Dfac,
  
   SMD_SMILEPARAM_SIZE  
}  
   SMILEPARAM_SMD; 


typedef enum  
{  
   E2Q_ArrayCount,
   E2Q_QLeft,  
   E2Q_QRight,  
   E2Q_FwdShift,
   E2Q_Amap,
   E2Q_Bmap,

   E2Q_SMILEPARAM_SIZE  
}  
   SMILEPARAM_E2Q; 
/* End of unified parameters -- include model choice */


/* Pack Zero-Curves */  

/* Task  : Packs up a single T_CURVE structure.  
*/  
int PackSingleTCurve( T_CURVE*      crv,  
                      IRDate        Today,
                      IRDate        ValueDate,   
                      char const*   SwapFreq, 
                      char const*   DCC, 
                      char const*   MMB, 
                      long          SpotDays,
                      IRDate const* ZeroDatesL,
                      double const* ZeroRatesL);

/* Task : Packs up all three zero curves.  
*/  
int PackTCurves(  T_CURVE    *Curves,  
                  IRDate     *BaseDatesL,
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  IRDate     *Zero0DatesL,
                  double     *Zero0RatesL,
                  IRDate     *Zero1DatesL,
                  double     *Zero1RatesL,
                  IRDate     *Zero2DatesL,
                  double     *Zero2RatesL );
  
/* Unpack Zero-Curves */  

/* Task  : Unpacks up a single T_CURVE structure.  
*/  
int UnPackSingleTCurve(T_CURVE const* crv,  
                      IRDate*         ValueDate,   
                      char*           SwapFreq, 
                      char*           DCC, 
                      char*           MMB, 
                      long*           SpotDays,
                      IRDate*         ZeroDatesL,
                      double*         ZeroRatesL );

/* Task : Unpacks up all three zero curves.
*/
int UnPackTCurves(
                  T_CURVE const* Curves,
                  IRDate      *BaseDatesL,
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  IRDate     *Zero0DatesL,
                  double     *Zero0RatesL,
                  IRDate     *Zero1DatesL,
                  double     *Zero1RatesL,
                  IRDate     *Zero2DatesL,
                  double     *Zero2RatesL );

#ifdef  __cplusplus
}
#endif


#endif  /* _dllutil_H */
