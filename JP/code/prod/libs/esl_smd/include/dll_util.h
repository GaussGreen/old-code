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
#include <esl_types.h>


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
    if ((int)arrayL[0] != (size)) {DR_Error("CHECK_ARRAY_LEN_L: \
    argument `%s' size (%d) is not %d.", #arrayL, (int)arrayL[0], size); return(FAILURE);}}

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

#define COPY_FROM_DATES_L(arrayL, MAX_SIZE, array)    {\
    int i;\
    if ((int)arrayL[0] >= (MAX_SIZE)) {DR_Error("COPY_FROM_ARRAY_L: \
    %s too large for target.", #arrayL); return(FAILURE);}\
    for (i=1;i<=(int)arrayL[0];i++) {array[i-1]=LDate2ADate(arrayL[i]);}}

/* Task: Convert to counted array of dates (no conversion)
 */
#define COPY_TO_DATES_L(arrayL, n, array)    {\
    int i;\
    if (n >= (MAXNBDATE)) {DR_Error("COPY_TO_ARRAY_L: \
    %s too large for target.", #array); return(FAILURE);}\
    arrayL[0]=n;\
    for (i=1;i<=n;i++) {arrayL[i]=ADate2LDate(array[i-1]);}}

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


typedef enum  
{  
   L_ArrayCount_TREEPARAM,
   L_Ppy,
   L_NbSigmaMax,  
   L_CvDiff,  
   L_CvDisc,  
  
   L_TREEPARAM_SIZE  
}  
   TREEPARAM;  

typedef enum
{
   L_ArrayCount_CALIBFLAG,
   L_VolSkipFlag,
   L_VolCalibFlag,
   L_VolFilterSpotVolFlag,
   L_VolSmoothingFlag,

   L_CALIBFLAG_SIZE
}
   CALIBFLAG;



typedef enum
{
   L_ArrayCount_ZCURVE,
   L_Zero0_curve,
   L_Zero1_curve,
   L_Zero2_curve,

   L_ZCURVE_SIZE
}
   ZCURVE;


/* Pack Zero-Curves */  

/* Task  : Packs up a single T_CURVE structure.  
*/  
int PackSingleTCurve( T_CURVE*      crv,  
                      long          Today,
                      long          ValueDate,   
                      char const*   SwapFreq, 
                      char const*   DCC, 
                      char const*   MMB, 
                      long          SpotDays,
                      long const*   ZeroDatesL,
                      double const* ZeroRatesL);

/* Task : Packs up all three zero curves.  
*/  
int PackTCurves(  T_CURVE    *Curves,  
                  long       *BaseDatesL,
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  long       *Zero0DatesL,
                  double     *Zero0RatesL,
                  long       *Zero1DatesL,
                  double     *Zero1RatesL,
                  long       *Zero2DatesL,
                  double     *Zero2RatesL );
  
/* Unpack Zero-Curves */  

/* Task  : Unpacks up a single T_CURVE structure.  
*/  
int UnPackSingleTCurve(T_CURVE const* crv,  
                      long*           ValueDate,   
                      char*           SwapFreq, 
                      char*           DCC, 
                      char*           MMB, 
                      long*           SpotDays,
                      long*           ZeroDatesL,
                      double*         ZeroRatesL );

/* Task : Unpacks up all three zero curves.
*/
int UnPackTCurves(
                  T_CURVE const* Curves,
                  long       *BaseDatesL,
                  char       **ZeroFrequenciesL,
                  char       **ZeroDCCsL,
                  char       **ZeroMMBsL,
                  long       *ZeroSpotDaysL,
                  long       *Zero0DatesL,
                  double     *Zero0RatesL,
                  long       *Zero1DatesL,
                  double     *Zero1RatesL,
                  long       *Zero2DatesL,
                  double     *Zero2RatesL );


#endif  /* _dllutil_H */
