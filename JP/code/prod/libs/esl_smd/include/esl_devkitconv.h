#ifndef ESL_DEVKITCONV_DOT_H
#define ESL_DEVKITCONV_DOT_H

/** NOTE: This file should be only included through 'esl_devkitconv.c'
 */


#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"


/* defining locally used macros */

#ifndef GTO_MONTHS_PER_YEAR
#define GTO_MONTHS_PER_YEAR     12
#endif


#ifndef GTO_TDATE_BASE_YEAR
#define GTO_TDATE_BASE_YEAR 1601
#endif

#ifndef DAYS_IN_1_YEAR
#define DAYS_IN_1_YEAR      365L
#endif

#ifndef DAYS_IN_4_YEARS
#define DAYS_IN_4_YEARS     1461L
#endif


#ifndef DAYS_IN_100_YEARS
#define DAYS_IN_100_YEARS  (DAYS_IN_4_YEARS * 25 - 1)
#endif

#ifndef DAYS_IN_400_YEARS
#define DAYS_IN_400_YEARS  (DAYS_IN_100_YEARS * 4 + 1)
#endif


#ifndef MULTIPLY_BY_1461
#define MULTIPLY_BY_1461(N) (((N) << 10) + \
                 ((N) << 8) + \
                 ((N) << 7) + \
                 ((N) << 5) + \
                 ((N) << 4) + \
                 ((N) << 2) + \
                 (N))
#endif

/*****  TDate2DRDate  ***************************************************/
/**
        Converts an alib TDate to a DR Date (YYYYMMDD)
 */
int TDate2DRDate(TDATE    date    /** (I) Alib date format */ 
                ,long    *DrDate  /** (O) YYYYMMDD         */
		);

/*****  DRDate2TDate  ***************************************************/
/**
        Converts a DR Date (YYYYMMDD) to an alib TDate
 */
int  DRDate2TDate(long     DrDate  /** (I) YYYYMMDD  */
                 ,TDATE   *odate   /** (O) Alib date */
		);

/*****  IsPtrNull  *******************************************************/
/**
        Returns: TRUE if ptr is NULL or if both 1) type is valid and 2) 
                 ptr points to the array {1,0} where both elements are 
                 represented in their corresponding type;
                 FALSE otherwise.
        WARNING: if pointer being tested intentionally holds {1,0}, or type
                 is invalid, then the function will return TRUE
*/
int    IsPtrNull (DEVKIT_TYPE   type    /** (I) Type    */
                 ,void         *ptr     /** (I) Pointer */
		);


#ifdef  __cplusplus
}
#endif


#endif



