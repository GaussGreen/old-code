#ifndef ESL_ZEROS_DOT_H
#define ESL_ZEROS_DOT_H


/** NOTE: This file should be only included through 'esl_zeros.c'
 */


#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_util.h"


/*****  ZeroPrice  ****************************************************/
/**
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      -- smooth-forward interp
 *      -- flat-forward interp
 */

double   ZeroPrice(long                 MatDate  /** (I) Mat date of the zero */
                  ,const T_CURVE  *zc      /** (I) Zero curve           */
		);


/*****	Get_Zero  *************************************************************/
/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *	   curve stored in t_curve style (deterministic).
 *
 */
int  Get_Zero(double   *OutZeroRate   /** (O) Zero rate	                  */
             ,double   *OutZeroPrice  /** (O) Zero price                  */
             ,const T_CURVE  *zc/** (I) Zero curve                  */
             ,long      Maturity      /** (I) Maturity of the zero        */
	     );




#ifdef  __cplusplus
}
#endif


#endif

