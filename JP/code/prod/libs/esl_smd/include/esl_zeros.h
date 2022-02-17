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

/** 0=Linear Zero Cpn; 1=Flat Fwd */ 
extern int  ZeroInterpTypeFlag;        
/** 0=Linear Stub; 1=Flat Stub    */ 
extern int  ZeroInterpTypeFlagStub;    

/*****  ZeroPrice  ****************************************************/
/**
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      (ZeroInterpTypeFlag = 0) Linear Zero Cpn:
 *      -- linear zero cpn interp
 *      -- flat zero cpn extrapolation on both sides
 *
 *      (ZeroInterpTypeFlag = 1) Flat Fwd:
 *      -- flat fwd interp and extrapolation
 */

double   ZeroPrice(long         MatDate    /** (I) Mat date of the zero       */
                  ,long         ValueDate  /** (I) Value date                 */
                  ,int          NbZero     /** (I) Number of zeros in curve   */
                  ,long const*  ZeroDates  /** (I) maturity dates in zero crv */
                  ,double const*ZeroRates  /** (I) Zero rates                 */
		);


/*****	Get_Zero  *************************************************************/
/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *	   curve stored in t_curve style (deterministic).
 *
 */
int  Get_Zero(double   *OutZeroRate   /** (O) Zero rate	                  */
             ,double   *OutZeroPrice  /** (O) Zero price                  */
             ,int       NbZero        /** (I) Number of zeros in the ZC   */
             ,double   *Zero          /** (I) Z rates on an ACT/365 basis */
             ,long     *ZeroDate      /** (I) Zero maturity dates         */
             ,long      CurrentDate   /** (I) Current date                */
             ,long      Maturity      /** (I) Maturity of the zero        */
	     );




#ifdef  __cplusplus
}
#endif


#endif



