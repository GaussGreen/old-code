#ifndef ESL_PARYIELD_NEW_DOT_H
#define ESL_PARYIELD_NEW_DOT_H

/** NOTE: This file should be only included through 'esl_paryield.c'
 */

#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_zeros.h"
#include "esl_util.h"
#include "esl_alloc.h"


/*****  Par_Yield_Plus ******************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int Par_Yield_Plus(
                   double*      ParYield      /** (O) Par forward yield          */
                  ,double*      Annuity       /** (O) Annuity                    */
                  ,const T_CURVE  *zc   /** (I) Zero curve                 */
                  ,long         StartDate     /** (I) Forward start date         */
                  ,char         IndexType     /** (I) Rate type (cash or swap)   */
                  ,int          IndexMat      /** (I) Index maturity             */
                  ,char         DayCount      /** (I) Index Day count convention */
                  ,char         IndexF        /** (I) Index payment frequency    */
		);


/*****  Par_Yield  **********************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int    Par_Yield (double*       ParYield    /** (O) Par forward yield         */
                 ,double*       Annuity     /** (O) Annuity                   */
                 ,const T_CURVE  *zc  /** (I) Zero curve                */
                 ,long          StartDate   /** (I) Forward start date        */
                 ,int           IndexMat    /** (I) Index maturity            */
                 ,char          DayCount    /** (I) Index Day count convention*/
                 ,char          IndexF      /** (I) Index payment frequency   */
		);

/** same - but the right way */
int    ParYld (double*              ParYield   /** (O) Par forward yield         */
              ,double*              Annuity    /** (O) Annuity                   */
              ,T_CURVE const* crv        /** (I) index zero curve          */
              ,long                 StartDate  /** (I) Forward start date        */
              ,int                  IndexMat   /** (I) Index maturity            */
              ,char                 DayCount   /** (I) Index Day count convention*/
              ,char                 IndexF     /** (I) Index payment frequency   */
		);

/*****  Par_Yield_From_Dates  ************************************************/
/**
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 */

 int    Par_Yield_From_Dates
                (double  *      ParYield  /** (O) Par forward yield          */
                ,double  *      Annuity   /** (O) Annuity                    */
                ,long           SwapSt    /** (I) Underlying swap start      */
                ,long           SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,const T_CURVE  *zc /** (I) Zero curve                 */
		);

/* like above but the right way */
int    ParYieldFromDates
                (double*                ParYield /** (O) Par forward yield          */
                ,double*                Annuity  /** (O) Annuity                    */
                ,long                   SwapSt   /** (I) Underlying swap start      */
                ,long                   SwapMat  /** (I) Underlying swap maturity   */
                ,char                   DCC      /** (I) Underlying day count conv. */
                ,char                   Freq     /** (I) Underlying frequency       */
                ,char                   StubConv /** (I) F, B or N(one allowed)     */
                ,T_CURVE const*   crv      /** zero curve */
		);

 /*****  ParYieldRatio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */

int  ParYieldRatio(double*      Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,long         StartDate1 /** (I) start date of yield 1      */
                  ,long         StartDate2 /** (I) start date of yield 2      */
                  ,double       Spread     /** (I) added to both yields       */
                  ,const T_CURVE  *zc/** (I) Zero curve                 */
                  ,int          IndexMat   /** (I) Index maturity             */
                  ,char         DayCount   /** (I) Index Day count convention */
                  ,char         IndexF     /** (I) Index payment frequency    */
		);

/*****  Par_Yield_Ratio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */
int  Par_Yield_Ratio(double *       Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                    ,long           StartDate1 /** (I) start date of yield 1      */
                    ,long           StartDate2 /** (I) start date of yield 2      */
                    ,double         Spread     /** (I) added to both yields       */
                    ,const T_CURVE  *zc  /** (I) Zero curve                 */
                    ,int            IndexMat   /** (I) Index maturity             */
                    ,char           DayCount   /** (I) Index Day count convention */
                    ,char           IndexF     /** (I) Index payment frequency    */
		);





int  ParYldRatio(double*        ratio       /** (O) Yield1+sprd / Yield2+sprd  */
                ,ESL_DATE       startDate1  /** (I) start date of yield 1      */
                ,ESL_DATE       startDate2  /** (I) start date of yield 2      */
                ,double         spread      /** (I) added to both yields       */
                ,const T_CURVE *crv   /** (I) index zero curve           */
                ,unsigned int   indexMat    /** (I) Index maturity             */
                ,ESL_DCC        dcc         /** (I) Index Day count convention */
                ,ESL_FREQ       freq        /** (I) Index payment frequency    */
		);


#ifdef  __cplusplus
}
#endif


#endif



