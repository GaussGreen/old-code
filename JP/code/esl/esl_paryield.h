#ifndef ESL_PARYIELD_DOT_H
#define ESL_PARYIELD_DOT_H

#include "esl_types.h"
#include "esl_date.h"
#include "esl_util.h"
#include "esl_alloc.h"

#ifdef  __cplusplus
extern "C" {
#endif

#ifndef ESL_NEW_CURVE
/* DEPRECATED OLD NON-IRX COMPATIBLE CODE.                                   */

/* DEPRECATED!!! Use ParYieldPlus instead.                                  */
int Par_Yield_Plus(
                   double*      ParYield      /** (O) Par forward yield          */
                  ,double*      Annuity       /** (O) Annuity                    */
                  ,int          NbZero        /** (I) Number of zeros            */
                  ,double const*Zero          /** (I) Zero rates                 */
                  ,long   const*ZeroDate      /** (I) Zero maturity dates        */
                  ,long         CurrentDate   /** (I) Current date               */
                  ,long         StartDate     /** (I) Forward start date         */
                  ,char         IndexType     /** (I) Rate type (cash or swap)   */
                  ,int          IndexMat      /** (I) Index maturity             */
                  ,char         DayCount      /** (I) Index Day count convention */
                  ,char         IndexF        /** (I) Index payment frequency    */
		);

/* DEPRECATED!!! Use ParYield instead.                                      */
int    Par_Yield (double*       ParYield   /** (O) Par forward yield         */
                 ,double*       Annuity    /** (O) Annuity                   */
                 ,int           NbZero      /** (I) Number of zeros           */
                 ,double const* Zero       /** (I) Zero rates                */
                 ,long   const* ZeroDate   /** (I) Zero maturity dates       */
                 ,long          CurrentDate /** (I) Current date              */
                 ,long          StartDate   /** (I) Forward start date        */
                 ,int           IndexMat    /** (I) Index maturity            */
                 ,char          DayCount    /** (I) Index Day count convention*/
                 ,char          IndexF      /** (I) Index payment frequency   */
		);

/* DEPRECATED!!! Use ParYieldRatio instead.                                 */
int  Par_Yield_Ratio(double*        Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,long           StartDate1 /** (I) start date of yield 1      */
                  ,long           StartDate2 /** (I) start date of yield 2      */
                  ,double         Spread     /** (I) added to both yields       */
                  ,int            NbZero     /** (I) Number of zeros            */
                  ,long    const* ZeroDates /** (I) Zero maturity dates        */
                  ,double  const* ZeroRates /** (I) Zero rates                 */
                  ,long           ValueDate  /** (I) Value date of zero curve   */
                  ,int            IndexMat   /** (I) Index maturity             */
                  ,char           DayCount   /** (I) Index Day count convention */
                  ,char           IndexF     /** (I) Index payment frequency    */
		);

/* DEPRECATED - Use ParYieldFromDates (above) instead.                       */
int    Par_Yield_From_Dates
                (double*        ParYield /** (O) Par forward yield           */
                ,double*        Annuity  /** (O) Annuity                     */
                ,long           SwapSt    /** (I) Underlying swap start      */
                ,long           SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,int            NbZero    /** (I) Number of zeros            */
                ,double  const* Zero     /** (I) Zero rates                  */
                ,long    const* ZeroDate /** (I) Zero maturity dates         */
                ,long           BaseDate  /** (I) Zero curve base date       */
		);

/*****  Swap_Yield_From_Dates  **********************************************/
/*
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 *       Uses the formula sum(FwdLibor*Z)/A, which is the correct formula to
 *       use in the presence of ccy basis.
 *       Please refer to Par_Yield_From_Dates for comparison
 */

int  Swap_Yield_From_Dates
             (double     *SwapYield,    /* (O) Par forward yield            */
              double     *Annuity,      /* (O) Annuity                      */
              long        SwapSt,       /* (I) Underlying swap start        */
              long        SwapMat,      /* (I) Underlying swap maturity     */
              char        DCC,          /* (I) Underlying day count conv.   */
              char        Freq,         /* (I) Underlying frequency         */
              char        StubConv,     /* (I) F, B or N(one allowed)       */
              const T_CURVE    *IdxZCurve, /* (I) index curve               */
              const T_CURVE    *DiscZCurve /* (I) discount curve            */
              );

#endif /* ifndef ESL_NEW_CURVE */


/*****  ParYieldPlus ******************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*       Supports cash or swap rates.
*/
int ParYieldPlus(double*        ParYield    /** (O) Par forward yield          */
                ,double*        Annuity     /** (O) Annuity                    */
                ,const T_CURVE *zc          /** (I) Zero curve                 */
                ,IRDate          StartDate   /** (I) Forward start date         */
                ,char           IndexType   /** (I) Rate type (cash or swap)   */
                ,int            IndexMat    /** (I) Index maturity             */
                ,char           DayCount    /** (I) Index Day count convention */
                ,char           IndexF      /** (I) Index payment frequency    */
		        );

/*****  ParYield  **********************************************************/
/**
*       Calculate par forward yield for a specific maturity from a zero 
*       coupon curve (deterministic). Does not allow for stub.
*/
int    ParYield (double*        ParYield    /** (O) Par forward yield         */
                ,double*        Annuity     /** (O) Annuity                   */
                ,const T_CURVE *zc          /** (I) Zero curve                */
                ,IRDate          StartDate   /** (I) Forward start date        */
                ,int            IndexMat    /** (I) Index maturity            */
                ,char           DayCount    /** (I) Index Day count convention*/
                ,char           IndexF      /** (I) Index payment frequency   */
		        );

/*****  ParYieldFromDates  ************************************************/
/**
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 */
 int    ParYieldFromDates
                (double*        ParYield  /** (O) Par forward yield          */
                ,double*        Annuity   /** (O) Annuity                    */
                ,IRDate          SwapSt    /** (I) Underlying swap start      */
                ,IRDate          SwapMat   /** (I) Underlying swap maturity   */
                ,char           DCC       /** (I) Underlying day count conv. */
                ,char           Freq      /** (I) Underlying frequency       */
                ,char           StubConv  /** (I) F, B or N(one allowed)     */
                ,const T_CURVE *zc        /** (I) Zero curve                 */
		        );

/*****  Swap_Yield_From_Dates  **********************************************/
/*
 *      Calculate par forward yield for a specific maturity from a zero 
 *      coupon curve (deterministic). Underlying swap is determined  by
 *      a swap start date and a swap end date, therefore allowing stubs.
 *
 *       Uses the formula sum(FwdLibor*Z)/A, which is the correct formula to
 *       use in the presence of ccy basis.
 *       Please refer to Par_Yield_From_Dates for comparison
 */
 int  SwapYieldFromDates
             (double     *SwapYield,    /* (O) Par forward yield            */
              double     *Annuity,      /* (O) Annuity                      */
              IRDate        SwapSt,       /* (I) Underlying swap start        */
              IRDate        SwapMat,      /* (I) Underlying swap maturity     */
              char        DCC,          /* (I) Underlying day count conv.   */
              char        Freq,         /* (I) Underlying frequency         */
              char        StubConv,     /* (I) F, B or N(one allowed)       */
              const T_CURVE    *IdxZCurve, /* (I) index curve               */
              const T_CURVE    *DiscZCurve /* (I) discount curve            */
              );



 /*****  ParYieldRatio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic par
 *       fwd yields
 *       
 */
int  ParYieldRatio(double*          Ratio      /** (O) Yield1+sprd / Yield2+sprd  */
                  ,IRDate            StartDate1 /** (I) start date of yield 1      */
                  ,IRDate            StartDate2 /** (I) start date of yield 2      */
                  ,double           Spread     /** (I) added to both yields       */
                  ,const T_CURVE   *zc         /** (I) Zero curve                 */
                  ,int              IndexMat   /** (I) Index maturity             */
                  ,char             DayCount   /** (I) Index Day count convention */
                  ,char             IndexF     /** (I) Index payment frequency    */
		          );



#ifdef  __cplusplus
}
#endif

#endif /* ESL_PARYIELD_DOT_H */
