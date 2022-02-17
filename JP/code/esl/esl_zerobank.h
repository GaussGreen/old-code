#ifndef ESL_ZEROBANK_DOT_H
#define ESL_ZEROBANK_DOT_H

/** NOTE: This file should be only included through 'esl_zerobank.c'
 */

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"



#ifdef  __cplusplus
extern "C" {
#endif

/*****  ZbkDLFromIdx  *****************************************************/
/**
 *      Given a list of reset dates and descriptions of an index
 *      generate a list of zero maturity dates for index estimation and their
 *      corresponding usage dates
 *      Returns SUCCESS or FAILURE
 */

int     ZbkDLFromIdx
            (int        NbResetDates    /** (I) size of reset date list  */
            ,long      *ResetDL         /** (I) reset date list          */
            ,long      *SwapStDL        /** (I) SwapRate start date list */
            ,int        IndexMat        /** (I) index tenor in months    */
            ,char       IndexFreq       /** (I) index freq               */
            ,int       *NbMatDates      /** (O) Nb of zero mat dates     */
            ,long     **MatDL           /** (O) zero mats                */
            ,int       *NbUseDates      /** (O) Nb zero use dates        */
            ,long     **UseDL           /** (O) zero use dates           */
	    );

/*****  ZbkBnd  *************************************************************/
/**
 *      Given the start date and the bounds date (the date the Bnds is
 *      calculated), returns a new date satisfying the following:
 *      if mode = +1:
 *          new date - start date = Bnd(bounds date)
 *      if mode = -1:
 *          start date - new date = Bnd(bounds date)
 */

long    ZbkBnd(long  Sdate      /** (I) date the offset is calculated from */
              ,long  Bdate      /** (I) date used to determine the bound   */
              ,long  ValueDate  /** (I) value date                         */
              ,int   mode       /** (I) +1 = Fwd; -1 = Backward            */
	      );

/*****  ZbkG  ***********************************************************/
/**
 *      Performs the first pass of the algorithm. This routine attempts
 *      to span a set of target maturities with the fewest number of 
 *      zero-bank maturity dates.
 *      Returns the zero-bank mats (in descending order) and their
 *      corresponding earliest usage dates, together with the index of 
 *      the smallest processed Mat.
 *      Returns -999 on FAILURE
 *
 *      Only called by ZbkOptDates, do NOT call directly
 */

int     ZbkG(int          NbMat      /** (I) size of target mats          */
            ,long        *Mat        /** (I) target zero mats             */
            ,long        *LaMat      /** (I) latest usage for zero mats   */
            ,long        *ErMat      /** (I) earliest usage for zero mats */
            ,long         ValueDate  /** (I) value date                   */
            ,int         *OutNbZ     /** (O) size of OutZ and ErOutZ      */
            ,long       **OutZ       /** (O) zerobank mat array           */
            ,long       **OutErZ     /** (O) zerobank earliest usage list */
            ,long       **OutFstMat  /** (O) 1st mat in each zero intval  */
	    );

/*****  ZbkH  ***********************************************************/
/**
 *      Performs the second pass of the algorithm. This routine attempts
 *      to shift the zero-bank mats from the first pass to the right, so 
 *      that they may coincide with the critical dates
 *      Returns the final zero-bank mats (in ascending order) and their
 *      corresponding earliest usage dates.
 *
 *      Only called by ZbkOptDates, do NOT call directly
 */

int     ZbkH(int       NbZ          /** (I) size of 1st pass zeros          */
            ,long     *Z            /** (I) 1st pass zero mats (desc order) */
            ,long     *ErZ          /** (I) earliest usage for zero mats    */
            ,long     *FirstMat     /** (I) 1st mat date in [Zi, Zi+1)      */
            ,int       NbC          /** (I) size of critical dates          */
            ,long     *C            /** (I) critical dates                  */
            ,long      ValueDate    /** (I) value date                      */
            ,int      *NbMatDates   /** (O) Nb of output zerobank mats      */
            ,long    **ZbkMats      /** (O) output zerobank mats            */
            ,int      *NbErDates    /** (O) Nb of output earliest usage     */
            ,long    **ZbkErs       /** (O) zerobank earliest usage list    */
	     );

/*****  ZbkProcessDL  **********************************************/
/**
 *      This routine performs the following to a raw list of zero mat
 *      dates and the corresponding usage dates:
 *      - checks that InpUseDL[i] is strictly < InpMatDL[i]
 *      - sort the dates in MatDates ascending order
 *      - remove duplicate MatDates 
 *      - set the earliest usage date as the minimum usage dates 
 *        of all duplicates
 *      - set the latest usage date as the max usage dates of all
 *        duplicates
 *
 *      Returns SUCCESS or FAILURE
 */

int     ZbkProcessDL(int    NbInpMat    /** (I) Nb of input mats          */
                    ,long  *InpMatDL    /** (I) zero mat datelist         */
                    ,long  *InpUseDL    /** (I) zero usage list           */
                    ,int   *NbOutMat    /** (O) Nb of processed zero mats */
                    ,long **Mat         /** (O) processed zero mats       */
                    ,long **LaMat       /** (O) latest usage date list    */ 
                    ,long **ErMat       /** (O) earliest usage date list  */ 
		);

/*****  ZbkOrdCritDates  ****************************************************/
/**
 *      Given a list of critical dates in the CRIT_DATE struct form,
 *      returns a list of unique critical dates sorted in ascending order
 */

int     ZbkOrdCritDates
            (int         NbCritDates    /** (I) size of CritDate struct     */
            ,CRIT_DATE  *CritDate       /** (I) CritDate struct             */
            ,int        *NbCritDL       /** (O) Nb of ordered crit dates    */
            ,long      **CritDL         /** (O) critical date list          */
	    );


/*****  ZbkOptDates  ********************************************************/
/**
 *      Given a set of zero maturities/usage dates that we want to eval,
 *      find an "optimal" set of zero-bank maturities according to the
 *      algorithm. The algorithm minimises the number of zero-bank dates
 *      subject to the Bnd(.) function, and attempts to place the final date
 *      on the critical dates
 */

int     ZbkOptDates
            (int        NbMatDates       /** (I) Nb of target mat dates      */
            ,long      *MatDL            /** (I) target mat date list        */
            ,int        NbUseDates       /** (I) Nb of target use dates      */
            ,long      *UseDL            /** (I) target use date list        */
            ,int        NbCrit           /** (I) Nb of critical dates        */
            ,CRIT_DATE *Crit             /** (I) critical date structs       */
            ,long       ValueDate        /** (I) value date                  */
            ,int       *NbZbkMats        /** (O) Nb of zerobank mat dates    */
            ,long     **ZbkMats          /** (O) zerobank mats               */
            ,int       *NbZbkErs         /** (O) Nb zbank earliest use dates */
            ,long     **ZbkErs           /** (O) zbank earliest use dates    */
	    );

/*****  CbkProcessDL  **********************************************/
/**
        This routine performs the following to a raw list of payoff known 
        dates (EvalDates) and the corresponding earliestUsage dates:
        - sort the dates by EvalDates ascending order
        - remove duplicate EvalDates and replace the earliest usage date by
          the minimum earliest usage dates of all duplicates
        - return the two lists and their size through the same arguments
        Returns SUCCESS or FAILURE
 */


int     CbkProcessDL
            (int       *NbEvDates,
             long     **EvDL,
             int       *NbErDates,
             long     **ErDL);

  




#ifdef  __cplusplus
}
#endif


#endif



