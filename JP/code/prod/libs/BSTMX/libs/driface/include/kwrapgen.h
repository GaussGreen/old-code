/*------------------------------------------------------------------
HEADER FILE:    kwrapgen.h

CREATED BY:     Leda Braga - April 1996
                David Fung
		Julia Chislenko

PURPOSE:        Data conversion functions to be used in the 
                Kapital DR wrapper in connection with  ALIB
                based executables.

$Header$
---------------------------------------------------------------------- */
#ifndef _KWRAPGEN_H
#define _KWRAPGEN_H


#include <stdio.h>
#include "bastypes.h"
/*#include "cfileio.h"*/

#define  BUFFER_SIZE   250



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetZC

CREATED BY:  Leda Braga -  April 1996

DESCRIPTION: Kapital wrapper function to read the zero.dat file (in DR
             format) into a TCurve zero curve
----------------------------------------------------------------------------*/
extern TCurve *  EXPORT DrKapWrapGetZC(char  *filename);      



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetVC

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper function to read the basevol.dat file (in 
             DR format) into a Library vol curve (i.e. TCurve format)
             
             Since a value date is not given in the DR format basevol.
             dat, this must be provided externally.   
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
extern TCurve *  EXPORT DrKapWrapGetVC
                     (char  *filename,      /*(I) Usually "basevol.dat"     */
                      TDate  vcValueDate);  /*(I) Not given in the file!    */
        

  
                
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetSwpvolGrid

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper function to read the swapvol.dat file(in 
             DR format) into an array of doubles.

             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/    

extern int EXPORT DrKapWrapGetSwpvolGrid
               (char           *filename,   /*(I) Usually "swapvol.dat"     */
                int            *numSwpvolExp,
                int            *numSwpvolMat,
                double        **swpvolExp,
                double        **swpvolMat,  
                double        **swpvolGrid);




/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetOneFacPars

CREATED BY:  Leda Braga -  October 1996

DESCRIPTION: Kapital wrapper function to read the interest rate one fac-
             tor model parameters.

             Restrictions are applied concerning the format of the data
             file.

             Returns: SUCCESS or FAILURE
----------------------------------------------------------------------------*/
int  EXPORT DrKapWrapGetOneFacPars
                (char     *filename,    /* (I) Usually "meanreversion.dat"  */
                 double   *meanRev1);   /* (O) Mean rev of 1st factor       */
                 




/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetTwoFacPars

CREATED BY:  Leda Braga -  October 1996

DESCRIPTION: Kapital wrapper function to read the interest rate two fac-
             tor model parameters.

             Restrictions are applied concerning the format of the data
             file.

             Returns: SUCCESS or FAILURE
----------------------------------------------------------------------------*/
int  EXPORT DrKapWrapGetTwoFacPars
                (char     *filename,    /* (I) "2factormeanreversion.dat"   */
                 double   *meanRev1,    /* (O) Mean rev of 1st factor       */
                 double   *meanRev2,    /* (O) Mean rev of 2nd factor       */
                 double   *correl,      /* (O) Correlation between factors  */
                 double   *weight);     /* (O) Relative weight of factors   */


/*----------------------------------------------------------------------------
FUNCTION:    DrDateToTDate

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper convenience function to convert a DR format
             date to a TDate.
----------------------------------------------------------------------------*/    
extern int EXPORT DrDateToTDate(long   DRDate,
                                TDate *libDate);      
   



/*----------------------------------------------------------------------------
FUNCTION:    DrFindAndSkipComLine

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper convenience function which finds and skips 
             a comment line  and flags an error if line is not found or 
             is not a comment.
 
                            FILE MUST BE OPEN!
----------------------------------------------------------------------------*/                                 
extern int EXPORT DrFindAndSkipComLine(FILE *fp,
                                        char  *filename);  
 



/*----------------------------------------------------------------------------
FUNCTION:    DrDRLabelToTDateInterval

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience to convert a DR date interval to
             a TDateInterval. It allows for "2Y" and "2A" as valid strings

----------------------------------------------------------------------------*/     
extern int EXPORT DrDRLabelToTDateInterval
                (char             *matLabel,
                 TDateInterval    *matInterval);   


                                                                                  
                                                                                  
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCheckTCurve

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience check a TCurve for dates in
             ascending order and for non-negative rates (either zeros
             or vols). 

----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapCheckTCurve(TCurve     *curve,
                                        char       *label);


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetRateConv

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to get standard rate
             conventions from a zero curve env file

----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapGetRateConv(
    char *fileName,     /* (I)  Kapital zero curve env file */
    long *stdMMDcc,     /* (O)  standard MM Rate Convention */
    long *stdSWDcc,     /* (O)  standard Swap Rate Convention */
    long *stdSWFreq);   /* (O)  standard Swap Rate Freq */    


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapDCCtoString

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to convert a clib day count
             convention to a string

----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapDCCtoString(
    long clibDCC,       /* (I) CLib Daycount conv identifier */
    char *DCCstring);   /* (O) string output                 */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapConvDCC

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to convert DR wrapper
             day count convention to clib convention

----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapConvDCC(
    char  drDCC,                /* (I) Daycount conv in DR wrapper format */
    long *clibDCC);             /* (I) CLib Daycount conv identifier      */



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapDatesInList

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to check if a list(2) of 
             dates is a subset of list (1)
             Both lists must be sorted in ascending order on entry 
             with no repeats
             list (1) = fullList
             list (2) = subList
----------------------------------------------------------------------------*/
extern int EXPORT DrKapWrapDatesInList(
    long   nbFullList,          /* (I) nb dates in full list */
    TDate *fullList,
    long   nbSubList,           /* (I) nb dates in sub list */
    TDate *subList);


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapExpandStepUpSchedule

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to calculate the current 
             rate level at each date in the target date list given the compact
             date/rate step-up schedule in a DR Wrapper format.


----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapExpandStepUpSchedule(
    long     nbSchDates,          /* (I) nb of dates in the schedule */
    TDate   *schDates,            /* (I) schedule dates */
    double  *schRates,            /* (I) schedule rate levels */
    TBoolean isOnOrAfter,         /* (I) T => apply rate on/after sch date 
                                         F => apply only after sch date    */
    long     nbTargetDates,       /* (I) nb of target dates */
    TDate   *targetDates,         /* (I) expanded list of dates */
    double **targetRates);        /* (O) expanded list of rates */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCalcAmortPrns

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to calculate the current
             notional at each date in a payment date list given the info
             from a Dr Wrapper Amortisation schedule

----------------------------------------------------------------------------*/    
extern int EXPORT DrKapWrapCalcAmortPrns(
    double   origNotional,    /* (I) original notional amount */
    double   origNotionalPct, /* (I) original notional %      */
    long     nbAmortDates,    /* (I) nb amort dates in the schedule */
    TDate   *amortDates,      /* (I) amort dates in the schedule */
    double  *amortPcts,       /* (I) amort % in the schedule */
    long     nbTargetDates,   /* (I) nb of target dates */
    TDate   *targetDates,     /* (I) target list of dates  */
    double **currNotionals);  /* (O) notional for each date in target list */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapParamInput

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to read in model parameters
             from either the 4 string lines in the DR wrapper text, or from
             the environment model parameters file.

----------------------------------------------------------------------------*/    
extern int DrKapWrapParamInput(                                                   
    long    nbFactor,                        /* (I) Number of factors       */
    char    overWriteString[3][BUFFER_SIZE], /* (I) string (alpha,beta,rho) */
    char   *fileName,                        /* (I) File name inc extension */
    double *alphas,                          /* (O) size 3 alloc'd on entry */
    double *betas,                           /* (O) size 3 alloc'd on entry */
    double *corrs);                          /* (O) size 3 alloc'd on entry */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapParamCheck

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to verify the model 
             parameters             

----------------------------------------------------------------------------*/    
extern int DrKapWrapParamCheck(
    long      nbFactor,     /* (I) Number of factors     */
    double   *alphas,       /* (I) alphas                */
    double   *betas,        /* (I) mean reversions       */
    double   *corrs);       /* (I) corrleations          */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCalibSpotVols

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to calculate spot vols.
             The routine takes a DR Vol calib label (e.g. 10yFix), the model 
             parameters, and returns a spot vol array.
----------------------------------------------------------------------------*/

extern int DrKapWrapCalibSpotVols(
    char     *DRVolLabel,   /* (I) e.g. 3m, 10yCms, 5yFix */
    TCurve   *zCurve,       /* (I) zero curve             */
    long      stdSWFreq,    /* (I) swap rate freq         */
    long      nbFactor,     /* (I) Number of factors      */
    double    cevPower,     /* (I) 0=Normal, 1=Lognormal  */
    double   *alphas,       /* (I) alphas                 */
    double   *betas,        /* (I) mean reversions        */
    double   *corrs,        /* (I) corrleations           */
    TDate    *volBaseDate,  /* (O) vol curve base date    */
    long     *nbSpotVols,   /* (O) nb of vol dates/rates  */
    TDate   **oVolDates,    /* (O) vol curve dates        */
    double  **oSpotVols);   /* (O) spot vols              */   


/*----------------------------------------------------------------------------
FUNCTION:    DrExtendTCurve

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to extend a curve (e.g. zero
             curve) to include dates implied by BaseDate and intervals
             NOTE: baseDate must be >= oriCurve->fBaseDate
----------------------------------------------------------------------------*/

extern int DrExtendTCurve(
    TCurve         *oriCurve,        /* (I) original curve        */
    TDate           baseDate,        /* (I) start date for intvls */
    long            numIntervals,    /* (I) nb intervals in list  */
    TDateInterval  *intervals,       /* (I) interval list         */
    TCurve        **newCurve);       /* (O) result curve          */

                                                                               
/*----------------------------------------------------------------------
  FUNCTION:       DRLWrapTDateIntvlToNum
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file maturity as a TDateInterval
                  and convert it into number of months.

*/
extern int DRLWrapTDateIntvlToNum(FILE *fp,
				  char  *fileName,
				  char *paramName,  /* for error message */
				  char intvlType,   /* 'M' or 'D' */
				  long *numMonths); /* (O) */

/*----------------------------------------------------------------------
  FUNCTION:       DRLWrapDayCountToString
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file London style day count conv
                  and convert it into a string.

*/   
extern int DRLWrapDayCountToString(FILE *fp,
				   char *fileName,
				   char *paramName,  /* for error message */
				   char *string);    /* (O) */

/*----------------------------------------------------------------------
  FUNCTION:       DRLWrapDayCountToString
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file London style frequency
                  and convert it into a num per year.

*/   
extern int DRLWrapFrequencyToLong(FILE *fp,
				  char *fileName,
				  char *paramName,     /* for error message */
				  long *frequency);    /* (O) */
                
#endif     
