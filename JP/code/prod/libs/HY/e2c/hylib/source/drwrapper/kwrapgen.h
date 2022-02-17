/*------------------------------------------------------------------
HEADER FILE:    kwrapgen.h

CREATED BY:     Leda Braga - April 1996
                David Fung
		        Julia Chislenko

MODIFIED :      Antonio Paras - Feb 1999 took away spot vol calibration
                                         added function to read FXVolatility File.

PURPOSE:        Data conversion functions to be used in the 
                Kapital DR wrapper in connection with  ALIB
                based executables.
---------------------------------------------------------------------- */
#ifndef _KWRAPGEN_H
#define _KWRAPGEN_H


#include <stdio.h>
#include "bastypes.h"
#include "tcurve.h"
/*#include "cfileio.h"*/

#define  BUFFER_SIZE   250

/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapReadEquitySta

CREATED BY:  Neil Yang -  Febuary 2000

DESCRIPTION: Kapital wrapper function to read the equaity1.sta file to get divident info  
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT (int) DrKapWrapReadEquitySta(char   *filename,       /*(I) Usually "zero.dat" */
										long   *numPoints,      /* (O) */
									    TDate  **dividentDates, /* (O) */
									    double **dividentRates); /* (O) */

/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapReadEquityDyn

CREATED BY:  Neil Yang -  Febuary 2000

DESCRIPTION: Kapital wrapper function to read the equaity1.dyn file to get value date,
			 price, and vol curve (i.e. TCurve format)  
             
             vol curve returned is checked for date order and for non-negative
             zero rates. 
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT (int) DrKapWrapReadEquityDyn(char   *filename,   /*(I) Usually "zero.dat" */
									  TDate  valueDate,  /* (O) */
									  double *stockPrice, /* (O) */
									  double *correlation, /* (O) */
									  TCurve **volCurve);    /* (O) */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetZC

CREATED BY:  Leda Braga -  April 1996

DESCRIPTION: Kapital wrapper function to read the zero.dat file (in DR
             format) into a TCurve zero curve
----------------------------------------------------------------------------*/
GTO_EXPORT(TCurve *)DrKapWrapGetZC(char  *filename);



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetFXVol

CREATED BY:  Antonio Paras -- Feb 1999 based on DrKapWrapGetVC

DESCRIPTION: Kapital wrapper function to read the fxvolatility.dat file (in 
             DR format) into a Library vol curve (i.e. TCurve format) and return
             the fxspot value.
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(TCurve *) DrKapWrapGetFXVol
                     (char      *filename,   /*(I) Usually "FXVolatility.dat"     */
                      double    *FXSpotRate);/*(O) FX Spot Rate  */



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
GTO_EXPORT(TCurve *) DrKapWrapGetVC
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

GTO_EXPORT(int ) DrKapWrapGetSwpvolGrid
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
GTO_EXPORT(int ) DrKapWrapGetOneFacPars
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
GTO_EXPORT(int ) DrKapWrapGetTwoFacPars
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
GTO_EXPORT(int ) DrDateToTDate(long   DRDate,
                                TDate *libDate);      
   



/*----------------------------------------------------------------------------
FUNCTION:    DrFindAndSkipComLine

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper convenience function which finds and skips 
             a comment line  and flags an error if line is not found or 
             is not a comment.
 
                            FILE MUST BE OPEN!
----------------------------------------------------------------------------*/                                 
GTO_EXPORT(int ) DrFindAndSkipComLine(FILE *fp,
                                        char  *filename);  
 



/*----------------------------------------------------------------------------
FUNCTION:    DrDRLabelToTDateInterval

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience to convert a DR date interval to
             a TDateInterval. It allows for "2Y" and "2A" as valid strings

----------------------------------------------------------------------------*/     
GTO_EXPORT(int ) DrDRLabelToTDateInterval
                (char             *matLabel,
                 TDateInterval    *matInterval);   


                                                                                  
                                                                                  
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCheckTCurve

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience check a TCurve for dates in
             ascending order and for non-negative rates (either zeros
             or vols). 

----------------------------------------------------------------------------*/    
GTO_EXPORT(int ) DrKapWrapCheckTCurve(TCurve     *curve,
                                        char       *label);


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetRateConv

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to get standard rate
             conventions from a zero curve env file

----------------------------------------------------------------------------*/    
GTO_EXPORT(int ) DrKapWrapGetRateConv(
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
GTO_EXPORT(int ) DrKapWrapDCCtoString(
    long clibDCC,       /* (I) CLib Daycount conv identifier */
    char *DCCstring);   /* (O) string output                 */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapConvDCC

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to convert DR wrapper
             day count convention to clib convention

----------------------------------------------------------------------------*/    
GTO_EXPORT(int ) DrKapWrapConvDCC(
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
GTO_EXPORT(int ) DrKapWrapDatesInList(
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
GTO_EXPORT(int ) DrKapWrapExpandStepUpSchedule(
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
GTO_EXPORT(int ) DrKapWrapCalcAmortPrns(
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
GTO_EXPORT(int ) DrKapWrapParamInput(                                                   
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
GTO_EXPORT(int ) DrKapWrapParamCheck(
    long      nbFactor,     /* (I) Number of factors     */
    double   *alphas,       /* (I) alphas                */
    double   *betas,        /* (I) mean reversions       */
    double   *corrs);       /* (I) corrleations          */


/*----------------------------------------------------------------------------
FUNCTION:    DrExtendTCurve

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to extend a curve (e.g. zero
             curve) to include dates implied by BaseDate and intervals
             NOTE: baseDate must be >= oriCurve->fBaseDate
----------------------------------------------------------------------------*/

GTO_EXPORT(int ) DrExtendTCurve(
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
GTO_EXPORT(int ) DRLWrapTDateIntvlToNum(FILE *fp,
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
GTO_EXPORT(int ) DRLWrapDayCountToString(FILE *fp,
				   char *fileName,
				   char *paramName,  /* for error message */
				   char *string);    /* (O) */

/*----------------------------------------------------------------------
  FUNCTION:       DRLWrapDayCountToString
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file London style frequency
                  and convert it into a num per year.

*/   
GTO_EXPORT(int ) DRLWrapFrequencyToLong(FILE *fp,
				  char *fileName,
				  char *paramName,     /* for error message */
				  long *frequency);    /* (O) */
                
#endif     
