#ifndef ESL_DATE_DOT_H
#define ESL_DATE_DOT_H

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_alloc.h"
#include "esl_util.h"

#ifdef  __cplusplus
extern "C" {
#endif

typedef long    ESL_DATE;

extern long noleap[13];
extern long leap[13];
extern long cumdays[12];

long ADate2LDate(long date);
long LDate2ADate(long date);


/** Returns month, day, year from an integer in                
*   the format YYY(Y)MMDD as in 840701 or 20011115 or 1120104  */
void Dsplit(long date_i, long *mm_o, long *dd_o, long *yy_o);

/** Returns zero if 4 digit argument is not a leap year  *
 *  otherwise returns 1.                                 */
int Isleap(long year);

/**  Returns TRUE if the date is an IMM date    */
int Isimm(long    date);

long   ThirdWed(long mth, long year);

/** This routine checks whether the date passed in is a holiday  * 
    or not according to the HOLIDAYS file. It returns (0) if     * 
    it is not and (-1) if it is. The format is YYYYMMDD.         */
int IsHoliday(long date_i);

/** Checks for validity of YYYYMMDD formatted dates  * 
 *  returns 0 if all is well                         * 
 *  returns 1 if bad year                            * 
 *  returns 2 if bad month                           * 
 *  returns 3 if bad day                             */

int Dateok(long date);

int Date_CheckAndReport(long	date); /**< (I) Date */

/** Converts a 2 digit year (84) to a 4 digit year (1984);       * 
 *  assumes that 2 digit years less than 50 are after year 2000  */
long Y2toy4(long year_i);

/** Packs mm_i,dd_i,yy_i into YYYYMMDD format  * 
 *  calls Y2toy4 on yy_i before packing        */
long Datepack(long mm_i, long dd_i, long yy_i);

/** Converts a 4 digit year to a 2 digit year;  * 
 *  years after 99 are 00, 01, etc.             */
long Y4toy2(long year_i);

/**  Convert YYYYMMDD to MM/DD/YY  */
void Y2date_str(long date, char *string);

/**  Convert MM/DD/YY to YYYYMMDD  */
long eval_date(char *datest);

/**  Convert DD-MMM-YYYY to YYYYMMDD  */
long eval_date2(char *datest);

/**  Calculates day of week (0-6) of YYYYMMDD formatted date  */
long Dayofwk(long date);

/*****************************************************************************/
/**                                                                          * 
 *  FUNCTION   Days360                                                       * 
 *                                                                           * 
 *     Counts days between two dates following the 30/360 covention as       * 
 *     specified in  the ISDA rule book. The method below agrees  with       * 
 *     the Analytics Library and with the routines implemented by  the       * 
 *     STIRT group (thanks to Neill Penney for useful discussions!)          * 
 *                                                                           * 
 *                                                                           * 
 *****************************************************************************/
long Days360(long date1_i, long date2_i);

long Months360(long date1_i, long date2_i);

long Daysact(long date1_i, long date2_i);

/** Returns a date in YYYYMMDD format based on moving        * 
 *  forward or backwards a number of calendar days.          */
long Nxtday(long date, long days);

/** Returns a date in YYYYMMDD format based on moving        * 
 *  forward or backwards (mths) from (date).                 * 
 *  if eom > 0 the returned does not spill over to the next  * 
 *  month. That is moving 1m from the 31rst of march gives   * 
 *  the 30th of april. If eom=0 it does spill over: 1m from  * 
 *  the 31rst of march is the 1rst of may.	             * 
 *  a positive entry for mths => move forward                * 
 *  a negative entry for mths => move backward.              */
long Nxtmth(long date, long mths, long eom);


/** a positive entry for nbPeriods => move forward     * 
 *  a negative entry for nbPeriods => move backward.   */
long  Nxtimm(long  date, long nbPeriods);

/**************************************************************************/
/**                                                                       * 
 * FUNCTION    DrDayCountFraction                                         * 
 *                                                                        * 
 * Calculates the daycount fraction between two dates according to one of * 
 * three methods: Actual/365Fixed, Actual/360 and 30/360.                 * 
 *                                                                        * 
 * RETURNS:  Success or Failure                                           * 
 *                                                                        */
/**************************************************************************/

int  DrDayCountFraction(long     Date1,      /**< (I) Start date             */
                        long     Date2,      /**< (I) End date               */
                        char     Conv,       /**< (I) Convention (A,0,3,5)   */
                        double  *Fraction);  /**< (O) Corresponding fraction */

double  DrDcf(ESL_DATE from,    /**< (I) Start date           */
              ESL_DATE to,      /**< (I) End date             */
              ESL_DCC  dcc);    /**< (I) Convention (0,3,5)   */

/** This routine returns a new date after advancing or going  * 
 *  backward a number of weekdays (Monday - Friday).          * 
 *  The format for date is YYYYMMDD.                          */
long Nxtwkday(long date, long advance);

/** This routine returns a new date after advancing or going  * 
 *  backward a number of business days (Monday - Friday &     * 
 *  holidays). It checks for holidays in HOLIDAY file.        * 
 *  The function returns FAILURE if it is not possible to     * 
 *  call the IsHoliday function successfully.                 * 
 *  The format for date is YYYYMMDD.                          */
long Nxtbusday(long date, long advance);

/*****************************************************************************/
/**                                                                          * 
 *   FUNCTION    DrNewEventListFromFreqWithInterpType                        * 
 *                                                                           * 
 * Better interface to DrNewEventListFromFreq where an interpolation         * 
 * method should be specified (instead of clients negating NbInpDates).      * 
 *                                                                           */
 EVENT_LIST * DrNewEventListFromFreqWithInterpType( 
	ESL_INTERP_TYPE interpType,  /**< (I) Interpolation type to use*/
        int         NbInpDates,  /**< (I) Nb of dates input directly   */
        long       *InpDates,    /**< (I) Dates given directly by user */
        char        Freq,        /**< (I) Frequency of event           */
        char        Stub,        /**< (I) Stub location Front or Back  */
        char        DatesIn,     /**< (I) Y=input dates must be in list*/
        double     *Curve0,      /**< (I) Set of values for event      */
        double     *Curve1,      /**< (I) Set of values for event      */
        double     *Curve2,      /**< (I) Set of values for event      */
        double     *Curve3,      /**< (I) Set of values for event      */
        double     *Curve4);     /**< (I) Set of values for event      */

/*****************************************************************************/
/**                                                                          *
 *   FUNCTION    DrNewEventListFromFreq                                      * 
 *                                                                           * 
 *   Added 8/04 -- allow staircase intepolation (use NEGATIVE NbInpDates)    * 
 *                                                                           */
 EVENT_LIST * DrNewEventListFromFreq
               (int    NbInpDates, /**< (I) Nb of dates input directly   */
                long   *InpDates,  /**< (I) Dates given directly by user */
                char    Freq,      /**< (I) Frequency of event           */
                char    Stub,      /**< (I) Stub location Front or Back  */
                char    DatesIn,   /**< (I) Y=input dates must be in list*/
                double *Curve0,    /**< (I) Set of values for event      */
                double *Curve1,    /**< (I) Set of values for event      */
                double *Curve2,    /**< (I) Set of values for event      */
                double *Curve3,    /**< (I) Set of values for event      */
                double *Curve4);   /**< (I) Set of values for event      */

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrDateFwdAny                                               * 
 *                                                                          * 
 *   Starting from a given date, this function calculates the date obtained * 
 *   by incrementing a number of periods.                                   * 
 *                                                                          */
int     DrDateFwdAny (  long    StartDate, /**< (I) Start date          */
                        int     NbPers,    /**< (I) Number of periods   */
                        char    Freq,      /**< (I) Frequency           */
                        char    FwdOrBwd,  /**< (I) Forward or backward */
                        long    *DateOut); /**< (O) Calculated date     */


/****************************************************************************/
/*                                                                          * 
 *   FUNCTION    DateListFromFreq                                           * 
 *                                                                          * 
 *   Generates a date list including the start and end dates given and in   * 
 *   accordance with the frequency and stub convention passed in.           * 
 *                                                                          */
 int DateListFromFreq
              (long       StartDate,     /**< (I) Start of date list       */
               long       EndDate,       /**< (I) End of date list         */
               char       Freq,          /**< (I) Frequency of list        */
               char       StubConv,      /**< (I) Stub location Front/back */
               int       *NbDates,       /**< (O) Number of dates in list  */
               long     **DateList);     /**< (O) List of dates asc order  */


/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrDatesInSchedule                                          * 
 *                                                                          * 
 *   Returns SUCCESS if all dates passed in are contained in the schedule   * 
 *   specified by a start date, an end date, a frequency and a stub conv.   * 
 *                                                                          * 
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       * 
 *                                                                          */
 int DrDatesInSchedule
              (int        NbDatesIn,     /**< (I) Number of dates to check   */
               long      *DatesIn,       /**< (I) User dates                 */
               long       StartDate,     /**< (I) Start of date list         */
               long       EndDate,       /**< (I) End of date list           */
               char       Freq,          /**< (I) Frequency of list          */
               char       StubConv);     /**< (I)                            */

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrSameDateSchedules                                        * 
 *                                                                          * 
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *  
 *   specifing a start date, an end date, a frequency, a stub conv,         * 
 *   arrears reset and a value date, which is used to cut the list before it* 
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       * 
 *                                                                          */
 int DrSameDateSchedules
              (int       *NbDatesIn,     /**< (I/0) Number of dates to check */
               long      *DatesIn,       /**< (I) User dates                 */
               long       StartDate,     /**< (I) Start of date list         */
               long       EndDate,       /**< (I) End of date list           */
               long       ValueDate,     /**< (I) Value date                 */
               char       Freq,          /**< (I) Frequency of list          */
               char       StubConv,      /**< (I)                            */
               char       Arrears);      /**< (I) 'Y' if reset-in-arrears    */


/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrSameDateSets                                             * 
 *                                                                          * 
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *  
 *   cuting a givne date list up to a value date.                           * 
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       * 
 *                                                                          */
int DrSameDateSets
              (int       *NbDates1,      /**< (I/0) Number of dates to check */
               long      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               long      *Dates2,        /**< (I) Cpm dates                  */
               long       ValueDate);    /**< (I) Value date                 */

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrSameDateSetsFlows                                        * 
 *                                                                          * 
 *   Returns SUCCESS if dates passed in are the same as dates generated by  *  
 *   cuting a givne date list up to a value date.                           * 
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       * 
 *                                                                          */
int DrSameDateSetsFlows
              (int       *NbDates1,      /**< (I/0) Number of dates to check */
               long      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               long      *Dates2,        /**< (I) Cpm dates                  */
               long       ValueDate);    /**< (I) Value date                 */


/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrDatesIn2FreqSchedule                                     * 
 *                                                                          * 
 *   Compares the passed-in date list to a date list generated in the       * 
 *   current (with respect to the value date) pmt period according to some  * 
 *   obs frequency. The current pmt period is found from the start date,    * 
 *   end date, pmt frequency, stub convention and a value date              * 
 *                                                                          */
int DrDatesIn2FreqSchedule
              (int        NbDatesIn,     /**< (I) Number of dates to check   */
               long      *DatesIn,       /**< (I) User dates                 */
               long       StartDate,     /**< (I) Start of date list         */
               long       EndDate,       /**< (I) End of date list           */
               long       ValueDate,     /**< (I) Value date                 */
               char       PmtFreq,       /**< (I) Frequency of pmts          */
               char       StubConv,      /**< (I)                            */
               char       ObsFreq,       /**< (I) Frequency of obs           */
               int       *NbDatesOut,    /**< (O) Number of found dates      */
               long      *DatesOut);     /**< (O) Found dates                */

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    DrDatesInSet                                               * 
 *                                                                          * 
 *   Returns SUCCESS if a first set of dates is contained in the second set * 
 *   of dates.                                                              * 
 *                                                                          * 
 *   DatesIn must be an arrays of ordered, strictly increasing dates.       * 
 *                                                                          */
int DrDatesInSet
              (int        NbDates1,      /**< (I) Number of dates to check   */
               long const*Dates1,        /**< (I) Dates to check             */
               int        NbDates2,      /**< (I) Number of org dates        */ 
               long const*Dates2);       /**< (I) Org dates                  */


/****************************************************************************/
/**                                                                         * 
 *   FUNCTION   DrFreeEventList                                             * 
 *                                                                          */
void DrFreeEventList(EVENT_LIST  *List);

/*****  AddDateToList  **********************************************/
/**
 *      Allocates memory and add a date to the end of a list of date
 *      Returns SUCCESS or FAILURE
 */

int     AddDateToList
            (int      *NbDates,
             long    **DateList,
             long      NewDate);

/*****  SortDateList  **********************************************/
/**
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order, 
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 */

int     SortDateList(int    NbDates, 
                     long  *DateList,
                     long  *SuppValue);


/*********** MergeDateLists ************************************************/
/**                                                                         * 
 *      1-merges two date lists,                                            * 
 *      2-sorts them                                                        * 
 *      3-removes possible duplicates                                       * 
 *                                                                          * 
 *                                                                          * 
 *      Note: that DateList (first input list) is preserved                 * 
 *      whereas MergedList (second input list) is modified by the function. * 
 *                                                                          * 
 *                                                                          * 
 ****************************************************************************/
int MergeDateLists(
    int  NbDatesList,    /**< (I)Nb dates in list to be merged with Mergedlist*/
    long *DateList,      /**< (I)List to be merged with MergedList            */
    int  *MergedListSize,/**< (I/O) MergedList size after duplicates removal  */
    long **MergedList);  /**< (I/O) Absorbing list                            */


/*****  GetDLOffset  ***************************************************/
/**
 *      Given a datelist and a targetDate,
 *      returns the offset of the element nearest to the targetdate
 *      if mode = 
 *          CbkEXACT:  returns -999 if no matching targetdate is found
 *
 *          CbkLOWER:  returns the nearest offset that is LOWER or equal;
 *                     or -999 if all dates > targetdate
 *          CbkHIGHER: returns the nearest offset that is HIGHER or equal;
 *                     or -999 if all dates < targetdate
 */

int     GetDLOffset(int         NbDates,
                    long const* DL,
                    long        targetDate,
                    int         mode);

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    hasStub                                                    * 
 *                                                                          * 
 *   Boolean function which returns 0 if the StartDate and EndDate define   *  
 *   an integer number (>=0) of periods (i.e. no stub), 1 if there exists   * 
 *   a stub.                                                                * 
 *   Assumes StartDate > EndDate and Freq = (A,S,Q,M,W)                     * 
 *                                                                          */

 int hasStub (long       StartDate,     /**< (I) Start date                   */
              long       EndDate,       /**< (I) End date list                */
              char       Freq);         /**< (I) Frequency                    */

/****************************************************************************/
/**                                                                         * 
 *   FUNCTION    hasStubLoc                                                 * 
 *                                                                          * 
 *   Boolean function which returns 0 if the StartDate and EndDate define   *  
 *   an integer number (>=0) of periods (i.e. no stub), 1 if there exists   * 
 *   a stub. The function accounts for stub location ('F' or 'B')           * 
 *   Assumes StartDate > EndDate and Freq = (A,S,Q,M,W)                     * 
 *                                                                          */

  int hasStubLoc (long       StartDate,     /* (I) Start date    */
                  long       EndDate,       /* (I) End date list */
                  char       Freq,          /* (I) Frequency     */
		  char       StubLoc);      /* (I) Stub location ('F' or 'B')  */

              
/** Generate output array by copying input dates and generating additional
 *  dates between first and last date. Output is sorted and has no duplicates.
 *  Input and output arrays can be the same or different, non-overlapping. 
 *  Returns the number of elements in the output array.
 */
unsigned 
ExpandDateSchedule(unsigned sSize,          /**< (I) size of the input array */
                   ESL_DATE const* sSched,  /**< (I) input array             */
                   unsigned tSize,          /**< (I) size of the output array*/
                   ESL_DATE* tSched,        /**< (I) output array            */
                   char frequency);         /**< (I) frequency-I,D,M,W,Q,S,A */

/*****  ExtendTCurve  ********************************************************/
/**
 *      Flat extend zero curve between fromDate and toDate dates 
 */
 int  ExtendTCurve(T_CURVE *tc, ESL_DATE fromDate, ESL_DATE toDate);

#ifdef  __cplusplus
}
#endif


#endif

