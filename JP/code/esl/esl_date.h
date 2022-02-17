/****************************************************************************/
/*      Date manipulation functions: conversion, day counting, settlement,  */
/*      calculation of relative dates, next business day, etc.              */
/****************************************************************************/

#ifndef ESL_DATE_DOT_H
#define ESL_DATE_DOT_H

#include "esl_types.h"

#ifdef  __cplusplus
extern "C" {
#endif

/*****************************************************************************
** The following typedefs are used to both improve code readability and
** insulate client code from changes in date representation (however, there
** are many places where "long" is still in use, and some code relies on
** implicit conversion to/from an integer type (i.e. long)).
**
** Note that the internal date representation (IRDate) depends upon whether or 
** not we are using IRX dates (i.e. -DESL_NEW_DATE).  Without IRX dates, the
** internal date representation (IRDate) is the same as YMDDate.
**
** IRX dates share the same internal representation as ALib.  There is a minor
** difference with QLib dates in that QLib uses "int" instead of "long".
** In reality, we (and the ALib) should use "int" (assuming that ints are at
** least 26-bits, since even this would take us until 31-Dec-3355).
**
** SEE ALSO:  Definition of IDATE, IDATE_PTR, and IDATE_D_PTR in esl_macros.h.
*****************************************************************************/
typedef long    IRDate;     /* Internal date representation  */
typedef long    YMDDate;    /* Dates represented as YYYYMMDD */

/*
extern long noleap[13];
extern long leap[13];
extern long cumdays[12];
*/

/** For consistency, we define a specific constant to represent an invalid  *
 *  date.  However, this is not the only possible invalid date!             *
 *  Any negative number is an invalid date when BOTH using and not using    *
 *  IRX dates.                                                              */
#define INVALID_DATE -999


/*****************************************************************************
** DATE CONSTRUCTION/CONVERSION
**
** The following functions are used to convert to and from the internal date
** representation (IRDate).
*****************************************************************************/

/** Answers an internal date (IRDate) from a given YMDDate (YYYYMMDD) */
IRDate IRDateFromYMDDate(YMDDate date); 

/** Answers a YMDDate from a given internal date (IRDate) */
YMDDate YMDDateFromIRDate(IRDate date);

/* Convert an ALib/IRX-style TDate into an internal IRDate.                  */
/* We use "long" instead of "TDate" to avoid dependencies on ALib headers.  */
/* We don't introduce our own typedef for TDate to avoid conflicts with     */
/* code that does use the ALib.                                             */
IRDate IRDateFromTDate(long tDate);

/**  Answers an internal date (IRDate) from a date string in "MM/DD/YY"      *
 *   format.                                                                */
IRDate eval_date(char *dateString);

/**  Answers an internal date (IRDate) from a date string in "DD-MON-YYYY"   *
 *   format (e.g. 12-Sep-2006), where MON is case-insensitive.              */
IRDate eval_date2(char *dateString);

/** Creates an IRDate from a given month, day, and year.
 *  Year (yy_i) is expanded to four digits via Y2toy4 */
IRDate Datepack(long mm_i, long dd_i, long yy_i);

/** Splits a date into its constituent month, day, and (4-digit) year.      */
void Dsplit(IRDate date_i, long *mm_o, long *dd_o, long *yy_o);

/**  Converts an IRDate to "DD-MON-YYYY" (e.g. 5-Sep-2006).                  *
 *   The "string" parameter is assumed to point to a pre-allocated memory   *
 *   location that is big enough to hold the resulting string (11 chars     *
 *   plus NULL terminator).                                                 */
void StringFromIRDate(IRDate date, char *string);

/**  Converts an IRDate to "MM/DD/YY".                                       * 
 *   The "string" parameter is assumed to point to a pre-allocated memory   *
 *   location that is big enough to hold the resulting string (8 chars      *
 *   plus NULL terminator).                                                 */
void Y2date_str(IRDate date, char *string);

/** Converts a 2 digit year (e.g. 84) to a 4 digit year (e.g. 1984);        * 
 *  assumes that 2 digit years less than 50 are after year 2000             */
long Y2toy4(long year_i);

/** Converts a 4 digit year to a 2 digit year;  * 
 *  years after 99 are 00, 01, etc.             */
long Y4toy2(long year_i);


/*****************************************************************************
** DATE MANIPULATION FUNCTIONS
*****************************************************************************/

/** Returns zero if 4 digit argument is not a leap year  *
 *  otherwise returns 1.                                 */
int Isleap(IRDate year);

/**  Returns TRUE if the date is an IMM date    */
int Isimm(IRDate    date);

IRDate   ThirdWed(long mth, long year);

/** This routine checks whether the date passed in is a holiday  * 
    or not according to the HOLIDAYS file. It returns (0) if     * 
    it is not and (-1) if it is.  The dates in the HOLIDAYS      *
    file are expected to be in YMDDate (YYYYMMDD) format.        */ 
int IsHoliday(IRDate date_i);

/** Checks for validity of an IRDate                  * 
 *  returns 0 if all is well                         * 
 *  returns 1 if bad year                            * 
 *  returns 2 if bad month                           * 
 *  returns 3 if bad day                             */

int Dateok(IRDate date);

int Date_CheckAndReport(long	date); /**< (I) Date */


/**  Calculates day of week (0-6) of IRDate */
long Dayofwk(IRDate date);

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
long Days360(IRDate date1_i, IRDate date2_i);

long Months360(IRDate date1_i, IRDate date2_i);

long Daysact(IRDate date1_i, IRDate date2_i);

/** Returns an IRDate based on moving forward or backwards a number of       *
 *  calendar days.                                                          */
IRDate Nxtday(IRDate date, long days);

/** Returns a date in IRDate format based on moving           * 
 *  forward or backwards (mths) from (date).                 * 
 *  if eom > 0 the returned does not spill over to the next  * 
 *  month. That is moving 1m from the 31rst of march gives   * 
 *  the 30th of april. If eom=0 it does spill over: 1m from  * 
 *  the 31rst of march is the 1rst of may.	             * 
 *  a positive entry for mths => move forward                * 
 *  a negative entry for mths => move backward.              */
IRDate Nxtmth(IRDate date, long mths, long eom);


/** a positive entry for nbPeriods => move forward     * 
 *  a negative entry for nbPeriods => move backward.   */
IRDate  Nxtimm(IRDate  date, long nbPeriods);

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

int  DrDayCountFraction(IRDate     Date1,      /**< (I) Start date             */
                        IRDate     Date2,      /**< (I) End date               */
                        char     Conv,       /**< (I) Convention (A,0,3,5)   */
                        double  *Fraction);  /**< (O) Corresponding fraction */

double  DrDcf(IRDate from,    /**< (I) Start date           */
              IRDate to,      /**< (I) End date             */
              ESL_DCC  dcc);    /**< (I) Convention (0,3,5)   */

/** This routine returns a new date after advancing or going  * 
 *  backward a number of weekdays (Monday - Friday).          */ 
IRDate Nxtwkday(IRDate date, long advance);

/** This routine returns a new date after advancing or going  * 
 *  backward a number of business days (Monday - Friday &     * 
 *  holidays). It checks for holidays in HOLIDAY file.        * 
 *  The function returns FAILURE if it is not possible to     * 
 *  call the IsHoliday function successfully.                 */
IRDate Nxtbusday(IRDate date, long advance);


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
        IRDate       *InpDates,    /**< (I) Dates given directly by user */
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
                IRDate   *InpDates,  /**< (I) Dates given directly by user */
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
int     DrDateFwdAny (  IRDate    StartDate, /**< (I) Start date          */
                        int     NbPers,    /**< (I) Number of periods   */
                        char    Freq,      /**< (I) Frequency           */
                        char    FwdOrBwd,  /**< (I) Forward or backward */
                        IRDate    *DateOut); /**< (O) Calculated date     */


/****************************************************************************/
/*                                                                          * 
 *   FUNCTION    DateListFromFreq                                           * 
 *                                                                          * 
 *   Generates a date list including the start and end dates given and in   * 
 *   accordance with the frequency and stub convention passed in.           * 
 *                                                                          */
 int DateListFromFreq
              (IRDate       StartDate,     /**< (I) Start of date list       */
               IRDate       EndDate,       /**< (I) End of date list         */
               char       Freq,          /**< (I) Frequency of list        */
               char       StubConv,      /**< (I) Stub location Front/back */
               int       *NbDates,       /**< (O) Number of dates in list  */
               IRDate     **DateList);     /**< (O) List of dates asc order  */


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
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
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
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
               IRDate       ValueDate,     /**< (I) Value date                 */
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
               IRDate      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               IRDate      *Dates2,        /**< (I) Cpm dates                  */
               IRDate       ValueDate);    /**< (I) Value date                 */

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
               IRDate      *Dates1,        /**< (I) User dates                 */
               int        NbDates2,      /**< (I) Number of dates to cmp     */
               IRDate      *Dates2,        /**< (I) Cpm dates                  */
               IRDate       ValueDate);    /**< (I) Value date                 */


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
               IRDate      *DatesIn,       /**< (I) User dates                 */
               IRDate       StartDate,     /**< (I) Start of date list         */
               IRDate       EndDate,       /**< (I) End of date list           */
               IRDate       ValueDate,     /**< (I) Value date                 */
               char       PmtFreq,       /**< (I) Frequency of pmts          */
               char       StubConv,      /**< (I)                            */
               char       ObsFreq,       /**< (I) Frequency of obs           */
               int       *NbDatesOut,    /**< (O) Number of found dates      */
               IRDate      *DatesOut);     /**< (O) Found dates                */

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
               IRDate const*Dates1,        /**< (I) Dates to check             */
               int        NbDates2,      /**< (I) Number of org dates        */ 
               IRDate const*Dates2);       /**< (I) Org dates                  */




/*****  AddDateToList  **********************************************/
/**
 *      Allocates memory and add a date to the end of a list of date
 *      Returns SUCCESS or FAILURE
 */

int     AddDateToList
            (int      *NbDates,
             IRDate    **DateList,
             IRDate      NewDate);

/*****  SortDateList  **********************************************/
/**
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order, 
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 */

int     SortDateList(int    NbDates, 
                     IRDate  *DateList,
                     IRDate  *SuppValue);


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
    IRDate *DateList,      /**< (I)List to be merged with MergedList            */
    int  *MergedListSize,/**< (I/O) MergedList size after duplicates removal  */
    IRDate **MergedList);  /**< (I/O) Absorbing list                            */


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
                    IRDate const* DL,
                    IRDate        targetDate,
                    int         mode);





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


/** Generate output array by copying input dates and generating additional
 *  dates between first and last date. Output is sorted and has no duplicates.
 *  Input and output arrays can be the same or different, non-overlapping. 
 *  Returns the number of elements in the output array.
 */
unsigned 
ExpandDateSchedule(unsigned sSize,          /**< (I) size of the input array */
                   IRDate const* sSched,  /**< (I) input array             */
                   unsigned tSize,          /**< (I) size of the output array*/
                   IRDate* tSched,        /**< (I) output array            */
                   char frequency);         /**< (I) frequency-I,D,M,W,Q,S,A */


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
 *   FUNCTION   DrFreeEventList                                             * 
 *                                                                          */
void DrFreeEventList(EVENT_LIST  *List);


/* other date functions */

/* Converts an Excel date as long (1L = 01-Jan-1900) to a date */
int     DrDateFromExcel(long xlDate, IRDate *date);

/* Converts a date to an Excel date as long */
int     DrDateToExcel(long date, IRDate *xlDate);

/**-------------------------------------------------------------
 * Returns Julian day number corresponding to mm, dd, and yyyy.
 */
int DrJulday(int mm, int id, int iyyy, long *julDate);


/**-------------------------------------------------------------
 * Returns Julian day number corresponding to "IRDate" 
 */
int DrDate2Julday(IRDate drDate, long *julDate);


/**-------------------------------------------------------------
 * Returns IRDate corresponding to a Julian day number "julDate".
 */
int DrJulday2Date(long julDate, IRDate *drDate);

/**-------------------------------------------------------------
 * Convert an Excel date to to Julian date.
 * The offset between Julian day and Excel day is supplied
 * by the API through global variable JULIAN_EXCEL_OFFSET.
 */
long DrJulday2Excel(long julDate);

long DrExcel2Julday(long xlDate);

#ifdef  __cplusplus
}
#endif

#endif  /* #ifndef ESL_DATE_DOT_H */
