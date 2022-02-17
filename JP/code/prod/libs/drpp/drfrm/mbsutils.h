/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsutils.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Robert Lenk
 *			   Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Version		:  1.64
 *	Extracted	:  6/18/97 at 13:06:12
 *	Last Updated	:  6/18/97 at 13:06:09
 ***************************************************************************
 *	Misc utilities for DR/MBS analytics base lib 
 *      (file i/o, memory, etc)
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/

#ifndef __mbsutils_h
#define __mbsutils_h

#ifdef __cplusplus
extern "C" {
#endif

#include "mdydate.h"
#include "macros.h"
#include "convert.h"
#include "dateconv.h"
#include "ldate.h"
#include "fltratea.h"
#include "zr2simp.h"
#include "zr2coup.h"
#include "cdate.h"
#include "cerror.h"
#include "bastypes.h"
#include "datelist.h"

/***************************************************************************
 *  CONSTANTS
 **************************************************************************/

/***************************************************************************
 *  PUBLIC FUNCTIONS
 **************************************************************************/

/***************************************************************************
 *  "Behavior control" funcs
 **************************************************************************/

/*****************************************************************
 *  set_behavior
 *  Turns specified behavior "on"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int set_behavior(char *behavior_label);	    /* (I) some label, e.g. "smooth oas" */

/*****************************************************************
 *  unset_behavior
 *  Turns specified behavior "off"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int unset_behavior(char *behavior_label);    /* (I) some label, e.g. "smooth oas" */

/*****************************************************************
 *  is_behavior
 *  BOOLEAN function: returns TRUE if behavior is "on"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int is_behavior(char *behavior_label);	    /* (I) some label, e.g. "smooth oas" */


/*****************************************************************************
 *  File I/O
 ****************************************************************************/

/****************************************************************************
 * MbsBuildFileName
 * File name (string) utility: prepends directory name to (raw) file
 * name, if that file name currently lacks a path
 * NOTE: calling code must have allocated fullFileName as string
 * of length MAXPATHLEN
 * NOTE: OK to use same file name string for input & output (raw/fullFileName)
 ****************************************************************************/
int MbsBuildFileName
   (char    *pathName,          /* (I) name of path (dir) to prepend */
    char    *rawFileName,       /* (I) raw name of file */
    char    *fullFileName);     /* (O) full name for file */

/****************************************************************************
 * MbsGetDefOutputPath
 * Determines the default output directory 
 * (at present, this is $HOME on UNIX, and c:\ on PC)
 * Calling code must have already allocated output string
 ****************************************************************************/
int MbsGetDefOutputPath
   (long    maxNumChars,        /* (I) max # chars to copy to outputPath */
    char   *outputPath);        /* (O) default output path */

/****************************************************************************
 *  MbsAddTrailSlash
 *  (DOS or UNIX) Adds trailing fwd (UNIX) or bkwd (DOS) slash to
 *  file/path spec string
 ****************************************************************************/
int MbsAddTrailSlash
   (char *pathSpec);            /* (I/O) path spec to alter */

/****************************************************************************
 *  MbsRemoveTrailSlash
 *  (DOS or UNIX) Removes trailing fwd (UNIX) or bkwd (DOS) slash 
 *  from file/path spec string; 
 *  Can handle case of root path ("\" or "/") in two ways
 ****************************************************************************/
int MbsRemoveTrailSlash
   (char     *pathSpec,             /* (I/O) path spec to alter */
    TBoolean  keepRootPath);        /* (I) if input path is root path
                                     * (i.e., just single slash), then if:
                                     *   TRUE: keep the trailing (only) slash
                                     *   FALSE: remove trailing slash (will
                                     *   then return empty string for path) */

/****************************************************************************
 *  mbs_filefound
 *  BOOLEAN func to check if file exists
 *    If rootname contains a path delim, is ordinary "file-find" function
 *    Otherwise, looks in all the "usual" data locations ($DATA/LATEST, etc);
 *  Returns boolean true if found
 ****************************************************************************/
int mbs_filefound(char *rootname);    /* (I) root name (w/o path) of file */

/****************************************************************************
 *  mbs_fopen
 *  Finds and opens data file for input (w/"smart" file-find funcs)
 *  Returns pointer to file (or NULL, if failed)
 *  NB: if rootname contains slash char, it's assumed to be full pathspec
 ****************************************************************************/
FILE *mbs_fopen(char *rootname,      /* (I) root name of file */
		char *mode);         /* (I) file mode ("r", "a", etc.) */

/****************************************************************************
 *  mbs_findfile
 *  Fundamental function to find data file
 *    If rootname contains a path delim, is ordinary "file-find" function
 *    Otherwise, looks in all the "usual" data locations ($DATA/LATEST, etc);
 *  Sets fullname to full name of file (if found), or empty string (if not found)
 *  "Usual" dirs: user-defined [if set], DRDATA, DATA/LATEST, 
 *                /opt/fims/1.0.0.0/data/LATEST
 *  NB: function itself does not fail (return FAILURE) if file not found;
 *  only fails if some fatal error occurs
 ****************************************************************************/
int mbs_findfile(char *rootname,       /* (I) root name of file */
		 long maxnamelen,       /* (I) max # chars to copy to fullname */
		 char *fullname);      /* (O) (if found) full name of file */

/***************************************************************
 * mbs_delfile
 * Returns SUCCESS/FAILURE
 ***************************************************************/
int mbs_delfile(char* fname);	      /* (I) name of file to delete */

/*****************************************************************************
 *  Arith. funcs
 ****************************************************************************/

/**************************************************
 *  mbs_pow
 *  Faster version of pow() (uses exp/ln)
 **************************************************/
double mbs_pow(double x, double y);

/**************************************************
 *  mbs_sgn
 *  Returns -1. for x<0, 0. for x=0, 1. for x>0
 **************************************************/
double mbs_sgn(double x);

/**************************************************
 *  mbs_round
 *  Returns rounded value 
 **************************************************/
long mbs_round(double x);

/**************************************************
 *  mbs_lininterp1
 *  1-d linear interpolation (assumes ascending xarr)
 *  error if not in ascending order (3/28/96 - DL) 
 *  Returns boolean true if successful
 *  NB: does NOT extrapolate
 **************************************************/
int mbs_lininterp1(long n,	      /* (I) # pts */
		   double xarr[],     /* (I) array of indep var */
		   double yarr[],     /* (I) array of dep var */
		   double x,	      /* (I) value at which to interpolate */
		   double *y);	      /* (O) interpolated value */


/**************************************************
 *  mbs_lininterp2
 *  1-d linear interpolation (assumes ascending xarr)
 *  error if not in ascending order (3/28/96 - DL) 
 *  Returns boolean true if successful
 *  NB: does extrapolate
 **************************************************/
int mbs_lininterp2(long n,	      /* (I) # pts */
		   double xarr[],     /* (I) array of indep var */
		   double yarr[],     /* (I) array of dep var */
		   double x,	      /* (I) value at which to interpolate */
		   double *y);	      /* (O) interpolated value */

/***************************************************************************
 *  Date funcs (supplementing GTO funcs)
 ***************************************************************************/

/************************************************************
 *  GetDayCountStr
 *  Converts GTO day count integer into string
 *  @@ someday, find GTO macro/func that does this
 *  and use it inside this func (couldn't find it today...)
 ************************************************************/
int EXPORT GetDayCountStr
   (long      dayCount,         /* (I) day count int */
    char     *dayCountStr);     /* (O) string rep */

/************************************************************
 * SimplifyDL()
 * Modifies date list (may actually shorten list):
 * sorts in order of increasing date,
 * and removes any repeated or zero dates 
 ************************************************************/
int EXPORT SimplifyDL
   (TDateList **DL);  /* (I/O) ptr to ptr for date list */

/************************************************************
 * Simply exchanges the values of the two date vars
 ************************************************************/
int EXPORT SwitchDates
   (TDate *date1,               /* (I/O) */
    TDate *date2);              /* (I/O) */

/***************************************************************************
 * MbsWala
 * Utility to compute WALA of MBS, using the 
 * (time-independent) eff orig term (WARM+WALA) and maturity date
 ***************************************************************************/
int EXPORT MbsWala
   (TDate   currDate,             /* (I) date for which WALA needed */
    TDate   mbsMatDate,           /* (I) maturity of MBS */
    long    mbsEffOrigTerm,       /* (I) WARM+WALA of mbs, in months */
    long   *wala);                /* (O) WALA, in months */

/***************************************************************************
 * MbsWarm
 * Utility to compute WARM of MBS, using the 
 * (time-independent) maturity date
 ***************************************************************************/
int EXPORT MbsWarm
   (TDate   currDate,             /* (I) date for which WARM needed */
    TDate   mbsMatDate,           /* (I) maturity of MBS */
    long   *warm);                /* (O) WARM, in months */


/***************************************************************************
 *  GetTime
 *  Copies string containing current time
 ***************************************************************************/
int EXPORT GetTime
   (char   *currTime);          /* (O) string w/current time */

/***********************************************************
 *  GetHolidayFileName
 *  Determines location/name of holiday file
 ***********************************************************/
int GetHolidayFileName
   (char   *holidayFile);       /* (O) full name (w/path) of holiday file */

/***********************************************************
 *  make_TDate
 *  Utility to set TDate to explicit month/day/year
 ***********************************************************/
int make_TDate(long month,	/* (I) 1-12 */
	       long day,		/* (I) 1-31 */
	       long year,	/* (I) should be 4-digit yyyy form */
	       TDate *date);	/* (O) set to corresponding TDate */

/***********************************************************
 *  TDateOf
 *  Converts YYYYMMDD date to GTO TDate
 ***********************************************************/
int TDateOf(long idate,         /* (I) YYYYMMDD date */
	    TDate *Date);       /* (O) TDate */

/***********************************************************
 *  IDateOf
 *  Converts GTO TDate to YYYYMMDD date
 ***********************************************************/
int IDateOf(TDate date,
	    long *iDate);

/***************************************************************************
 *  FindDateInDatelist()
 *  Determines if specified date is in date list, and if so,
 *  returns index in date list
 ***************************************************************************/
int FindDateInDatelist
   (TDate       date,           /* (I) date to find in date list */
    TDateList  *dateList,       /* (I) datelist to search */
    TBoolean   *dateFound,      /* (O) TRUE if date found */
    long       *dlIdx           /* (O) (if found) index of this date
                                 * in datelist; -1 if not found */
    );

/************************************************************************
 * AdvDate()
 * Advances given date by nIntervals
 ************************************************************************/
int AdvDate(TDate *date,                  /* (I/O) date to alter */
            long nIntervals,               /* (I) # intervals */
            TDateInterval *dateInterval); /* (I) (single) interval */


/************************************************************************
 * NxtMth()
 * Advances date by nMonths months;
 * Note: OK to use on same date (e.g., NxtMth(date,&date))
 * Note: if day-of-month beyond last valid day in new month,
 * function sets day-of-month to last valid day;
 * Note also that function tries to preserve the day-of-month
 * of startDate--if we start at Jan 31, and advance 2 months
 * (i.e., past Feb), the result would be Mar 31
 ************************************************************************/
int NxtMth(TDate startDate,
           long nMonths,
           TDate * endDate
           );

/****************************************************************************
 * IDateNxtMth
 * Computes new YYYYMMDD date nMonths fwd/back from start date
 ****************************************************************************/
int EXPORT IDateNxtMth
   (long     startIDate,        /* (I) start date (YYYYMMDD) */
    long     nMonths,           /* (I) # months to advance (or bkwd) */
    TBoolean preserveEndOfMon,  /* (I) TRUE: if start date is last day of mon
                                 * (e.g., 2/28), then ensure that newIDate
                                 * is also last day of new month (e.g., 4/31);
                                 * else, if FALSE, always preserve day-of-month*/
    long    *newIDate);         /* (O) new date (YYYYMMDD)  */

/************************************************************************
 * MonthNumOfTDate
 * Computes "monthnum" from TDate: 12*year+month
 * (useful in some monthly date arithmetic)
 ************************************************************************/
int MonthNumOfTDate(TDate    date,      /* (I) input date */
                    long    *monthNum); /* (O) monthnum: 12*yr+month */


/************************************************************************
 * TDateOfMonthNum
 * Computes TDate from "monthnum" (12*year+month),
 * where lost day-of-month information is supplied
 * as additional input
 ************************************************************************/
int TDateOfMonthNum(long     monthNum, /* (I) monthnum (12*yr+mon) */
                    long     dom,      /* (I) day-of-mon (1-31) to
                                        * which to set output date */
                    TDate   *date);    /* (O) output date */


/************************************************************************
 * MonsDiff()
 * Computes the number of months between two dates (date2 - date1)
 * (day-of-month for dates ignored)
 ************************************************************************/
int MonsDiff(TDate date1,       /* (I) date */
             TDate date2,       /* (I) date */
             long *nmons);      /* (O) difference, in months */


/************************************************************************
 * GetDOM()
 * Returns the day-of-month of a TDate
 ************************************************************************/
int GetDOM(TDate  date,         /* (I) date */
           long  *dom);         /* (O) day-of-month (1-31)*/


/************************************************************************
 * SetDOM()
 * Sets the day-of-month of a TDate to specified value
 ************************************************************************/
int SetDOM(long dom,            /* (I) day-of-month (1-31)*/
           TBoolean endOfMonth, /* (I) determines handling of a day-of-month
                                 * later than last day in mon (e.g., 
                                 * dom=31 for February):
                                 * TRUE=set dom to last valid day,
                                 * FALSE=declare error */
           TDate * date);       /* (I/O) date to alter */

/************************************************************************
 * ForceToClosestDOM()
 * Forces the date to the nearest date having the desired day-of-month
 ************************************************************************/
int ForceToClosestDOM
   (long desiredDom,            /* (I) desired day-of-month (1-31)*/
    TBoolean endOfMonth,        /* (I) determines handling of a day-of-month
                                 * later than last day in mon (e.g., 
                                 * dom=31 for February):
                                 * TRUE=set dom to last valid day,
                                 * FALSE=declare error */
    TDate * date);              /* (I/O) date to alter */

/************************************************************************
 * GetMOY()
 * Returns the month-of-year of a TDate
 ************************************************************************/
int GetMOY(TDate  date,         /* (I) date */
           long  *moy);         /* (O) month-of-year (1-12)*/


/***************************************************************************
 *  Misc utility functions
 ***************************************************************************/

/***************************************************************************
 *  CopyDouble
 *  Copies array of doubles from one array to another
 ***************************************************************************/
int EXPORT CopyDouble
   (long      nValues,          /* (I) # values to copy */
    double   *inpArray,         /* (I) input array */
    double   *outArray);        /* (O) output array */


/***************************************************************************
 *  Misc GTO-interface utility functions
 ***************************************************************************/

/***************************************************************************
 *  SetSingleFloatArray
 *  Constructor for TFloatRateArray to define a single
 *  floating rate (i.e., weight=1.0)
 ***************************************************************************/
int EXPORT SetSingleFloatArray
   (long     cpnsPerYear,           /* (I) cpns per year */
    long     matInMonths,           /* (I) maturity in months */
    long     dayCountConv,          /* (I) day count convention */
    long     curveIndex,            /* (I) which curve (e.g.,
                                       GTO_CURVE_DISCOUNT) */
    long     numSettleDays,         /* (I) (aka spotOffsetDays) # bus. days
                                     * between trade & settle for this rate;
                                     * equiv. to # bus. days from observation of 
                                     * rate to day on which it is true spot */
    TFloatRateArray **floatIndex);  /* (O) ptr to new TFloatRateArray */



/***************************************************************************
 *  Rate funcs
 ***************************************************************************/

/********************************************************************************
 *  ResetArmCpn
 *  Function to reset ARM cpn (one reset)
 ********************************************************************************/
int EXPORT ResetArmCpn
   (double    uncappedNewCpn,     /* (I) rate+sprd, before caps */
    double    prevCpn,            /* (I) previous cpn */
    double    pdCapSprd,          /* (I) sticky cap spread */
    double    pdFlrSprd,          /* (I) sticky floor spread */
    double    lifeCap,            /* (I) life cap on cpn */
    double    lifeFlr,            /* (I) life floor on cpn */
    double    roundingMultiple,   /* (I) rounding multiple (e.g., 0.00125) */
    TBoolean  useSmoothing,       /* (I) TRUE to use smoothing of MAX/MIN,
                                   * FALSE to use ordinary (discontin) MAX/MIN */
    double   *maxDiffs,           /* (I) Max diff from other nodes */
    double   *newCpn);            /* (O) new (capped/rounded) cpn */

/***********************************************************************
 * GetIndexRate
 * Given rate definition, reset date and zero curves, computes
 * index rate (incl wgt/spread, if any)
 ***********************************************************************/
int EXPORT GetIndexRate
   (TDate             effResetDate,     /* (I) eff. reset date */
    TFloatRateArray  *indexRateDefn,    /* (I) defn of index rate */
    TCurve           *discZeroCurve,    /* (I) disc zero curve */
    TCurve          **indexZeroCurves,  /* (I) index zero curves */
    char             *holidayFile,      /* (I) holiday file */
    long              interpType,       /* (I) interp for zero curve;
                                         * e.g., GTO_LINEAR_INTERP */
    long              stubType,         /* (I) for coupon-bearing;
                                         * e.g., GTO_STUB_BOND */
    TBoolean          stubAtEnd,        /* (I) for coupon-bearing;
                                         * eg., FALSE */
    double           *indexRate);       /* (O) index rate (w/wgt,sprd) */


/***************************************************************************
 *  Misc funcs
 ***************************************************************************/

/***************************************************************************
 *  MbsCopyDArray
 *  Creates copy of double array
 *  Sets copyArray to NULL if no values or orig is NULL
 *  Special behavior: if origArray is NULL, but nValues > 0
 *  (user wants array copy, but no input array exists),
 *  create output array of size nValues and fill w/def value
 ***************************************************************************/
int EXPORT MbsCopyDArray
   (double  *origArray,         /* (I) original array to copy from */
    long     nValues,           /* (I) # elements to fill in copyArray */
    double   defValue,          /* (I) default value (used if origArray=NULL)*/
    double **copyArray);        /* (O) copy */

/***************************************************************************
 *  MbsCopyLArray
 *  Creates copy of long int array
 *  Sets copyArray to NULL if no values or orig is NULL
 ***************************************************************************/
int EXPORT MbsCopyLArray
   (long     nValues,           /* (I) # elements in origArray */
    long    *origArray,         /* (I) original */
    long   **copyArray);        /* (O) copy */



/***************************************************************************
 *  Fundamental memory alloc funcs
 ***************************************************************************/

/***********************************************************
 *  mbs_calloc
 *  Central "calloc" function
 *  Examples:
 *       XONFAIL( mbs_calloc(2,1000,(char**) &somevec) );
 *    XONFAILMSG( mbs_calloc(2,1000,(char**) &somevec), "failed to alloc somevec" );
 ***********************************************************/
int mbs_calloc(long nvars,	     /* (I) # vars/objs to alloc */
	       long varsize, 	     /* (I) size of each var/obj */
	       char** ptr);	     /* (O) ptr to allocated var/obj */

/**************************************************
 *  mbs_malloc
 *  Central "malloc" function
 *  Examples:
 *       XONFAIL( mbs_malloc(1000,(char**) &somevec) );
 *    XONFAILMSG( mbs_malloc(1000,(char**) &somevec), "failed to alloc somevec" );
 **************************************************/
int mbs_malloc(long varsize,	     /* (I) size of var/obj */
	       char** ptr);	     /* (O) ptr to allocated var/obj */

/**************************************************
 *  mbs_free
 *  Central memory de-allocation function
 **************************************************/
void mbs_free(char* ptr);

/***************************************************************************
 *  Port of Numerical Recipes memory functions
 ***************************************************************************/
int mbs_vector(long nl, long nh, float** p);
int mbs_ivector(long nl, long nh, long **p);
int mbs_dvector(long nl, long nh, double **p);
int mbs_matrix(long nrl, long nrh, long ncl, long nch, float ***p);
int mbs_dmatrix(long nrl, long nrh, long ncl, long nch, double ***p);
int mbs_imatrix(long nrl, long nrh, long ncl, long nch, long ***p);
void mbs_free_vector(float *v, long nl, long nh);
void mbs_free_ivector(long *v, long nl, long nh);
void mbs_free_dvector(double *v, long nl, long nh);
void mbs_free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void mbs_free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void mbs_free_imatrix(long **m, long nrl, long nrh, long ncl, long nch);


/***************************************************************************
 *  Misc. string functions
 ***************************************************************************/

/******************************************************************************
 * remove_blanks
 * Removes all blanks (or tabs) from string
 ******************************************************************************/
int remove_blanks(char *s);

/******************************************************************************
 * remove_lead_blanks
 * Removes leading blanks (or tabs) from string
 ******************************************************************************/
int remove_lead_blanks(char *s);

/******************************************************************************
 * remove_trail_blanks
 * Removes trailing blanks (or tabs) from string
 ******************************************************************************/
int remove_trail_blanks(char *s);

/******************************************************************************
 * remove_end_blanks
 * Removes leading and trailing blanks (or tabs) from string
 ******************************************************************************/
int remove_end_blanks(char *s);

/************************************************************************************
 * remove_trail_cr
 * Removes trailing carriage return from string;
 ************************************************************************************/
int remove_trail_cr(char *s);

/******************************************************************************
 * dncase
 * Simple function to convert string to lowercase
 * (created this local routine in case we need special handling)
 ******************************************************************************/
int dncase(char *s);

/******************************************************************************
 * upcase
 * Simple function to convert string to uppercase
 * (created this local routine in case we need special handling)
 ******************************************************************************/
int upcase(char *s);

/******************************************************************************
 * str_npos
 * Finds position of nth occurence of char c in string
 * (or -1 if not found)
******************************************************************************/
int str_npos(char *s,		/* (I) string to search */
	     char c,		/* (I) char to find */
             long n,             /* (I) occurence to search for */
	     long *ipos);	/* (O) position of char */

/******************************************************************************
 * str_pos
 * Finds position of 1st occurence of char c in string
 * (or -1 if not found)
******************************************************************************/
int str_pos(char *s,		/* (I) string to search */
	    char c,		/* (I) char to find */
	    long *ipos);		/* (O) position of char */

/******************************************************************************
 * str_match
 * Determines if strings match (case-sensitive)
 * NB: differs from strcmp() in that it can handle strings of
 * different length--compares shorter string to the same-length
 * starting piece of longer string
******************************************************************************/
int str_match(char *s1,		/* (I) string 1 */
	      char *s2,		/* (I) string 2 */
	      long *match);	/* (O) boolean true if matching */

/******************************************************************************
 * parse_bond_price
 * Converts string representation of bond price to floating point #;
 * works with either "nn-mm+" (tick format) or "nn.mm" (decimal) forms
 * Returns boolean (true if successful)
******************************************************************************/
int parse_bond_price(char *sz,	 /* (I) string to parse */
		     double *x); /* (O) converted decimal value */


/***************************************************************************
 *  "Perl" data-file parsing functions
 ***************************************************************************/

/*********************************************************************************
 *  read_perl_long
 *  Gets long int data from matching line in "perl" data file
 *********************************************************************************/
int read_perl_long(FILE *f,           /* (I) file (already open) to read from */
		   char *perl_lbl,    /* (I) label of data line to get */
		   long nvalues,       /* (I) # of values to get */
		   long *outarr);     /* (O) values copied to this array */

/*********************************************************************************
 *  read_perl_dbl
 *  Gets double data from matching line in "perl" data file
 *********************************************************************************/
int read_perl_dbl(FILE *f,           /* (I) file (already open) to read from */
	          char *perl_lbl,    /* (I) label of data line to get */
		  long nvalues,       /* (I) # of values to get */
		  double *outarr);   /* (O) values copied to this array */

/*********************************************************************************
 *  read_perl_str
 *  Gets strings from matching line in "perl" data file
 *  NB: strings are copied consecutively to outstr, 
 *  separated by null chars
 *********************************************************************************/
int read_perl_str(FILE *f,           /* (I) file (already open) to read from */
		  char *perl_lbl,    /* (I) label of data to get */
		  long nvalues,       /* (I) # of strings to get */
                  long maxlen,        /* (I) # chars to copy to outstr */
                  long trim_ends,     /* (I) boolean: if true, removes
                                        lead/trail blanks from strings */
		  char *outstr);     /* (O) values copied to this string */

#ifdef __cplusplus
}
#endif

#endif
