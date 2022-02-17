/*********************************************************************************
 * L_DATE.H 
 * defines London dates
 *
 ********************************************************************************/

#ifndef __L_DATE_H__
#define __L_DATE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <alib/bastypes.h>
#include <alib/dateconv.h>
#include <alib/tcurve.h>
#include <alib/ldate.h>

#define  MAX_DATE           500           /* Max nb of elemts in input date arr */
#define  MAX_VOL            200           /* Max nb of elemts in input date arr */
#define  MAXBUFF            250
#define  MAX_DISC_DIV       2000          /* Max nb of elements in discreet
                                             dividend array                     */
#define  MAX_CONT_DIV       200           /* Max nb of elements in continuous
                                             dividend array                     */
#define  MAX_SETTLE         2000          /* Max nb of elements in settlement
                                             array                              */
#define  MAX_CURVE          5             /* Max nb of curves                   */

/* different mode choices for GetDLOffset function in date.c */
#define  CbkEXACT   0
#define  CbkLOWER   -1
#define  CbkHIGHER  1

/*********************************************************************************
 *  Checks for validity of YYYYMMDD formatted dates  
 *  returns 0 if all is well                        
 *  returns 1 if bad year                           
 *  returns 2 if bad month                          
 *  returns 3 if bad day                            
 *
 ********************************************************************************/
int Dateok(long date);

/*********************************************************************************
 * split date into month, day, year
 *
 ********************************************************************************/

void Dsplit(long   date_i,                /* (I) london date                    */
            long   *mm_o,                 /* (O) month                          */
            long   *dd_o,                 /* (O) day                            */
            long   *yy_o);                /* (O) year                           */
    
/*********************************************************************************
 * check whether year is leap year
 *
 ********************************************************************************/
int Isleap(long year);

/*********************************************************************************
 * calculate the number of days between date1_i and date2_i
 *
 ********************************************************************************/
long Daysact(long date1_i, long date2_i);

/*********************************************************************************
 * Given a DateList and its size
 * sort its contents in 1st: date ascending order, 
 *                      2nd: value ascending order. (if not NULL)
 * Returns SUCCESS or FAILURE
 ********************************************************************************/
int   CrxSortDateList(int    NbDates,      /* (I)   num of dates                 */
                     long  *DateList,     /* (I/O) date list to be sorted       */
                     long  *SuppValue);   /* (I)   what to compare              */

/*********************************************************************************
 * merge two date lists and sort it in ascending order
 * redundant entries are also removed
 *
 ********************************************************************************/
int CrxSortMergeDateList(                    
	int   *NbMerge,                       /* (O) Num of points in merged list   */
	long  *DLMerge,                       /* (O) Merged & sorted date list      */
	int    Nb1,                           /* (I) Num of points in 1st list      */
	long  *DL1,                           /* (I) 1st date list                  */
	int    Nb2,                           /* (I) num of points in 2nd list      */
	long  *DL2);                          /* (I) 2nd date list                  */

/*****  GetDLOffset  ***************************************************/
/*
 *      Given a datelist and a targetDate,
 *      returns the offset of the element nearest to the targetdate
 *
 *      if mode = 
 *          CbkEXACT:  returns -999 if no matching targetdate is found
 *
 *          CbkLOWER:  returns the nearest offset that is LOWER or equal;
 *                     or -999 if all dates > targetdate
 *
 *          CbkHIGHER: returns the nearest offset that is HIGHER or equal;
 *                     or -999 if all dates < targetdate
 *
 *      NOTE:   datelist has to be increasing! 
 *
 *              In case some of the dates in the datelist are the same, 
 *
 *              (i) if targetDate = repeated date (this case is not well defined!)
 *                  then any of the offsets of the repeated date is returned.
 *
 *              (ii) CbkHIGHER: if repeated date is the date to be returned (i.e. just above
 *                   targetDate), then the lowest offset of the repeated date is returned.
 *
 *              (iii) CbkLOWER: if repeated date is the date to be returned (i.e. just below
 *                    targetDate), then the highest offset of the repeated date is returned.
 *
 */
int     GetDLOffset(int         NbDates,
                    long       *DL,
                    long        targetDate,
                    int         mode);


    
/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif

