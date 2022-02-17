#ifndef ESL_TIME_DOT_H
#define ESL_TIME_DOT_H

/** NOTE: This file should be only included through 'esl_time.c'
 */


#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"
#include "esl_util.h"

#ifdef  __cplusplus
extern "C" {
#endif

/*****  Add_To_DateList  ******************************************************/
/**
*       Add a set of dates to the critical date list.
*/
int     Add_To_DateList (int        *NbCritDate  /** (I/O) Number of dates    */
                        ,CRIT_DATE  **CritDate 	 /** (I/O) Critical date list */
                        ,long       Date         /** (I) New date             */
                        ,int        Type         /** (I) New date type        */
                        ,double	    Value0       /** (I) New date values      */
                        ,double	    Value1       /** (I)                      */
                        ,double	    Value2       /** (I)                      */
                        ,double	    Value3       /** (I)                      */
                        ,double	    Value4       /** (I)                      */
                        ,long       SuppDate0    /** (I) Supplementary dates  */
                        ,long       SuppDate1    
                        ,long       SuppDate2);

/*****  TimeInterp  *********************************************************/
/**
*       Utility routine. Eliminates dates falling before the value date and 
*       interpolate the corresponding value.
*       It assumes dates are entered in ascending order.
*/
int     TimeInterp (long    ValueDate  /** (I) Value date                     */
                   ,char    *Name      /** (I) Array name                     */
                   ,char    Type       /** (I) Type of interpolation          */
                   ,int     *NbDate    /** (I/O) Number of dates in the array */
                   ,long    *Date      /** (I/O) Dates in ascending order     */
                   ,double  *Curve0    /** (I/O) Set of associated curves     */
                   ,double  *Curve1 
                   ,double  *Curve2 
                   ,double  *Curve3 
                   ,double  *Curve4);

/*****  Sort_CritDate  ******************************************************/
/**
*       Sort routine of critical dates. Used in time.c.
*/
int     Sort_CritDate (	int	        n       /** Number of dates to sort */
                       ,CRIT_DATE	*Date   /** Dates to be sorted      */
		);

/*****  DrExtendAmerSch  ****************************************************/
/**
 *       Utility routine to extend an input American date schedule to a full
 *       European style schedule using the 'Increasing time-step' method.
 *       It also does the following:
 *       - eliminates dates STRICTLY BEFORE value date
 *       - interpolates supplementary values (if available)
 *
 *       NOTE: the extension is performed on the NotifDates if they are given,
 *             the SchDates are then interpolated. Otherwise, the extension is
 *             done on the SchDates themselves.
 */

int DrExtendAmerSch(long      ValueDate    /** (I)   Value date          */
                   ,int       Ppy          /** (I)   ppy in tree         */
                   ,int      *NbSchDates   /** (I/O) Nb schedule dates   */
                   ,long    **SchDates     /** (I/O) Schedule dates      */
                   ,long    **NotifDates   /** (I/O) Notif dates         */
                   ,char      InterpStyle  /** (I)   Curve interp style  */
                   ,double  **Curve0       /** (I/O) Curve 0 values      */
                   ,double  **Curve1       /** (I/O) Curve 1 values      */
                   ,double  **Curve2       /** (I/O) Curve 2 values      */
                   ,double  **Curve3       /** (I/O) Curve 3 values      */
                   ,double  **Curve4       /** (I/O) Curve 4 values      */
		);



#ifdef  __cplusplus
}
#endif


#endif



