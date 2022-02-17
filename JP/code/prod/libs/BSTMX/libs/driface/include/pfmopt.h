/****************************************************************
 * Module:	driface
 * File:	pfmopt.h
 * Function:	Out performance options on indices 
 * Author:	David Liu, Nov. 1998
 *****************************************************************/
#ifndef	_pfmopt_H
#define	_pfmopt_H

#include <stdio.h>

#include "bastypes.h"
#include "dritkwrp.h"		/* TDrWrapperData */
#include "drieq.h"		/* TEqStatData */


/*------------------------------------------------------------
  FUNCTION:       DriOutPerformanceIdxOption
 
  CREATED BY:     David Liu  Sept, 1998
 
  PURPOSE:        Outperformance Option on two indices.
 
		  Max, Min, Diff on two indices.
*/
DLL_EXPORT(int)
DriOutPerformanceIdxOption(   
TDate       	today,  	/* (I) spot date */
TDrWrapperData 	*drWrap,        /* (I) DR Wrapper data  */
int         	optType,   	/* (I) 1 = Max 
				 *     2 = Min 
				 *     3 = Diff  
				 *     4 = Min-Min
                                 *     5 = Spread   */
TCurve	   	**indxVolCurve,	/* (I) Index vol curves  */
double	   	corr12,		/* (I) Correlation between two indices */
TEqStatData 	**eqStatData,	/* (I) data in equity.sta files	*/
char       	*holidayFile,   /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
long	   	busDayConv,	/* (I) Business day convention  */
double     	*currIdxPrice,  /* (I) Spot indx level */
double		*lastIdxPrice,	/* (I) Initial idx level */
long        	numResetDates,  /* (I) Number of index reset dates */
TDate      	*resetDate,	/* (I) Observation date > today */
TDate		*payDate,	/* (I) Payment date >= Observation date */
double     	strike,         /* (I) strike  */
double     	spread,         /* (I) spread added to fwdIdx1  */
double     	*pv); 		/* (O) present value */


/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriOutPerformanceIdxOption}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "dritridx_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriOutPerformanceIdxOptionW(char *dataFnam);

#endif	/* _pfmopt_H */


