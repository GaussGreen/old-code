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


DLL_EXPORT(int)

DriOutSpreadStrikeOption( 
		TDate          today,		/* (I) spot date */
		TDrWrapperData *drWrap,         /* (I) DR Wrapper Data  */
		int            optType,   	/* (I) Max, Min, Diff */
		TDate		   optExpiry,
		TDate	       payDate,
		double		   vol,
		char		   distType,
		double		   spot1,
		TDate	       matDate1,
		TDateInterval  freq1,
		double         coupon1,
		long	       dcc1,
		double         spot2,
		TDate	       matDate2,
		TDateInterval  freq2,
		double         coupon2,
		long	       dcc2,
		TEqStatData    **eqStatData,	/* (I) data in equity.sta files	*/
		double         strike,
		double	       *pv);

/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriOutPerformanceIdxOption}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "dritridx_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)

DriOutSpreadStrikeptionW(char *dataFnam);

#endif	/* _pfmopt_H */


