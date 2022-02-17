/****************************************************************
 * Module:	driface
 * File:	dritridx.h
 * Function:	Total return index swap and cap/floor
 * Author:	Julia Chislenko May 1998
 * $Header$
 *****************************************************************/
#ifndef	_dritridx_H
#define	_dritridx_H

#include <stdio.h>

#include "bastypes.h"
#include "dritkwrp.h"		/* TDrWrapperData */


/*---------------------------------------------------------------------
 * Pricing routine for total return index swap/cap/floor
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriTotalRetIndxSwapCapFloor(   
TDate       today,              /* (I) spot date */
TCurve     *discZC,             /* (I) discount zero curve */
TCurve     *indxZC,             /* (I) index zero curve */
char        instrType,   	/* (I) 'S'wap, 'C'ap, 'F'loor */
TSwaptionMatrix2D *swMat,       /* (I) CMS swaption matrix
				 *     can be NULL if swap only */
TCurve     *indxVolCurve,       /* (I) index volatility curve
				 *     can be NULL if swap only  */
double      correlation,        /* (I) between indx and interest rate */    
long        settleDays,         /* (I) Num bus days btw indx reset and paymt*/
long	    busDayConv,		/* (I) Business day conv      */
char       *holidayFile,        /* (I) "NONE" for weekends only
				 *     "No_Weekends" for no adjustments */
double      spotPrice,     	/* (I) Spot indx price 		*/
long        numIdxFixing,       /* (I) Number of past refixings */
TDate      *refixingDates,      /* (I) Refixing dates           */
TDate      *refixingEffDates,   /* (I) Refixing effective dates */
double     *indexFixings,       /* (I) Past refixings           */
long        numResetDates,      /* (I) Number of index reset dates */
TDate      *resetStDates,	    /* (I) Reset st dates > today */
TDate      *resetEndDates,	    /* (I) Reset end dates      */
TDate      *payDates,	        /* (I) Payment dates        */
double     *strikes,            /* (I) [numResetDates] for cap/floor
				 *     NULL if swap only */
double     *pv);		/* (O) present value */


/*---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriIdxForwardSwapCapFloorW}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "dritridt.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriIdxForwardSwapCapFloorW(char *dataFnam);

#endif	/* _dritridx_H */


