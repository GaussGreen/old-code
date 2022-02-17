/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	driwbsp.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_driwbsp_H
#define	_driwbsp_H
#include "drlstd.h"

#include <stdio.h>
#include "drlsmat.h"
#include "drlts.h"

extern	int	DriSpreadBarrierBinary(
	TDate expDate,		/* (I) expiration date */
	TDate resetDate,	/* (I) rate reset date */

	TDateInterval rateMat,	/* (I) rate maturity interval */
	int rateFreq,		/* (I) rate frequency */
	TDayCount rateDcc,	/* (I) rate day count convention */

	double strike,		/* (I) used if asset */
	double barrier,		/* (I) barrier level */

	char *callPutS,		/* (I) "C"all, "P"ut */
	char *upDownS,		/* (I) "U"p, "D"own */
	char *inOutS,		/* (I) "I"n, "O"ut */
	char *assetCashS,	/* (I) "A"sset, "C"ash */
	char *hitExpS,		/* (I) "H"it, "E"xp */

	double spdVol,		/* (I) spread volatility in bp/day */
	char *distTypeS,	/* (I) "L"ognormal, "N"ormal */

	TDate today,		/* (I) todays's date */
	TCurve *discZcCurve,	/* (I) discount zero curve */
	TCurve *indxZcCurve,	/* (I) cmt zero curve */
	TSwaptionMatrix2D *swMat,/* (I) CMS swaption matrix */
	double *pv);		/* (O) present value */


extern	int	DriSpreadBarrierBinaryW(char *dataFnam);


#endif	/* _driwbsp_H */

