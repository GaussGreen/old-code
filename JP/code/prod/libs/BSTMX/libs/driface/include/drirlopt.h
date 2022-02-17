/****************************************************************
 * Module:	DRI
 * Submodule:	
 * File:	drirlopt.h
 * Function:	DR Wrapper for Rate Lock Option.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drirlopt_H
#define	_drirlopt_H
#include "drlstd.h"

#include <stdio.h>
#include "drlsmat.h"		/* TSwaptionMatrix2D */
#include "drlts.h"		/* TCurve */
#include "dritkwrp.h"		/* TSmile3Data */

extern	int	DriRLockOption(
	TDate rateObsDate,		/* (I) notification date */
	TDate rateEffDate,		/* (I) start date */

	TDateInterval rateMat1,	/* (I) rate maturity interval */
	int rateFreq1,			/* (I) rate frequency */
	TDayCount rateDcc1,		/* (I) day count convention */
	double rateWeight1,		/* (I) rate weight */

	TDateInterval rateMat2,	/* (I) rate maturity interval */
	int rateFreq2,			/* (I) rate frequency */
	TDayCount rateDcc2,		/* (I) day count convention */
	double rateWeight2,		/* (I) rate weight */

	TDate payDate,			/* (I) payment date */

	TDateInterval annMat,		/* (I) payment maturity interval */
	int annFreq,			/* (I) payment frequency */

	double strikeRate,		/* (I) strike rate for options */
	char *callPutType,		/* (I) "C"all, "P"ut, "F"orward */

	TDate todayDate,		/* (I) todayDates's date */
	TCurve *discZcCurve,		/* (I) discount zero curve */
	TCurve *indZcCurve1,		/* (I) index zero curve rate 1 */
	TCurve *indZcCurve2,		/* (I) index zero curve rate 2 */
	TSwaptionMatrix2D *swMat,	/* (I) CMS swaption matrix */
	TSmile3Data *smlData,		/* (I) smile data (or NULL) */
	VnfmData *vnfmData,		/* (I) calibration info */
	int noAdj,			/* (I) disable adjustment */

	double *pv);			/* (O) present value */



extern	int	DriRLockOptionW(char *dataFnam);




#endif	/* _drirlopt_H */

