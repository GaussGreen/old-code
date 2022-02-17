/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	drlkwrap.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlkwrap_H
#define	_drlkwrap_H
#include "drlstd.h"

#include <stdio.h>
#ifdef	DRL_CLIB
#include "bastypes.h"
#endif
#include "drlsmat.h"
#include "drlts.h"



/*t---------------------------------------------------------------
 * A data structure that holds the market information of
 * type II DR wrapper.
 */

typedef	struct	{
	DDate		fToday;			/* today's date */
	DCurve		*fDiscZcCurve;		/* Discount zero curve */
	DCurve		*fZcCurve;		/* Index zero curve */
	DCurve		*fRiskZcCurve;		/* Risky zero curve */
	DCurve		*fBasisZcCurve;		/* Basis zero curve */
	DCurve		*fBvCurve;		/* Base volatility curve */
	DCurve		*fBSVolCurve;		/* Base volatility curve */
	DSwopMat	 *fCmsSwMat;		/* Swaption matrix */
        long            fMMDenom;               /* 360 or 365 */
        DDayCount       fSwDcc;                 /* Swap day count conv */

	double		f1Beta;			/* 1F parameters */
	double		f1Weight;
	int		f1Ppy;

	double		f2Beta1;		/* 2F parameters */
	double		f2Beta2;
	double		f2Weight1;
	double		f2Weight2;
	double		f2Corr12;
	int		f2Ppy;

	double		f3Beta1;		/* 3F parameters */
	double		f3Beta2;
	double		f3Beta3;
	double		f3Weight1;
	double		f3Weight2;
	double		f3Weight3;
	double		f3Corr12;
	double		f3Corr13;
	double		f3Corr23;
	int		f3Ppy;

} DDrWrapData;
/*e*/


#define	DRL_DRW_TYPE2_2CURVES		(0x0001L)
#define	DRL_DRW_TYPE2_3CURVES		(0x0002L)
#define	DRL_DRW_BASIS			(0x0003L)

extern	int	DrlDDrWrapDataImport(
	char *pathdir,
	long options,
	DDrWrapData **that);

extern	int	DrlDDrWrapDataImportType2(
	char *pathdir,
	DDrWrapData **that);

extern	int	DrlDDrWrapDataFpWrite(
	DDrWrapData *drwData,
	FILE *fp);


extern	int	DrlDDrWrapDataExport(
	char *pathdir,
	long options,
	DDrWrapData *drwData);


extern	int	DrlDDrWrapDataFree(DDrWrapData *that);


extern	int	DrlDDrWrapDataPutPrice(double value);

extern	int	DrlDDrWrapDataGetInterpVol(
	DDrWrapData *drWrapper,	/* (I) wrapper data */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	DDate calibMatDate,		/* (I) final mat (used if final) */
	DInterval calibMatInt,	/* (I) fwd mat (used if cms) */
	int *numVols,			/* (O) # of vol pnts (or NULL) */
	DDate **volDates,		/* (O) vol exp dates (or NULL) */
	double **volMat,		/* (O) vol fwd mat (or NULL) */
	int **volFreq,			/* (O) vol freq (or NULL) */
	double **volRates);		/* (O) vol (can be NULL) */


extern	DCurve*	DrlDDrWrapDataGetDCurve(
	DDrWrapData *drWrapData,	/* (I) wrapper data */
	char zcType);			/* (I) (D)isc,(Z)idx,(R)isky,(B)asis */


#endif	/* _drlkwrap_H */

