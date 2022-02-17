/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	drlts.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlts_H
#define	_drlts_H
#include "drlstd.h"

#include <stdio.h>

/*----------------------------------------------------------------------
 * Term Structure (zero coupons, base vol, etc ...)
 */

#if !defined(DRL_CLIB)
/*
 * Using FIX3-like curve
 */

#define  TCMAXNBDATE     500    /* Max nb of elements in input date array */

/*t---------------------------------------------------------------------
 * Data type for zero coupon curve (but can also
 * be used for e.g. volatilities, etc.). 
 */
typedef struct
{
    /* Base date information */
    DDate   Today;                  /* Today's date */
    int     SpotDays;               /* Spot days    */ 
    DDate   ValueDate;              /* Value date   */
    
    /* Zero curve */
    char    CurveDCC;		    /* DCC for interp */
    char    CurveFreq;              /* rate compounding */
    int     NbZero;                 /* Number of zeros                   */ 
    DDate   ZeroDate[TCMAXNBDATE];  /* Zero maturity dates               */ 
    double  Zero[TCMAXNBDATE];      /* Zero rates (Annual ACT/365 basis) */ 
    
} DCurve;
/*e*/

#define	DRLTCURVE_BASEDATE(object)	((object)->ValueDate)
#define	DRLTCURVE_LASTDATE(object)	((object)->ZeroDate[(object)->NbZero-1])
/*#define	DRLTCURVE_BASIS(object)		((object)->CurveFreq)*/
#define	DRLTCURVE_DCC(object)		((object)->CurveDCC)
#define	DRLTCURVE_NUMITEMS(object)	((object)->NbZero)
#define	DRLTCURVE_DATE(object, index)	((object)->ZeroDate[(index)])
#define	DRLTCURVE_RATE(object, index)	((object)->Zero[(index)])

#else
/*
 * Analytics library curve
 */
# include "tcurve.h"	/* C analytics library */
typedef	TCurve	DCurve;

#define	DRLTCURVE_BASEDATE(object)	((object)->fBaseDate)
#define	DRLTCURVE_LASTDATE(object)	((object)->fArray[(object)->fNumItems-1].fDate)
/*#define	DRLTCURVE_BASIS(object)		((object)->fBasis)*/
#define	DRLTCURVE_DCC(object)		((object)->fDayCountConv)
#define	DRLTCURVE_NUMITEMS(object)	((object)->fNumItems)
#define	DRLTCURVE_DATE(object, index)	((object)->fArray[(index)].fDate)
#define	DRLTCURVE_RATE(object, index)	((object)->fArray[(index)].fRate)

#endif



/*
 * I/O
 */

#define	DRL_TCURVE_FMT_STD	(1L)
#define	DRL_TCURVE_FMT_CUR	(2L)
#define	DRL_TCURVE_FMT_PERCENT	(3L)
#define	DRL_TCURVE_FMT_LON	(4L)
#define	DRL_TCURVE_FMT_WRAPZC	(4L)	/* wrapper zero.dat */
#define	DRL_TCURVE_FMT_WRAPBV	(6L)	/* wrapper basevol.dat */
#define	DRL_TCURVE_FMT_REGTEST	(5L)


extern	DCurve*	DrlDCurveNew(DDate valueDate, int numItems,
				double freq, DDayCount dayCountConv);
extern	DCurve*	DrlDCurveNewCopy(DCurve* that);
extern	int	DrlDCurveFree(DCurve* that);	

extern	int	DrlDCurveFileRead(DCurve **that, char *, long fmt);
extern	int	DrlDCurveFileWrite(DCurve *that, char *, long fmt);
extern	int	DrlDCurveFpRead(DCurve **that, FILE *fp, long fmt);
extern	int	DrlDCurveFpWrite(DCurve *that, FILE *fp, long fmt);

extern	int	DrlDCurveLondonBaseVolFpRead(DCurve **that,
				DDate baseDate, FILE *fp);
extern	int	DrlDCurveLondonBaseVolFileRead(DCurve **that,
				DDate baseDate, char *fnam);
extern	char*	DrlDCurvePrintYields(DCurve *that, char *s);

extern	int	DrlDCurveWrapRead(DCurve **that, long *refDate,
				long *zDate, double *zRate);
extern	int	DrlDCurveWrapWrite(DCurve *that, long *refDate,
				long *zDate, double *zRate);


extern	int DrlDCurveFreqSet(
	DCurve *that,		/* (O) zero curve */
	int freq);		/* (I) basis (0, 1, 2, 4, 12) */

extern	int DrlDCurveFreq(
	DCurve *that,		/* (I) zero curve */
	int *freq);		/* (O) basis (0, 1, 2, 4, 12) */

extern	int DrlDCurveInterp(
	DCurve* that,		/* (I) zero curve */
	DDate matDate,		/* (I) reset date */
	double *discRate);	/* (O) forward yield */

extern	int	DrlDCurveDiscFact(
	DCurve* that,		/* (I) zero curve */
	DDate matDate,		/* (I) reset date */
	double *discFact);	/* (O) forward yield */

/* Forward rate calculation */

extern	int	DrlDCurveForwardRate(
	DCurve* that,		/* (I) zero curve */
	DDate startDate,	/* (I) start date */
	DDate maturityDate,	/* (I) maturity date */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield);		/* (I) forward yield */


extern	int	DrlDCurveForwardRate2(
	DCurve* that,		/* (I) zero curve */
	DDate startDate,	/* (I) reset date */
	DInterval maturity,	/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield);		/* (I) forward yield */

extern	int	DrlDCurveForwardRate3(
	DCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	DDayCount dayCountConv,	/* (I) day count convention */
	double *yield);		/* (O) forward rate */


/*
 *
 */
extern	DDate	DrlDCurveRefDate(DCurve *that);
extern	int	DrlDCurveDatesArray(DCurve *that, int *nDates, DDate **dates);




#endif	/* _drlts_H */

