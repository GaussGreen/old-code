/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	dritkwrp.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_dritkwrp_H
#define	_dritkwrp_H
#include "drlstd.h"

#include <stdio.h>
#include "bastypes.h"
#include "drlsmat.h"
#include "drlts.h"
#include "vnfmanly.h"



/*t-@CDOC(idxn=TDrWrapperData,catn=structdef)---------------------
 * A data structure that holds the market information of
 * type II DR wrapper.
 */

typedef struct  {
    TDate       fToday;             /* today's date */
    TCurve      *fDiscZcCurve;      /* Discount zero curve */
    TCurve      *fZcCurve;          /* Index zero curve */
    TCurve      *fRiskZcCurve;      /* Risky zero curve */
    TCurve      *fBasisZcCurve;     /* Basis zero curve */
    TCurve      *fBvCurve;          /* Base volatility curve */
    TCurve      *fBSVolCurve;       /* Base volatility curve */
    TSwaptionMatrix2D *fCmsSwMat;   /* Swaption matrix */
    long        fMMDenom;           /* 360 or 365 */
    TDayCount   fSwDcc;             /* Swap day count conv */

    double      f1Beta;             /* 1F parameters */
    double      f1Weight;
    int         f1Ppy;

    double      f2Beta1;            /* 2F parameters */
    double      f2Beta2;
    double      f2Weight1;
    double      f2Weight2;
    double      f2Corr12;
    int         f2Ppy;

    double      f3Beta1;            /* 3F parameters */
    double      f3Beta2;
    double      f3Beta3;
    double      f3Weight1;
    double      f3Weight2;
    double      f3Weight3;
    double      f3Corr12;
    double      f3Corr13;
    double      f3Corr23;
    int         f3Ppy;

    double      LQ;                 /* default 2Q parameters*/
    double      RQ;
    double      FwdShift;
    int         nbIter;

} TDrWrapperData;
/*e*/


#define	DRI_DRW_TYPE2_2CURVES		(0x0001L)
#define	DRI_DRW_TYPE2_3CURVES		(0x0002L)
#define	DRI_DRW_BASIS			(0x0003L)

extern	DLL_EXPORT(int)	DriTDrWrapperDataGetFull(
	char *pathdir,
	long options,
	TDrWrapperData **that);

extern	DLL_EXPORT(int)	DriTDrWrapperDataGet(
	char *pathdir,
	TDrWrapperData **that);

extern	DLL_EXPORT(int)	DriTDrWrapperDataGetLiverate(
	const char *pathdir,		/* (I) tmp direcory name (or NULL) */
	const char *curCode,		/* (I) currency code */
	int mmDen,			/* (I) money market denominator */
	int bvFreq,			/* (I) base vol frequency */
	int swapFreq,			/* (I) swap frequency */
	TDayCount swapDcc,		/* (I) swap day count convention */
	const char* yldcrvBnam,		/* (I) yield curve file basename */
	const char* baseAtmBnam,	/* (I) ATM base vol file basename */
	const char* swoAtmBnam,		/* (I) ATM swaption vol file basename */
	TDrWrapperData **that);		/* (O) dr wrapper data */


extern	DLL_EXPORT(int)	DriTDrWrapperDataFree(TDrWrapperData *that);


extern	DLL_EXPORT(int)	DriTDrWrapperDataFpWrite(
	TDrWrapperData *drwData,	/* (I) dr wrapper data */
	FILE *fp);			/* (I) FILE to write to (or NULL) */

extern	DLL_EXPORT(int)	DriTDrWrapperDataPutPrice(double value);

extern	DLL_EXPORT(VnfmData*)	DriTDrWrapperDataGetVnfmData(
	TDrWrapperData *drWrapper,	/* (I) wrapper data */
	int numFact,			/* (I) number of factors */
	double backBoneQ,		/* (I) back bone q */
	int datesFrom);			/* (I) dates: 0=flat, 1=bv, 2=swvol */

extern	DLL_EXPORT(int)	DriTDrWrapperDataGetInterpVol(
	TDrWrapperData *drWrapper,	/* (I) wrapper data */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	TDate calibMatDate,		/* (I) final mat (used if final) */
	TDateInterval calibMatInt,	/* (I) fwd mat (used if cms) */
	int *numVols,			/* (O) # of vol points (can be NULL) */
	TDate **volDates,		/* (O) vol exp dates (can be NULL) */
	double **volMat,		/* (O) vol fwd mat (can be NULL) */
	int **volFreq,			/* (O) vol freq (can be NULL) */
	double **volRates);		/* (O) vol (can be NULL) */


extern	DLL_EXPORT(TCurve*)	DriTDrWrapperDataGetTCurve(
	TDrWrapperData *drWrapData,	/* (I) wrapper data */
	char zcType);			/* (I) (D)isc,(Z)idx,(R)isky,(B)asis */




/*t-@CDOC(idxn=TSmile3Data,catn=structdef)---------------------
 * Data structure for 3-pts swaption volatility cube.
 * The "smileType" field determines the type of smile function used:
 * \begin{itemize}
 * \item if "smileType" is 0, no smile function is computed.
 * \item if "smileType" is 1, a 2-points smile power function is
 * computed using $V_{+}$, $h$, $q$ and $\alpha$.
 * \item if "smileType" is 2, a 3-points smile power/parabolic function is
 * computed using $V_{-}$, $V_{+}$, $h$, $q$ and $\alpha$.
 * \end{itemize}
 */
typedef	struct	{
	int			fSmileType;	/* 0=no,1=2pts,2=3pts */
	TSwaptionMatrix2D	*fQcoeff;	/* q coeff matrix */
	TSwaptionMatrix2D	*fVhicoeff;	/* v_hi coeff matrix */
	double			fHcoeff;	/* h coefficient */
	double			fAlpha;		/* alpha coefficient */
	TSwaptionMatrix2D	*fBackBoneQ;	/* back bone coeff matrix */
	TSwaptionMatrix2D	*fAtmRate;	/* atm rate (for back bone) */
} TSmile3Data;

extern	TSmile3Data*	DrlTSmile3DataNewEmpty(void);
extern	void		DrlTSmile3DataDelete(TSmile3Data *that);
extern	int		DrlTSmile3DataCleanup(TSmile3Data *that);
extern	int		DrlTSmile3DataFpWrite(TSmile3Data *that, FILE *fp);

extern	int		DrlTSmile3DataReadFromExport(
	TSmile3Data *that,		/* (B) smile data */
	char *curCode,			/* (I) currency code (USD,etc)  */
	char *dirname,			/* (I) directory  */
        TDate baseDate);		/* (I) reference date */

extern	int	DrlTSmile3DataInterp(
	TCurve *zcCurve,		/* (I) zero curve */
	TSmile3Data *that,		/* (I) smile data */
	TDate baseDate,			/* (I) base date */
	double strikeRate,		/* (I) strike rate */
	int freq,			/* (I) rate freq (0,1,2,4,12) */
	long dayCountConv,		/* (I) rate day count conv */
	TDate volDate,			/* (I) vol exp date */
	double volMat,			/* (I) vol fwd maturity */
	int volFreq,			/* (I) vol frequency */
	double volRate,			/* (I) vol */
	double *volAdj);		/* (O) vol smile adjustment */

extern	int	DrlTSmile3DataInterpVolCurve(
	TSmile3Data *that,		/* (I) smile data */
	TCurve *zcCurve,		/* (I) zero curve */

	int freq,			/* (I) rate freq (0,1,2,4,12) */
	long dayCountConv,		/* (I) day count convention */

	int numVolDates,		/* (I) # of volatility dates */
	TDate *volDates,		/* (I) array of vol dates */
	double *strikeRates,		/* (I) array of strike rate */
	double *volMat,			/* (I) array of vol fwd maturities */
	int *volFreq,			/* (I) array of vol frequency */
	double *volRates,		/* (I/B) array of vol */
	double *volAdj);		/* (O) array of vol adj (or NULL)*/



#endif	/* _dritkwrp_H */

