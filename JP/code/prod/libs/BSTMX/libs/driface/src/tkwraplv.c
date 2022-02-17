/************************************************************************
 * Module:	DRIIFACE
 * Function:	DRWrapper Utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>

#include "date_sup.h"
#include "gtonpi.h"

#include "drlio.h"
#include "drlsmat.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlptable.h"		/* TCSheet */
#include "drlvtype.h"
#include "drltime.h"		/* DrlTDateScanYMD */

#include "dritkwrp.h"	/* prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__

static	int	generateZcFromYldcrv(
	TCSheet *yldcrvlc,	/* (I) "yldcrv liverate format */
	TDate valueDate,	/* (I) value date */
	int mmDen,		/* (I) money market denominator */
	int swapFreq,		/* (I) swap frequency */
	TDayCount swapDcc,	/* (I) swap day count convention */
	TCurve **zcCurveO);	/* (O) zero curve */

static	int	getSwaptionMatrix(
	TCSheet *swo_vol,		/* (I) "swo_vol" liverate format */
	int swapFreq,			/* (I) swap frequency */
	TSwaptionMatrix2D **swmatO);	/* (O) zero curve */

static	int	getBaseVol(
	TCSheet *base_atm,	/* (I) "base_atm" liverate format */
	TDate baseDate,		/* (I) base date */
	int bvFreq,		/* (I) base vol frequency */
	TCurve **bvCurveO);	/* (O) zero curve */


/*f--------------------------------------------------------------
 * Reads the market data from a liverate directory
 * and creates a TDrWrapperData structure.
 * The following files are read:
 * The live yield curve "yldcrvlc.cur",
 * the atm base volatility "base_atm.cur" and the 
 * swaption volatility matrix "swo_atm.cur".
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * The 3-letter currency code "curCode" is appened to all file names.
 * The arguments
 *	"mmDen" (money market denominator), 
 *	"bvFreq" (base volatility frequency), 
 *	"swapFreq" (swap frequency),
 *	"swapDcc" (swap day count convention),
 * 	"yldcrvBnam" (yield curve file basename),
 * 	"baseAtmBnam" (ATM base volatility file basename),
 * 	"swoAtmBnam, (ATM swaption vol file basename)
 * are necessary for the curve generation.
 * Remark: because the data does not contain the value date information,
 * but just the S/N date in the base volatility matrix "base_atm" file,
 * the value date is set equal to the today date which is assumed to be 
 * the S/N less one day.
 * If succeed, then SUCCESS is returned and a new that is allocated.
 * Otherwise, returns FAILURE.
 */

DLL_EXPORT(int)
DriTDrWrapperDataGetLiverate(
	const char *pathdir,		/* (I) tmp direcory name (or NULL) */
	const char *curCode,		/* (I) currency code */
	int mmDen,			/* (I) money market denominator */
	int bvFreq,			/* (I) base vol frequency */
	int swapFreq,			/* (I) swap frequency */
	TDayCount swapDcc,		/* (I) swap day count convention */
	const char* yldcrvBnam,		/* (I) yield curve file basename */
	const char* baseAtmBnam,	/* (I) ATM base vol file basename */
	const char* swoAtmBnam,		/* (I) ATM swaption vol file basename */
	TDrWrapperData **that)	/* (O) dr wrapper data */
{
static	char	routine[] = "DriTDrWrapperDataGet";
	int	status = FAILURE;

	char	fnam[256];

	TCSheet	*swo_atm = NULL;		/* swaption matrix */
	TCSheet	*base_atm = NULL;		/* base vol matrix */
	TCSheet	*yldcrvlc = NULL;		/* yield curve */

	TDate	baseDate;

#undef	PATHNAME
#define	PATHNAME(basename)	\
	(sprintf(fnam, "%s%s%s.%s", (pathdir&&pathdir[0] ?pathdir:""),\
		(pathdir&&pathdir[0]?"/":""), \
		(basename), curCode), fnam)

	/* Initialize */
	if ((*that = NEW(TDrWrapperData)) == NULL) goto done;
	(*that)->fDiscZcCurve = NULL;
	(*that)->fZcCurve = NULL;
	(*that)->fRiskZcCurve = NULL;
	(*that)->fBasisZcCurve = NULL;
	(*that)->fBvCurve = NULL;
	(*that)->fBSVolCurve = NULL;
	(*that)->fCmsSwMat = NULL;


	/*
	 * Read inputs and put in table 
	 */
	IF_FAILED_DONE( DrlTCSheetFileRead(&swo_atm,  PATHNAME(swoAtmBnam)));
	IF_FAILED_DONE( DrlTCSheetFileRead(&base_atm, PATHNAME(baseAtmBnam)));
	IF_FAILED_DONE( DrlTCSheetFileRead(&yldcrvlc, PATHNAME(yldcrvBnam)));


	GTO_IF_LOGGING( \
		DrlFPrintf(NULL, "swo_atm:\n"); \
		IF_FAILED_DONE( DrlTCSheetFpWrite(swo_atm,  NULL)); \
		DrlFPrintf(NULL, "base_atm:\n"); \
		IF_FAILED_DONE( DrlTCSheetFpWrite(base_atm, NULL)); \
		DrlFPrintf(NULL, "yldcrvlc:\n"); \
		IF_FAILED_DONE( DrlTCSheetFpWrite(yldcrvlc, NULL)));
	

	/*
	 * Get base Date: REMARK base date = value date.
	 */
	/* read spot next */
	IF_FAILED_DONE(DrlTCSheetGetVType(
		base_atm, 1, 1, DRL_TDATE_T, (void*) &baseDate));
	baseDate -= 1L;
	(*that)->fToday = baseDate;


	/*
	 * Build zero curve
	 */

	IF_FAILED_DONE( generateZcFromYldcrv(
		yldcrvlc,
		baseDate,
		mmDen,
		swapFreq,
		swapDcc,
		&((*that)->fDiscZcCurve)));



	/*
	 * Get swaption volatility
	 */
	IF_FAILED_DONE( getSwaptionMatrix(
		swo_atm,
		swapFreq,
		&(*that)->fCmsSwMat));



	/*
	 * Base volatility
	 */
	IF_FAILED_DONE( getBaseVol(
		base_atm,
		baseDate,
		bvFreq,
		&((*that)->fBvCurve)));

	GTO_IF_LOGGING( \
		DrlFPrintf(NULL, "baseDate: %10s\n", \
			DrlTDatePrint(NULL, baseDate)); \
		DrlTCurveFpWrite((*that)->fDiscZcCurve, \
			NULL,  DRL_TCURVE_FMT_PERCENT); \
		DrlTSwaptionMatrix2DFpWrite((*that)->fCmsSwMat, \
			NULL, TSWAPTION_MATRIX_FMT_STD); \
		DrlTCurveFpWrite((*that)->fBvCurve, \
			NULL,  DRL_TCURVE_FMT_PERCENT));



	/* OK */
	status = SUCCESS;
done:
	DrlTCSheetFree(swo_atm);
	DrlTCSheetFree(base_atm);
	DrlTCSheetFree(yldcrvlc);

	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


/*---------------------------------------------------------------
 * Convenience routine to read "yldcrv" format
 */

static	int
generateZcFromYldcrv(
	TCSheet *yldcrvlc,	/* (I) "yldcrv" liverate format */
	TDate valueDate,	/* (I) value date */
	int mmDen,		/* (I) money market denominator */
	int swapFreq,		/* (I) swap frequency */
	TDayCount swapDcc,	/* (I) swap day count convention */
	TCurve **zcCurveO)	/* (O) zero curve */
{
static	char	routine[] = "generateZcFromYldcrv";
	int	status = FAILURE;

#define	DEF_SWAPRATE_MAX	64
	int	i;

	double	rateArr[DEF_SWAPRATE_MAX];
	TDate	dateArr[DEF_SWAPRATE_MAX];
	double	priceArr[DEF_SWAPRATE_MAX];
	char	nameArr[DEF_SWAPRATE_MAX];
	int	useDefault = GtoSETINSTR;
	int	numRate;
	double	*futArr = NULL;
	TDate 	*futDateArr = NULL;
	int	numFuture = 0;
	double	yieldVolatility = 0;
	int	stubMethod = GtoSTUB3M;
	double	*stubData = NULL;
	int	stubDataSize = 0;
	int	couponInterpMethod = GtoLINEARINTERP;
	int	zeroInterpMethod = GtoLINEARINTERP;
	char  	fwdLength = 'Q';

	TCurve	*zcCurve;	/* zero curve */
	char	buf[32], *p;
	TDateInterval	intv;
	double	value;

	/*
	 *
	 */
	numRate = yldcrvlc->numRows;
	for (i=0; i< numRate; i++) {
		/* Read interval and mm type */
		strcpy(buf, DrlTCSheetGet(yldcrvlc, i, 0));

		if ((p = strstr(buf, "MMY")) != NULL) {
			nameArr[i] = GtoMONEYNAME;
		} else if ((p = strstr(buf, "BEY")) != NULL) {
			nameArr[i] = GtoSWAPNAME;
		} else {
			GtoErrMsg("%s: can read `%s'.\n", routine, buf);
			goto done;
		}
		*p = '\0';
		if (DrlTDateIntervalScan(buf, &intv) != SUCCESS) {
			GtoErrMsg("%s: can read `%s'.\n", routine, buf);
			goto done;
		}

		/* calculate maturity date */
		IF_FAILED_DONE(GtoDtFwdAny(
			valueDate,
			&intv, 
			&dateArr[i]));

		/* par bonds */
		priceArr[i] = 1.0;

		/* get rate */
		rateArr[i] = 0e0;

		IF_FAILED_DONE( DrlTCSheetGetVType(
			yldcrvlc,
			i, 1,
			DRL_DOUBLE_T,
			(void*) &value));

		rateArr[i] += value;

		IF_FAILED_DONE( DrlTCSheetGetVType(
			yldcrvlc,
			i, 4,
			DRL_DOUBLE_T,
			(void*) &value));

		rateArr[i] += value*1e-4;

	}

	/*
	 * Generate zero curve
	 */
	zcCurve = GtoNPiZC(
		valueDate,
		mmDen,
		swapFreq,
		swapDcc,
		rateArr,
		dateArr,
		priceArr,
		nameArr,
		useDefault,
		numRate,
		futArr,
		futDateArr,
		numFuture,
		0, NULL,
		yieldVolatility,
		stubMethod,
		stubData,
		stubDataSize,
		couponInterpMethod,
		zeroInterpMethod,
		fwdLength,
		0
		);

	if (zcCurve == NULL) goto done;

	*zcCurveO = zcCurve;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*---------------------------------------------------------------
 * Convenience routine to read "yldcrv" format
 */

static	int
getSwaptionMatrix(
	TCSheet *swo_vol,		/* (I) "swo_vol" liverate format */
	int swapFreq,			/* (I) swap frequency */
	TSwaptionMatrix2D **swmatO)	/* (O) zero curve */
{
static	char	routine[] = "getSwaptionMatrix";
	int	status = FAILURE;
	int	numExp, numMat, idxE, idxM;
	double	value;
	TSwaptionMatrix2D	*swmat = NULL;

	/*
	 *
	 */
	numExp = swo_vol->numRows-1;
	numMat = swo_vol->numCols-1;

 
	/* malloc memory */
	swmat = DrlTSwaptionMatrix2DNew(
                FALSE,  /* vertical */
                (long) swapFreq,
                numExp,
		NULL,
		numMat,
		NULL);
	if (swmat == NULL) goto done;
 
 
	/* Read data */
	for (idxE=0; idxE<numExp; idxE++) {
		IF_FAILED_DONE( DrlTCSheetGetVType(
			swo_vol,
			idxE+1, 0,
			DRL_DOUBLE_T,
			(void*) &value));

		TSWAPTION_MATRIX2D_EXP(swmat, idxE) = value;
	}


	for (idxM=0; idxM<numMat; idxM++) {
		IF_FAILED_DONE( DrlTCSheetGetVType(
			swo_vol,
			0, idxM+1,
			DRL_DOUBLE_T,
			(void*) &value));

		TSWAPTION_MATRIX2D_MAT(swmat, idxM) = value;
	}


	for (idxE=0; idxE<numExp; idxE++)
	for (idxM=0; idxM<numMat; idxM++) {
		IF_FAILED_DONE( DrlTCSheetGetVType(
			swo_vol,
			idxE+1, idxM+1,
			DRL_DOUBLE_T,
			(void*) &value));

		TSWAPTION_MATRIX2D_VOL(swmat, idxE, idxM) = value;
	}


	*swmatO = swmat;

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlTSwaptionMatrix2DFree(swmat);
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*---------------------------------------------------------------
 * Convenience routine to read "yldcrv" format
 */

static	int
getBaseVol(
	TCSheet *base_atm,	/* (I) "base_atm" liverate format */
	TDate baseDate,		/* (I) base date */
	int bvFreq,		/* (I) base vol frequency */
	TCurve **bvCurveO)	/* (O) zero curve */
{
static	char	routine[] = "getBaseVol";
	int	status = FAILURE;
	TCurve	*bvCurve = NULL;
	double	bvMat;
	int	numExp, numMat, idxE, idxM;
	long	dayCountConv = GTO_ACT_365F;
	TDate	vDate;
	double	vDouble;

	/* */
	bvMat = 1e0 / ((double) bvFreq);

	numMat = base_atm->numCols-2;
	numExp = base_atm->numRows-1;


	for (idxM=0; idxM<numMat; idxM++) {
		IF_FAILED_DONE( DrlTCSheetGetVType(
			base_atm,
			0, idxM+2,
			DRL_DOUBLE_T,
			(void*) &vDouble));

		if (IS_ALMOST_ZERO(vDouble - bvMat)) {
			break;
		}
	}
	if (idxM == numMat) {
		GtoErrMsg("%s: did not find maturity %lf in matrix.\n",
			routine, bvMat);
		goto done;
	}

	/*
	 * Create curve.
	 */
	bvCurve = GtoNewTCurve(
		baseDate,
		numExp,
		(double) bvFreq,
		dayCountConv);
	if (bvCurve == NULL) goto done;

	for (idxE=0; idxE<numExp; idxE++) {
		IF_FAILED_DONE( DrlTCSheetGetVType(
			base_atm,
			idxE+1, 1,
			DRL_TDATE_T,
			(void*) &vDate));

		IF_FAILED_DONE( DrlTCSheetGetVType(
			base_atm,
			idxE+1, idxM+2,
			DRL_DOUBLE_T,
			(void*) &vDouble));

		bvCurve->fArray[idxE].fDate = vDate;
		bvCurve->fArray[idxE].fRate = vDouble;
	}


	*bvCurveO = bvCurve;

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoFreeTCurve(bvCurve);
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}





