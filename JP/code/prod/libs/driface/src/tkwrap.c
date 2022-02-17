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

#include "drlio.h"
#include "drlsmat.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlvtype.h"
#include "drltime.h"		/* DrlTDateScanYMD */

#include "dritkwrp.h"	/* prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__



static	char	todayFnam[] = "today.dat";
static	char	discZcFnam[] = "disczero.dat";
static	char	riskZcFnam[] = "riskzero.dat";
static	char	zcFnam[] = "zero.dat";
static	char	basisZcFnam[] = "basiszero.dat";
static	char	bvFnam[] = "basevol.dat";
static	char	swvolFnam[] = "swapvol.dat";
static	char	bsvolFnam[] = "basisvol.dat";
static	char	mparamFnam[] = "modelParameters.dat";




/*f--------------------------------------------------------------
 * Reads the market data from a DR Wrapper type 2
 * and creates a TDrWrapperData.
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * Convenience routine that calls 
 * "DriTDrWrapperDataGetFull" with DRI_DRW_TYPE2_2CURVES.
 */

DLL_EXPORT(int)
DriTDrWrapperDataGet(
	char *pathdir,			/* (I) tmp direcory name (or NULL) */
	TDrWrapperData **drWrapper)	/* (O) dr wrapper data */
{
	return DriTDrWrapperDataGetFull(
		pathdir,
		DRI_DRW_TYPE2_2CURVES,
		drWrapper);
}


/*f--------------------------------------------------------------
 * Reads the market data from a DR wrapper
 * and creates a "TDrWrapperData".
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * The "options" argument dteremines the type of the wrapper
 * and the version. Possible values are:\\
 * DRI_DRW_TYPE2_2CURVES: Reads a type 2 wrapper
 * with two zero curves (disczero.dat and zero.dat)
 * DRI_DRW_TYPE2_3CURVES: Reads a type 2 wrapper
 * with three curves (disczero.dat, zero.dat and riskzero.dat).
 * DRI_DRW_BASIS: Reads a basis type wrapper
 * with four curves (disczero.dat, zero.dat, riskzero.dat, basiszero.dat),
 * plus basis spot vol curve.
 * Return NULL TCurve pointer if file not available.
 */

DLL_EXPORT(int)
DriTDrWrapperDataGetFull(
	char *pathdir,			/* (I) tmp direcory name (or NULL) */
	long options,			/* (I) see description */
	TDrWrapperData **drWrapper)	/* (O) dr wrapper data */
{
static	char	routine[] = "DriTDrWrapperDataGet";
	int	status = FAILURE;
	char	fnam[256];

        FILE   *fp = NULL;
        char    buf1[256], *buf = buf1;
/*
 NOW STATIC (see above)
static	char	todayFnam[] = "today.dat";
static	char	discZcFnam[] = "disczero.dat";
static	char	riskZcFnam[] = "riskzero.dat";
static	char	zcFnam[] = "zero.dat";
static	char	basisZcFnam[] = "basiszero.dat";
static	char	bvFnam[] = "basevol.dat";
static	char	swvolFnam[] = "swapvol.dat";
static	char	bsvolFnam[] = "basisvol.dat";
static	char	mparamFnam[] = "modelParameters.dat";
*/

	/* */
	if ((*drWrapper = NEW(TDrWrapperData)) == NULL) {
		goto done;
	}

	(*drWrapper)->fDiscZcCurve = NULL;
	(*drWrapper)->fZcCurve = NULL;
	(*drWrapper)->fRiskZcCurve = NULL;
	(*drWrapper)->fBasisZcCurve = NULL;
	(*drWrapper)->fBvCurve = NULL;
	(*drWrapper)->fBSVolCurve = NULL;
	(*drWrapper)->fCmsSwMat = NULL;

    /* Swaption vol Curve */
    if (pathdir != NULL)
        sprintf(fnam, "%s/%s", pathdir, swvolFnam);
    else
        strcpy(fnam, swvolFnam);
    if (fopen(fnam, "r") == NULL)
        (*drWrapper)->fCmsSwMat = NULL;
    else if (DrlTSwaptionMatrix2DFileRead(&(*drWrapper)->fCmsSwMat, fnam,
        TSWAPTION_MATRIX_FMT_LON) != SUCCESS) goto done;


	/* Disc Zero Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, discZcFnam);
	} else {
		strcpy(fnam, discZcFnam);
	}
	if (DrlTCurveFileRead(&((*drWrapper)->fDiscZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;


	/* Zero Curve (index) */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, zcFnam);
	} else {
		strcpy(fnam, zcFnam);
	}
	if (DrlTCurveFileRead(&((*drWrapper)->fZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;

        /* Read additional parameters from zero.dat - day count conv's
	 * and swap frequency */
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    GtoErrMsg("%s: can't open `%s' (%s).\n",
		      routine, fnam, strerror(errno));
	    goto done;
	}

	/* Zero Curve (risky) */
	if (options >= DRI_DRW_TYPE2_3CURVES) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, riskZcFnam);
	    } else {
		strcpy(fnam, riskZcFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fRiskZcCurve = NULL;
	    else if (DrlTCurveFileRead(&((*drWrapper)->fRiskZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;
	} else {
	    (*drWrapper)->fRiskZcCurve = NULL;
	}

	/* Bais Zero Curve (basiszero) */
	if (options == DRI_DRW_BASIS) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, basisZcFnam);
	    } else {
		strcpy(fnam, basisZcFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fBasisZcCurve = NULL;
	    else if (DrlTCurveFileRead(&((*drWrapper)->fBasisZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;

	} else {
	    (*drWrapper)->fBasisZcCurve = NULL;
	}


#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
    
    READ_DATA(DRL_STRING_T, &buf,	        "start date");
    READ_DATA(DRL_LONG_T, &((*drWrapper)->fMMDenom),  "mm basis");

    READ_DATA(DRL_STRING_T, &buf,	        "swap frequency");
    if ((*drWrapper)->fCmsSwMat != NULL)
    {
        if (DrlStrLongValueScan(buf, "swap frequency", 
            &((*drWrapper)->fCmsSwMat->swapPayFreq),
            "A", 1L,
            "S", 2L,
            NULL) != SUCCESS) goto done;
    }

    READ_DATA(DRL_STRING_T, &buf,	        "swap basis");
    if (DrlStrLongValueScan(buf, "swap basis", &((*drWrapper)->fSwDcc),
		"ACT", (TDayCount) GTO_B30_360,
		"360", (TDayCount) GTO_ACT_360,
		"365", (TDayCount) GTO_ACT_365F,
		NULL) != SUCCESS) goto done;

        fclose(fp); fp = NULL;
        
    /* Base vol Curve */
    if (pathdir != NULL) 
        sprintf(fnam, "%s/%s", pathdir, bvFnam);
    else
        strcpy(fnam, bvFnam);
    if (fopen(fnam, "r") == NULL)
        (*drWrapper)->fBvCurve = NULL;
    else if (DrlTCurveLondonBaseVolFileRead(&((*drWrapper)->fBvCurve), fnam)
        != SUCCESS) goto done;

	/* Basis vol Curve */
	if (options == DRI_DRW_BASIS) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, bsvolFnam);
	    } else {
		strcpy(fnam, bsvolFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fBSVolCurve = NULL;
	    else if (DrlTCurveLondonBaseVolFileRead(
			&((*drWrapper)->fBSVolCurve), 
			fnam) != SUCCESS) 
		goto done;
	} else {
	    (*drWrapper)->fBSVolCurve = NULL;
	}


	/*
	 * Set today date to value date (Kapital)
	 */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, todayFnam);
	} else {
		strcpy(fnam, todayFnam);
	}
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    GtoErrMsg("%s: can't open `%s' (%s).\n",
		      routine, fnam, strerror(errno));
	    goto done;
	}
	READ_DATA(DRL_STRING_T, &buf,	        "today");
	if (DrlTDateScanYMD(buf, &((*drWrapper)->fToday)) != SUCCESS) {
		goto done;
	}
        fclose(fp); fp = NULL;

#ifdef	_SKIP
	/* We don't have today, so take base date */
       	(*drWrapper)->fToday = (*drWrapper)->fBvCurve->fBaseDate;
#endif

	/*
	 * Read model parameters
	 */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, mparamFnam);
	} else {
		strcpy(fnam, mparamFnam);
	}
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    GtoErrMsg("%s: can't open `%s' (%s).\n",
		      routine, fnam, strerror(errno));
	    goto done;
	}

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f1Beta),
		"oneFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f1Weight),
		"oneFactorVolatility1");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f1Ppy),
		"oneFactorPPY");

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Beta1),
		"twoFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Beta2),
		"twoFactorMeanReversion2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Weight1),
		"twoFactorVolatility1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Weight2),
		"twoFactorVolatility2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f2Corr12),
		"twoFactorCorrelation1and2");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f2Ppy),
		"twoFactorPPY");

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta1),
		"threeFactorMeanReversion1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta2),
		"threeFactorMeanReversion2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Beta3),
		"threeFactorMeanReversion3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight1),
		"threeFactorVolatility1");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight2),
		"threeFactorVolatility2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Weight3),
		"threeFactorVolatility3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr12),
		"threeFactorCorrelation1and2");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr13),
		"threeFactorCorrelation1and3");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->f3Corr23),
		"threeFactorCorrelation2and3");
	READ_DATA(DRL_INT_T, &((*drWrapper)->f3Ppy),
		"threeFactorPPY");

	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->LQ),
		"LHQ");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->RQ),
		"RHQ");
	READ_DATA(DRL_DOUBLE_T, &((*drWrapper)->FwdShift),
		"FWD shift");
	READ_DATA(DRL_INT_T, &((*drWrapper)->nbIter),
		"Nb of Iterations");


        fclose(fp); fp = NULL;



	/* OK */
	status = SUCCESS;
done:

        if (fp != NULL) fclose(fp);
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



/*f--------------------------------------------------------------
 * Frees a DRW allocated by "DriTDrWrapperDataGet".
 */

DLL_EXPORT(int)
DriTDrWrapperDataFree(TDrWrapperData *drWrapper)
{
	if (drWrapper != NULL) {
		GtoFreeTCurve(drWrapper->fDiscZcCurve);
		GtoFreeTCurve(drWrapper->fZcCurve);
		GtoFreeTCurve(drWrapper->fRiskZcCurve);
		GtoFreeTCurve(drWrapper->fBasisZcCurve);
		GtoFreeTCurve(drWrapper->fBvCurve);
		GtoFreeTCurve(drWrapper->fBSVolCurve);
		DrlTSwaptionMatrix2DFree(drWrapper->fCmsSwMat);
		FREE(drWrapper);
	}
	return(SUCCESS);
}

/*f--------------------------------------------------------------
 * Prints a wrapper data on a file pointer "fp" (for debugging).
 */


DLL_EXPORT(int)
DriTDrWrapperDataFpWrite(
	TDrWrapperData *drwData,	/* (I) dr wrapper data */
	FILE *fp)			/* (I) FILE to write to (or NULL) */
{
static	char	routine[] = "DriTDrWrapperDataFpWrite";
	int	status = FAILURE;

	DrlFPrintf(fp, "TODAY:\n\t%s\n",
		DrlTDatePrint(NULL, drwData->fToday));

	DrlFPrintf(fp, "DISC_ZC_CURVE:\n");
	if (drwData->fDiscZcCurve) {
		DrlTCurveFpWrite(drwData->fDiscZcCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}

	DrlFPrintf(fp, "ZC_CURVE:\n");
	if (drwData->fZcCurve) {
		DrlTCurveFpWrite(drwData->fZcCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}
#ifdef	_TO_BE_DONE
	TCurve		*fRiskZcCurve;		/* Risky zero curve */
	TCurve		*fBasisZcCurve;		/* Basis zero curve */
#endif

	DrlFPrintf(fp, "BV_CURVE:\n");
	if (drwData->fBvCurve) {
		DrlTCurveFpWrite(drwData->fBvCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}

	DrlFPrintf(fp, "SW_MAT:\n");
	if (drwData->fCmsSwMat) {
		DrlTSwaptionMatrix2DFpWrite(drwData->fCmsSwMat, fp,
			TSWAPTION_MATRIX_FMT_STD);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}



#ifdef	_TO_BE_DONE
	TCurve		*fBvCurve;		/* Base volatility curve */
	TCurve		*fBSVolCurve;		/* Base volatility curve */

        long            fMMDenom;               /* 360 or 365 */
        TDayCount       fSwDcc;                 /* Swap day count conv */

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
#endif


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Creates the output price file for a DRW type II.
 */

DLL_EXPORT(int)
DriTDrWrapperDataPutPrice(double value)
{
static	char	routine[] = "DriTDrWrapperDataPutPrice";
	FILE	*fp = NULL;
static	char	fnam[] = "price";

	if ((fp = fopen(fnam, "w")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, strerror(errno));
		return(FAILURE);
	}
	fprintf(fp, "%.8f", value);
	fclose(fp);
	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * Returns a VnfmData set with the parameters passed from the wrapper,
 * with number of factor equal to "numFact".
 * The dates used in the timeline are selected using 
 * the "datesFrom" argument: 0 for flat volatility structure (no timeline),
 * 1 for the base volatiity dates, 2 for the swaption matrix dates.\\
 */

DLL_EXPORT(VnfmData*)
DriTDrWrapperDataGetVnfmData(
	TDrWrapperData *drWrapper,	/* (I) wrapper data */
	int numFact,			/* (I) number of factors */
	double backBoneQ,		/* (I) back bone q */
	int datesFrom)			/* (I) 0=flat, 1=bv dates, 2=swvol dates */
{
static	char	routine[] = "DriTDrWrapperDataGetVnfmData";
	int	status = FAILURE;

	VnfmData *vnfmData = NULL;
	TDate	*dates = NULL;
	int	idxD, nDates;


	/* Create timeline */
	switch (datesFrom) {
	case 0:
		nDates = 2;
		if ((dates = NEW_ARRAY(TDate, nDates)) == NULL) goto done;
		dates[0] = drWrapper->fToday;
		if (GtoTDateAdvanceYears(drWrapper->fToday,
			10e0, &dates[1]) != SUCCESS) goto done;
		break;
	case 1:
	case 2:
		GtoErrMsg("%s: bad argument datesFrom %d (not implemented).\n",
			routine, datesFrom);
		goto done;
	default:
		GtoErrMsg("%s: bad argument datesFrom %d.\n", routine, datesFrom);
		goto done;
	}


	/* Create vnfm structure */
	vnfmData = VnfmNewTimeLine(
		numFact,
		backBoneQ,
		dates[0],
		nDates,
		dates,
		NULL,
		drWrapper->fDiscZcCurve);
	if (vnfmData == NULL) goto done;


	/* Fill parameters */
	switch (numFact) {
	case 1:
		vnfmData->fBeta[0]  = drWrapper->f1Beta;
		vnfmData->fAlpha[0] = drWrapper->f1Weight;
		for(idxD=0; idxD<=nDates-1; idxD++) {
			vnfmData->fSigma[0][idxD] = 1e0;
		}
		break;
	case 2:
		vnfmData->fBeta[0]  = drWrapper->f2Beta1;
		vnfmData->fBeta[1]  = drWrapper->f2Beta2;
		vnfmData->fAlpha[0] = drWrapper->f2Weight1;
		vnfmData->fAlpha[1] = drWrapper->f2Weight2;
		for(idxD=0; idxD<=nDates-1; idxD++) {
			vnfmData->fSigma[0][idxD] = 1e0;
			vnfmData->fSigma[2][idxD] = 1e0;
			vnfmData->fRho[0][idxD] = drWrapper->f2Corr12;
		}
		break;
	case 3:
		vnfmData->fBeta[0]  = drWrapper->f3Beta1;
		vnfmData->fBeta[1]  = drWrapper->f3Beta2;
		vnfmData->fBeta[2]  = drWrapper->f3Beta3;
		vnfmData->fAlpha[0] = drWrapper->f3Weight1;
		vnfmData->fAlpha[1] = drWrapper->f3Weight2;
		vnfmData->fAlpha[2] = drWrapper->f3Weight3;
		for(idxD=0; idxD<=nDates-1; idxD++) {
			vnfmData->fSigma[0][idxD] = 1e0;
			vnfmData->fSigma[1][idxD] = 1e0;
			vnfmData->fSigma[2][idxD] = 1e0;
			vnfmData->fRho[0][idxD] = drWrapper->f3Corr12;
			vnfmData->fRho[1][idxD] = drWrapper->f3Corr13;
			vnfmData->fRho[2][idxD] = drWrapper->f3Corr23;
		}
		break;
	default:
		GtoErrMsg("%s: bad number of factors.\n", routine, numFact);
		goto done;
	}

	/* Recalc */
	if (VnfmComputeCoeff(vnfmData) != SUCCESS)
		goto done;



	/* OK */
	status = SUCCESS;
done:
	FREE(dates);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
		VnfmFree(vnfmData);
		return(NULL);
	}
	return(vnfmData);
}







/*f--------------------------------------------------------------
 * Given DR wrapper market data "drWrapper" and a calibration
 * index specification,
 * interpolates and returns the swaptions or base volatilities
 * The calibration index is specified by the arguments
 * "calibFinal", "calibMatDate" and "calibMatInt". \\
 * If "calibFinal" is TRUE, the calibration is done along a diagonal
 * and the swaption volatilities are used. If "calibMatDate" is
 * not 0, it is used to determine the diagonal, otherwise
 * "calibMatInt" is used to generate a final maturity date for the diagonal.\\
 * If "calibFinal" is FALSE, the calibration is done along a vertical
 * and the maturity used is obtained from "calibMatInt".
 * If "calibMatInt" is less or equal to the base volatility maturity,
 * the base volatility is returned. Otherwise, the swaption volatilities
 * of the corresponding maturities are used.\\
 * On exit, "numVols" contains the number of interpolated points,
 * and "volDates", "volMat", "volFreq", "volRates" are allocated.
 * They must be freed by the calling routine.
 */

DLL_EXPORT(int)
DriTDrWrapperDataGetInterpVol(
	TDrWrapperData *drWrapper,	/* (I) wrapper data */

	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	TDate calibMatDate,		/* (I) final mat (used if final) */
	TDateInterval calibMatInt,	/* (I) fwd mat (used if cms) */

	int *numVols,			/* (O) # of vol points (can be NULL) */
	TDate **volDates,		/* (O) vol exp dates (can be NULL) */
	double **volMat,		/* (O) vol fwd mat (can be NULL) */
	int **volFreq,			/* (O) vol freq (can be NULL) */
	double **volRates)		/* (O) vol (can be NULL) */
{
static	char	routine[] = "DriTDrWrapperDataGetInterpVol";
	int	status = FAILURE;

	TDate		baseDate = drWrapper->fToday;
	double		bvMat = 0.25,
			calibMatYrs,
			minMat = 0.25;


    /* base vol maturity */
    if ( drWrapper->fBvCurve != NULL)
        bvMat = 1e0 / (double) drWrapper->fBvCurve->fBasis;

	if ((!calibFinal) || (calibMatDate == 0L)) {
		IF_FAILED_DONE( GtoDateIntervalToYears(
			&calibMatInt, &calibMatYrs));
	}


	/*
	 * If final calibration or if CMS calib with
	 * mat larger than base vol mat, use swaptions
	 * Otherwise use base vol.
	 */
	if ((calibFinal) ||
	    (calibMatYrs > bvMat + 1e-4)) {
		/*
		 * Interp from swaption matrix
		 */

        if (drWrapper->fCmsSwMat == NULL ){
            GtoErrMsg("%s failed: CMS vol curve is NULL.\n",
                      routine);
            goto done;
        }
        minMat = drWrapper->fCmsSwMat->table->dim2Values[0];

	
		/* If final, but not date, use interval
		 */
		if ((calibFinal) && (calibMatDate == 0L)) {
			IF_FAILED_DONE( GtoDtFwdAny(
				baseDate,
				&calibMatInt,
				&calibMatDate));
		}


		/* Perform interp
		 */
		IF_FAILED_DONE( DrlTSwaptionMatrix2DInterpVolCurve(
			drWrapper->fCmsSwMat,
			baseDate,
			calibFinal,
			calibMatDate,
			calibMatInt,
			minMat,	
			FALSE,	/* TRUE=rebucket, FALSE=interp */
			numVols,
			volDates,
			NULL,
			volMat,
			volFreq,
			volRates,
			NULL));




	} else {
		/*
		 * Interp from base volatility
		 */
		TCurve	*bvCurve = drWrapper->fBvCurve;
		int	idx;

        if (bvCurve == NULL ){
            GtoErrMsg("%s failed: base vol curve is NULL.\n",
                      routine);
            goto done;
        }

		*numVols = bvCurve->fNumItems;

		ASSERT_OR_DONE(
			(*volDates = NEW_ARRAY(TDate, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volMat = NEW_ARRAY(double, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volFreq = NEW_ARRAY(int, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volRates = NEW_ARRAY(double, *numVols)) != NULL);


		for (idx=0; idx<=*numVols-1; idx++) {
			(*volDates)[idx] = bvCurve->fArray[idx].fDate;
			(*volMat)[idx]   = bvMat;
			(*volFreq)[idx]  = (int) (1e0 / bvMat);
			(*volRates)[idx] = bvCurve->fArray[idx].fRate;
		}


	}
	
	
	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f--------------------------------------------------------------
 * Returns a pointer to the zero curve specified by the field zcType:
 * "D" for discount,
 * "Z" for the index curve,
 * "R" for the risky,
 * "B" basis.
 * Returns NULL if failed.
 */

DLL_EXPORT(TCurve*)
DriTDrWrapperDataGetTCurve(
	TDrWrapperData *drWrapData,	/* (I) wrapper data */
	char zcType)			/* (I) (D)isc,(Z)idx,(R)isky,(B)asis */
{
static  char    routine[] = "DriTDrWrapperDataGetTCurve";
	TCurve  *zc = NULL;
 
	zcType = toupper(zcType);
 
	switch (zcType) {
	case 'D':
		zc = drWrapData->fDiscZcCurve;
		break;
	case 'Z':
		zc = drWrapData->fZcCurve;
		break;
	case 'R':
		if (zc == NULL) {
			GtoErrMsg("%s: no risky curve (%s) in environment.\n",
				routine, riskZcFnam);
		}
		zc = drWrapData->fRiskZcCurve;
		break;
	case 'B':
		zc = drWrapData->fBasisZcCurve;
		if (zc == NULL) {
			GtoErrMsg("%s: no basis curve (%s) in environment.\n",
				routine, basisZcFnam);
		}
		break;
	default:
		GtoErrMsg("%s: invalid zero curve type (%c).\n", 
				routine, zcType);
		return(NULL);
	}
 
	return(zc);
}



