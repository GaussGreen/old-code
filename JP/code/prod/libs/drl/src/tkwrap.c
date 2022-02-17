/************************************************************************
 * Module:	DRLIFACE
 * Function:	DRWrapper Utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>


#include "drlio.h"
#include "drlsmat.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlvtype.h"
#include "drltime.h"		/* DrlDDateScanYMD */
#ifdef	DRL_CLIB
#include "date_sup.h"
#endif

#include "drlkwrap.h"	/* prototype consistency */


/* Standard wrapper filenames */

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
 * and creates a DDrWrapData.
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * Convenience routine that calls 
 * "DrlDDrWrapDataGetFull" with DRL_DRW_TYPE2_2CURVES.
 */

int
DrlDDrWrapDataImportType2(
	char *pathdir,			/* (I) tmp dir name (or NULL) */
	DDrWrapData **drWrapper)	/* (O) dr wrapper data */
{
	return DrlDDrWrapDataImport(
		pathdir,
		DRL_DRW_TYPE2_2CURVES,
		drWrapper);
}


/*f--------------------------------------------------------------
 * Reads the market data from a DR wrapper
 * and creates a "DDrWrapData".
 * If "pathdir" is not "NULL", the files are read
 * in the corresponding directory (current directory otherwise).
 * The "options" argument dteremines the type of the wrapper
 * and the version. Possible values are:\\
 * DRL_DRW_TYPE2_2CURVES: Reads a type 2 wrapper
 * with two zero curves (disczero.dat and zero.dat)
 * DRL_DRW_TYPE2_3CURVES: Reads a type 2 wrapper
 * with three curves (disczero.dat, zero.dat and riskzero.dat).
 * DRL_DRW_BASIS: Reads a basis type wrapper
 * with four curves (disczero.dat, zero.dat, riskzero.dat, basiszero.dat),
 * plus basis spot vol curve.
 * Return NULL DCurve pointer if file not available.
 */

int
DrlDDrWrapDataImport(
	char *pathdir,			/* (I) tmp direcory name (or NULL) */
	long options,			/* (I) see description */
	DDrWrapData **drWrapper)	/* (O) dr wrapper data */
{
static	char	routine[] = "DrlDDrWrapDataImport";
	int	status = FAILURE;
	char	fnam[256];

        FILE   *fp = NULL;
        char    buf1[256], *buf = buf1;
	DDate	todayDate;
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
	if ((*drWrapper = NEW(DDrWrapData)) == NULL) {
		goto done;
	}

	(*drWrapper)->fDiscZcCurve = NULL;
	(*drWrapper)->fZcCurve = NULL;
	(*drWrapper)->fRiskZcCurve = NULL;
	(*drWrapper)->fBasisZcCurve = NULL;
	(*drWrapper)->fBvCurve = NULL;
	(*drWrapper)->fBSVolCurve = NULL;
	(*drWrapper)->fCmsSwMat = NULL;

#define	READ_DATA(type,ptr,str)	\
	{ if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	{ DrlErrMsg("%s: can't read %s.\n", routine, str); \
	goto done;}}


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
	    DrlErrMsg("%s: can't open `%s' (%s).\n",
		      routine, fnam, strerror(errno));
	    goto done;
	}
	READ_DATA(DRL_STRING_T, &buf,	        "today");
	if (DrlDDateScanYMD(buf, &((*drWrapper)->fToday)) != SUCCESS) {
		goto done;
	}
	todayDate = (*drWrapper)->fToday;
        fclose(fp); fp = NULL;

	/*
	 * Swaption vol Curve
	 */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, swvolFnam);
	} else {
		strcpy(fnam, swvolFnam);
	}

	if (DrlDSwopMatFileRead(&(*drWrapper)->fCmsSwMat, fnam,
		TSWAPTION_MATRIX_FMT_LON) != SUCCESS) goto done;


	/* Disc Zero Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, discZcFnam);
	} else {
		strcpy(fnam, discZcFnam);
	}
	if (DrlDCurveFileRead(&((*drWrapper)->fDiscZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;


	/* Zero Curve (index) */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, zcFnam);
	} else {
		strcpy(fnam, zcFnam);
	}
	if (DrlDCurveFileRead(&((*drWrapper)->fZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;

        /* Read additional parameters from zero.dat - day count conv's
	 * and swap frequency */
        if ((fp = fopen(fnam, "r")) == NULL) 
	{
	    DrlErrMsg("%s: can't open `%s' (%s).\n",
		      routine, fnam, strerror(errno));
	    goto done;
	}

	/* Zero Curve (risky) */
	if (options >= DRL_DRW_TYPE2_3CURVES) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, riskZcFnam);
	    } else {
		strcpy(fnam, riskZcFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fRiskZcCurve = NULL;
	    else if (DrlDCurveFileRead(&((*drWrapper)->fRiskZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;
	} else {
	    (*drWrapper)->fRiskZcCurve = NULL;
	}

	/* Bais Zero Curve (basiszero) */
	if (options == DRL_DRW_BASIS) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, basisZcFnam);
	    } else {
		strcpy(fnam, basisZcFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fBasisZcCurve = NULL;
	    else if (DrlDCurveFileRead(&((*drWrapper)->fBasisZcCurve), fnam,
		DRL_TCURVE_FMT_PERCENT) != SUCCESS)
		goto done;

	} else {
	    (*drWrapper)->fBasisZcCurve = NULL;
	}


    
    READ_DATA(DRL_STRING_T, &buf,	        "start date");
    READ_DATA(DRL_LONG_T, &((*drWrapper)->fMMDenom),  "mm basis");

    READ_DATA(DRL_STRING_T, &buf,	        "swap frequency");
    if (DrlStrLongValueScan(buf, "swap frequency", 
			    &((*drWrapper)->fCmsSwMat->swapPayFreq),
		"A", 1L,
		"S", 2L,
		NULL) != SUCCESS) goto done;


	READ_DATA(DRL_STRING_T, &buf,	        "swap basis");
	if (!strcmp("ACT", buf)) {
		((*drWrapper)->fSwDcc) = (DDayCount) DRL_B30_360;
	} else if (!strcmp("360", buf)) {
		((*drWrapper)->fSwDcc) = (DDayCount) DRL_ACT_360;
	} else if (!strcmp("365", buf)) {
		((*drWrapper)->fSwDcc) = (DDayCount) DRL_ACT_365F;
	} else {
		DrlErrMsg("%s: can't read %s.\n", routine, buf);
		goto done;
	}
	/* if (DrlStrLongValueScan(buf, "swap basis", &((*drWrapper)->fSwDcc),
		"ACT", (DDayCount) DRL_B30_360,
		"360", (DDayCount) DRL_ACT_360,
		"365", (DDayCount) DRL_ACT_365F,
		NULL) != SUCCESS) goto done;*/



        fclose(fp); fp = NULL;
        
	/* Base vol Curve */
	if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, bvFnam);
	} else {
		strcpy(fnam, bvFnam);
	}
	if (DrlDCurveLondonBaseVolFileRead(
		&((*drWrapper)->fBvCurve),
		todayDate,
		fnam)
			!= SUCCESS) goto done;

	/* Basis vol Curve */
	if (options == DRL_DRW_BASIS) {
	    if (pathdir != NULL) {
		sprintf(fnam, "%s/%s", pathdir, bsvolFnam);
	    } else {
		strcpy(fnam, bsvolFnam);
	    }

	    if (fopen(fnam, "r") == NULL)
		(*drWrapper)->fBSVolCurve = NULL;
	    else if (DrlDCurveLondonBaseVolFileRead(
		&((*drWrapper)->fBSVolCurve), 
		todayDate,
		fnam) != SUCCESS) 
		goto done;
	} else {
	    (*drWrapper)->fBSVolCurve = NULL;
	}



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
	    DrlErrMsg("%s: can't open `%s' (%s).\n",
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



        fclose(fp); fp = NULL;



	/* OK */
	status = SUCCESS;
done:

        if (fp != NULL) fclose(fp);
	if (status != SUCCESS)
		DrlErrMsg("%s: failed.\n", routine);
	return(status);
}



/*f--------------------------------------------------------------
 * Frees a DRW allocated by "DrlDDrWrapDataGet".
 */

int
DrlDDrWrapDataFree(DDrWrapData *drWrapper)
{
	if (drWrapper != NULL) {
		DrlDCurveFree(drWrapper->fDiscZcCurve);
		DrlDCurveFree(drWrapper->fZcCurve);
		DrlDCurveFree(drWrapper->fRiskZcCurve);
		DrlDCurveFree(drWrapper->fBasisZcCurve);
		DrlDCurveFree(drWrapper->fBvCurve);
		DrlDCurveFree(drWrapper->fBSVolCurve);
		DrlDSwopMatFree(drWrapper->fCmsSwMat);
		FREE(drWrapper);
	}
	return(SUCCESS);
}

/*f--------------------------------------------------------------
 * Prints a wrapper data on a file pointer "fp" (for debugging).
 */


int
DrlDDrWrapDataFpWrite(
	DDrWrapData *drwData,	/* (I) dr wrapper data */
	FILE *fp)			/* (I) FILE to write to (or NULL) */
{
static	char	routine[] = "DrlDDrWrapDataFpWrite";
	int	status = FAILURE;

	DrlFPrintf(fp, "TODAY:\n\t%s\n",
		DrlDDatePrint(NULL, drwData->fToday));

	DrlFPrintf(fp, "DISC_ZC_CURVE:\n");
	if (drwData->fDiscZcCurve) {
		DrlDCurveFpWrite(drwData->fDiscZcCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}

	DrlFPrintf(fp, "ZC_CURVE:\n");
	if (drwData->fZcCurve) {
		DrlDCurveFpWrite(drwData->fZcCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}
#ifdef	_TO_BE_DONE
	DCurve		*fRiskZcCurve;		/* Risky zero curve */
	DCurve		*fBasisZcCurve;		/* Basis zero curve */
#endif

	DrlFPrintf(fp, "BV_CURVE:\n");
	if (drwData->fBvCurve) {
		DrlDCurveFpWrite(drwData->fBvCurve, fp,
			DRL_TCURVE_FMT_PERCENT);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}

	DrlFPrintf(fp, "SW_MAT:\n");
	if (drwData->fCmsSwMat) {
		DrlDSwopMatFpWrite(drwData->fCmsSwMat, fp,
			TSWAPTION_MATRIX_FMT_STD);
	} else {
		DrlFPrintf(fp, "\t(nil)\n");
	}



#ifdef	_TO_BE_DONE
	DCurve		*fBvCurve;		/* Base volatility curve */
	DCurve		*fBSVolCurve;		/* Base volatility curve */

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
#endif


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Creates the output price file for a DRW type II.
 */

int
DrlDDrWrapDataPutPrice(double value)
{
static	char	routine[] = "DrlDDrWrapDataPutPrice";
	FILE	*fp = NULL;
static	char	fnam[] = "price";

	if ((fp = fopen(fnam, "w")) == NULL) {
		DrlErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, strerror(errno));
		return(FAILURE);
	}
	fprintf(fp, "%.8f", value);
	fclose(fp);
	return(SUCCESS);
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

int
DrlDDrWrapDataGetInterpVol(
	DDrWrapData *drWrapper,	/* (I) wrapper data */

	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	DDate calibMatDate,		/* (I) final mat (used if final) */
	DInterval calibMatInt,	/* (I) fwd mat (used if cms) */

	int *numVols,			/* (O) # of vol points (can be NULL) */
	DDate **volDates,		/* (O) vol exp dates (can be NULL) */
	double **volMat,		/* (O) vol fwd mat (can be NULL) */
	int **volFreq,			/* (O) vol freq (can be NULL) */
	double **volRates)		/* (O) vol (can be NULL) */
{
static	char	routine[] = "DrlDDrWrapDataGetInterpVol";
	int	status = FAILURE;

	DDate		baseDate = drWrapper->fToday;
	int		bvFreq;
	double		bvMat,
			calibMatYrs,
			minMat = 0.25;


	/* base vol maturity is given by the curve frequency */
	IF_FAILED_DONE( DrlDCurveFreq(drWrapper->fBvCurve, &bvFreq));
	bvMat = 1e0 / (double) bvFreq;
	minMat = bvMat;

	if (!calibFinal) {
		IF_FAILED_DONE( DrlDIntervalToYears(
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
	
		/* Perform interp
		 */
		IF_FAILED_DONE( DrlDSwopMatInterpVolCurve(
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
		DCurve	*bvCurve = drWrapper->fBvCurve;
		int	idx;

		*numVols = DRLTCURVE_NUMITEMS(bvCurve);

		ASSERT_OR_DONE(
			(*volDates = NEW_ARRAY(DDate, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volMat = NEW_ARRAY(double, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volFreq = NEW_ARRAY(int, *numVols)) != NULL);
		ASSERT_OR_DONE(
			(*volRates = NEW_ARRAY(double, *numVols)) != NULL);


		for (idx=0; idx<=*numVols-1; idx++) {
			(*volDates)[idx] = DRLTCURVE_DATE(bvCurve, idx);
			(*volMat)[idx]   = bvMat;
			(*volFreq)[idx]  = (int) (1e0 / bvMat);
			(*volRates)[idx] = DRLTCURVE_RATE(bvCurve, idx);
		}


	}
	
	
	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
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

DCurve*
DrlDDrWrapDataGetDCurve(
	DDrWrapData *drWrapData,	/* (I) wrapper data */
	char zcType)		/* (I) (D)isc,(Z)idx,(R)isky,(B)asis */
{
static  char    routine[] = "DrlDDrWrapDataGetDCurve";
	DCurve  *zc = NULL;
 
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
			DrlErrMsg("%s: no risky curve (%s) in environment.\n",
				routine, riskZcFnam);
		}
		zc = drWrapData->fRiskZcCurve;
		break;
	case 'B':
		zc = drWrapData->fBasisZcCurve;
		if (zc == NULL) {
			DrlErrMsg("%s: no basis curve (%s) in environment.\n",
				routine, basisZcFnam);
		}
		break;
	default:
		DrlErrMsg("%s: invalid zero curve type (%c).\n", 
				routine, zcType);
		return(NULL);
	}
 
	return(zc);
}


/*f--------------------------------------------------------------
 * Exports a wrapper data to a directory.
 * See DrlDDrWrapDataImport.
 */


int
DrlDDrWrapDataExport(
	char *pathdir,		/* (I) path to director (or NULL) */
	long options,			/* (I) see description */
	DDrWrapData *drwData)	/* (I) dr wrapper data */
{
static	char	routine[] = "DrlDDrWrapDataExportDrw";
	int	status = FAILURE;
	char	fnam[2048];

	/* Export today's date */
	sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"today.dat");
	IF_FAILED_DONE ( DrlFilePrintf(
		fnam, "%ld\n", drwData->fToday));

	/* Export zero curve */
	if (drwData->fDiscZcCurve) {
		sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"disczero.dat");
		IF_FAILED_DONE( DrlDCurveFileWrite(
			drwData->fDiscZcCurve,
			fnam,
			DRL_TCURVE_FMT_WRAPZC));
	}

	if (drwData->fZcCurve) {
		sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"zero.dat");
		IF_FAILED_DONE( DrlDCurveFileWrite(
			drwData->fZcCurve,
			fnam,
			DRL_TCURVE_FMT_WRAPZC));
	}

	if (drwData->fZcCurve) {
		sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"riskzero.dat");
		IF_FAILED_DONE( DrlDCurveFileWrite(
			drwData->fRiskZcCurve,
			fnam,
			DRL_TCURVE_FMT_WRAPZC));
	}


	/* Export base volatility curve */
	if (drwData->fBvCurve) {
		sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"basevol.dat");
		IF_FAILED_DONE( DrlDCurveFileWrite(
			drwData->fBvCurve,
			fnam,
			DRL_TCURVE_FMT_WRAPBV));
	}


	/* Export swaption matrix */
	if (drwData->fCmsSwMat) {
		sprintf(fnam, "%s/%s", (pathdir != NULL ? pathdir : "."),
			"swapvol.dat");
		IF_FAILED_DONE( DrlDSwopMatFileWrite(
			drwData->fCmsSwMat,
			fnam,
			TSWAPTION_MATRIX_FMT_LON));
	}

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




