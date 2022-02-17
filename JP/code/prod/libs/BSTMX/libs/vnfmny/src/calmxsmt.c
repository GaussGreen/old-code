/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "drlmem.h"
#include "drllno.h"
#include "drlsort.h"	/* DrlDoubleArrayFloorIdx */
#include "drlstr.h"	/* DrlFloatPrint() */
#include "drlinter.h"	/* DrlLinearInterp1dWeights() */
#include "drlio.h"	/* DrlFPrintf() */

#define	_vnfm_SOURCE
#include "vnfmwrap.h"
#include "vnfmopca.h"

#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif


#undef	SWNEXP
#undef	SWNMAT
#undef	SWFREQ
#undef	SWTEXP
#undef	SWTMAT
#undef	SWVOL

#define	SWNEXP(mtx)		(mtx->table->matrix->numDim1)
#define	SWNMAT(mtx)		(mtx->table->matrix->numDim2)
#define	SWFREQ(mtx)		(mtx->swapPayFreq)
#define	SWTEXP(mtx, idxExp)	(mtx->table->dim1Values[idxExp])
#define	SWTMAT(mtx, idxMat)	(mtx->table->dim2Values[idxMat])
#define	SWVOL(mtx, idxExp, idxMat) \
			(mtx->table->matrix->data[idxExp][idxMat])


/*
 *
 */

typedef	struct	{
	VnfmData	**tfData;
	int		idxE,
			nMat,
			*indexCnln;
	TSwaptionMatrix2D *midMkt,
			*matWeigth,
			*expWeigth;
	double		*blcnln,
			*bucnln;
	int		nIterF,
			nIterC;
} TOptimData;	

static	TOptimData	optimData;

static	int	objectiveFunction(int *MODE, int *N, double *X,
			void *userData,
			double *OBJF, double *OBJGRD);
static	int	constraintFunction(int *NCNLN, int *N, int *MODE, int *NEEDC,
			double *X,
			void *userData,
			double *C, double *cjac1);




/*--------------------------------------------------------------
 * DO NOT USE
 */

DLL_EXPORT(int)
VnfmSmoothSwaptionMatrix(
	TCurve *zcCurve,		/* (I) zero curve */
	TSwaptionMatrix2D *midMkt,	/* (I) */
	TSwaptionMatrix2D *bidToMid,	/* (I) */
	TSwaptionMatrix2D *matWeigth,	/* (I) */
	TSwaptionMatrix2D *expWeigth,	/* (I) */
	int optimType,			/* (I) */
	int normType,			/* (I) */
	double smoothParam,		/* (I) */
	double *nfParamsInL,		/* (I) LIL array of parameters */
	TSwaptionMatrix2D *origMkt,	/* (I) used for timeline only */
	TSwaptionMatrix2D **outputMat,	/* (O) output matrix */
	TSwaptionMatrix2D **spvolMat,	/* (O) output spot vol */
	double *nfParamsOut)		/* (O) */
{
static	char	routine[] = "VnfmSmoothSwaptionMatrix";
	int	status = FAILURE;

	TDate	baseDate,
		*expDatesL = NULL, *expDates = NULL;
	int	idxE, nExp,
		idxM, nMat,
		idxF, idxMlo, idxMhi, idxEv,
		freq;
	double	floatScalarsL[2];
	double	swvol, bmid, tExp, tMat, wlo, whi;

	int	NVAR = 0;		/* number of variables */
	double	*blsc=NULL,
		*busc=NULL;		/* constraints on variables */

	int	NCLIN = 0;		/* # linear constraints */
	double	*blclin=NULL,
		*buclin=NULL,
		**a=NULL;

	int	NCNLN = 0;
	double	*blcnln=NULL,
		*bucnln=NULL;

	int	ITERMAX;
	double	*C=NULL;
	double	OBJF;
	double	*X=NULL;

	VnfmData	*outTfData = NULL;

	/**
	 ** Optimization
	 **/

	/* */
	baseDate = zcCurve->fBaseDate;
	nMat = SWNMAT(midMkt);

	/* create time line */
	if (DrlTSwaptionMatrix2DCreateTimeLine(
		midMkt, 
		baseDate,
		&nExp,
		&expDates) != SUCCESS)
			goto done;

	/* create LIL array of dates for routine VnfmWrapReadSimple */
	if ((expDatesL = NEW_ARRAY(TDate, nExp+2)) == NULL) goto done;
	expDatesL[0] = nExp+1;
	expDatesL[1] = zcCurve->fBaseDate;
	for (idxE=0; idxE<= nExp-1; idxE++) expDatesL[idxE+2] = expDates[idxE];
	floatScalarsL[0] = 1e0;
	floatScalarsL[1] = VNFM_LOGNORMAL_DIST;


	/* create tf parameter data structure */
	if ((optimData.tfData = NEW_ARRAY(VnfmData*, nMat)) == NULL)
		goto done;
	for (idxM=0; idxM<=nMat-1; idxM++) optimData.tfData[idxM] = NULL;

	for (idxM=0; idxM<=nMat-1; idxM++) {
	    /* create tf calib data */

	    if (VnfmWrapReadSimple(
		&optimData.tfData[idxM],
		floatScalarsL,
		nfParamsInL,
		expDatesL,
		zcCurve) != SUCCESS)
				goto done;

#ifdef	_SKIP
	    optimData.tfData[idxM] = VnfmNew2FactSimple(
		zcCurve->fBaseDate,
		nExp,
		expDates,
		NULL,
		zcCurve,
		VNFM_LOGNORMAL_DIST,
		nfParamsIn[0],
		nfParamsIn[1],
		1e0,
		nfParamsIn[2],
		0.20e0,
		0.20e0,
		nfParamsIn[3]);
#endif
	    if (optimData.tfData[idxM] == NULL)
		goto done;

	    if (VnfmComputeCoeff(optimData.tfData[idxM]) != SUCCESS)
		goto done;
	}




	/* allocate memory for constraints */
	if ((X      = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;
	if ((blsc   = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;
	if ((busc   = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;

	if ((blcnln = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;
	if ((bucnln = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;
	if ((C      = DrlDoubleVectAlloc(0, nMat-1)) == NULL) goto done;
	if ((optimData.indexCnln = DrlIntVectAlloc(0, nMat-1)) == NULL)
		goto done;


	/**
	 ** Loop; step idxE, adjust Vnfm volatilities of index idxE-1
	 ** to match swaption volatilities od index idxE-1.
	 **/
	for (idxE=1; idxE<= nExp; idxE++) {
#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "start: ROW: %2d \n", idxE));
#endif

		/**
		 ** Optimization
		 **/

		/*
		 * BOUNDS ON VARIABLES
		 */
		NVAR = nMat;
		for (idxM=0; idxM<=NVAR-1; idxM++) {
			blsc[idxM] = 0e0;
			busc[idxM] = 100e0;
		}

		/*
		 * LINEAR CONSTRAINTS
		 */
		NCLIN = 0;

		/*
		 * NONLINEAR CONSTRAINTS:
		 * (vol-bidtomid)^2 <= vol^2 <= (vol+bidtomid)^2
		 */
		NCNLN = 0;
		for (idxM=0; idxM<=nMat-1; idxM++) {
		    bmid  =SWVOL(bidToMid, idxE-1, idxM);

		    if (bmid >= 0e0) {
			swvol = SWVOL(midMkt, idxE-1, idxM);
			bmid  = MAX(bmid, 1e-4);

			blcnln[NCNLN] = SQR(MAX(swvol - bmid, 0));
			bucnln[NCNLN] = SQR(swvol + bmid);

			optimData.indexCnln[NCNLN] = idxM;
			NCNLN++;
		    }
		}
	
		/*
		 * INITIAL VALUES
		 */
		ITERMAX = 1000;
		OBJF = 0e0;

		for (idxM=0; idxM<=NVAR-1; idxM++) {
			X[idxM] = 20.0e-2;
		}

		/* store data in structure */
		optimData.idxE = idxE;
		optimData.nMat = nMat;
		optimData.midMkt = midMkt;
		optimData.matWeigth = matWeigth;
		optimData.expWeigth = expWeigth;
		optimData.blcnln = blcnln;
		optimData.bucnln = bucnln;
		optimData.nIterF = 0;
		optimData.nIterC = 0;



		/*
		 * Do The Minimization
		 */
		if (DrlNLinProg(
			NVAR, blsc, busc,
			NCLIN, blclin, buclin, a,
			NCNLN, blcnln, bucnln, constraintFunction,
			objectiveFunction,
			ITERMAX, C, &OBJF, X,
			(void*)NULL,	/* user data */
			DRL_LNO_IMSL_NLN,	/* IMSL non linear optim */
			0)
		    != 0) {
			GtoErrMsg("%s: optimization failed.\n", routine);
			goto done;
		}


#ifndef	NO_LOGGING
		GTO_IF_LOGGING(for (idxM=0; idxM<=optimData.nMat-1; idxM++) {\
			DrlFPrintf(vnfmFpLog, "done: ROW: %2d COLUMN %2d:\n",\
				idxE, idxM);\
			VnfmPrintCoeff(optimData.tfData[idxM], vnfmFpLog);\
		});
#endif

	}

	/**
	 **  Create output matrix 
	 **/

	*outputMat = DrlTSwaptionMatrix2DNewCopy(origMkt);
	if (*outputMat == NULL) goto done;

	*spvolMat = DrlTSwaptionMatrix2DNewCopy(origMkt);
	if (*spvolMat == NULL) goto done;


	/* create output time line */
	if (DrlTSwaptionMatrix2DCreateTimeLine(
		midMkt, 
		baseDate,
		&nExp,
		&expDates) != SUCCESS)
			goto done;

	/* create tf calib data */

	if (VnfmWrapReadSimple(
	    	&outTfData,
		floatScalarsL,
		nfParamsInL,
		expDatesL,
		zcCurve) != SUCCESS)
				goto done;

#ifdef	_SKIP
	outTfData = VnfmNew2FactSimple(
		zcCurve->fBaseDate,
		nExp,
		expDates,
		NULL,
		zcCurve,
		VNFM_LOGNORMAL_DIST,
		nfParamsIn[0],
		nfParamsIn[1],
		1e0,
		nfParamsIn[2],
		0.20e0,
		0.20e0,
		nfParamsIn[3]);
	if (outTfData == NULL)
		goto done;
#endif



	freq = SWFREQ((*outputMat));

	for (idxM=0; idxM<=SWNMAT((*outputMat))-1; idxM++) {

	    tMat = SWTMAT((*outputMat), idxM);

	    if (DrlLinearInterp1dWeights(
		&SWTMAT(midMkt, 0), SWNMAT(midMkt),
		tMat,
		&idxMlo, &idxMhi, &wlo, &whi) != SUCCESS)
			goto done;

	    /* fill with interp vol */
	    for (idxE=0; idxE<=outTfData->fNDates-1; idxE++) {
		for (idxF=0; idxF<=outTfData->fNf-1; idxF++) {
		    outTfData->fSigma[idxF][idxE] =
		    wlo * optimData.tfData[idxMlo]->fSigma[idxF][idxE] + 
		    whi * optimData.tfData[idxMhi]->fSigma[idxF][idxE];
		}
	    }
	    if (VnfmComputeCoeff(outTfData) != 0)
		return(FAILURE);
#if !defined(NO_LOGGING) && defined(__DEBUG__)
	    GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "Interpolated VnfmData idxM=%d tMat=%lf\n",\
		idxM, tMat);\
	    VnfmFpWrite(outTfData, vnfmFpLog));
#endif


	    /* fill the putput matrices */
	    for (idxE=0; idxE<=SWNEXP((*outputMat))-1; idxE++) {
		tExp = SWTEXP((*outputMat), idxE);

		/* WARNING: index offset betwen Vnfm and swaption matrix */
		DrlDoubleArrayFloorIdx(outTfData->fTime, outTfData->fNDates,
			MAX(tExp-1e-3, 0e0), &idxEv);
		ASSERT_OR_DONE(idxEv >= 0);

		/* store spot vol */
		SWVOL((*spvolMat), idxE, idxM) = 
		    outTfData->fSigma[0][idxEv];

		/* swaption vol */
		if (VnfmAvgQBVol(
			outTfData,
			0e0, tExp, tExp, tMat, freq,
			LOGVOL,
			&SWVOL((*outputMat), idxE, idxM)) != SUCCESS)
				goto done;
#if !defined(NO_LOGGING) && defined(__DEBUG__)
		GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "Fill Matrix: idxE=%2d idxEv=%d idxM=%2d"\
			" tExp=%6.4f tMat=%6.4f freq=%d vol=%lf spvol=%lf\n",\
			idxE, idxEv, idxM, tExp, tMat, freq,\
			SWVOL((*outputMat), idxE, idxM),\
			SWVOL((*spvolMat), idxE, idxM)));
#endif

	    }
	}





#ifdef	_SKIP
	/* create output matrix: swaption and spot ON vol */
	modelMat = DrlTSwaptionMatrix2DNewCopy(midMkt);
	if (modelMat == NULL) goto done;

	spvolModelMat = DrlTSwaptionMatrix2DNewCopy(midMkt);
	if (spvolModelMat == NULL) goto done;

	for (idxE=0; idxE<=SWNEXP(modelMat)-1; idxE++)
	for (idxM=0; idxM<=SWNMAT(modelMat)-1; idxM++) {

		double	tExp = SWTEXP(modelMat, idxE);
		double	tMat = SWTMAT(modelMat, idxM);
		int	freq = SWFREQ(modelMat);

		/* swaption vol */
		if (VnfmAvgQBVol(
			optimData.tfData[idxM],
			0e0, tExp, tExp, tMat, freq,
			LOGVOL,
			&SWVOL(modelMat, idxE, idxM)) != SUCCESS)
				goto done;

		/* O/N vol */
		if (VnfmAvgQBVol(
			optimData.tfData[idxM],
			tExp, tExp+1e-2, tExp+1e-2, 0e0, 0,
			LOGVOL,
			&SWVOL(spvolModelMat, idxE, idxM)) != SUCCESS)
				goto done;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
		/*DrlFPrintf(vnfmFpLog, "QBVol: idxE=%2d idxM=%2d tExp=%6.4f "
			"tMat=%6.4f freq=%d vol=%lf\n",
			idxE, idxM, tExp, tMat, freq,
			SWVOL((*modelMat), idxE, idxM));*/
#endif
	}
	*outputMat = DrlTSwaptionMatrix2DNewCopy(origMkt);
	if (*outputMat == NULL) goto done;

	*spvolMat = DrlTSwaptionMatrix2DNewCopy(origMkt);
	if (*spvolMat == NULL) goto done;

	if (DrlTSwaptionMatrix2DInterpMatrix(
		modelMat,
		*outputMat,
		baseDate,
		FALSE) != SUCCESS)	/* TRUE=rebucket, FALSE=interp */
			goto done;

	if (DrlTSwaptionMatrix2DInterpMatrix(
		spvolModelMat,
		*spvolMat,
		baseDate,
		FALSE) != SUCCESS)	/* TRUE=rebucket, FALSE=interp */
			goto done;
#endif


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "OUTPUT_MAT:\n");\
		DrlTSwaptionMatrix2DFpWrite(*outputMat, vnfmFpLog,\
			TSWAPTION_MATRIX_FMT_STD);\
		DrlFPrintf(vnfmFpLog, "SPOTVOL_MAT:\n");\
		DrlTSwaptionMatrix2DFpWrite(*spvolMat,\
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD));
#endif




	/* made it through OK */
	status = SUCCESS;
done:
	/* free memory */
	DrlDoubleVectFree(X,      0, nMat-1);
	DrlDoubleVectFree(blsc,   0, nMat-1);
	DrlDoubleVectFree(busc,   0, nMat-1);
	DrlDoubleVectFree(blcnln, 0, nMat-1);
	DrlDoubleVectFree(bucnln, 0, nMat-1);
	DrlDoubleVectFree(C,      0, nMat-1);
	DrlIntVectFree(optimData.indexCnln, 0, nMat-1);

	if (expDates)  FREE(expDates);
	if (expDatesL) FREE(expDatesL);

	VnfmFree(outTfData);

	for (idxM=0; idxM<=nMat-1; idxM++)
		VnfmFree(optimData.tfData[idxM]);
	FREE(optimData.tfData);


	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}






/*--------------------------------------------------------------
 *
 */

static	int
bvNorm(double *value)
{
	int	im;

	*value = 0e0;

	for (im=0; im<=optimData.nMat-2; im++) {
		*value += SQR(optimData.tfData[im]->fSigma[0][optimData.idxE-1]
			- optimData.tfData[im+1]->fSigma[0][optimData.idxE-1])
			* SWVOL(optimData.matWeigth, optimData.idxE-1, im);
	}

	/* BV with respect to previous */
	if (optimData.idxE != 1) {
	     for (im=0; im<=optimData.nMat-1; im++) {
		*value += SQR(optimData.tfData[im]->fSigma[0][optimData.idxE-1]
			- optimData.tfData[im]->fSigma[0][optimData.idxE-2])
			* SWVOL(optimData.expWeigth, optimData.idxE-1, im);
	     }
	}



	return(SUCCESS);
}


/*--------------------------------------------------------------
 * Fills the spot vols in the 2F parameter data structures
 * and computa all internal coefficients.
 */

static	int
fillData(double *x)
{
	int	im, ie, idxF, iemax, idxFMax;

	iemax = optimData.tfData[0]->fNDates;
	idxFMax = optimData.tfData[0]->fNf;

	for (im=0; im<=optimData.nMat-1; im++) {
	    for (ie=optimData.idxE-1; ie<=iemax-1; ie++) {
		for (idxF=0; idxF<=idxFMax-1; idxF++)
		    optimData.tfData[im]->fSigma[idxF][ie] = x[im];
	    }

/*
	    if (VnfmComputeCoeff(optimData.tfData[im]) != 0)
		return(FAILURE);
	    if (VnfmComputeJ(optimData.tfData[im]) != 0)
		return(FAILURE);
*/
	    if (VnfmComputeJPartial(
			optimData.tfData[im],
	    		optimData.idxE-1,
	    		optimData.idxE+1) != 0)
		return(FAILURE);

	}

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	/*for (im=0; im<=optimData.nMat-1; im++) {
		DrlFPrintf(vnfmFpLog, "fillData: COLUMN %d:\n", im);
		VnfmPrintCoeff(optimData.tfData[im], vnfmFpLog);
	}*/
#endif
	return(SUCCESS);
}






/*ARGSUSED*/
/*--------------------------------------------------------------
 * Objective function for the minimization routine.
 */

static	int
objectiveFunction(
	int	*MODE,		/* (I) not used */
	int	*N,		/* (I) number of variables */
	double	*X,		/* (I) point where obective valued */
	void	*userData,	/* (I) user data */
	double	*OBJF,		/* (O) value of the objective */
	double	*OBJGRD) 	/* (O) not used */
{
	int	status = FAILURE;

	/*
	 * 1: new values of the model parameters
	 */
	if (fillData(X) != SUCCESS) goto done;

	/*
	 * 2: compute the objective function
	 */
	if (bvNorm(OBJF) != SUCCESS)
		goto done;

	++optimData.nIterF;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING({	int	i;\
	    DrlFPrintf(vnfmFpLog, "XOPT:");\
	    for(i=0; i<=*N-1; i++)\
		DrlFPrintf(vnfmFpLog, " %s", DrlFloatPrint(NULL, X[i],8));\
	    DrlFPrintf(vnfmFpLog, " o=%.4e\n", *OBJF) ;\
	});
#endif
	status = SUCCESS;
done:
	return(status);
}


/*ARGSUSED*/
/*--------------------------------------------------------------
 * Nonlinear constraints for minimization.
 */

static	int
constraintFunction(
	int	*NCNLN,		/* (I) number of nonlinear constraints */
	int	*N,		/* (I) number of variables */
	int	*MODE,		/* (I) not used */
	int	*NEEDC,		/* (I) array indicate if constraint needed */
	double	*X,		/* (I) point wher constraints valued */
	void	*userData,	/* (I) user data */
	double	*C,		/* (O) array of values of the constraints */
	double	*cjac1) 	/* (O) not used */
{
	int	status = FAILURE;
	int	idxC, idxM;
	int	freq = SWFREQ(optimData.midMkt);
	double	tMat,
		tExp;


	if (fillData(X) != SUCCESS) goto done;


	for(idxC=0; idxC<=*NCNLN-1; idxC++) {
	    /* Check if constraint need to be computed */
	    if (NEEDC[idxC] > 0) {

		tExp = SWTEXP(optimData.midMkt, optimData.idxE-1);
		idxM = optimData.indexCnln[idxC];
		tMat = SWTMAT(optimData.midMkt, idxM);

		if (VnfmAvgQBVol(
			optimData.tfData[idxM],
			0e0, tExp, tExp, tMat, freq,
			LOGVOL,
			&C[idxC]) != SUCCESS)
				goto done;

		C[idxC] = SQR(C[idxC]);


#if !defined(NO_LOGGING) && defined(__DEBUG__)
		/*DrlFPrintf(vnfmFpLog, "\tidxC=%2d idxM=%2d tExp=%6.4f "
			"tMat=%6.4f vol=%lf\n",
			idxC, idxM, tExp, tMat, C[idxC]);*/
		GTO_IF_LOGGING({int	i;\
		    DrlFPrintf(vnfmFpLog, "XOPT:");\
		    for(i=0; i<=*N-1; i++)\
		      DrlFPrintf(vnfmFpLog, " %s", DrlFloatPrint(NULL, X[i],8));\
		      DrlFPrintf(vnfmFpLog, " C[%2d]=%8.4f %c (%4.2fx%4.2f:%8.4f)\n",\
			    idxC, C[idxC],\
			    (C[idxC]< optimData.blcnln[idxC]-DBL_EPSILON ?'-':\
			    (C[idxC]>=optimData.bucnln[idxC]+DBL_EPSILON ?'+':\
			    's')),\
			    tExp, tMat, sqrt(C[idxC]*1e2)\
			);\
		});
#endif
	    }
	}


	status = SUCCESS;
done:
	return(status);
}


