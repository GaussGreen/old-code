/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher - Arnon Levy
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "drlmem.h"
#include "drllno.h"
#include "drlproc.h"			/* DrlProc */

#define	_vnfm_SOURCE
#include "vnfmopca.h"

#include "drlsmat.h"
#include "drlts.h"			/* DrlTCurveWrapRead() */
#include "drlvtype.h"			/* DrlLilVectLoggingFile() */
#include "drlio.h"


#include "bastypes.h"
#include "macros.h"

#include "vnfmwrap.h"	/* Prototype Consistency */


#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif


#define	GtoSwaptionMatrix2DNewEmpty(diag, freq, nExp, nMat)	\
		DrlTSwaptionMatrix2DNew(diag, freq, nExp, NULL, nMat, NULL)
#define	GtoSwaptionMatrix2DFree	DrlTSwaptionMatrix2DFree

static	long shrinkMatrix( TSwaptionMatrix2D *midMkt,
                   TSwaptionMatrix2D *bidToMid,
                   TSwaptionMatrix2D *matWeight,
                   TSwaptionMatrix2D *expWeight,
                   TSwaptionMatrix2D **shrunkMidMkt,
                   TSwaptionMatrix2D **shrunkBidToMid,
                   TSwaptionMatrix2D **shrunkMatWeight,
                   TSwaptionMatrix2D **shrunkExpWeight
                 );

/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_SMOOTH_SWMAT.
 *                                                             
 * <br><br>
 * Wrapper for swaption matrix smoothing.
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
 * <br>[swTypeL] array of 2 elements:\\
 * 	(1) swaption matrix output type (0 for vertical, 1 for diagonal)
 *		for the output <i> SwVolL</i>.\\
 * 	(2) swaption matrix volatility frequency (0,1,2,4,12)
 *		for the output <i> SwVolL</i>.
 * <br>[swMatL] array of maturity intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[swExpL] array of expration intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[midMktL] matrix of mid-market swaption volatilities.
 * <br>[bidToMidL] matrix of bid to mid tolerance
 * (a point for which no tolerance is provided should contain $-1$).
 * <br>[matWeightL] matrix of maturity weights for objective function.
 * <br>[expWeightL] matrix of expiration weights for objective function.

 * <br>[integerScalarsL] array of 2 numeric scalars:\\
 *	(1) NOT CURRENTLY USED (pass 0)\\ %optimType 
 *	(2) NOT CURRENTLY USED (pass 0)\\ %normType 
 * <br>[doubleScalarsL] array of 1 numeric scalar:\\
 *	(1) NOT CURRENTLY USED (pass 0) %smoothParam 
 * <br>[nfParamsInL] array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * that has ${n(n+3)/ 2}$ elements,
 * or a ``short'' array format where the 1st weight,
 * assumed to be 1.0, is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * that has ${n(n+3)/ 2}-1$ elements,
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * that has $n(n+2)$ elements,
 * {\bf Remark: because of possible ambiguity in the format
 * (e.g. 8 elements can correspond to $n=3$ of the short array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[nfParamsOutL] NOT CURRENTLY USED (pass dummay array of size 4).
 * <br>[outputMktL] output smoothed swaption matrix.
 * <br>[spvolMatL] output matrix of spot volatilities  matrix.
 * <br>
 */

DLL_EXPORT(int)
VnfmSmoothSwaptionMatrixL(
	TDateL *refDateL,	/* 01 'D' (I) reference date */
	TDateL *zcDateL,	/* 02 'D' (I) zero coupondates */
	FloatL *zcRateL,	/* 03 'F' (I) zero coupon rates */

	IntL *swTypeL,		/* 04 'L' (I) array of matrix param [2]: */
				/*        [0] matr type (0=vertical, 1=diag) */
				/*        [1] vol frequency (1,2,4,12) */
	double *swMatL,		/* 05 'F' (I) array of mat intervals */
	double *swExpL,		/* 06 'F' (I) array of exp intervals */
	double *midMktL,	/* 07 'F' (I) mid market matrix */
	double *bidToMidL,	/* 08 'F' (I) bid to mid matrix (>=0) */
	double *matWeightL,	/* 09 'F' (I) maturity weights */
	double *expWeightL,	/* 10 'F' (I) expiration weights */

	long *integerScalarsL,	/* 11 'L' (I) numeric scalars */
				/*        [1] optimType */
				/*        [2] normType */
	double *doubleScalarsL,	/* 12 'F' (I) float scalars */
				/*        [1] smoothParam */
	double *nfParamsInL,	/* 13 'F' (I) array of parameters */

	double *nfParamsOutL,	/* 14 'F' (O) output parameters */
	double *outputMktL,	/* 15 'F' (O) output swaption matrix */
	double *spvolMatL)	/* 16 'F' (O) output spot vol matrix */
{
static	char		routine[] = "VnfmSmoothSwaptionMatrixL";
	int		status = FAILURE;

	TCurve		*zcCurve = NULL;
	TSwaptionMatrix2D *midMkt = NULL;
	TSwaptionMatrix2D *bidToMid = NULL;
	TSwaptionMatrix2D *matWeight = NULL;
	TSwaptionMatrix2D *expWeight = NULL;
	TSwaptionMatrix2D *shrunkMidMkt= NULL;
    TSwaptionMatrix2D *shrunkBidToMid= NULL;
    TSwaptionMatrix2D *shrunkMatWeight= NULL;
    TSwaptionMatrix2D *shrunkExpWeight= NULL;
	int		optimType;
	int		normType;
	double		smoothParam;
	/*double		*nfParamsIn = NULL;*/
	TSwaptionMatrix2D *outputMkt = NULL;
	TSwaptionMatrix2D *spvolMat = NULL;
	double		*nfParamsOut = NULL;

	TSwaptionMatrix2D	*diffMat = NULL;	/* testing */

	VNFM_LOGOPEN

	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "w", "TF_SMOOTH_SWMAT",
		DRL_TDATE_L,  refDateL, "ZC_REFDATE",
		DRL_TDATE_L,  zcDateL, "ZC_DATES",
		DRL_FLOAT_L,  zcRateL, "ZC_RATES",
		DRL_LONG_L,   swTypeL, "SW_TYPE",
		DRL_FLOAT_L,  swMatL, "SW_MAT",
		DRL_FLOAT_L,  swExpL, "SW_EXP",
		DRL_FLOAT_L,  midMktL, "MID_MKT",
		DRL_FLOAT_L,  bidToMidL, "BIFMID",
		DRL_FLOAT_L,  matWeightL, "EXPWEI",
		DRL_FLOAT_L,  expWeightL, "MATWEI",
		DRL_LONG_L,   integerScalarsL, "INT_SCALARS",
		DRL_FLOAT_L,  doubleScalarsL, "DOUBLE_SCALARS",
		DRL_FLOAT_L,  nfParamsInL, "PARAMS_IN",
		DRL_FLOAT_L,  nfParamsOutL, "PARAMS_OUT",
		DRL_FLOAT_L,  outputMktL, "OUTMKT",
		DRL_FLOAT_L,  spvolMatL, "SPVOL",
		0);
	}



	/*
	 * Wrap arguments
	 */

	/* read the zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* get matrices */
	if (DrlTSwaptionMatrix2DWrapRead(&midMkt ,
		swTypeL, swMatL, swExpL, midMktL) != SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTSwaptionMatrix2DFpWrite(midMkt, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));
#endif

	if (DrlTSwaptionMatrix2DWrapRead(&bidToMid ,
		swTypeL, swMatL, swExpL, bidToMidL) != SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTSwaptionMatrix2DFpWrite(bidToMid, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));
#endif

	if (DrlTSwaptionMatrix2DWrapRead(&matWeight ,
		swTypeL, swMatL, swExpL, matWeightL) != SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "matWeight matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(matWeight, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));
#endif

	if (DrlTSwaptionMatrix2DWrapRead(&expWeight ,
		swTypeL, swMatL, swExpL, expWeightL) != SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "expWeight matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(expWeight, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));
#endif

	WRAP_CHECK_VECTOR_LEN(integerScalarsL, 2);
	WRAP_CHECK_VECTOR_LEN(doubleScalarsL, 1);

	WRAP_CHECK_VECTOR(nfParamsInL);

    /* shrink */

    if (shrinkMatrix( midMkt,
                      bidToMid,
                      matWeight,
                      expWeight,
                      &shrunkMidMkt,
		      &shrunkBidToMid,
                      &shrunkMatWeight,
                      &shrunkExpWeight
                 ) ISNT SUCCESS) goto done;


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "Shrunk midMkt matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(shrunkMidMkt,\
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD);\
	DrlFPrintf(vnfmFpLog, "Shrunk bidToMid matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(shrunkBidToMid,\
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD);
	DrlFPrintf(vnfmFpLog, "Shrunk matWeight matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(shrunkMatWeight,\
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD);\
	DrlFPrintf(vnfmFpLog, "Shrunk expWeight matrix:\n");\
	DrlTSwaptionMatrix2DFpWrite(shrunkExpWeight,\
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD));
#endif


	/*
	 * Do it
	 */
	if (VnfmSmoothSwaptionMatrix(
		zcCurve,
		shrunkMidMkt,
		shrunkBidToMid,
		shrunkMatWeight,
		shrunkExpWeight,
		integerScalarsL[1], 	/* optimType */
		integerScalarsL[2],	/* normType */
		doubleScalarsL[1],	/* smoothParam */
		nfParamsInL,		/* LIL array of parameters */
		midMkt,
		&outputMkt,
		&spvolMat,
		&nfParamsOutL[1]) != SUCCESS)
			goto done;


	/* test differnce */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	diffMat = DrlTSwaptionMatrix2DNewCopy(midMkt);\
	if (diffMat == NULL) goto done;\
	DrlTSwaptionMatrix2DOperScalar(diffMat, "=", 0e0);\
	DrlTSwaptionMatrix2DOperMatrix(diffMat, "+=", outputMkt);\
	DrlTSwaptionMatrix2DOperMatrix(diffMat, "-=", midMkt));

	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "INPUT_MAT:\n");\
	DrlTSwaptionMatrix2DFpWrite(midMkt, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));

	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "OUTPUT_MAT:\n");\
	DrlTSwaptionMatrix2DFpWrite(outputMkt, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));

	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "DIFF_MAT:\n");\
	DrlTSwaptionMatrix2DFpWrite(diffMat, vnfmFpLog,\
		TSWAPTION_MATRIX_FMT_STD));
#endif


	/* wrap out */
	if (DrlTSwaptionMatrix2DWrapWrite(outputMkt, NULL, NULL,
		NULL, outputMktL) != SUCCESS)
			goto done;

	if (DrlTSwaptionMatrix2DWrapWrite(spvolMat, NULL, NULL,
		NULL, spvolMatL) != SUCCESS)
			goto done;



	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	DrlTSwaptionMatrix2DFree(midMkt);
	DrlTSwaptionMatrix2DFree(bidToMid);
	DrlTSwaptionMatrix2DFree(matWeight);
	DrlTSwaptionMatrix2DFree(expWeight);
	DrlTSwaptionMatrix2DFree(outputMkt);
	DrlTSwaptionMatrix2DFree(spvolMat);
	DrlTSwaptionMatrix2DFree(diffMat);

	DrlTSwaptionMatrix2DFree(shrunkMidMkt);
	DrlTSwaptionMatrix2DFree(shrunkBidToMid);
	DrlTSwaptionMatrix2DFree(shrunkExpWeight);
	DrlTSwaptionMatrix2DFree(shrunkMatWeight);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*---------------------------------------------------------------
 * The matrix shrinker.
 */

static
long shrinkMatrix( TSwaptionMatrix2D *midMkt,
                   TSwaptionMatrix2D *bidToMid,
                   TSwaptionMatrix2D *matWeight,
                   TSwaptionMatrix2D *expWeight,
                   TSwaptionMatrix2D **shrunkMidMktO,
                   TSwaptionMatrix2D **shrunkBidToMidO,
                   TSwaptionMatrix2D **shrunkMatWeightO,
                   TSwaptionMatrix2D **shrunkExpWeightO
                 )
{
 static char routine[] = "shrinkMatrix";
	int	status = FAILURE;

 TTable2D *bidToMidTable = bidToMid->table;
 double   **bidToMidMtrx = bidToMidTable->matrix->data;
/*
 double   *resets =  bidToMidTable->dim1Values;
 double   *maturities = bidToMidTable->dim2Values;
*/
 long      numResets = bidToMidTable->matrix->numDim1;
 long      numMat = bidToMidTable->matrix->numDim2;
 long     *includeReset = NULL;
 long     *includeMat = NULL;
 long      resetIdx,matIdx,shrunkResIdx,shrunkMatIdx;
 long      newNumMat, newNumResets;
 TSwaptionMatrix2D *shrunkMidMkt = NULL;
 TSwaptionMatrix2D *shrunkBidToMid = NULL;
 TSwaptionMatrix2D *shrunkMatWeight = NULL;
 TSwaptionMatrix2D *shrunkExpWeight = NULL;


 shrunkMidMkt = NULL;
 shrunkBidToMid = NULL;

 includeReset = NEW_ARRAY(long, numResets);
 if (includeReset IS NULL)
  {
    GtoErrMsg("%s: could not allocate memory for includeReset",routine);
    goto done;
  }
 includeMat = NEW_ARRAY(long, numMat);
 if (includeMat IS NULL)
  {
    GtoErrMsg("%s: could not allocate memory for includeMat",routine);
    goto done;
  }
 /* columns or rows that have all negative bidToMids are discarded from 
   the matrix. Here we find them 
 */

 for (resetIdx = 0; resetIdx < numResets; resetIdx++)
  for ( matIdx = 0; matIdx < numMat; matIdx++)
     if ( bidToMidMtrx[resetIdx][matIdx]>=0 ) 
      {
       includeReset[resetIdx]=1;
       includeMat[matIdx]=1; 
      }      

 for (newNumResets=0, resetIdx = 0; resetIdx < numResets; resetIdx++)
   if ( includeReset[resetIdx] IS 1) newNumResets++; 

 for (newNumMat=0, matIdx = 0; matIdx < numMat; matIdx++)
   if ( includeMat[matIdx] IS 1) newNumMat++;

 shrunkMidMkt = GtoSwaptionMatrix2DNewEmpty( FALSE,
                                             midMkt->swapPayFreq,
                                             newNumResets,
                                             newNumMat);
 if ( shrunkMidMkt IS NULL)
 {
    GtoErrMsg("%s: could not create srunkMidMkt \n",routine);
    goto done;
  }

  shrunkBidToMid = GtoSwaptionMatrix2DNewEmpty( FALSE,
                                              midMkt->swapPayFreq,
                                              newNumResets,
                                              newNumMat);
 if ( shrunkBidToMid IS NULL)
 {
    GtoErrMsg("%s: could not create shrunkBidToMid\n",routine);
    goto done;
  }
  shrunkExpWeight = GtoSwaptionMatrix2DNewEmpty( FALSE,
                                              midMkt->swapPayFreq,
                                              newNumResets,
                                              newNumMat);
 if ( shrunkExpWeight IS NULL)
 {
    GtoErrMsg("%s: could not create shrunkExpWeight\n",routine);
    goto done;
  }
  shrunkMatWeight = GtoSwaptionMatrix2DNewEmpty( FALSE,
                                                 midMkt->swapPayFreq,
                                                 newNumResets,
                                                 newNumMat);
 if ( shrunkMatWeight IS NULL)
 {
    GtoErrMsg("%s: could not create shrunkMatWeight\n",routine);
    goto done;
  }


 for (resetIdx = 0, shrunkResIdx=0; resetIdx < numResets; resetIdx++)
 {
  if ( includeReset[resetIdx] IS 1)
  {
   shrunkMidMkt->table->dim1Values[shrunkResIdx] =  
                                           bidToMidTable->dim1Values[resetIdx];
   shrunkBidToMid->table->dim1Values[shrunkResIdx] =  
                                           bidToMidTable->dim1Values[resetIdx];
   shrunkExpWeight->table->dim1Values[shrunkResIdx] =  
                                           bidToMidTable->dim1Values[resetIdx];     shrunkMatWeight->table->dim1Values[shrunkResIdx] =  
                                           bidToMidTable->dim1Values[resetIdx];
   shrunkResIdx++;
  }
 }

 for (matIdx = 0, shrunkMatIdx=0; matIdx < numMat; matIdx++)
 {
  if ( includeMat[matIdx] IS 1)
  {
   shrunkMidMkt->table->dim2Values[shrunkMatIdx] =  
                                           bidToMidTable->dim2Values[matIdx];
   shrunkBidToMid->table->dim2Values[shrunkMatIdx] =  
                                           bidToMidTable->dim2Values[matIdx];
   shrunkExpWeight->table->dim2Values[shrunkMatIdx] =  
                                           bidToMidTable->dim2Values[matIdx];
   shrunkMatWeight->table->dim2Values[shrunkMatIdx] =  
                                           bidToMidTable->dim2Values[matIdx];
   shrunkMatIdx++;
  }
 }

 for (resetIdx = 0, shrunkResIdx=0 ; resetIdx < numResets; resetIdx++)
 {
   if (includeReset[resetIdx] IS 1)
   {
     for ( matIdx = 0, shrunkMatIdx = 0 ; matIdx < numMat; matIdx++)
     {
       if (includeMat[matIdx] IS 1)
       {
        shrunkMidMkt->table->matrix->data[shrunkResIdx][shrunkMatIdx] = 
                           midMkt->table->matrix->data[resetIdx][matIdx];
        shrunkBidToMid->table->matrix->data[shrunkResIdx][shrunkMatIdx] = 
                           bidToMid->table->matrix->data[resetIdx][matIdx];
        shrunkExpWeight->table->matrix->data[shrunkResIdx][shrunkMatIdx] = 
                           expWeight->table->matrix->data[resetIdx][matIdx];
        shrunkMatWeight->table->matrix->data[shrunkResIdx][shrunkMatIdx] = 
                           matWeight->table->matrix->data[resetIdx][matIdx];
        shrunkMatIdx++; 
       }           
     }
    shrunkResIdx++;
   }
 }



	*shrunkMidMktO = shrunkMidMkt;
	*shrunkBidToMidO = shrunkBidToMid;
	*shrunkMatWeightO = shrunkMatWeight;
	*shrunkExpWeightO = shrunkExpWeight;


  status = SUCCESS;

  done:

  FREE(includeReset);
  FREE(includeMat);

  if (status ISNT SUCCESS)
  {
   GtoErrMsg("%s: failed.\n", routine);
   GtoSwaptionMatrix2DFree(shrunkMidMkt);
   GtoSwaptionMatrix2DFree(shrunkBidToMid);
   GtoSwaptionMatrix2DFree(shrunkExpWeight);
   GtoSwaptionMatrix2DFree(shrunkMatWeight);
  }

 return status;
}

