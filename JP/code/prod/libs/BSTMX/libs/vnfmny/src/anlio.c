/****************************************************************
 * Module:	VNFM
 * Submodule:	ANLY
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>
#include <string.h>
#ifndef	_WIN32
# include <errno.h>
#endif

#include "date_sup.h"
#include "convert.h"

#include "drlio.h"
#include "drlstr.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"

#define xDrKey

/*
 * Global File Pointer for Logging
 */
#if defined(CDDEV)
	FILE	*vnfmFpLog = stdout;
#else
	FILE	*vnfmFpLog = NULL;		/* default to GtoErrMsg */
#endif


/*f-------------------------------------------------------------
 * Read a VnfmData from a file.
 *                                                             
 * <br><br>
 * Reads a set of parameters in a file "fnam" and put the 
 * result in the data structure "that".
 * Returns 0 iff OK.
 */

int
VnfmFileRead(VnfmData **that, char *fnam)
{
static	char	routine[] = "VnfmFileRead";
	FILE	*fp;
	int	status;

	if ((fp = fopen(fnam, "r")) == NULL) {
#ifdef	_WIN32
		GtoErrMsg("%s: can't open `%s'.\n",
			routine, fnam);
#else
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, strerror(errno));
#endif
		return(-1);
	} else {
		status = VnfmFpRead(that, fp);
		fclose(fp);
		return(status);
	}
}


/*f-------------------------------------------------------------
 * Write a VnfmData to a file.
 *                                                             
 * <br><br>
 * Writes a set of parameters from data structure "that"
 * in a file "fnam". \\
 * WARNING: if the file already exists, it is truncated to zero.
 * Returns 0 iff OK.
 */

int
VnfmFileWrite(VnfmData *that, char *fnam)
{
static	char	routine[] = "VnfmFileWrite";
	FILE	*fp;
	int	status;

	if ((fp = fopen(fnam, "w")) == NULL) {
#ifdef	_WIN32
		GtoErrMsg("%s: can't open `%s'.\n",
			routine, fnam);
#else
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, strerror(errno));
#endif
		return(FAILURE);
	} else {
		status = VnfmFpWrite(that, fp) ;
		fclose(fp);
		return(status);
	}
}


/*f-------------------------------------------------------------
 * Read a nfmData from a file pointer.
 *                                                             
 * <br><br>
 * Reads a set of parameters in a file pointer "fp" and put the 
 * result in the data structure "that".
 * Returns 0 iff OK.
 */

int
VnfmFpRead(VnfmData **that, FILE *fp)
{
static	char	routine[] = "VnfmFpRead";
	char	buf[2048];
	int	status = FAILURE,
		n, i, j, k, nDim, nDates,
		line = 0;
	double	beta[VNFM_NDIMMAX],
		alpha[VNFM_NDIMMAX];
	double	backboneq;



	/* */
	status = DrlFScanStruct(fp, &line, " \t",
		DRL_CVAR_T, "backbone-q", DRL_DOUBLE_T, (void*) &backboneq,
		DRL_CVAR_T, "nDim",     DRL_INT_T,  (void*) &nDim,
		DRL_CARRAY_T, "beta", &n, (int) 2,
			DRL_DOUBLE_T, FALSE, (void*) beta,
			DRL_DOUBLE_T, FALSE, (void*) alpha,
		DRL_CVAR_T, "nDates",   DRL_INT_T,  (void*) &nDates,
		0);
	if (status != 0) goto done;

	/* check distribution */

	if (n != nDim) {
		GtoErrMsg("%s: number of MR (%d) != number of factors (%d)\n",
			routine, n, nDim);
		goto done;
	}

	/* allocate memory */
	if (nDates < 0) {
		GtoErrMsg("%s: bad number of dates (%d)\n", routine, nDates);
		goto done;
	}
	if ((nDim < 0) || (nDim >= VNFM_NDIMMAX)) {
		GtoErrMsg("%s: bad number of dim (%d)\n", routine, nDim);
		goto done;
	}

	/* allocate memory */
	if ((*that = VnfmNew(nDates, nDim, 0)) == NULL) goto done;
	nDim = (*that)->fNf;
	(*that)->fBackBoneQ = 1e0 - 2e0 * backboneq;

	/* mr coeffs */
	for (j=0; j<=(*that)->fNf-1; j++) {
		(*that)->fBeta[j] = beta[j];
		(*that)->fAlpha[j] = alpha[j];
	}


	/* get time dependent volatilities and correlations */
	for (i=0; i<=(*that)->fNDates-1; i++) {
		if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL)
			goto done;

		if (DrlStrTokScan(buf, DRL_TDATE_T, &(*that)->fDate[i]) != 0)
			goto done;
		for (j=0; j<=nDim-1; j++) {
			if (DrlStrTokScan(NULL, DRL_DOUBLE_T,
				&(*that)->fSigma[j][i]) != 0)
				goto done;
		}

		for (j=0; j<=nDim-1; j++) 
		for (k=j+1; k<=nDim-1; k++) {
			if (DrlStrTokScan(NULL, DRL_DOUBLE_T,
				&(*that)->fRho[RHOIDX(j,k)][i]) != 0)
				goto done;
		}

	}

	/* compute day count fraction */
	for (i=0; i<=(*that)->fNDates-1; i++) {
		GtoDayCountFraction((*that)->fDate[0], (*that)->fDate[i],
			GTO_ACT_365F, &(*that)->fTime[i]);
	}


	/*
	 *
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: error line %d, code %d\n",
			routine, line, status);
		VnfmFree(*that);
		*that = NULL;
	}
	return(status);
#undef	CHECK
}


/*f-------------------------------------------------------------
 * Write a VnfmData to a file pointer.
 *                                                             
 * <br><br>
 * Writes a set of parameters from data structure "that"
 * in a file pointer "fp".
 * Returns 0 iff OK.
 */

int
VnfmFpWrite(VnfmData *that, FILE *fp)
{
static	char	routine[] = "VnfmFpWrite";
	char	buf[2048];
	int	status = FAILURE,
		i, j, k, nDim = that->fNf,
		line = 0;
	double	backboneq = (1e0 - that->fBackBoneQ) * 0.5e0;


	/* */
	DrlFPrintStruct(fp, &line, '\t',
		DRL_CVAR_T, "backbone-q", DRL_DOUBLE_T, (void*) &backboneq,
		DRL_CVAR_T, "nDim",     DRL_INT_T,  (void*) &NDIM,
		DRL_CARRAY_T, "beta", NDIM, (int) 2,
			DRL_DOUBLE_T, (void*) that->fBeta,
			DRL_DOUBLE_T, (void*) that->fAlpha,
		DRL_CVAR_T, "nDates",   DRL_INT_T,  (void*) &that->fNDates,
		0);

	/* print time dependent volatilities and correlations */
	for (i=0; i<=that->fNDates-1; i++) {

		strcpy(buf, "");
		DrlStrTokPrint(buf, DRL_TDATE_T, &that->fDate[i]);
		for (j=0; j<=NDIM-1; j++) {
			DrlStrTokPrint(buf, DRL_DOUBLE_T, &that->fSigma[j][i]);
		}

		for (j=0; j<=NDIM-1; j++) 
		for (k=j+1; k<=NDIM-1; k++) {
			DrlStrTokPrint(buf, DRL_DOUBLE_T, 
				&that->fRho[RHOIDX(j,k)][i]);
		}
		strcat(buf, "\n");

		DrlFPrintf(fp, "%s", buf);

	}


	/* we are done */
	status = SUCCESS;
/*done:*/
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);

}




/*f-------------------------------------------------------------
 * Prints all cached coefficients to a file pointer.
 *                                                             
 * <br><br>
 * Prints all parameters and internally stored coefficients on the
 * file pointer "fp" (for debugging purposes).
 * Returns 0 iff OK.
 */
int
VnfmPrintCoeff(
	VnfmData *that,	/* (I) model parameters */
	FILE *fp)
{
static	char	routine[] = "VnfmPrintCoeff";
	int	i, j, m,
		rIdx, jIdx,
		nDim = NDIM;

	DrlFPrintf(fp, "%s: \n", routine);
	DrlFPrintf(fp, "Bbq: %lf\n", that->fBackBoneQ);
	DrlFPrintf(fp, "VolNorm: %lf\n", that->fVolNorm);
	DrlFPrintf(fp, "VolLogn: %lf\n", that->fVolLogn);
	DrlFPrintf(fp, "QLeft:   %lf\n", that->fQLeft);
	DrlFPrintf(fp, "QRight:  %lf\n", that->fQRight);
	DrlFPrintf(fp, "QFShift: %lf\n", that->fQFShift);


	DrlFPrintf(fp, "Beta:\t");
	for (j=0; j<=nDim-1; j++) 
		DrlFPrintf(fp, "%7.4f ", BETA[j]);
	DrlFPrintf(fp, "\n");

	DrlFPrintf(fp, "Alpha:\t");
	for (j=0; j<=nDim-1; j++) 
		DrlFPrintf(fp, "%7.4f ", ALPHA[j]);
	DrlFPrintf(fp, "\n");

	DrlFPrintf(fp, "Volatilities: nDates = %d\n", NDATES);
	for(m=0; m<=NDATES-1; m++) {
		DrlFPrintf(fp, "[%3d,%12.8f,%10s]",
			m,
			that->fTime[m],
			GtoFormatDate(that->fDate[m]));

		for (j=0; j<=nDim-1; j++) 
			DrlFPrintf(fp, "%7.4f ", SIGMA[j][m]);
		DrlFPrintf(fp, "|");

		for (i=0; i<=nDim-1; i++)
		for (j=i+1; j<=nDim-1; j++) {
			rIdx = RHOIDX(i, j);
			DrlFPrintf(fp, "%7.4f ", RHO[rIdx][m]);
		}

		DrlFPrintf(fp, "\n");
	}



	/* print out fZTime , fRate and fZero */
	DrlFPrintf(fp, "Zero Curve Timeline Info:\n");
	DrlFPrintf(fp, "[idx,time ,date      ] Rate    Zero    \n");
	DrlFPrintf(fp, "[%3d,%5.2f,%10s]",
		0, ZTT[0], GtoFormatDate(REFDATE));
	DrlFPrintf(fp, "%7.4f, %7.4f \n",
		that->fRate[0], that->fZero[0]);

	for (m=1; m<=NZDATES-1;m++) {
		DrlFPrintf(fp, "[%3d,%5.2f,%10s]",
			m,
			ZTT[m],
#ifndef	VNFM_V5X
			GtoFormatDate((that->fZcCurve)->fArray[m-1].fDate));
#else
			GtoFormatDate(that->fDate[m]));
#endif	/*VNFM_V5X*/
		DrlFPrintf(fp, "%7.4f ", that->fRate[m]);
		DrlFPrintf(fp, "%7.4f ", that->fZero[m]);
		DrlFPrintf(fp, "\n");
	}




	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		DrlFPrintf(fp, " J[%d,%d]=%2d ",
			i, j, jIdx);
	}
	DrlFPrintf(fp, "\n");

	DrlFPrintf(fp, "Coefficients: nDates = %d\n", NDATES);
	for(m=0; m<=NDATES-1; m++) {
		DrlFPrintf(fp, "[%3d,%5.2f,%10s]",
			m,
			that->fTime[m],
			GtoFormatDate(that->fDate[m]));

		for (i=0; i<=nDim-1; i++)
		for (j=i; j<=nDim-1; j++) {
			jIdx = JIDX(i, j);
			DrlFPrintf(fp, " %12.8f ",
				that->fJ[jIdx][m]);
		}

		DrlFPrintf(fp, "\n");
	}



	return(0);
}



/*f-------------------------------------------------------------
 * Read a VnfmData from a LIL interface.
 *                                                             
 * <br><br>
 * Reads a set of parameters from a wrapper interface
 * and creates a new parameter data structure "that".
 * <br>
 * <br>[backboneqL] LIL array of length 1, containting the
 * backbone $q$ coefficient (0 for lognormal, 0.5 for normal)
 * <br>[betaL] LIL-array of length equal to the number of dimensions
 *  containing the mean-reversions,
 * <br>[alphaL] LIL-array of length equal to the number of dimensions
 *  containing the factor weights,
 * <br>[dateL] LIL-array of length $M$ (arbitrary) containing 
 * the timeline dates,
 * <br>[sigmaL] LIL-array containing the spot volatilities: should
 * be a range of dimension $Nx M$ with the factor index in
 * X-coordinate (horizontal) and timeline index in Y-coordinate (vertical),
 * <br>[rhoL] LIL-array containing the spot volatilities: should
 * be a range of dimension $Nx {M(M-1)/ 2}$ with the correlation 
 * index $k$ in the  X-coordinate,
 * where the correlation between factors $i$ and $j$  has index 
 * $$ k = {N(N-1)/ 2} - {(N-i)*(N-i-1)/ 2} + (j-i-1),$$
 * and timeline index in Y-coordinate (vertical),
 * <br>[zcCurve] the zero curve to be used in the parameter
 * computations.
 * <br>
 * Returns 0 iff OK.
 */

int
VnfmWrapRead(
	VnfmData **that,	/*     (O) parameters */
	FloatL *backboneqL,	/* 'F' (I) back bone Q */
	FloatL *betaL,		/* 'F' (I) array of mr coeff */
	FloatL *alphaL,		/* 'F' (I) array of weight coeff */
	TDateL *dateL,		/* 'D' (I) array of dates */
	FloatL *sigmaL,		/* 'F' (I) volatility arrays */
	FloatL *rhoL,		/* 'F' (I) correlation arrays */
	TCurve *zcCurve)	/*     (I) zero curve */
{
static	char	routine[] = "VnfmWrapRead";
	int	status = FAILURE,
		i, j, k, rIdx, nDim, nDates;


extern	int	DrSecurity(unsigned int c);

#ifdef DrKey
	/*DR_SECURITY_CHECK();*/
	if (DrSecurity((unsigned int) 2) != SUCCESS)
		goto done;
#endif


#undef	CHECK
#define	CHECK(cond)	{if (!(cond)) {GtoErrMsg("%s: assertion failed (%s)\n",\
			routine, #cond); goto done;}}

	*that = (VnfmData*)NULL;

	CHECK(ARGSIZE(backboneqL) == 1);

	/* get # of dimensions */
	nDim = ARGSIZE(betaL);
	if ((nDim <= 0) || (nDim >= VNFM_NDIMMAX)) {
		GtoErrMsg("%s: bad dimension (%d)\n", routine, nDim);
		goto done;
	}


	/* check array of arrays */
	CHECK(ARGSIZE(betaL) == (unsigned)nDim);
	CHECK(ARGSIZE(alphaL) == (unsigned)nDim);
	CHECK(ARGSIZE(sigmaL) == (unsigned)(nDim*ARGSIZE(dateL)));
	CHECK(zcCurve ISNT NULL);
	CHECK(zcCurve->fNumItems >= 1);

	if (nDim > 1)  {
	    CHECK(ARGSIZE(rhoL) == (unsigned)(nDim*(nDim-1)*ARGSIZE(dateL)/2));
	}

	/* get number of nonzero dates */
	for (i=0; (unsigned)i<=ARGSIZE(dateL)-1 && dateL[i+1] > 1; i++)
		nDates = i+1;
	if (nDates <= 0) {
		GtoErrMsg("%s: bad dates array\n", routine, nDates);
		goto done;
	}


	/* now allocate memory */
	*that = VnfmNew(nDates, nDim, zcCurve->fNumItems+1);
	if (*that == NULL) goto done;
	(*that)->fBackBoneQ = 1e0 - 2e0 * backboneqL[1];

	/* fill the data structure */
	for(j=0; j<=(*that)->fNf-1; j++) {
		(*that)->fBeta[j] = betaL[j+1];
		(*that)->fAlpha[j] = alphaL[j+1];
	}

	for (i=0; i<=(*that)->fNDates-1; i++) {
		(*that)->fDate[i] = dateL[i+1];
		GtoDayCountFraction((*that)->fDate[0], (*that)->fDate[i],
			GTO_ACT_365F, &(*that)->fTime[i]);

		for(j=0; j<=(*that)->fNf-1; j++) {
			(*that)->fSigma[j][i] = sigmaL[1+j+nDim*i];
		}

		for(j=0; j<=(*that)->fNf-1; j++)
		for(k=j+1; k<=(*that)->fNf-1; k++) {
			rIdx = RHOIDX(j,k);
			(*that)->fRho[rIdx][i] =
				rhoL[1+rIdx+nDim*(nDim-1)*i/2];
		}
	}

	/* set the zero curve */
	(*that)->fZcCurve = zcCurve;

#ifndef	VNFM_V5X
	/* set fZTime */
	(*that)->fZTime[0] = 0.0;    /* equals REFDATE or (*that)->fDate[0] */
	for (i=1;i<=zcCurve->fNumItems;i++) {
	    GtoDayCountFraction((*that)->fDate[0], zcCurve->fArray[i-1].fDate,
		                      GTO_ACT_365F, &((*that)->fZTime[i]));
	}    
#endif	/*VNFM_V5X*/


	if (VnfmCheckValid(*that) != SUCCESS) goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != 0) {
	    if (*that != (VnfmData*)NULL) VnfmFree(*that);
	    *that = (VnfmData*)NULL;
	    GtoErrMsg("%s: failed (code %d)\n", routine, status);
	}
	return(status) ;
#undef	CHECK
}


/*f-------------------------------------------------------------
 * Write a VnfmData to a LIL interface.
 *
 *                                                             
 * <br><br>
 * Writes a set of parameters from the data structure "that"
 * to a wrapper interface.
 * For a description of the output arguments, refer to the 
 * description of the function <i> VnfmWrapRead</i>.
 * Returns 0 iff OK.
 */

int
VnfmWrapWrite(
	VnfmData *that,		/* (I)     parameters */
	FloatL *backboneqL,	/* (O) 'F' back bone Q (or NULL) */
	FloatL *betaL,		/* (O) 'F' array of mt coeff (or NULL) */
	FloatL *alphaL,		/* (O) 'F' array of weight coeff (or NULL) */
	TDateL *dateL,		/* (O) 'D' array of dates */
	FloatL *sigmaL,		/* (O) 'F' volatility arrays */
	FloatL *rhoL)		/* (O) 'F' correlation arrays */
{
static	char	routine[] = "VnfmWrapWrite";
	int	status = FAILURE,
		i, j, k, rIdx, nDim = that->fNf;

#undef	CHECK
#define	CHECK(cond)	{if (!(cond)) {GtoErrMsg("%s: assertion failed "\
			"(%s)\n", routine, #cond); goto done;}}

	if (backboneqL != NULL)     { CHECK(ARGSIZE(backboneqL)  == 1); }
	if (betaL != NULL)     { CHECK(ARGSIZE(betaL)  >= 1); }
	if (alphaL != NULL)    { CHECK(ARGSIZE(alphaL) >= 1); }
	CHECK(ARGSIZE(dateL)  >= 1);
	CHECK(ARGSIZE(sigmaL) >= 1);
	if (nDim > 1) {
		CHECK(ARGSIZE(rhoL)   >= 1);
	}

	/* get # of dimensions */

	/* check array of arrays */
	if (ARGSIZE(dateL) < (unsigned)that->fNDates) {
		GtoErrMsg("%s: not enough dates in array (expect %d, got %d)\n",
			routine, that->fNDates, ARGSIZE(dateL));
		goto done;
	}
	if (betaL  != NULL) { CHECK(ARGSIZE(betaL) == (unsigned)nDim); }
	if (alphaL != NULL) { CHECK(ARGSIZE(alphaL) == (unsigned)nDim); }
	CHECK(ARGSIZE(sigmaL) == nDim*ARGSIZE(dateL));

	if (nDim > 1) {
	    CHECK(ARGSIZE(rhoL) == nDim*(nDim-1)*ARGSIZE(dateL)/2);
	}




	/* fill the output arguments */
 	if (backboneqL != NULL) {
 	    backboneqL[1] = (1e0 - that->fBackBoneQ)*0.5e0;
	}

	if (betaL != NULL) {
	    for(j=0; j<=that->fNf-1; j++)
		betaL[j+1]  = that->fBeta[j];
	}
	if (alphaL != NULL) {
	    for(j=0; j<=that->fNf-1; j++)
		alphaL[j+1] = that->fAlpha[j];
	}

	for (i=0; i<=that->fNDates-1; i++) {
		dateL[i+1] = that->fDate[i];

		for(j=0; j<=that->fNf-1; j++) {
			sigmaL[1+j+nDim*i] = that->fSigma[j][i];
		}

		for(j=0; j<=that->fNf-1; j++)
		for(k=j+1; k<=that->fNf-1; k++) {
			rIdx = RHOIDX(j,k);
			rhoL[1+rIdx+nDim*(nDim-1)*i/2] =
				that->fRho[rIdx][i];
		}
	}
	for (i=that->fNDates; (unsigned)i<=ARGSIZE(dateL)-1; i++) {
		dateL[i+1] = 1L;
	}


	/* made it through */
	status = SUCCESS;
done:
	if (status != 0) {
		GtoErrMsg("%s: failed (code %d)\n", routine, status);
	}
	return(status);
#undef	CHECK
}

/*f-------------------------------------------------------------
 * Read a VnfmData from a LIL interface (TCurve given).
 *                                                             
 * <br><br>
 * Reads a set of parameters from a wrapper interface
 * and creates a new parameter data structure <i> thatp</i>.
 * A simplified version of <i> VnfmWrapRead</i> that does not
 * allow time dependent correlations and volatilities to be entered.
 * <br>
 * <br>[that] on successful exit, points to a newly allocated structure.
 * <br>[backboneqL] LIL array of length 1, containting the
 * backbone $q$ coefficient (0 for lognormal, 0.5 for normal)
 * <br>[paramsL] LIL array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * or an array format where the 1st weight assumed to be 1.0 is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * <br>[dateL] LIL-array of length $M$ (arbitrary) containing 
 * the timeline dates.
 * <br>[zcCurve] the zero curve to be used in the parameter
 * computations.
 * <br>
 * Returns 0 iff OK.
 */

int
VnfmWrapReadSimple(
	VnfmData **thatp,	/*     (O) parameters */
	FloatL *backboneqL,	/* 'F' (I) back bone Q */
	FloatL *paramsL,	/* 'F' (I) array of mr, weights, corrs */
	TDateL *datesL,		/* 'D' (I) array of dates */
	TCurve *zcCurve)	/*     (I) zero curve */
{
static	char	routine[] = "VnfmWrapReadSimple";
	int	status = FAILURE;

	int		nF, nDim, nDates, idxF1, idxF2, idxT, idxR;
	VnfmData	*that = NULL;
extern	int		DrSecurity(unsigned int c);

#define	NFA(nParams)	(sqrt(2e0*((double)(nParams))+2.25e0)-1.5e0)
#define	FNFA(nParams)	IS_ALMOST_ZERO((nF = (int) NFA(nParams)) - NFA(nParams))
#define	NFM(nParams)	(sqrt(((double)(nParams))+1e0)-1e0)
#define	FNFM(nParams)	IS_ALMOST_ZERO((nF = (int) NFM(nParams)) - NFM(nParams))


	WRAP_CHECK_VECTOR(backboneqL);
	WRAP_CHECK_VECTOR(paramsL);
	WRAP_CHECK_VECTOR(datesL);


#ifdef DrKey
	/* security check */
	if (DrSecurity((unsigned int) 2) != SUCCESS)
		goto done;
#endif

	/* get number of nonzero dates */
	for (idxT=0;
		(unsigned)idxT <= ARGSIZE(datesL)-1 && datesL[idxT+1] > 1;
		idxT++)
			nDates = idxT+1;
	if (nDates <= 0) {
		GtoErrMsg("%s: empty dates array.\n", routine, nDates);
		goto done;
	}


	/* */
	if (FNFA(ARGSIZE(paramsL))) {
		/**
		 ** array (b_1,..,b_n, alpha_1,...,alpha_n,
		 **        rho_{1,2},..,rho_{1,n},rho_{2,1},...,rho_{n-1,n}).
		 **/

		if ((that = VnfmNew(nDates, nF, zcCurve->fNumItems+1)) == NULL)
			goto done;
		nDim = that->fNf;

		/* copy parameters */
		for (idxF1=0; idxF1<=nF-1; idxF1++) {
			that->fBeta[idxF1] = (IS_ALMOST_ZERO(paramsL[1+idxF1]) ?
						1e-8 : paramsL[1+idxF1]);
			that->fAlpha[idxF1] = paramsL[1+idxF1+nF];
		}


		for (idxF1=0;       idxF1<=nF-1; idxF1++)
		for (idxF2=idxF1+1; idxF2<=nF-1; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			for (idxT=0; idxT<=that->fNDates-1; idxT++) {
				that->fRho[idxR][idxT] = paramsL[1+idxR+nF+nF];
			}
		}

	} else if (FNFA(ARGSIZE(paramsL)+1)) {
		/** ARRAY FORM WITH NO WEIGHT 1
		 ** array (b_1,..,b_n, alpha_2,...,alpha_n,
		 **        rho_{1,2},..,rho_{1,n},rho_{2,1},...,rho_{n-1,n}).
		 **/

		if ((that = VnfmNew(nDates, nF, zcCurve->fNumItems+1)) == NULL)
			goto done;
		nDim = that->fNf;

		/* copy parameters */
		for (idxF1=0; idxF1<=nF-1; idxF1++) {
			that->fBeta[idxF1] = (IS_ALMOST_ZERO(paramsL[1+idxF1]) ?
						1e-8 : paramsL[1+idxF1]);
		}
		that->fAlpha[0] = 1e0;
		for (idxF1=1; idxF1<=nF-1; idxF1++) {
			that->fAlpha[idxF1] = paramsL[1+idxF1+nF-1];
		}


		for (idxF1=0;       idxF1<=nF-1; idxF1++)
		for (idxF2=idxF1+1; idxF2<=nF-1; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			for (idxT=0; idxT<=that->fNDates-1; idxT++) {
			    that->fRho[idxR][idxT] = paramsL[1+idxR+nF+nF-1];
			}
		}


	} else if (FNFM(ARGSIZE(paramsL))) {
		/** TABLE FORM
		 ** b_1,..,b_n, alpha_1,...,alpha_n,
		 ** rho_{1,1},..,rho_{1,n},rho_{2,1},...,rho_{n,n}.
		 **/
		if ((that = VnfmNew(nDates, nF, zcCurve->fNumItems+1)) == NULL)
			goto done;
		nDim = that->fNf;

		/* copy parameters */
		for (idxF1=0; idxF1<=nF-1; idxF1++) {
			that->fBeta[idxF1] = (IS_ALMOST_ZERO(paramsL[1+idxF1]) ?
						1e-8 : paramsL[1+idxF1]);
		}
		for (idxF1=0; idxF1<=nF-1; idxF1++) {
			that->fAlpha[idxF1] = paramsL[1+idxF1+nF];
		}


		for (idxF1=0;       idxF1<=nF-1; idxF1++)
		for (idxF2=idxF1+1; idxF2<=nF-1; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			for (idxT=0; idxT<=that->fNDates-1; idxT++) {
				that->fRho[idxR][idxT] = paramsL[1+2*nF+
					nF*idxF1+idxF2];
			}
		}


	} else {
		GtoErrMsg("%s: Parameters array has bad size (%d).\n",
			routine);
		GtoErrMsg("%s: Must be n(n+3)/2, n(n+3)/2-1 or n(n+2).\n",
			routine);
		goto done;
	}



	/* Set time dependent volatility */
	for (idxF1=0;       idxF1<=nF-1; idxF1++)
	for (idxT=0; idxT<=that->fNDates-1; idxT++) {
		that->fSigma[idxF1][idxT] = 1e0;
	}

	/* Set other parameters */
	that->fBackBoneQ = 1e0 - 2e0 * backboneqL[1];

	/* set dates in timeline */
	for (idxT=0; idxT<=that->fNDates-1; idxT++) {
		that->fDate[idxT] = datesL[idxT+1];
		if (GtoDayCountFraction(REFDATE, that->fDate[idxT],
			GTO_ACT_365F, &that->fTime[idxT]) != SUCCESS)
				goto done;
	}
	for (idxT=0; idxT<=that->fNDates-2; idxT++) {
		that->fDt[idxT] = that->fTime[idxT+1] - that->fTime[idxT];
	}
	that->fDt[that->fNDates-1] = 0.;


	/* Set zero curve */
	that->fZcCurve = zcCurve;

#ifndef	VNFM_V5X
	/* Set fZTime */
	that->fZTime[0] = 0.0;    /* equals REFDATE */
	for (idxT=1; idxT<=NZDATES-1; idxT++) {
		GtoDayCountFraction(REFDATE, zcCurve->fArray[idxT-1].fDate,
			GTO_ACT_365F, &that->fZTime[idxT]);
	}
#endif	/*VNFM_V5X*/


	/* Check valid paramerers */
	if (VnfmCheckValid(that) != SUCCESS) goto done;

	/* set output pointer */
	*thatp = that;

	/* made it through */
	status = SUCCESS;
done:
	if (status != 0) {
	    if (that != (VnfmData*)NULL) VnfmFree(that);
	    *thatp = (VnfmData*)NULL;
	    GtoErrMsg("%s: failed.\n", routine);
	}
	return(status) ;
}



/*f-------------------------------------------------------------
 * ALIB array format parameter conversion.
 *                                                             
 * <br><br>
 * This routine extracts the parameters
 * (mean-reversion coefficients, spot volatilities and correlations)
 * from the calibration data <i> that</i> and puts it in arrays
 * compatible with the ALIB conventions
 * (i.e. the spot volatility $i$ applies between date $i-1$ and date $i$)
 * It is to be used with the ALIB routine <i> GtoVolDefIRNew</i>.\\
 * <br>
 * <br>[meanRevers] on exit, contains mean reversion coefficients
 *		(should be allocated before routine call).
 * <br>[backboneq] on exit, 
 *      the back bone coefficient $q$ (0 for log-normal, 0.5 for normal).
 * %<br>[numFact] size of arrays passed to routines
 * %	(should be larger or equal to actual number of factors
 * %	in <i> that</i>).
 * <br>[today] on exit, contains the volatility start date.
 * <br>[dates] on exit, contains dates where spot volatility and
 *	correlations are changing,
 *	constistent with the ALIB conventions,
 *	i.e. the spot volatility $i$ applies between date $i-1$ and date $i$.
 * 	(should be larger or equal to actual number of dates in <i> that</i>).
 * <br>[spotVols] on exit, 2D array spot volatilities,
 * 	indexed  [0..numFact-1][0..numDates-1].
 * 	(should be larger or equal to actual number of dates in <i> that</i>).
 * <br>[correlations] on exit, 2D array of correlations,
 * 	indexed [0..nF*(nF-1)/2-1][0..nDates-1].
 * 	(should be larger or equal to actual number of dates in <i> that</i>).
 * %<br>[numDates] size of arrays passed to routine
 * %	(should be larger or equal to actual number of dates in <i> that</i>).
 * <br>
 * The number of factors can be retrieved from the structure
 * <i> VtfnData</i> in the field named <i> fNf</i>.
 * The number of dates in the field named <i> fNDates</i>.
 * Returns SUCCESS/FAILURE.
 */


DLL_EXPORT(int)
VnfmToGtoParameters(
	VnfmData *that,		/* (I) model parameters */
	double *backboneq,	/* (O) back bone Q */
	double *meanRevers,	/* (O) array of mean revers [0..numFact-1] */
	TDate *today,		/* (O) where volatility starts */
	TDate *dates,		/* (O) dates up to which spot vols apply */
	double **spotVols,	/* (O) sp vols[0..numFact-1][0..numDates-1] */
	double **correlations)	/* (O) corr[0..nF*(nF-1)/2-1][0..nDates-1] */
{
static	char	routine[] = "VnfmToGtoParameters";
	int	status = FAILURE;

	int	idxF, numSpvol, idxS;

#ifdef	_SKIP
	/*
	 * Check sizes 
	 */
	if (numFact < that->fNf) {
		GtoErrMsg("%s: array not large enough: "
			"numFact (%d) < actual number of factors (%d).\n",
			routine, ((int) numFact), that->fNf);
		goto done;
	}
	if (numDates < that->fNDates) {
		GtoErrMsg("%s: array not large enough: "
			"numDates (%d) < actual number of dates (%d).\n",
			routine, ((int) numDates), that->fNDates);
		goto done;
	}
#endif


	/*
	 * Convert to Gto vol array conventions
	 */
	if (meanRevers != NULL) {
	    for (idxF=0; idxF<=that->fNf-1; idxF++) {
		meanRevers[idxF] = that->fBeta[idxF];
	    }
	}

	if (backboneq != NULL)
		*backboneq = (1e0 - that->fBackBoneQ) * 0.5e0;

	/*
	 * Convert to Gto vol array conventions
	 */
	if (today != NULL) {
		*today = that->fDate[0];
	}
	numSpvol = that->fNDates;
	for (idxS=0; idxS<=numSpvol-1; idxS++) {
	    if (dates != NULL) {
	    	if ((idxS != numSpvol-1) && (numSpvol != 1)) {
		    dates[idxS] = that->fDate[idxS+1];
	    	} else {
		    dates[idxS] = that->fDate[idxS] + 3650L;
	    	}
	    }

	    if (spotVols != NULL) {
	    	for (idxF=0; idxF<=that->fNf-1; idxF++) {
		    spotVols[idxF][idxS] = that->fSigma[idxF][idxS] *
						that->fAlpha[idxF];
		}
	    }

	    if (correlations != NULL) {
	    	for (idxF=0; idxF<=that->fNf*(that->fNf-1)/2-1; idxF++) {
		    correlations[idxF][idxS] = that->fRho[idxF][idxS];
		}
	    }


	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: today %15s\n",\
		routine, (today != NULL ? GtoFormatDate(*today) : "N/A"));\
	for (idxS=0; idxS<=numSpvol-1; idxS++) {\
	    DrlFPrintf(vnfmFpLog, " [%2d] %12s", \
		idxS, (dates != NULL ? GtoFormatDate(dates[idxS]) : "N/A"));\
	    if (spotVols != NULL) {\
	    for (idxF=0; idxF<=that->fNf-1; idxF++) {\
		DrlFPrintf(vnfmFpLog, " %6.4f", spotVols[idxF][idxS]);\
	    }}\
	    if (correlations != NULL) {\
	    DrlFPrintf(vnfmFpLog, " |");\
	    for (idxF=0; idxF<=that->fNf*(that->fNf-1)/2-1; idxF++) {\
		DrlFPrintf(vnfmFpLog, " %6.4f", correlations[idxF][idxS]);\
	    }}\
	    DrlFPrintf(vnfmFpLog, "\n");\
	});
#endif



	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



