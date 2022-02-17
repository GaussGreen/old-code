/****************************************************************
 * Module:	DRL
 * Submodule:	SMAT
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>	
#include <float.h>	
#include <stdarg.h>	

#include "gtomat.h"
#include "date_sup.h"
#include "convert.h"
#include "tcurve.h"
#include "check.h"		/* checkFrequency */

#include "drlinter.h"		/* DrlLinInterp2d */
#include "drlstr.h"		/* DrlFloatPrint */
#include "drlio.h"		/* DrlFScanDouble */
#include "drltime.h"		/* DrlTDateIntervalToYears */
#include "drlproc.h"		/* DrlStrerror */


#include "drlsmat.h"

#undef	SWNEXP
#undef	SWNMAT
#undef	SWFREQ
#undef	SWTEXP
#undef	SWTMAT
#undef	SWVOL

#define	SWNEXP		(that->table->matrix->numDim1)
#define	SWNMAT		(that->table->matrix->numDim2)
#define	SWFREQ		(that->swapPayFreq)
#define	SWTEXP(idxExp)	(that->table->dim1Values[idxExp])
#define	SWTMAT(idxMat)	(that->table->dim2Values[idxMat])
#define	SWVOL(idxExp, idxMat) \
			(that->table->matrix->data[idxExp][idxMat])


#define	SMATNEXP(sw)		((sw)->table->matrix->numDim1)
#define	SMATNMAT(sw)		((sw)->table->matrix->numDim2)
#define	SMATFREQ(sw)		((sw)->swapPayFreq)
#define	SMATTEXP(sw, idxExp)	((sw)->table->dim1Values[idxExp])
#define	SMATTMAT(sw, idxMat)	((sw)->table->dim2Values[idxMat])
#define	SMATVOL(sw, idxExp, idxMat) \
				((sw)->table->matrix->data[idxExp][idxMat])

#define	SW2NEXP		(swMat->table->matrix->numDim1)
#define	SW2NMAT		(swMat->table->matrix->numDim2)
#define	SW2FREQ		(swMat->swapPayFreq)
#define	SW2TEXP(idxExp)	(swMat->table->dim1Values[idxExp])
#define	SW2TMAT(idxMat)	(swMat->table->dim2Values[idxMat])
#define	SW2VOL(idxExp, idxMat) \
			(swMat->table->matrix->data[idxExp][idxMat])


	/* check frequency */
/* until it is fixed */
#define	GtoErrMsgSet(x)	(!(x))


/*f--------------------------------------------------------------
 * Creates and returns a pointer to a new <i> TSwaptionMatrix2D</i>.
 * Returns NULL of failed.
 */

DLL_EXPORT(TSwaptionMatrix2D*)
DrlTSwaptionMatrix2DNew(
	TBoolean diagonal,	/* (I) TRUE=diagonal, FALSE=vertical */
	long freq,		/* (I) rate frequency (1,2,4,12) */
	int nExp,		/* (I) # of expirations */
	double *tExp,		/* (I) array of expirations (or NULL) */
	int nMat,		/* (I) # of maturities */
	double *tMat)		/* (I) array of maturities (or NULL) */
{
static	char		routine[] = "DrlTSwaptionMatrix2DNew";
	TSwaptionMatrix2D *that = NULL;
	int		status = FAILURE,
			i, j;

#undef	ASSERT
#define	ASSERT(cond)	{if (!(cond)) {GtoErrMsg("%s: `%s' failed\n",routine,\
			#cond); goto done;}}

	/* alloc memory */
	ASSERT((that = NEW(TSwaptionMatrix2D)) != NULL);

	/* check type */
	ASSERT((diagonal == 0) || (diagonal == 1));
	that->diagonal = diagonal;

	/* check frequency */
	switch (freq) {
	case 0:
	case 1:
	case 2:
	case 4:
	case 12:
		break;
	default:
		GtoErrMsg("%s: bad frequency %d.\n", routine, freq);
		goto done;
	}

	/*if (checkFrequency(routine, (long) freq) != TRUE) {
		GtoErrMsg("%s: bad frequency %d\n", routine, freq);
		goto done;
	}*/
	that->swapPayFreq = freq;

	/* malloc memory */
	ASSERT((that->table = NEW(TTable2D)) != NULL)
	ASSERT((that->table->dim1Values = NEW_ARRAY(double, nExp)) != NULL)
	ASSERT((that->table->dim2Values = NEW_ARRAY(double, nMat)) != NULL)

	ASSERT((that->table->matrix = GtoMatrixNewEmpty(nExp, nMat)) != NULL)

	/* copy */
	if (tExp != NULL) {
	    for (i=0; i<=nExp-1; i++) {
		that->table->dim1Values[i] = tExp[i];
	    }
	}
	if (tMat != NULL) {
	    for (j=0; j<=nMat-1; j++) {
		that->table->dim2Values[j] = tMat[j];
	    }
	}

	for (j=0; j<=nExp-1; j++)
	for (i=0; i<=nMat-1; i++)
		that->table->matrix->data[j][i] = 0e0;


	/*
	 *
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    if (that != NULL) FREE((char*) that);
	    GtoErrMsg("%s: failed (# mat %d x # exp %d)\n",
		routine, nMat, nExp);
	    return(NULL);
	}
	return(that);
#undef	ASSERT
}


/*f--------------------------------------------------------------
 * Deletes a pointer to a <i> TSwaptionMatrix2D</i>
 * crated by <i> DrlTSwaptionMatrix2DNew</i>.
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DFree(TSwaptionMatrix2D *that)
{
	if (that != NULL) {

	    GtoMatrixFree(that->table->matrix);

	    FREE((void*) that->table->dim1Values);
	    FREE((void*) that->table->dim2Values);
	    FREE((void*) that->table);
	    FREE((void*) that);
	}
	return(0);
}

/*f--------------------------------------------------------------
 * Checks that a <i> TSwaptionMatrix2D</i> contains valid
 * parameters.
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DCheck(TSwaptionMatrix2D *that)
{
static	char	routine[] = "DrlTSwaptionMatrix2DCheck";
	int	status = FAILURE;
	int	n;

	/* check frequency */
	if (checkFrequency(routine, (long) SWFREQ) != TRUE) {
		GtoErrMsg("%s: bad frequency %d\n", routine, SWFREQ);
		goto done;
	}

	if (SWTEXP(0) < 0e0) {
	    GtoErrMsg("%s: exp # %d (%lf) < 0\n",
		routine, 1, SWTEXP(0));
	    goto done;
	}
	for (n=1; n<=SWNEXP-1; n++) {
	    if (SWTEXP(n-1) >= SWTEXP(n)) {
		GtoErrMsg("%s: exp #%d (%lf) > exp #%d (%lf).\n",
		   routine, n, SWTEXP(n-1), n+1, SWTEXP(n));
		goto done;
	    }
	}

	if (SWTMAT(0) < 0e0) {
	    GtoErrMsg("%s: mat # %d (%lf) < 0\n",
		routine, 1, SWTMAT(0));
	    goto done;
	}
	for (n=1; n<=SWNMAT-1; n++) {
	    if (SWTMAT(n-1) >= SWTMAT(n)) {
		GtoErrMsg("%s: mat #%d (%lf) > mat #%d (%lf).\n",
		   routine, n, SWTMAT(n-1), n+1, SWTMAT(n));
		goto done;
	    }
	}


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f--------------------------------------------------------------
 * Creates and returns a pointer to a copy of the 
 * <i> TSwaptionMatrix2D</i> "fromMatrix".
 * Returns NULL if failed.
 */

DLL_EXPORT(TSwaptionMatrix2D*)
DrlTSwaptionMatrix2DNewCopy(TSwaptionMatrix2D *fromMatrix)
{
static	char	routine[] = "DrlTSwaptionMatrix2DNewCopy";
	int	status = FAILURE;
	TSwaptionMatrix2D *that = NULL;
	int	i,j;

	/* Check if ptr is null */
	if (fromMatrix == NULL) {
		GtoErrMsg("%s: null pointer received.\n", routine);
		return(NULL);
	}

	if ((that = DrlTSwaptionMatrix2DNew(
		fromMatrix->diagonal,
		fromMatrix->swapPayFreq,
		fromMatrix->table->matrix->numDim1,
	    	fromMatrix->table->dim1Values,
		fromMatrix->table->matrix->numDim2,
	    	fromMatrix->table->dim2Values)) == NULL)
			goto done;

	for (j=0; j<=that->table->matrix->numDim1-1; j++)
	for (i=0; i<=that->table->matrix->numDim2-1; i++) {
		that->table->matrix->data[j][i] =
			fromMatrix->table->matrix->data[j][i];
	}

	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlTSwaptionMatrix2DFree(that);
	    GtoErrMsg("%s: failed.\n", routine);
	    return(NULL);
	}
	return(that);
}


/*f--------------------------------------------------------------
 * Sums the values of a <i> TSwaptionMatrix2D</i> to "value".
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DSum(TSwaptionMatrix2D *that, double *value)
{
register int	i, j;
	*value = 0e0;
	for (i=0; i<=that->table->matrix->numDim2-1; i++)
	for (j=0; j<=that->table->matrix->numDim1-1; j++) {
		*value += that->table->matrix->data[j][i];
	}
	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * Performs a matrix operation on the swaption matrix <i> that</i>
 * (changed on exit).
 * <i> operation</i> can be one of the following:
 * $=$, $+$, $-$, $*$, $/$
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DOperScalar(
	TSwaptionMatrix2D *that,	/* (B) changed on exit */
	char *operation,		/* (I) operation see below */
	double value)			/* (I) input argument */
{
static	char	routine[] = "DrlTSwaptionMatrix2DOperScalar";
register int	idxM, idxE;

#undef	DOLOOP
#define	DOLOOP(statement)	{for (idxM=0; idxM<=SWNMAT-1; idxM++)\
				for (idxE=0; idxE<=SWNEXP-1; idxE++) {\
				statement;}}

	switch (operation[0]) {
	case '=':
		DOLOOP(SWVOL(idxE,idxM) = value);
		return(SUCCESS);
	case '+':
		DOLOOP(SWVOL(idxE,idxM) += value);
		return(SUCCESS);
	case '-':
		DOLOOP(SWVOL(idxE,idxM) -= value);
		return(SUCCESS);
	case '*':
		DOLOOP(SWVOL(idxE,idxM) *= value);
		return(SUCCESS);
	case '/':
		DOLOOP(SWVOL(idxE,idxM) /= value);
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad operation `%c'.\n", routine);
		return(FAILURE);
	}
#undef	DOLOOP
}


/*f--------------------------------------------------------------
 * Performs a matrix operation on the swaption matrix <i> that</i>
 * (changed on exit).
 * <i> operation</i> can be one of the following:
 * $=$, $+$, $-$, $*$, $/$
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DOperMatrix(
	TSwaptionMatrix2D *that,	/* (B) changed on exit */
	char *operation,		/* (I) operation see below */
	TSwaptionMatrix2D *swMat)	/* (I) matrix argument (unchanged) */
{
static	char	routine[] = "DrlTSwaptionMatrix2DOperMatrix";
register int	idxM, idxE;

#undef	DOLOOP
#define	DOLOOP(statement)	{for (idxM=0; idxM<=SWNMAT-1; idxM++)\
				for (idxE=0; idxE<=SWNEXP-1; idxE++) {\
				statement;}}

	switch (operation[0]) {
	case '=':
		DOLOOP(SWVOL(idxE,idxM) = SW2VOL(idxE,idxM))
		return(SUCCESS);
	case '+':
		DOLOOP(SWVOL(idxE,idxM) += SW2VOL(idxE,idxM))
		return(SUCCESS);
	case '-':
		DOLOOP(SWVOL(idxE,idxM) -= SW2VOL(idxE,idxM))
		return(SUCCESS);
	case '*':
		DOLOOP(SWVOL(idxE,idxM) *= SW2VOL(idxE,idxM))
		return(SUCCESS);
	case '/':
		DOLOOP(SWVOL(idxE,idxM) /= SW2VOL(idxE,idxM))
		return(SUCCESS);
	default:
		GtoErrMsg("%s: bad operation `%c'.\n", routine);
		return(FAILURE);
	}
#undef	DOLOOP
}

/*f--------------------------------------------------------------
 * Computes a numerical functional of the matrix <i> that</i>.
 * <i> operation</i> can be one of the following:
 * <br>
 * <br>[<i> "MIN"</i>] minimum value of all buckets,
 * <br>[<i> "MAX"</i>] maximum value of all buckets,
 * <br>[<i> "SUM"</i>] sum of all buckets,
 * <br>
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DNumVal(
	TSwaptionMatrix2D *that,	/* (I) input matrix */
	char *operation,		/* (I) operation see below */
	double *value)			/* (I) output value */
{
static	char	routine[] = "DrlTSwaptionMatrix2DOperScalar";
register int	idxM, idxE;

#undef	DOLOOP
#define	DOLOOP(statement)	{for (idxM=0; idxM<=SWNMAT-1; idxM++)\
				for (idxE=0; idxE<=SWNEXP-1; idxE++) {\
				statement;}}

	if (!strcmp(operation, "MIN")) {
		*value = 1e32;
		DOLOOP(*value = MIN(*value, SWVOL(idxE,idxM)));
		return(SUCCESS);
	} if (!strcmp(operation, "MAX")) {
		*value = -1e32;
		DOLOOP(*value = MAX(*value, SWVOL(idxE,idxM)));
		return(SUCCESS);
	} if (!strcmp(operation, "SUM")) {
		*value = 0e0;
		DOLOOP(*value += SWVOL(idxE,idxM));
		return(SUCCESS);
	} else {
		GtoErrMsg("%s: bad operation `%s'.\n", routine, operation);
		return(FAILURE);
	}
#undef	DOLOOP
}


/*f--------------------------------------------------------------
 * Checks that two matrices are similar.
 * Returns TRUE/FALSE iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DIsSameType(
	TSwaptionMatrix2D *that,	/* (I) matrix #1 */
	TSwaptionMatrix2D *swMat)	/* (I) matrix #2 */
{
static	char	routine[] = "DrlTSwaptionMatrix2DIsSameType";
register int	n;

	/* check frequency */
	if (that->diagonal != swMat->diagonal) {
	    GtoErrMsg("%s: diag matrix1 %d != matrix2 %d.\n",
		routine, that->diagonal, swMat->diagonal);
	    return(FALSE);
	}
	if (SWFREQ != SW2FREQ) {
	    GtoErrMsg("%s: freq matrix1 %d != matrix2 %d.\n",
		routine, SWFREQ ,SW2FREQ);
	    return(FALSE);
	}
	if (SWNEXP != SW2NEXP) {
	    GtoErrMsg("%s: # exp matrix1 %d != matrix2 %d.\n",
		routine, SWNEXP, SW2NEXP);
	    return(FALSE);
	}
	if (SWNMAT != SW2NMAT) {
	    GtoErrMsg("%s: # mat matrix1 %d != matrix2 %d.\n",
		routine, SWNMAT, SW2NMAT);
	    return(FALSE);
	}

	for (n=1; n<=SWNEXP-1; n++)  {
	    /*if (!IS_ALMOST_ZERO(SWTEXP(n) - SW2TEXP(n))) {*/
	    if (fabs(SWTEXP(n) - SW2TEXP(n)) >= 1.36e-3 /* 1/2 day */) {
		GtoErrMsg("%s: exp %d matrix1 %lf != matrix2 %lf.\n",
			routine, n, SWTEXP(n), SW2TEXP(n));
		return (FALSE);
	    }
	}

	for (n=1; n<=SWNMAT-1; n++) {
	    /*if (!IS_ALMOST_ZERO(SWTMAT(n) - SW2TMAT(n))) {*/
	    if (fabs(SWTMAT(n) - SW2TMAT(n)) >= 1.36e-3 /* 1/2 day */) {
		GtoErrMsg("%s: mat %d matrix1 %lf != matrix2 %lf.\n",
			routine, n, SWTMAT(n), SW2TMAT(n));
		return (FALSE);
	    }
	}
	return(TRUE);
}


/*f--------------------------------------------------------------
 * Reads a <i> TSwaptionMatrix2D</i> on a file "fnam" and, if successful,
 * creates a pointer to a new <i> TSwaptionMatrix2D</i> "that".
 * Calls <i> DrlTSwaptionMatrix2DFpRead</i>.
 * Returns 0 iff success.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DFileRead(
	TSwaptionMatrix2D **that,	/* (O) swaption matrix */
	char *fnam,			/* (I) file name */
	int format)			/* (I) format definition */
{
static	char	routine[] = "DrlTSwaptionMatrix2DFileRead";
	int	status = FAILURE;
	FILE	*fp = NULL;

	if ((fp = fopen(fnam, "r")) == NULL) {
	    GtoErrMsg("%s: can't open `%s' (%s)\n",
		routine, fnam, DrlStrerror());
	     goto done;
	} else {
	     if (DrlTSwaptionMatrix2DFpRead(that, fp, format) != SUCCESS)
		goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (fp != NULL) fclose(fp);
	if (status != SUCCESS) {
	    GtoErrMsg("%s: reading file `%s' failed.\n",
			routine, fnam);
	    *that = NULL;
	}
	return(status);
}

/*f--------------------------------------------------------------
 * Reads a <i> TSwaptionMatrix2D</i> on a file "fnam" and, if successful,
 * creates a pointer to a new <i> TSwaptionMatrix2D</i> "that".
 * Calls <i> DrlTSwaptionMatrix2DFpRead_Mod</i>.
 * Returns 0 iff success.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DFileRead_Mod(
	TSwaptionMatrix2D ***that,	/* (O) swaption matrix */
	char *fnam,			/* (I) file name */
	int format,
    int nbMatrix)			/* (I) format definition */
{
static	char	routine[] = "DrlTSwaptionMatrix2DFileRead_Mod";
	int	status = FAILURE;
	FILE	*fp = NULL;

	if ((fp = fopen(fnam, "r")) == NULL) {
	    GtoErrMsg("%s: can't open `%s' (%s)\n",
		routine, fnam, DrlStrerror());
	     goto done;
	} else {
	     if (DrlTSwaptionMatrix2DFpRead_Mod(that, fp, format, nbMatrix) != SUCCESS)
		goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (fp != NULL) fclose(fp);
	if (status != SUCCESS) {
	    GtoErrMsg("%s: reading file `%s' failed.\n",
			routine, fnam);
	    *that = NULL;
	}
	return(status);
}

/*f--------------------------------------------------------------
 * Writes a <i> TSwaptionMatrix2D</i> to a file "fnam".
 * Returns 0 iff successful.\\
 * <b> Remark:</b> automatically recognizes the London
 * swaption matrix format (DR Wrapper).
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DFileWrite(
	TSwaptionMatrix2D *that,	/* (I) */
	char *fnam,			/* (I) file name */
	int format)			/* (I) TSWAPTION_MATRIX_FMT_XXX */
{
static	char	routine[] = "DrlTSwaptionMatrix2DFileWrite";
	int	status = FAILURE;
	FILE	*fp = NULL;

	if ((fp = fopen(fnam, "w")) == NULL) {
	    GtoErrMsg("%s: can't open `%s' (%s)\n",
		routine, fnam, DrlStrerror());
	    goto done;
	} else {
	    if (DrlTSwaptionMatrix2DFpWrite(that, fp, format) != SUCCESS)
		goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (fp != NULL) fclose(fp);
	if (status != SUCCESS) {
	    GtoErrMsg("%s: reading file `%s' failed.\n", routine, fnam);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Reads a <i> TSwaptionMatrix2D</i> on a file pointer "fp"
 * and, if successful,
 * creates a pointer to a new <i> TSwaptionMatrix2D</i> "that".
 * Returns 0 iff successful. \\
 * <b> Remark:</b> recognizes the London swaption matrix format (DR Wrapper).
 */


DLL_EXPORT(int)
DrlTSwaptionMatrix2DFpRead(
	TSwaptionMatrix2D **that,	/* (O) */
	FILE *fp,			/* (I) file pointer */
	int format)			/* (I) see description */
{
static	char	routine[] = "DrlTSwaptionMatrix2DFpRead";
	char	buf[256];
	int	status = FAILURE,
		c, i, j, idxE, idxM, nMat, nExp, nPts,
		line=0;
	int	*iExp = NULL,
		*iMat = NULL;
	double	*vol = NULL;
	/*int		errMsgStatus;*/
	/*TDateInterval	dateInterval;*/
	double		val;

#define	SCAN_NEXT_LINE	{if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) {\
			GtoErrMsg("%s: %s\n", routine, DrlStrerror()); \
			goto done;}}

	/* get err msg status */
	/*errMsgStatus = GtoErrMsgSet(errMsgStatus = GtoErrMsgSet(FALSE));*/
	*that = NULL;

	/*
	 * Identify format
	 */
	switch (format) {
	case TSWAPTION_MATRIX_FMT_LON:
	    /*
	     * (1) London format
	     */
	    if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) goto done;
	    if (sscanf(buf, "%d", &nExp) != 1) goto done;
	    if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) goto done;
	    if (sscanf(buf, "%d", &nMat) != 1) goto done;

	    /* malloc memory */
	    if ((*that = DrlTSwaptionMatrix2DNew(
		    FALSE, 	/* vertical */
		    2,
		    nExp, NULL, nMat, NULL)) == NULL) goto done;

	    for (idxM=0; idxM<=SMATNMAT(*that)-1; idxM++) {
		if (DrlFScanDouble(fp, &SMATTMAT(*that, idxM)) != SUCCESS)
			goto done;
	    }

	    for (idxE=0; idxE<=SMATNEXP(*that)-1; idxE++) {
		if (DrlFScanDouble(fp, &SMATTEXP(*that, idxE)) != SUCCESS)
			goto done;
		SMATTEXP(*that, idxE) /= 12e0;

		for (idxM=0; idxM<=SMATNMAT(*that)-1; idxM++) {
		    if (DrlFScanDouble(fp, &SMATVOL(*that, idxE, idxM))
				!= SUCCESS) goto done;
		    SMATVOL(*that, idxE, idxM) *= 1e-2;
		}
	    }

    
	    if (DrlTSwaptionMatrix2DCheck(*that) != SUCCESS)
		    goto done;

	    break;

	case TSWAPTION_MATRIX_FMT_STD:
	case TSWAPTION_MATRIX_FMT_NUM:
	case TSWAPTION_MATRIX_FMT_CUR:
	    /*
	     * (2) CD format
	     */
	    /* get size of swaption matrix */
	    SCAN_NEXT_LINE;
	    if (sscanf(buf, "NUMBER_OF_MAT: %d", &nMat) != 1)
		goto done;

	    SCAN_NEXT_LINE;
	    if (sscanf(buf, "NUMBER_OF_EXP: %d", &nExp) != 1)
		goto done;


	    /* malloc memory */
	    if ((*that = DrlTSwaptionMatrix2DNew(
		FALSE, 	/* vertical */
		2,
		nExp, NULL, nMat, NULL)) == NULL) goto done;


	    /* read data */
	    SCAN_NEXT_LINE;
	    if (sscanf(buf, "FREQ: %d", &c) != 1)
		goto done;
	    (*that)->swapPayFreq = (long) c;

	    SCAN_NEXT_LINE;
	    if (sscanf(buf, "DIAG: %d", &c) != 1)
		goto done;
	    (*that)->diagonal = (c == 0 ? FALSE : TRUE);



	    line++;
	    for (i=0; i<=nMat-1;i++) {
		if (DrlFScanString(fp, buf) != SUCCESS) goto done;
		if (DrlDoubleIntervalScan(buf, &(*that)->table->dim2Values[i])
			!= SUCCESS) goto done;
	    }

	    for (j=0; j<=nExp-1; j++) {
		line++;
		if (DrlFScanString(fp, buf) != SUCCESS) goto done;
		if (DrlDoubleIntervalScan(buf, &(*that)->table->dim1Values[j])
			!= SUCCESS) goto done;
	        for (i=0; i<=nMat-1;i++) {
		    if (DrlFScanString(fp, buf) != SUCCESS) goto done;

		    /*c = sscanf(buf, "%lf", &val);*/
		    c = (DrlCurScan(buf, &val) == SUCCESS ? 1 : 0);

		    if ((*that)->diagonal == TRUE)
		    {
		        if (c == 1)
		        {
			    (*that)->table->matrix->data[j][i] = val;
		        } else {
			    (*that)->table->matrix->data[j][i] = 0e0;
		        }
		    } else {
		        if (c != 1) goto done;
		        (*that)->table->matrix->data[j][i] = val;
		    }

		    switch (format) {
		    case TSWAPTION_MATRIX_FMT_STD:
			(*that)->table->matrix->data[j][i] *= 1e-2;
			break;
		    case TSWAPTION_MATRIX_FMT_NUM:
		    case TSWAPTION_MATRIX_FMT_CUR:
			break;
		    }

	        }
	    }

	    if (DrlTSwaptionMatrix2DCheck(*that) != SUCCESS)
		goto done;

	    break;

	default:
	    goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	/*GtoErrMsgSet(errMsgStatus);*/
	if (iExp != NULL) FREE((void*) iExp);
	if (iMat != NULL) FREE((void*) iMat);
	if (vol != NULL)  FREE((void*) vol);

	if (status != SUCCESS) {
	    if (*that != NULL) {
		DrlTSwaptionMatrix2DFree(*that);
	    }
	    GtoErrMsg("%s: failed (line %d)\n", routine, line);
	    *that = NULL;
	}
	return(status);
#undef	SCAN_NEXT_LINE
}


/*f--------------------------------------------------------------
 * Reads a <i> TSwaptionMatrix2D</i> on a file pointer "fp"
 * and, if successful,
 * creates a pointer to a new <i> TSwaptionMatrix2D</i> "that".
 * Returns 0 iff successful. \\
 * <b> Remark:</b> recognizes the London swaption matrix format (DR Wrapper).
 */


DLL_EXPORT(int)
DrlTSwaptionMatrix2DFpRead_Mod(
    TSwaptionMatrix2D ***that,  /* (O) */
    FILE *fp,                   /* (I) file pointer */
    int format,                 /* (I) see description */
    int nbMatrix)               /* (I) nb of matrices  */
{
    static  char    routine[] = "DrlTSwaptionMatrix2DFpRead_Mod";
    char    buf[256];
    int     status = FAILURE,
            i, idxE, idxM, nMat, nExp,
            line=0;
    int     *iExp = NULL,
            *iMat = NULL;
    double  *vol = NULL;
    TSwaptionMatrix2D *CurrMatrix = NULL;

#define	SCAN_NEXT_LINE	{if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) {\
            GtoErrMsg("%s: %s\n", routine, DrlStrerror()); \
            goto done;}}*/
	/* get err msg status */
	/*errMsgStatus = GtoErrMsgSet(errMsgStatus = GtoErrMsgSet(FALSE));*/

    *that = NULL;

    /*
     * Identify format
     */
    switch (format) 
    {
    case TSWAPTION_MATRIX_FMT_LON:
        /*
         * (1) London format
         */
        if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) goto done;
        if (sscanf(buf, "%d", &nExp) != 1) goto done;
        if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL) goto done;
        if (sscanf(buf, "%d", &nMat) != 1) goto done;

        /* malloc memory */
        if ( (*that = NEW_ARRAY(
                      TSwaptionMatrix2D*,
                      nbMatrix)) == NULL) goto done;

        for (i = 0; i < nbMatrix; i++)
        {
            if ( (CurrMatrix = DrlTSwaptionMatrix2DNew(
                               FALSE,
                               2,
                               nExp,
                               NULL,
                               nMat,
                               NULL)) == NULL) goto done;
            (*that)[i] = CurrMatrix;
            CurrMatrix = NULL;
        }


        for (i = 0; i < nbMatrix; i++)
        {
            CurrMatrix = *((*that) + i);
            for (idxM=0; idxM<=SMATNMAT(CurrMatrix)-1; idxM++) 
            {
                if (DrlFScanDouble(fp, &SMATTMAT(CurrMatrix, idxM)) != SUCCESS)
                goto done;
            }

            for (idxE=0; idxE<=SMATNEXP(CurrMatrix)-1; idxE++) 
            {
                if (DrlFScanDouble(fp, &SMATTEXP(CurrMatrix, idxE)) != SUCCESS)
                    goto done;
                SMATTEXP(CurrMatrix, idxE) /= 12e0;


                for (idxM=0; idxM<=SMATNMAT(CurrMatrix)-1; idxM++) 
                {
                    if (DrlFScanDouble(fp, &SMATVOL(CurrMatrix, idxE, idxM))
                        != SUCCESS) goto done;
                }
            }
            if (DrlTSwaptionMatrix2DCheck(CurrMatrix) != SUCCESS)
                goto done;
        }
        break;

    case TSWAPTION_MATRIX_FMT_STD:
    case TSWAPTION_MATRIX_FMT_NUM:
    case TSWAPTION_MATRIX_FMT_CUR:
        goto done;

    default:
        goto done;
    }

    /* made it through OK */
    status = SUCCESS;
done:
    /*GtoErrMsgSet(errMsgStatus);*/
    if (iExp != NULL) FREE((void*) iExp);
    if (iMat != NULL) FREE((void*) iMat);
    if (vol != NULL)  FREE((void*) vol);

    if (status != SUCCESS) 
    {
        for (i = 0; i < nbMatrix; i++)
        {
            CurrMatrix = *((*that)+i);
            if (CurrMatrix != NULL) 
            {
                DrlTSwaptionMatrix2DFree(CurrMatrix);
            }
        }
        GtoErrMsg("%s: failed (line %d)\n", routine, line);
        if ( *that != NULL) FREE((void*) *that);
    }
    return(status);
#undef	SCAN_NEXT_LINE
}



/*
 * Convenoence routine
 */

static	char*
printInterval(double dblValue)
{
static	char		buf[] = "ERR";
	TDateInterval	interval;
	if (GtoYearsToDateInterval(dblValue, &interval) != SUCCESS) {
		return buf;
	} else {
		return GtoFormatDateInterval(&interval);
	}
}






/*f--------------------------------------------------------------
 * Writes a <i> TSwaptionMatrix2D</i> to a file pointer "fp".
 * The fommowing options are available for the "format"
 * argument:
 * <br>
 * <br>[TSWAPTION\_MATRIX\_FMT\_STD] standard format,
 * <br>[TSWAPTION\_MATRIX\_FMT\_LON] London swaption matrix format
 * <br>[TSWAPTION\_MATRIX\_FMT\_TXT] test file format for import
 * in spreadsheet: intervals as numbers, text enclosed in quotes.
 * <br>
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DFpWrite(
	TSwaptionMatrix2D *that,	/* (I) */
	FILE *fp,			/* (I) file pointer */
	int format,			/* (I) TSWAPTION_MATRIX_FMT_XXX */
	...)
{
static	char		routine[] = "DrlTSwaptionMatrix2DFpWrite";
	int		i,j;
	int		idxE, idxM;
	TDateInterval	interval;
	char		*tagName;
	va_list		ap;

	switch (format) {
	case TSWAPTION_MATRIX_FMT_STD:
	case TSWAPTION_MATRIX_FMT_NUM:
	case TSWAPTION_MATRIX_FMT_CURM:
	case TSWAPTION_MATRIX_FMT_CURK:
	    /*
	     * Standard format
	     */
	    DrlFPrintf(fp, "NUMBER_OF_MAT: %d\n", that->table->matrix->numDim2);
	    DrlFPrintf(fp, "NUMBER_OF_EXP: %d\n", that->table->matrix->numDim1);
	    DrlFPrintf(fp, "FREQ: %d\n", (int) that->swapPayFreq);
	    DrlFPrintf(fp, "DIAG: %d\n", (that->diagonal == TRUE));

	    DrlFPrintf(fp, "        ");
	    for (i=0; i<=that->table->matrix->numDim2-1; i++) {
		if (GtoYearsToDateInterval(
		    that->table->dim2Values[i], &interval) != SUCCESS) {
			DrlFPrintf(fp, "\t%4s", "ERR");
	    	} else {
		    DrlFPrintf(fp, "\t%4s", GtoFormatDateInterval(&interval));
	        }
	    }
	    DrlFPrintf(fp, "\n");

	    for (j=0; j<=that->table->matrix->numDim1-1;j++) {
	        if (GtoYearsToDateInterval(
	    	    that->table->dim1Values[j], &interval) != SUCCESS) {
	    	    DrlFPrintf(fp, "%8s", "ERR");
	        } else {
	    	    DrlFPrintf(fp, "%8s", GtoFormatDateInterval(&interval));
	        }

	        for (i=0; i<=that->table->matrix->numDim2-1; i++) {
		    switch (format) {
		    case TSWAPTION_MATRIX_FMT_STD:
			DrlFPrintf(fp, "\t%-.2f",
				that->table->matrix->data[j][i] * 1e2) ;
			break;
		    case TSWAPTION_MATRIX_FMT_NUM:
			DrlFPrintf(fp, " %s", DrlFloatPrint(NULL, SWVOL(j, i), 8));
			break;
		    case TSWAPTION_MATRIX_FMT_CUR:
			DrlFPrintf(fp, " %15s", DrlCurPrint(NULL, SWVOL(j, i), 0));
			break;
		    case TSWAPTION_MATRIX_FMT_CURK:
			DrlFPrintf(fp, " %11s",
				DrlCurPrint(NULL, SWVOL(j, i)*1e-3, 0));
			break;
		    case TSWAPTION_MATRIX_FMT_CURM:
			DrlFPrintf(fp, " %7s",
				DrlCurPrint(NULL, SWVOL(j, i)*1e-6, 0));
			break;
		    }
	        }
	        DrlFPrintf(fp, "\n");
	    }
	    break;

	case TSWAPTION_MATRIX_FMT_CUR:
	    /*
	     * Standard format
	     */
	    DrlFPrintf(fp, "NUMBER_OF_MAT: %d\n", SWNMAT);
	    DrlFPrintf(fp, "NUMBER_OF_EXP: %d\n", SWNEXP);
	    DrlFPrintf(fp, "FREQ: %d\n", (int) SWFREQ);
	    DrlFPrintf(fp, "DIAG: %d\n", (that->diagonal == TRUE));

	    DrlFPrintf(fp, "        ");
	    for (i=0; i<=SWNMAT-1; i++) {
		DrlFPrintf(fp, " %15s", printInterval(SWTMAT(i)));
	    }
	    DrlFPrintf(fp, "\n");

	    for (j=0; j<=SWNEXP-1;j++) {
	    	DrlFPrintf(fp, "%8s", printInterval(SWTEXP(j)));
	        for (i=0; i<=SWNMAT-1; i++) {
			DrlFPrintf(fp, " %15s", DrlCurPrint(NULL, SWVOL(j, i), 0));
	        }
	        DrlFPrintf(fp, "\n");
	    }
	    break;

	case TSWAPTION_MATRIX_FMT_TXT:
	case TSWAPTION_MATRIX_FMT_PRN:
	    /*
	     * Format compatible with import is in spreadsheets
	     */
	    DrlFPrintf(fp, "\"NUMBER_OF_MAT\"\t%d\n",
		that->table->matrix->numDim2);
	    DrlFPrintf(fp, "\"NUMBER_OF_EXP\"\t%d\n",
		that->table->matrix->numDim1);
	    DrlFPrintf(fp, "\"FREQ\"\t%d\n", (int) that->swapPayFreq);
	    DrlFPrintf(fp, "\"DIAG\"\t%d\n", (that->diagonal == TRUE));

	    DrlFPrintf(fp, "\"VOLS\"");
	    for (i=0; i<=that->table->matrix->numDim2-1; i++) {
		DrlFPrintf(fp, "%8.4f", that->table->dim2Values[i]);
	    }
	    DrlFPrintf(fp, "\n");

	    for (j=0; j<=that->table->matrix->numDim1-1;j++) {
		DrlFPrintf(fp, "\t%8.4f", that->table->dim1Values[j]);
	        for (i=0; i<=that->table->matrix->numDim2-1; i++) {
		    DrlFPrintf(fp, "\t%12.8f",
			that->table->matrix->data[j][i]) ;
	        }
	        DrlFPrintf(fp, "\n");
	    }
	    break;
	case TSWAPTION_MATRIX_FMT_LON:
	    /*
	     * Export for Dr Wrapper Type 2
	     */
            DrlFPrintf(fp, "# Number of expiries\n %d\n", SWNEXP);
            DrlFPrintf(fp, "# No of forward maturities\n %d \n", SWNMAT);
            DrlFPrintf(fp, "# Forward maturities (years), "
                "then expiries (months) and volatilities in %%\n");

	    DrlFPrintf(fp, "        ");
	    for (i=0; i<=SWNMAT-1; i++) {
		DrlFPrintf(fp, "\t%4.0f", SWTMAT(i));
	    }
	    DrlFPrintf(fp, "\n");

	    for (j=0; j<=SWNEXP-1;j++) {
	    	DrlFPrintf(fp, "  %4.0f  ", SWTEXP(j)*12e0);

	        for (i=0; i<=SWNMAT-1; i++) {
			DrlFPrintf(fp, "\t%-.2f", SWVOL(j, i) * 1e2);
		}
	        DrlFPrintf(fp, "\n");
	    }
	    break;

	case TSWAPTION_MATRIX_FMT_LIST:
	    /*
	     * List format '<tagName> <exp> <mat> <value>'
	     */
	    va_start(ap, format);
	    tagName = (char*) va_arg(ap, char*);
	    va_end(ap);

	    for (idxM=0; idxM<=SMATNMAT(that)-1; idxM++)
	    for (idxE=0; idxE<=SMATNEXP(that)-1; idxE++) {
		DrlFPrintf(fp, " %15s", tagName);
	    	DrlFPrintf(fp, " %8.4f", SMATTEXP(that, idxE));
	    	DrlFPrintf(fp, " %4s", printInterval(SMATTMAT(that, idxM)));
		DrlFPrintf(fp, " %15s", DrlCurPrint(NULL,
				SMATVOL(that, idxE, idxM), 0));
		DrlFPrintf(fp, "\n");
	    }

	    break;

	default:
	    GtoErrMsg("%s: unknowm format %d.\n", routine, format);
	    return(FAILURE);
	}


	return(SUCCESS);
}

/*f--------------------------------------------------------------
 * Reads data <i> from</i> the wrapper interface and create a new C structure.
 * <br>
 * <br>[typeL]	LIL array of length 2.\\
 *              typeL[0]: type of the matrix (0=diagonal, 1=vertical),\\
 * 		typeL[1]: frequency of rate.
 * <br>[tMatL] LIL array of float containing the swaption maturities,
 * <br>[tExpL] LIL array of float containing the swaption expiration,
 * <br>[volL]  LIL array of float (range) containing the swaption volatilities.
 * <br>
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DWrapRead(
	TSwaptionMatrix2D **that,	/* (O) new swaption data structure */
	long *typeL,			/* (I) 'L' [0]=diag, [1]=freq */
	double *tMatL,			/* (I) 'F' */
	double *tExpL,			/* (I) 'F' */
	double *volL)			/* (I) 'F' volatilities (can be NULL) */
{
static	char	routine[] = "DrlTSwaptionMatrix2DWrapRead";
	int	i, j,
		nMat, nExp,
		status = FAILURE;

	WRAP_CHECK_VECTOR_LEN(typeL, 2);
	WRAP_CHECK_VECTOR(tMatL);
	WRAP_CHECK_VECTOR(tExpL);

	nMat = (int) tMatL[0];
	nExp = (int) tExpL[0];

	if (volL != NULL) {
		WRAP_CHECK_VECTOR(volL);
		if (nExp * nMat != ((int) volL[0])) {
		    GtoErrMsg("%s: bad volatility  range size "
			"(#exp. (%d) x #mat (%d) != #vols %d)\n",
			routine, nMat, nExp, (int) volL[0]);
	    	    goto done;
		}
	}


	*that = DrlTSwaptionMatrix2DNew(
		(int) typeL[1],
		(int) typeL[2],
		(int) tExpL[0], tExpL+1,
		(int) tMatL[0], tMatL+1);
	if (*that == NULL) goto done;

	/* put vols */
	if (volL != NULL) {
	    for (i=0; i<=nMat-1; i++)
	    for (j=0; j<=nExp-1; j++) {
		(*that)->table->matrix->data[j][i] =
			volL[i + j*nMat + 1];
	    }
	} else {
	    for (i=0; i<=nMat-1; i++)
	    for (j=0; j<=nExp-1; j++) {
		(*that)->table->matrix->data[j][i] = 0e0;
	    }
	}

#if defined(__DEBUG__)
	DrlFPrintf(stdout, "%s:\n", routine);
	DrlTSwaptionMatrix2DFpWrite(*that, stdout, TSWAPTION_MATRIX_FMT_STD);
#endif

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed \n", routine);
	    *that = NULL;
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Writes a C data structure to the wrapper interface.
 * For a description of the arguments, see <i> DrlTSwaptionMatrix2DWrapRead</i>.
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DWrapWrite(
	TSwaptionMatrix2D *that,	/* (I) swaption matrix */
	long *typeL,			/* (O) 'L' (can be NULL) */
	double *tMatL,			/* (O) 'F' (can be NULL) */
	double *tExpL,			/* (O) 'F' (can be NULL) */
	double *volL)			/* (O) 'F' */
{
static	char	routine[] = "DrlTSwaptionMatrix2DWrapWrite";
	int	status = FAILURE;
	int	i, j,
		nMat, nExp;


	if (typeL != NULL) WRAP_CHECK_VECTOR_LEN(typeL, 2);
	if (tMatL != NULL) WRAP_CHECK_VECTOR(tMatL);
	if (tExpL != NULL) WRAP_CHECK_VECTOR(tExpL);

	/* check matrix size consistent */
	if ((tMatL != NULL) && (tExpL != NULL)) {
		nMat = (int) tMatL[0];
		nExp = (int) tExpL[0];
	} else {
		nMat = that->table->matrix->numDim2;
		nExp = that->table->matrix->numDim1;
	}

	if (nExp*nMat != ((int) volL[0])) {
		GtoErrMsg("%s: bad array length (%dx%d != %d)\n",
			routine, nExp, nMat, (int) volL[0]);
		goto done;
	}

	/* put mat and exp times */
	if (tMatL != NULL) {
	    for (i=0; i<=MIN(that->table->matrix->numDim2, nMat)-1; i++) 
		tMatL[i+1] = that->table->dim2Values[i];
	}
	if (tExpL != NULL) {
	    for (j=0; j<=MIN(that->table->matrix->numDim1, nExp)-1; j++) 
		tExpL[j+1] = that->table->dim1Values[i];
	}

	/* put vols */
	for (j=0; j<=MIN(that->table->matrix->numDim1, nExp)-1; j++) 
	for (i=0; i<=MIN(that->table->matrix->numDim2, nMat)-1; i++) 
	{
	    volL[i + j*nMat + 1] = that->table->matrix->data[j][i];
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Reads a <i> TSwaptionMatrix2D</i> on a file pointer <i> fp</i>
 * assuming the base volatility matrix format (possibly IMM dates)
 * of the form
 * \begin{verbatim}
 * "Maturity"   "Exp\Und"  0.08330 0.25000 0.50000 ...
 * "SN"      "20-May-1998" 0.43000 0.45000 0.41000 ...
 * "IMM1"    "17-Jun-1998" 0.45000 0.42000 0.42000 ...
 * "IMM2"    "16-Sep-1998" 0.42000 0.43000 0.46000 ...
 * "IMM3"    "16-Dec-1998" 0.43000 0.44000 0.48000 ...
 * ...        ...          ...     ...     ...     ...
 * \end{verbatim}
 * The date <i> baseDate</i> is used to convert intervals 
 * to years (assuming a 30/360 day count).
 * If <i> baseDate</i> is 0, the reference date is assumed
 * to be the first date (SN).
 * Returns 0 iff successful. \\
 */


DLL_EXPORT(int)
DrlTSwaptionMatrix2DFpReadBVFmt(
	TSwaptionMatrix2D **that,	/* (O) new matrix */
	FILE *fp,			/* (I) file pointer */
	TDate baseDate,			/* (I) reference date */
	int freq)			/* (I) frequency */
{
static	char	routine[] = "DrlTSwaptionMatrix2DFpReadBVFmt";
	int	status = FAILURE;
	char	buf[1024];
	int	line=0;
	char	s0[1024], *s1, *p;
static	char	tok[] = " \t";
	TDateInterval	interval;
	TDate	expDate;
	int	nExp, nMat, idxE, idxM;
#define		nExpMax	128
#define		nMatMax	128
	double	tExp[nExpMax], tMat[nMatMax], values[nExpMax][nMatMax];

        TDayCount	dayCount = GTO_B30_360; /* dcc for year fraction */


	/* read first line */
	if (DrlFGetLine(buf, sizeof(buf), fp, &line) == NULL)
		goto done;
	if ((p = DrlStrToken(buf, tok, s0, &s1)) == NULL)
		goto done;
	if ((p = DrlStrToken(NULL, tok, s0, &s1)) == NULL)
		goto done;
	nMat = 0;
	while ((p = DrlStrToken(NULL, tok, s0, &s1)) != NULL) {
		if (sscanf(p, "%lf", &tMat[nMat]) != 1)
			goto done;
		nMat++;
		if (nMat >= nMatMax) {
			GtoErrMsg("%s: too many maturities.\n", routine);
			goto done;
		}
	}
	if (nMat <= 0) goto done;


	/* read lines */
	nExp = 0;
	while (DrlFGetLine(buf, sizeof(buf), fp, &line) != NULL) {
		/* read interval */
		if ((p = DrlStrToken(buf, tok, s0, &s1)) == NULL)
			goto done;

#ifdef	_SKIP
		if (DrlTDateIntervalScan(p, &interval) != SUCCESS)
			goto done;

		if (DrlTDateIntervalToYears(
			baseDate,
 			interval,
			dayCount,
			&tExp[nExp]) != SUCCESS)
				goto done;
#endif

		/* read date */
		if ((p = DrlStrToken(NULL, tok, s0, &s1)) == NULL)
			goto done;
		if (DrlTDateScan(p, &expDate) != SUCCESS)
			goto done;


		if ((nExp == 0) && (baseDate <= 0L))
			baseDate = expDate;


		if (GtoDayCountFraction(baseDate, expDate, dayCount,
			&tExp[nExp]) != SUCCESS)
				goto done;


		/* read values */
		for (idxM=0; idxM<=nMat-1; idxM++) {
			if ((p = DrlStrToken(NULL, tok, s0, &s1)) == NULL)
				goto done;
			if (sscanf(p, "%lf", &values[nExp][idxM]) != 1)
				goto done;
		}

		nExp++;
		if (nExp >= nExpMax) {
			GtoErrMsg("%s: too many maturities.\n", routine);
			goto done;
		}
	}


	/* Create Matrix */

	*that = DrlTSwaptionMatrix2DNew(
		FALSE,		/* TRUE=diagonal, FALSE=vertical */
		(long) freq,
		nExp,
		tExp,
		nMat,
		tMat);
	if (*that == NULL)
		goto done;

	for (idxM=0; idxM<=SMATNMAT(*that)-1; idxM++)
	for (idxE=0; idxE<=SMATNEXP(*that)-1; idxE++) {
		SMATVOL(*that, idxE, idxM) = values[idxE][idxM];
	}



	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed (line %d).\n", routine, line);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Allocates an array "dates" of length "nDates" corresponding
 * to the expirations of the swaption matrix "that".
 * The dates are offset from "todayDate".
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlTSwaptionMatrix2DCreateTimeLine(
	TSwaptionMatrix2D *that,	/* (I) swaption matrix */
	TDate todayDate,		/* (I) reference date */
	int *nDates,			/* (O) # dates */
	TDate **dates)			/* (O) array of dates */
{
static	char	routine[] = "DrlTSwaptionMatrix2DCreateTimeLine";
	int	status = FAILURE;
	int	j;
	TDateInterval	interval;

	*nDates = 0;
	*dates = NULL;
	if ((*dates = NEW_ARRAY(TDate, that->table->matrix->numDim1))
	    == NULL) goto done;

	for (j=0; j<=that->table->matrix->numDim1-1;j++) {
	    if (GtoYearsToDateInterval(that->table->dim1Values[j],
		&interval) != SUCCESS) goto done;
	    if (GtoDtFwdAny(todayDate, &interval, (*dates) +j)
		!= SUCCESS) goto done;
	}
	*nDates = that->table->matrix->numDim1;

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    if (*dates != NULL) FREE(*dates);
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

