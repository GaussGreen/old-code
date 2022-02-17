/************************************************************************
 * Module:	DRL
 * Submodule:	VTYPE
 * File:	
 * Function:	Variable Type
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>		/* IEEE finite() */
#include <float.h>		/* DBL_MIN, etc. */

#include "ldate.h"
#include "date_sup.h"
#include "convert.h"
#include "yearfrac.h"

#include "drltime.h"		/* DrlTDateScan() */
#include "drlstr.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlvtype.h"		/* DrlVType...() */
#include "drlproc.h"		/* DrlStrerror() */

/* for WINDLL*/
#if defined(CLIB) && defined(_WINDLL) && !(defined(WIN32) || defined(_WIN32))
# include "cfileio.h"
# define	sscanf		GtoSscanf
# define	sprintf		sprintf
#endif

#include "drlptable.h"		/* Prototype Consistency */


#define	MAXSTR		256		/* maximum string length */

/*#define	__DEBUG__*/

/*f-------------------------------------------------------------
 * Creates a new (empty) table with <i> numFields</i> fields
 * with types given the array {\t fieldType} (of length <i> numFields</i>.
 * Returns NULL if failed.
 */

TVTable*
DrlTVTableNew(
	int numFields,		/* (I) number of fields */
	TVType *fieldType)	/* (I) field type [numFields] */
{
static	char	routine[] = "DrlTVTableNew";
	int	status= FAILURE;

	int	idxF;
	TVTable	*that = NULL;

	/*
	 *
	 */
	if ((that = NEW(TVTable)) == NULL) goto done;

	that->numFields =  numFields;
	that->numData = 0;
	that->numDataMax = 0;

	if ((that->fieldType = NEW_ARRAY(TVType, numFields)) == NULL)
		goto done;
	if ((that->fieldName = NEW_ARRAY(char*, numFields)) == NULL)
		goto done;
	if ((that->dataValue = NEW_ARRAY(void*, numFields)) == NULL)
		goto done;
	for (idxF=0; idxF<=numFields-1; idxF++) {
		that->fieldType[idxF] = fieldType[idxF];
		if ((that->fieldName[idxF] = NEW_ARRAY(char, MAXSTR)) == NULL)
			goto done;
		that->dataValue[idxF] = NULL;
	}

	if ((that->numFieldValues = NEW_ARRAY(int, numFields)) == NULL)
		goto done;
	if ((that->numFieldValuesMax = NEW_ARRAY(int, numFields)) == NULL)
		goto done;
	if ((that->fieldValues = NEW_ARRAY(void*, numFields)) == NULL)
		goto done;
	for (idxF=0; idxF<=numFields-1; idxF++) {
		that->numFieldValues[idxF] = 0;
		that->numFieldValuesMax[idxF] = 0;
		that->fieldValues[idxF] = NULL;
	}

	if ((that->tmpData = NEW_ARRAY(void*, numFields)) == NULL)
		goto done;
	for (idxF=0; idxF<=numFields-1; idxF++) {
		if ((that->tmpData[idxF] = DrlVTypeVectAlloc(1,
			that->fieldType[idxF])) == NULL) goto done;
	}

#ifdef	__DEBUG__
	GtoErrMsg("%s:\n", routine);
	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		GtoErrMsg("\t%1d/%1d\t%s\n", idxF+1, that->numFields,
			DrlVTypeName(that->fieldType[idxF]));
	}
#endif


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
		return(NULL);
	}
	return(that);
}


/*f-------------------------------------------------------------
 * Frees a table.
 */

int
DrlTVTableFree(TVTable* that)
{
	int	idxF;

	if (that == NULL) return(SUCCESS);

	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		DrlVTypeVectFree(
			that->dataValue[idxF],
			that->numDataMax,
			that->fieldType[idxF]);
	}


	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		DrlVTypeVectFree(
			that->fieldValues[idxF],
			that->numFieldValuesMax[idxF],
			that->fieldType[idxF]);
	}

	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		DrlVTypeVectFree(
			that->tmpData[idxF],
			1,
			that->fieldType[idxF]);

		FREE(that->fieldName[idxF]);
	}

	FREE(that->tmpData);

	FREE(that->numFieldValues);
	FREE(that->numFieldValuesMax);
	FREE(that->fieldValues);
	FREE(that->fieldType);
	FREE(that->fieldName);
	FREE(that->dataValue);

	FREE(that);

	return (SUCCESS);
}



/*f-------------------------------------------------------------
 * Reads data from a file <i> fp</i>
 * and creates a new table with <i> numFields</i> fields
 * with types given the array {\t fieldType} (of length <i> numFields</i>.
 * The format for text file is (N stands for <i> numField</i>):
 * \begin{verbatim}
 * <field name 1>  ..... <field name N>
 * <data 1 field 1>  ..... <data 1 field N>
 * ...
 * <data K field 1>  ..... <data K field N>
 * ...
 * \end{verbatim}
 * The file scanning stops when end EOF is encountered.
 * Returns NULL if failed. 
 */

TVTable*
DrlTVTableNewFpRead(
	FILE *fp,		/* (I) file pointer */
	int numFields,		/* (I) number of fields (or -1) */
	TVType *fieldTypes)	/* (I) field type [numFields] */
{
static	char	routine[] = "DrlTVTableNewFpRead";
	int	status= FAILURE;

	TVTable* that = NULL;
	TVType	fieldTypes1[32];
	int	line = 1;
	char	bufline[2048], buf[2048], *p;

	int	idxF, n, nMax;

	/* If field info is not provided, */
	if (numFields <= 0) {
		fieldTypes = fieldTypes1;

		/* Read size  & types */
		if (DrlFScanInt(fp, &numFields) != SUCCESS) {
			goto done;
		}
		for (idxF=0; idxF<=numFields-1; idxF++) {
			if (DrlFScanString(fp, buf) != SUCCESS)
				goto done;
			if (DrlVTypeNameScan(&fieldTypes[idxF], buf) != SUCCESS)
				goto done;
		}
	}

	/* Create table */
	if ((that = DrlTVTableNew(numFields, fieldTypes)) == NULL)
			goto done;


	/* Read field names */
	if (DrlFGetLine(bufline, sizeof(bufline), fp, &line) == NULL) {
		GtoErrMsg("%s: can't read line containing field names.\n",
			routine);
		goto done;
	}

	p = bufline;
	for (idxF=0; idxF<=numFields-1; idxF++) {
		if (DrlStrScanString(&p, buf) == NULL) {
			GtoErrMsg("%s: can't read field name %d/%d.\n",
				routine, idxF+1, numFields);
			goto done;
		}
		strncpy(that->fieldName[idxF], buf, MAXSTR);
	}


	/* Read data */
	while (DrlFGetLine(bufline, sizeof(bufline), fp, &line) != NULL) {
		/* Read fields of data point */
		p = bufline;
		for (idxF=0; idxF<=numFields-1; idxF++) {
			if (DrlStrScanString(&p, buf) == NULL) {
				GtoErrMsg("%s: element %d can't read field "
					"name %d/%d", routine, that->numData,
					idxF+1, numFields);
				goto done;
			}
			if (DrlVTypeScan(
				buf,
				that->fieldType[idxF],
				that->tmpData[idxF]) != SUCCESS)
					goto done;
#ifdef	__DEBUG__
			GtoErrMsg("%s: scanned field #%d type %-10s: "
				"`%s' -> `%s'\n", routine,
				idxF, DrlVTypeName(that->fieldType[idxF]),
				buf,
				DrlVTypePrint(NULL, that->fieldType[idxF],
					that->tmpData[idxF]));
#endif
		}

		/* Add element to data set */
		for (idxF=0; idxF<=numFields-1; idxF++) {
			n    = that->numData;
			nMax = that->numDataMax;
			if (DrlVTypeVectAdd(
				&that->dataValue[idxF],
				&n,
				&nMax,
				1024,
				that->tmpData[idxF],
				FALSE,	/* do not remove double items */
				that->fieldType[idxF]) != SUCCESS)
					goto done;
		}
		that->numDataMax = nMax;
		that->numData += 1;

		/* Add field values */
		for (idxF=0; idxF<=numFields-1; idxF++) {
			if (DrlVTypeVectAdd(
				&that->fieldValues[idxF],
				&that->numFieldValues[idxF],
				&that->numFieldValuesMax[idxF],
				1024,
				that->tmpData[idxF],
				TRUE,
				that->fieldType[idxF]) != SUCCESS)
					goto done;
		}
#ifdef	__DEBUG__
		/*fprintf(stdout, "\b\b\b\b%4d", that->numData);
		fflush(stdout);*/
#endif
	}

	/**
	 ** Sort the field values
	 **/
	for (idxF=0; idxF<=numFields-1; idxF++) {
#ifdef	__DEBUG__
		DrlVTypeVectPrint(
			that->fieldValues[idxF],
			that->numFieldValues[idxF],
			that->fieldType[idxF],
			NULL);
#endif

		if (DrlVTypeVectSort(
			that->fieldValues[idxF],
			&that->numFieldValues[idxF],
			that->fieldType[idxF],
			FALSE) != SUCCESS)
				goto done;

#ifdef	__DEBUG__
		DrlVTypeVectPrint(
			that->fieldValues[idxF],
			that->numFieldValues[idxF],
			that->fieldType[idxF],
			NULL);
#endif


	}

#ifdef	_SKIP
	/**
	 ** Test consistency
	 **/
	{
		int	nCheck = 1;

		for (idxF=0; idxF<=numFields-1; idxF++)
			nCheck *= that->fieldValues[idxF];

		if (nCheck != that->numData) {
			GtoErrMsg("%s: inconsistent number of elements (%d) "
				" to form an array:\n",
				routine, that->numData);
			for (idxF=0; idxF<=numFields-1; idxF++) {
				GtoErrMsg("\t Field %2d (type %-20s) has %2d"
			}

		}


	}
#endif





	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure (line %d).\n", routine, line);
		return(NULL);
	}
	return(that);
}


/*f-------------------------------------------------------------
 * Reads data from a text file <i> fnam</i>.
 * Convenience routine for <i> DrlTVTableNewFpRead</i>.
 * Returns NULL if failed. 
 */

TVTable*
DrlTVTableNewFileRead(
	char *fnam,		/* (I) file name */
	int numFields,		/* (I) number of fields (or -1) */
	TVType *fieldTypes)	/* (I) field type [numFields] */
{
static	char	routine[] = "DrlTVTableNewFileRead";
	int	status= FAILURE;
	TVTable* that = NULL;
	FILE	*fp = NULL;

	if ((fp = fopen(fnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, DrlStrerror());
		return (NULL);
        }
	that = DrlTVTableNewFpRead(fp, numFields, fieldTypes);
	fclose(fp);
	if (that == NULL) {
		GtoErrMsg("%s: failed reading `%s'.\n", routine, fnam);
		return (NULL);
	}
	return(that);
}

/*f-------------------------------------------------------------
 * Writes the table to a file pointer <i> fp</i>.
 */

int
DrlTVTableFpWrite(
	TVTable *that,		/* (I) table */
	FILE *fp)		/* (I) if NULL, write to err log */
{
static	char	routine[] = "DrlTVTableFpWrite";
	int	status= FAILURE;

	int	idxF, idxD;

	DrlFPrintf(fp, "DATA:\n");
	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		DrlFPrintf(fp, "\t'%s'", that->fieldName[idxF]);
	}
	DrlFPrintf(fp, "\n");

	for (idxD=0; idxD<=that->numData-1; idxD++) {
	    for (idxF=0; idxF<=that->numFields-1; idxF++) {
		DrlFPrintf(fp, "\t%s", 
		    DrlVTypePrint(
			NULL,
			that->fieldType[idxF],
			DrlVTypeOffsetVect(
				that->dataValue[idxF],
				idxD,
				that->fieldType[idxF])));
	    }
	    DrlFPrintf(fp, "\n");
	}

	for (idxF=0; idxF<=that->numFields-1; idxF++) {
	    DrlFPrintf(fp, "FIELD %2d '%s' (type %s):\n",
		idxF+1,that->fieldName[idxF],
		DrlVTypeName(that->fieldType[idxF]));

	    for (idxD=0; idxD<=that->numFieldValues[idxF]-1; idxD++) {
		DrlFPrintf(fp, "\t%2d/%2d -> `%s'\n", 
		    idxD+1, that->numFieldValues[idxF],
		    DrlVTypePrint(
			NULL,
			that->fieldType[idxF],
			DrlVTypeOffsetVect(
				that->fieldValues[idxF],
				idxD,
				that->fieldType[idxF])));
	    }
	}


	/* OK */
	status = SUCCESS;
/*done:*/
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
	}
	return(status);
}





/*f-------------------------------------------------------------
 * Extracts a pivot table.
 */

int
DrlTVTableExtractMatrix(
	TVTable *that,		/* (I) table */
	int indexFieldX,	/* (I) field X index (output row) */
	int indexFieldY,	/* (I) field Y index (output col) */
	int indexFieldP,	/* (I) field page index */
	void *valFieldP,	/* (I) field page value (or NULL) */
	int indexFieldXY,	/* (I) field XY index */
	int *numItemsX,		/* (O) number of elem in X dim (output row) */
	int *numItemsY,		/* (O) number of elem in Y dim (output col) */
	TDate baseDate,		/* (I) base date for interval conversion */
	double ***matrix)	/* (O) XY matrix of values */
{
static	char	routine[] = "DrlTVTableExtractMatrix";
	int	status = FAILURE;

	double	value, **vmat = NULL;
	int	nX, nY, neqX, neqY, neqP,
		idxX, idxY, idxD;

	/*
	 *
	 */
	if ((indexFieldX < 0)  || (indexFieldX > that->numFields-1)) {
		GtoErrMsg("%s: bad field index X (%d), not in [0,%d].\n",
			routine, indexFieldX, that->numFields-1);
		goto done;
	}

	if ((indexFieldY < 0)  || (indexFieldY > that->numFields-1)) {
		GtoErrMsg("%s: bad field index Y (%d), not in [0,%d].\n",
			routine, indexFieldY, that->numFields-1);
		goto done;
	}

	if ((indexFieldP < 0)  || (indexFieldP > that->numFields-1)) {
		GtoErrMsg("%s: bad field index P (%d), not in [0,%d].\n",
			routine, indexFieldP, that->numFields-1);
		goto done;
	}

	if ((indexFieldXY < 0)  || (indexFieldXY > that->numFields-1)) {
		GtoErrMsg("%s: bad field index Z (%d), not in [0,%d].\n",
			routine, indexFieldXY, that->numFields-1);
		goto done;
	}


	nX = that->numFieldValues[indexFieldX];
	nY = that->numFieldValues[indexFieldY];

	if ((vmat = DrlDoubleMatrAlloc(0, nX-1, 0, nY-1)) == NULL)
		goto done;
	for (idxX=0; idxX<=nX-1; idxX++)
	for (idxY=0; idxY<=nY-1; idxY++)
		vmat[idxX][idxY] = 0e0;



	for (idxD=0; idxD<=that->numData-1; idxD++) {

	    if (valFieldP != NULL) {
		neqP = DrlVTypeCompare(
			DrlVTypeOffsetVect(
				that->dataValue[indexFieldP],
				idxD,
				that->fieldType[indexFieldP]),
			valFieldP,
			that->fieldType[indexFieldP]);
		if (neqP) continue;
	    }


	    for (idxX=0; idxX<=nX-1; idxX++)
	    for (idxY=0; idxY<=nY-1; idxY++) {

		neqX = DrlVTypeCompare(
			DrlVTypeOffsetVect(
				that->dataValue[indexFieldX],
				idxD,
				that->fieldType[indexFieldX]),
			DrlVTypeOffsetVect(
				that->fieldValues[indexFieldX],
				idxX,
				that->fieldType[indexFieldX]),
			that->fieldType[indexFieldX]);

		neqY = DrlVTypeCompare(
			DrlVTypeOffsetVect(
				that->dataValue[indexFieldY],
				idxD,
				that->fieldType[indexFieldY]),
			DrlVTypeOffsetVect(
				that->fieldValues[indexFieldY],
				idxY,
				that->fieldType[indexFieldY]),
			that->fieldType[indexFieldY]);
		if ((!neqX) && (!neqY)) {
			if (DrlVTypeNumValue(
				DrlVTypeOffsetVect(
					that->dataValue[indexFieldXY],
					idxD,
					that->fieldType[indexFieldXY]),
				that->fieldType[indexFieldXY],
				((void*) &baseDate),
				&value) != SUCCESS)
					goto done;
			vmat[idxX][idxY] += value;
		}
	    }
	}


#ifdef	__DEBUG__
	/* Logging */
	GtoErrMsg("%s: output table:\n", routine);
	GtoErrMsg("X : field index %d `%s' (type %s), nX=%d.\n", indexFieldX,
			that->fieldName[indexFieldX],
			DrlVTypeName(that->fieldType[indexFieldX]), nX);
	GtoErrMsg("Y : field index %d `%s' (type %s), nY=%d\n", indexFieldY,
			that->fieldName[indexFieldY],
			DrlVTypeName(that->fieldType[indexFieldY]), nY);
	if (valFieldP != NULL) {
	    GtoErrMsg("P : field index %d `%s' (type %s), value=%s\n",
			indexFieldP,
			that->fieldName[indexFieldP],
			DrlVTypeName(that->fieldType[indexFieldY]),
			DrlVTypePrint(NULL, that->fieldType[indexFieldY],
				valFieldP));
	}

	GtoErrMsg("XY: field index %d `%s' (type %s).\n", indexFieldXY,
			that->fieldName[indexFieldXY],
			DrlVTypeName(that->fieldType[indexFieldXY]));
	GtoErrMsg("\tTBL");
	for (idxY=0; idxY<=nY-1; idxY++) {
		GtoErrMsg("\t%s",
			DrlVTypePrint(
			    NULL,
			    that->fieldType[indexFieldY],
			    DrlVTypeOffsetVect(
				that->fieldValues[indexFieldY],
				idxY,
				that->fieldType[indexFieldY])));
	}
	GtoErrMsg("\n");

	for (idxX=0; idxX<=nX-1; idxX++) {
		GtoErrMsg("\t%s",
			DrlVTypePrint(
			    NULL,
			    that->fieldType[indexFieldX],
			    DrlVTypeOffsetVect(
				that->fieldValues[indexFieldX],
				idxX,
				that->fieldType[indexFieldX])));

		for (idxY=0; idxY<=nY-1; idxY++) {
			GtoErrMsg("\t%g", vmat[idxX][idxY]);
		}
		GtoErrMsg("\n");
	}
#endif


	*numItemsX = nX;
	*numItemsY = nY;
	*matrix = vmat;

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
	}
	return(status);
}

/*f-------------------------------------------------------------
 * Extracts a pivot table by field names.
 * Convenience routine for <i> DrlTVTableExtractMatrix</i>.
 */

int
DrlTVTableExtractMatrixN(
	TVTable *that,		/* (I) table */
	char *nameFieldX,	/* (I) field X name (output row) */
	char *nameFieldY,	/* (I) field Y name (output col) */
	char *nameFieldP,	/* (I) field page index */
	void *valFieldP,	/* (I) field page value (or NULL) */
	char *nameFieldXY,	/* (I) field XY index */
	int *numItemsX,		/* (O) number of elem in X dim (output row) */
	int *numItemsY,		/* (O) number of elem in Y dim (output col) */
	TDate baseDate,		/* (I) base date for interval conversion */
	double ***matrix)	/* (O) XY matrix of values */
{
static	char	routine[] = "DrlTVTableExtractMatrixN";
	int	status = FAILURE;

	int	indexFieldX;
	int	indexFieldY;
	int	indexFieldP;
	int	indexFieldXY;


	if ((indexFieldX = DrlTVTableFindFieldName(that, nameFieldX)) == -1)
		goto done;
	if ((indexFieldY = DrlTVTableFindFieldName(that, nameFieldY)) == -1)
		goto done;
	if ((indexFieldP = DrlTVTableFindFieldName(that, nameFieldP)) == -1)
		goto done;
	if ((indexFieldXY = DrlTVTableFindFieldName(that, nameFieldXY)) == -1)
		goto done;



	if (DrlTVTableExtractMatrix(
		that,
		indexFieldX,
		indexFieldY,
		indexFieldP,
		valFieldP,
		indexFieldXY,
		numItemsX,
		numItemsY,
		baseDate,
		matrix) != SUCCESS)
			goto done;

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
	}
	return(status);
}


/*--------------------------------------------------------------
 * Finds a field name
 */

int
DrlTVTableFindFieldName(
	TVTable *that,		/* (I) table */
	char *nameField)	/* (I) field X name */
{
	int	idxF;

	for (idxF=0; idxF<=that->numFields-1; idxF++) {
		if (!strcmp(nameField, that->fieldName[idxF]))
			return(idxF);
	}
	GtoErrMsg("drlTVTableFindFieldName: can't find name `%s':\n",
		nameField);
	for (idxF=0; idxF<=that->numFields-1; idxF++) {
	    GtoErrMsg("\tFIELD %2d/%2d '%s' (type %s)\n",
		idxF+1, that->numFields,
		that->fieldName[idxF],
		DrlVTypeName(that->fieldType[idxF]));
	}

	return(-1);
}


/*f-------------------------------------------------------------
 * Extracts a pivot table by field names
 * as a <i> TSwaptionMatrix2D</i> data structure.
 * Convenience routine for <i> DrlTVTableExtractMatrix</i>.
 */

int
DrlTVTableExtractTSwaptionMatrix2D(
	TVTable *that,		/* (I) table */
	char *nameFieldX,	/* (I) field X index (row) */
	char *nameFieldY,	/* (I) field Y index (column) */
	char *nameFieldP,	/* (I) field page index */
	void *valFieldP,	/* (I) field page value (or NULL) */
	char *nameFieldXY,	/* (I) field XY index */
	TDate baseDate,		/* (I) base date for interval conversion */
        long freq,              /* (I) swmat frequency (1,2,4,12) */
	TSwaptionMatrix2D **mat)/* (O) */
{
static	char	routine[] = "drlTVTableExtractSwMat";
	int	status = FAILURE;

	int	numItemsX;
	int	numItemsY;
	double	**matrix = NULL;
	TSwaptionMatrix2D *swMat = NULL;
	int	idxE, idxM;

        int	indexFieldX, indexFieldY;
	double	value;


        if ((indexFieldX = DrlTVTableFindFieldName(that, nameFieldX)) == -1)
                goto done;
        if ((indexFieldY = DrlTVTableFindFieldName(that, nameFieldY)) == -1)
                goto done;



	if (DrlTVTableExtractMatrixN(
		that,
		nameFieldX,
		nameFieldY,
		nameFieldP,
		valFieldP,
		nameFieldXY,
		&numItemsX,
		&numItemsY,
		baseDate,
		&matrix) != SUCCESS)
			goto done;



	if ((swMat = DrlTSwaptionMatrix2DNew(
        	FALSE,
        	freq,
        	numItemsX, NULL,
        	numItemsY, NULL)) == NULL)
			goto done;


	for (idxE=0; idxE<=numItemsX-1; idxE++) {
		if (DrlVTypeNumValue(
			DrlVTypeOffsetVect(
				that->fieldValues[indexFieldX],
				idxE,
				that->fieldType[indexFieldX]),
			that->fieldType[indexFieldX],
			((void*) &baseDate),
			&value) != SUCCESS)
				goto done;

		TSWAPTION_MATRIX2D_EXP(swMat, idxE) = value;
	}


	for (idxM=0; idxM<=numItemsY-1; idxM++) {
		if (DrlVTypeNumValue(
			DrlVTypeOffsetVect(
				that->fieldValues[indexFieldY],
				idxM,
				that->fieldType[indexFieldY]),
			that->fieldType[indexFieldY],
			((void*) &baseDate),
			&value) != SUCCESS)
				goto done;

		TSWAPTION_MATRIX2D_MAT(swMat, idxM) = value;
	}

	for (idxE=0; idxE<=numItemsX-1; idxE++) 
	for (idxM=0; idxM<=numItemsY-1; idxM++) {
		TSWAPTION_MATRIX2D_VOL(swMat, idxE, idxM) = matrix[idxE][idxM];
	}



	*mat = swMat;

	/* OK */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(matrix, 0, numItemsX-1, 0, numItemsY-1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
	}
	return(status);
}






