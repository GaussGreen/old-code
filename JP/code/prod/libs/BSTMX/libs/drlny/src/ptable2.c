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

TCSheet*
DrlTCSheetNewEmpty(
	int numRowsMax,		/* (I) max number of rows */
	int numColsMax)		/* (I) max number of cols */
{
static	char	routine[] = "DrlTCSheetNewEmpty";
	int	status= FAILURE;

	TCSheet	*that = NULL;
	int	idxR, idxC;

	that = NEW(TCSheet);
	ASSERT_OR_DONE(that != NULL);

	that->cells = DrlCharPMatrAlloc(
		0, numRowsMax-1,
		0, numColsMax-1);
	ASSERT_OR_DONE(that->cells != NULL);
	that->numRowsMax = numRowsMax;
	that->numColsMax = numColsMax;

	for (idxR=0; idxR<=that->numRowsMax-1; idxR++)
	for (idxC=0; idxC<=that->numColsMax-1; idxC++) {
		that->cells[idxR][idxC] = NULL;
	}
	that->numRows = 0;
	that->numCols = 0;



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
DrlTCSheetFree(TCSheet* that)
{
	int	idxR, idxC;

	if (that == NULL) return(SUCCESS);

	for (idxR=0; idxR<=that->numRowsMax-1; idxR++)
	for (idxC=0; idxC<=that->numColsMax-1; idxC++) {
		if (that->cells[idxR][idxC])
			FREE(that->cells[idxR][idxC]);
	}

	DrlCharPMatrFree(that->cells,
		0, that->numRowsMax-1,
		0, that->numColsMax-1);

	FREE(that);

	return (SUCCESS);
}

/*f-------------------------------------------------------------
 * Sets an element.
 */

int
DrlTCSheetSet(TCSheet* that, int row, int col, const char* value)
{
static	char	routine[] = "DrlTCSheetSet";

	if ((row < 0) || (row > that->numRowsMax)) {
		GtoErrMsg("%s: bad row %d (min 0, max %d).\n",
			routine, 0, that->numColsMax);
		return(NULL);
	}
	if ((col < 0) || (col > that->numColsMax)) {
		GtoErrMsg("%s: bad col %d (min 0, max %d).\n",
			routine, 0, that->numColsMax);
		return(NULL);
	}
	if (that->cells[row][col] != NULL) {
		FREE(that->cells[row][col]);
		that->cells[row][col] = NULL;
	}
	that->cells[row][col] = NEW_ARRAY(char, strlen(value)+1);
	if (that->cells[row][col] == NULL) {
		return(FAILURE);
	}
	strcpy(that->cells[row][col], value);
	return(SUCCESS);
}


/*f-------------------------------------------------------------
 * Gets an element a string.
 */

char*
DrlTCSheetGet(TCSheet* that, int row, int col)
{
static	char	routine[] = "DrlTCSheetGet";
static	char	emptyString = '\0';

	if ((row < 0) || (row > that->numRowsMax)) {
		GtoErrMsg("%s: bad row %d (min 0, max %d).\n",
			routine, 0, that->numColsMax);
		return(NULL);
	}
	if ((col < 0) || (col > that->numColsMax)) {
		GtoErrMsg("%s: bad col %d (min 0, max %d).\n",
			routine, 0, that->numColsMax);
		return(NULL);
	}
	return (that->cells[row][col] != NULL ?
		that->cells[row][col] : &emptyString);
}


/*f-------------------------------------------------------------
 * Gets an element of a specific type "fieldType" and
 * puts it in the pointer "ptr" (that must be of corrseponding type).
 */


int
DrlTCSheetGetVType(
	TCSheet* that,		/* (I) table */
	int row,		/* (I) row number */
	int col, 		/* (I) column number */
	TVType fieldType,	/* (I) element type */
	void *ptr)		/* (O) */
{
static	char	routine[] = "DrlTCSheetGetVType";
	int	status = FAILURE;
	char	*p;

	/* get string element */
	if ((p = DrlTCSheetGet(that, row, col)) == NULL)
		goto done;

	IF_FAILED_DONE( DrlVTypeScan(p, fieldType, ptr));

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure (element %d x %d).\n",
			routine, row, col);
	}
	return(status);
}



/*f-------------------------------------------------------------
 * Reads data from a file <i> fp</i>
 * and creates a new range 
 * The format for text file is
 * \begin{verbatim}
 * <data 1  1>  ..... <data 1  Nc>
 * ...
 * <data Nr 1>  ..... <data Nr Nc>
 * ...
 * \end{verbatim}
 * The file scanning stops when an empty line or EOF is encountered.
 * Returns NULL if failed. 
 */

int
DrlTCSheetFpRead(
	TCSheet **that,		/* (O) output sheet */
	FILE *fp)		/* (I) file pointer */
{
static	char	routine[] = "DrlTCSheetNewFpRead";
	int	status= FAILURE;

	TCSheet	*csheet = NULL;
	int	line = 1;
	char	bufline[2048], buf[2048], *p;
	int	idxR, idxC;


	/*$$$ preallocate: need to be changed !!! */
	csheet = DrlTCSheetNewEmpty(64, 64);


	/* Read data */
	idxR = 0;
	while (DrlFGetLine(bufline, sizeof(bufline), fp, &line) != NULL) {
	    /* Empty line is end of the table */
	    if (strlen(bufline) == 0)
		break;

	    /* Read fields of data point */
	    p = bufline;
	    idxC = 0;
	    while (DrlStrScanString(&p, buf) != NULL) {
		IF_FAILED_DONE( DrlTCSheetSet(csheet, idxR, idxC, buf));
		idxC++;
	    }
	    idxR++;
	}
	csheet->numRows = idxR;
	csheet->numCols = idxC;

	*that = csheet;

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure (line %d).\n", routine, line);
	}
	return(status);
}


/*f-------------------------------------------------------------
 * Reads data from a text file <i> fnam</i>.
 * Convenience routine for <i> DrlTCSheetNewFpRead</i>.
 * Returns NULL if failed. 
 */

int
DrlTCSheetFileRead(
	TCSheet **that,		/* (O) output sheet */
	const char *fnam)	/* (I) file name */
{
static	char	routine[] = "DrlTCSheetNewFileRead";
	int	status= FAILURE;
	FILE	*fp = NULL;

	if ((fp = fopen(fnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, fnam, DrlStrerror());
		goto done;
        }
	status = DrlTCSheetFpRead(that, fp);
	fclose(fp);
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed reading `%s'.\n", routine, fnam);
	}
	return(status);
}

/*f-------------------------------------------------------------
 * Writes the table to a file pointer <i> fp</i>.
 */

int
DrlTCSheetFpWrite(
	TCSheet *that,		/* (I) table */
	FILE *fp)		/* (I) if NULL, write to err log */
{
static	char	routine[] = "DrlTCSheetFpWrite";
	int	status= FAILURE;

	int	idxR, idxC;


	for (idxR=0; idxR<=that->numRows-1; idxR++) {
	    for (idxC=0; idxC<=that->numCols-1; idxC++) {
		DrlFPrintf(fp,"\"%s\"\t",
			(that->cells[idxR][idxC] ? 
			that->cells[idxR][idxC] : ""));
	    }
	    DrlFPrintf(fp,"\n");
	}

	/* OK */
	status = SUCCESS;
/*done:*/
	if (status != SUCCESS) {
		GtoErrMsg("%s: failure.\n", routine);
	}
	return(status);
}



