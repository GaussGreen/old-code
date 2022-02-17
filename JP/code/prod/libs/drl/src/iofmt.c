/************************************************************************
 * Module:	DRL - IO
 * Function:	I/O utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#define DRL_SYSCALL
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#ifdef	DRL_CLIB
# include "convert.h"
# include "date_sup.h"
# include "yearfrac.h"
# include "tcurve.h"
#endif

#include "drlmem.h"	/* XXXVectAlloc */
#include "drlstr.h"
#include "drlvtype.h"

#include "drlio.h"	/* prototype consistency */

#define	MAXARG		32	/* maximum # of args per line */
#define	MAXLEN		256	/* maximum # of elements in array */

#define	__DEBUG__
#undef	__DEBUG__


/*f-------------------------------------------------------------
 * I/O: read next line (skip comments).
 *
 * <br><br>
 * Acts like the standard ANSI routine <i> fgets</i>, except
 * that is discards all lines starting by '\#' (commented)
 * and counts number of lines: the counter <i> line</i>
 * is incremented every time a line is read (useful for
 * error reporting).
 */

char*
DrlFGetLine(char *s, int n, FILE *fp, int *line)
{
	int	acceptFlag = 0;
	char	*p;

	do {
		/* read next line */
		if (fgets(s, n, fp) == NULL) return(NULL);
		if (line != NULL) (*line)++;

		/* check for empty line (comment # ;) */
		p = s;
		while ((isspace(*p)) && (*p != '\0')) p++;
		if ((*p != '#') && (*p != ';') && (*p != '\0'))
			acceptFlag = 1;
	} while (acceptFlag == 0);
	return(s);
}


/*f-------------------------------------------------------------
 * I/O: advance to next char (skip comments).
 *
 * <br><br>
 * Advances in the file to the next non space character
 * (without reading it).
 * Skips comments (i.e. discards all lines starting by the pound sign).
 * Returns SUCCESS or FAILURE if EOF is encountered.
 */

int
DrlFAdvanceToNextChar(FILE *fp)
{
	int	c;

#ifdef	_OLD
	while (((c = getc(fp)) != EOF) && (isspace(c)));
#else
	/* skip over '#' commented lines */
	while (((c = getc(fp)) == '#') || (isspace(c))) {
	    if (c == '#') {
		while (((c = getc(fp)) != EOF) && (c != '\n'));
		if (c == EOF) return(FAILURE);
	    }
	}
#endif
	if (c == EOF) {
		return(FAILURE);
	}

	ungetc(c, fp);
	return(SUCCESS);
}


/*f-------------------------------------------------------------
 * I/O: advance to specific token (skip comments).
 *
 * <br><br>
 * Advances in the file util the string <i> token</i> is encountered.
 * (a block of characters either enclosed
 * in double quotes or not containing any space).
 * Returns SUCCESS/FAILURE.
 */

int
DrlFAdvanceToToken(FILE *fp, char *token)
{
static	char	buf[64];
	for (;;) {
	    if (DrlFGetToken(fp, NULL, buf, sizeof(buf))!= SUCCESS)
		return (FAILURE);
	    if (!strcmp(buf, token)) {
	 	return(SUCCESS);
	    }
	}
	return(FAILURE);
}


/*f-------------------------------------------------------------
 * I/O: read next token.
 *
 * <br><br>
 * Get the next token (a block of characters either enclosed
 * in double brackets or not containing any space)
 * and puts it in the string <i> s</i> of length <i> sz</i>
 * (removes any double quotes).\\
 * <i> WARNING: does not skip comments '\#'.</i> 
 * Returns SUCCESS/FAILURE.
 */

int
DrlFGetToken(FILE *fp, char *sep, char *s, int sz)
{
	char	*q;
	int	c;
static	char	buf[256];

	q = (s != NULL ? s : buf);

	while (((c = getc(fp)) != EOF) && (isspace(c)));
	if (c == EOF)
		return(FAILURE);

	if (c != '"') {
		*q++ = c;
	        while (((c = getc(fp)) != EOF) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	} else {
	        while (((c = getc(fp)) != EOF) && (c != '"'))
			*q++ = c;
		*q = '\0';
	}

	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * I/O: read double (skip comments).
 *
 * <br><br>
 * Scans a double value in the file <i> fp</i> and puts
 * the result in <i> val</i>.
 * Supports the currency format (e.g. 100,000,000.00).
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanDouble(FILE *fp, double *val)
{
	char	buf[64];
	if (DrlFScanString(fp, buf) != SUCCESS) return(FAILURE);
	return (DrlCurScan(buf , val));
}


/*f--------------------------------------------------------------
 * I/O: read long (skip comments).
 *
 * <br><br>
 * Scans a long value in the file <i> fp</i> and puts
 * the result in <i> val</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanLong(FILE *fp, long *val)
{
	char	buf[64];
	if (DrlFScanString(fp, buf) != SUCCESS) return(FAILURE);
	return (sscanf(buf, "%ld", val) != 1);
}




/*f---------------------------------------------------------------------
 * I/O: read long in predefined token values (skip comments).
 *
 * <br><br>
 * Scans a long value in a file <i> fp</i> by trying to match
 * the string with different values. On successful return, the
 * value is put in <i> value</i>. See <i> StringLongValueScan</i>.
 */

int
DrlFScanStringLongValue(
	FILE *fp,		/* (I) input file */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...)
{
static	char		routine[] = "DrlFScanStringLongValue";
	int		status = FAILURE;
	char		buf[64];
	va_list		ap;

	va_start(ap, value);
	if (DrlFScanString(fp, buf) != SUCCESS)
		goto done;
	if (DrlStrLongValueScanV(buf, "[?]", value, ap) != SUCCESS)
		goto done;
	status = SUCCESS;
done:
	va_end(ap);
	if (status != SUCCESS) {
	    DrlErrMsg("%s: scanning `%s' failed.\n", routine, buf);
	}
	return(status);
}

/*f---------------------------------------------------------------------
 * I/O: print long in predefined token values.
 *
 * <br><br>
 * Prints a long value in a string <i> s</i> by trying to match
 * it with different possible values. See <i> DrlString\-Long\-Value\-Print</i>
 * Returns <i> s</i>.
 */

int
DrlFPrintStringLongValue(
	FILE *fp,		/* (O) output file */
	long value,		/* (I) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...)
{
	va_list		ap;
	va_start(ap, value);
	DrlFPrintf(fp, "%s", 
		DrlStrLongValuePrintV(NULL, "[?]", value, ap));
	va_end(ap);
	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * I/O: read integer (skip comments).
 *
 * <br><br>
 * Scans a integer value in the file <i> fp</i> and puts
 * the result in <i> val</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanInt(FILE *fp, int *val)
{
	char	buf[64];
	if (DrlFScanString(fp, buf) != SUCCESS) return(FAILURE);
	return (sscanf(buf, "%d", val) != 1);
}


/*f--------------------------------------------------------------
 * I/O: read char string (skip comments).
 *
 * <br><br>
 * Scans a string value in the file pointer <i> fp</i> and puts
 * the result in <i> s</i>
 * (a string being a block of characters either enclosed
 * in double brackets or not containing any space).
 * Skips all end of lines after '\#' is encountered.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanString(FILE *fp, char *s)
{
	char	*q;
	int	c;

	q = s;


#ifdef	_OLD
	while (((c = getc(fp)) != EOF) && (isspace(c)));
#else
	/* skip over '#' commented lines */
	while (((c = getc(fp)) == '#') || (isspace(c))) {
	    if (c == '#') {
		while (((c = getc(fp)) != EOF) && (c != '\n'));
		if (c == EOF) return(FAILURE);
	    }
	}
#endif
	if (c == EOF)
		return(FAILURE);

	if (c != '"') {
		/*return (fscanf(fp, "%s", s) != 1);*/
		*q++ = c;
	        while (((c = getc(fp)) != EOF) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	} else {
	        while (((c = getc(fp)) != EOF) && (c != '"'))
			*q++ = c;
		*q = '\0';
	}

	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * I/O: read variable type (skip comments).
 *
 * <br><br>
 * Scans a variable type value of type <i> type</i>
 * in the file <i> fp</i> and puts the result in <i> valPtr</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanVType(FILE *fp, DVType type, void *valPtr)
{
static	char	buf[1024];
	if (DrlFScanString(fp, buf) != SUCCESS)
		return(FAILURE);
	/*if (DrlFGetToken(fp, NULL, buf, sizeof(buf)) != SUCCESS)
		return(FAILURE);*/

	if (DrlVTypeScan(buf, type, valPtr) != SUCCESS)
		return(FAILURE);

	return(SUCCESS);
}



/*f--------------------------------------------------------------
 * I/O: print variable type.
 *
 * <br><br>
 * Prints a variable type value <i> valPtr</i> of type <i> type</i>
 * in the file <i> fp</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFPrintVType(FILE *fp, DVType type, void *valPtr)
{
static	char	buf[256];

	strcpy(buf, DrlVTypePrint(NULL, type, valPtr));

	if (fp == NULL) {
		DrlErrMsg(buf);
	} else {
		fputs(buf, fp);
	}

	return(SUCCESS);
}


/*f---------------------------------------------------------------------
 * I/O: read structured input (skip comments).
 *
 * <br><br>
 * Performs a formatted scan of variables and arrays of variables
 * on a file pointer.
 * The routine has variable number of arguments
 * and  <b> the last argument should always be </b><i> 0</i>.
 * The three first arguments are mandatory:
 * <br>
 * <br>[fp] is the file pointer on which the scan is to be done.
 * <br>[line] is a counter that is incremented every time a line
 * is read (useful for debugging).
 * <br>[tokens] is a string defining the characters that are admitted
 * as separators when reading an array.
 * <br>
 * The remaining arguments define the type of formatted scan
 * to be performed.
 * \vspace{2mm}\par\noindent
 * %
 * %
 * <b> Single variable scan:</b> the recognized format is
 * of the form
 * \begin{verbatim}
 * <tagName> = <variable>
 * \end{verbatim}
 * The arguments of the function are
 * \begin{verbatim}
 *    ..., VAR_T, char *tagName, int varType, void *var, ...
 * \end{verbatim}
 * where
 * <br>
 * <br>[DRL\_CVAR\_T] is a key (defined in the header).
 * <br>[tagName] is the identifier in the text file,
 * <br>[varType] the type of the variable: use one of the predefined
 *  key listed below,
 * <br>[var] the address of variable in 
 *     which the data is to be stored casted to <i> void \*</i>.
 * <br>
 * <b> C Example:</b> to read the line
 * \begin{verbatim}
 * # volatility (in %)
 * vol = 20.45
 * \end{verbatim}
 * call the routine as in 
 * \begin{verbatim}
 *     double vol;
 *     int    line=0;
 *     ...
 *     VarStructScan(stdin, &line, NULL,
 *         VAR_T, "vol", DRL_DOUBLE_T, (void*) &vol,
       0);
 * \end{verbatim}
 * %
 * 
 * \vspace{5mm}\par\noindent
 * <b> Array scan:</b> the recoignized format is
 * of the form
 * \begin{verbatim}
 *   <tagname>
 *   <var#1 of array #1>  ...  <var#1 of array #NVAR> 
 *   ...
 *   <var#K of array #1>  ...  <var#K of array #NVAR> 
 *   end
 * \end{verbatim}
 * The end of the array is specified by the keyword <i> end</i>
 * (with no space before it).
 * The arguments of the function are
 * \begin{verbatim}
 *      ARRAY_T, char *tagName, int *nItems, int nVar, 
 *      int varType1, int varAllocType1, void* var1,
 *      ...
 *      int varTypeN, int varAllocTypeN, void* varN, 
 * \end{verbatim}
 * where
 * <br>
 * <br>[DRL\_CARRAY\_T] is a key (defined in the header).
 * <br>[tagName] is the identifier of the array in the text file,
 * <br>[nItems] on exit, contains the length of the array read,
 * <br>[nVar] on entry, the number of arrays expected,
 * <br>[varType1..N] the type of the variable: use one of the predefined
 *   key listed below.
 * <br>[varAllocType1..N] if TRUE, expects the adress of a pointer
 *  and performs a malloc. If FALSE, expects a pointer to a vector
 *  already allocated (should be large enough).
 * <br>[var1..N] the adress of a pointer or a pointer in 
 *   which the data is to be stored casted to <i> void \*</i>
 *   (see previous argument).
 * <br>
 */

int
DrlFScanStruct(
	FILE *fp,
	int *line,
	char *toks,
	/* DRL_CVAR_T, char *tagName, DVType varType, void *var,
	 * DRL_CARRAY_T, char *tagName, int *nItems, int nVar,
	 *	DVType varType1, int varAllocType1, void* var1,
	 *	...
	 *	DVType varTypeN, int varAllocTypeN, void* varN,
	 * DRL_CVECTOR_T, char *tagName, int *nItems,
	 *	DVType varType, void* var,
	 * DRL_CDVECTOR_T, char *tagName,
	 *	int *nItems, double **v,
	 * DRL_CDMATRIX_T, char *tagName,
	 *	int *nx, int *ny, double ***v,
	 * DRL_LILARRAY_L, char *tagName, DVType varType, void **var,
	 * 0)
	 */ ...)
{
static	char		routine[] = "DrlFScanStruct";
	int		status = FAILURE;

static	char		tokDef[] = " \t;";
	va_list		ap;
	int		i, j, m,
			numItems,
			*nItems, *nx, *ny,
			lineDef = 0,
			nVar,
			varAllocType[MAXARG];
	DVType		varType[MAXARG],
			strType;
	char		*tagName = NULL,
			nextLine[256],
			*tokens, 
			*p, *q,
			t0[128], *t1;			/* tmp variables */
	void		*var[MAXARG],
			*vtv[MAXARG],
			*vp;

#undef	ERROR
#define	ERROR(str)	{DrlErrMsg("%s: %s (name `%s', line %d)\n",\
			 routine, (str), (tagName != NULL ? tagName : "N/A"),\
			*line); status = FAILURE; goto done;}
#undef	GETNEXTLINE
#define	GETNEXTLINE	{if (DrlFGetLine(nextLine, sizeof(nextLine), \
			fp, line) == NULL) ERROR("EOF");}



	/* initialize variable */
	for (i=0; i<=MAXARG-1; i++) vtv[i] = NULL;
	tokens = (toks == NULL ? tokDef : toks);
	line = (line == NULL ? &lineDef : line);

	/*
	 *
	 */
	va_start(ap, toks);


	while ((strType = (DVType) va_arg(ap, DVType)) != 0) {
	    switch (strType) {
	    /*
	     * Single variable
	     */
	    case DRL_CVAR_T:
		/* get arguments */
		tagName    = (char*)    va_arg(ap, char*);
		varType[0] = (DVType) va_arg(ap, DVType);
		var[0]     = (void*)    va_arg(ap, void*);

		/* read the variable */
		GETNEXTLINE;

		if (tagName != NULL) {
		    /* format is "<tagName>=<value>" */
		    if ((p = strchr(nextLine, '=')) == NULL)
			ERROR("bad format");
		    *(p++) = '\0';

		    /* scan tagName */
		    if ((sscanf(nextLine, "%s", t0) != 1) ||
			    (strncmp(t0, tagName, strlen(tagName)) != 0))
			    ERROR("bad tag name");
		} else {
		    /* format is "<value>" */
		    p = nextLine;
		}

		/* scan value */
		if (DrlVTypeScan(p, varType[0], var[0]) != SUCCESS)
				ERROR("bad format");
		break;

	    /*
	     * Arrays of variables
	     */
	    case DRL_CARRAY_T:
		/* get arguments */
		tagName = (char*) va_arg(ap, char*);
		nItems  = (int*) va_arg(ap, int*);
		nVar    = (int) va_arg(ap, int);
		for (i=0; i <= nVar-1; i++) {
		    varType[i]      = (DVType) va_arg(ap, DVType);
		    varAllocType[i] = (int) va_arg(ap, int);
		    var[i]          = (void*) va_arg(ap, void*);
		    if (DrlVTypeCheckValidCType(varType[i]) != 0) 
			ERROR("bad C type in ARRAY_T");
		}

		/* Get number of elements */
		if (DrlFScanInt(fp, nItems) != SUCCESS) {
			ERROR("can't read num items in CARRAY_T");
		}
		/* If length is zero, nothing else to do */
		if (*nItems <= 0) {
		    for (i=0; i<=nVar-1; i++)
			if (varAllocType[i] != 0)
			    *((void**) var[i]) = NULL;
		    status = SUCCESS;
		    goto done;
		}

		/* allocate memory */
		for (i=0; i<=nVar-1; i++) {
		    if (varAllocType[i] != 0) {
			/* malloc needed */
			*((void**) var[i]) =
				DrlVTypeVectAlloc(*nItems, varType[i]);
			if (*((void**) var[i]) == NULL)
				ERROR("malloc failed");
			var[i] = *((void**) var[i]);
		    }
		}

		/* read data */
		for (m=0; m<=*nItems-1; m++) {
		    for (i=0; i<=nVar-1; i++) {
			if ((DrlFScanString(fp, t0) != SUCCESS) ||
			    (DrlVTypeScan(t0, varType[i], 
				DrlVTypeOffsetVect(var[i],
				m, varType[i])) != SUCCESS)) {
			    DrlErrMsg("%s: CARRAY_T `%s': can't read "
				"item %d in line %d (type %s).\n",
				routine, 
				(tagName != NULL ? tagName : "N/A"),
				i+1, m+1,
				DrlVTypeName(varType[i]));
			    goto done;
			}
		    }
		}

		break;


	    /*
	     * Vector of arbitary type
	     */
	    case DRL_CVECTOR_T:
		tagName    = (char*) va_arg(ap, char*);
		nItems     = (int*) va_arg(ap, int*);
		varType[0] = (DVType) va_arg(ap, DVType);
		var[0]     = (void*) va_arg(ap, void*);

		if (DrlFScanInt(fp, nItems) != SUCCESS)
			ERROR("can't read vector length");

		vtv[0] = DrlVTypeVectAlloc(*nItems, varType[0]);
		if (vtv[0] == NULL) goto done;

		for (i=0; i <= *nItems-1; i++) {
			if (DrlFScanVType(fp, varType[0], 
				DrlVTypeOffsetVect(vtv[0], i, varType[0]))
				!= SUCCESS) {
				DrlErrMsg("%s: can' read DRL_CVECT_T element "
					" #%d/%d (name `%s', line %d)\n",
					routine, i+1, *nItems,
					(tagName != NULL ? tagName : "N/A"),
					*line);
				goto done;
			}
		}

		*((void**) var[0]) = vtv[0];
		break;


	    /*
	     * Vector of doubles
	     */
	    case DRL_CDVECTOR_T:
		tagName = (char*) va_arg(ap, char*);
		nItems  = (int*) va_arg(ap, int*);
		var[0]  = (void*) va_arg(ap, double**);

		if (DrlFScanInt(fp, nItems) != SUCCESS)
			ERROR("bad format");
		vp = (void*) DrlVTypeVectAlloc(*nItems, DRL_DOUBLE_T);
		if (vp == NULL) ERROR("malloc failed");

		for (i=0; i<=*nItems-1; i++) {
			if (DrlFScanDouble(fp, ((double*)vp)+i) != SUCCESS) 
				ERROR("bad format");
		}
		*((double**) var[0]) = (double*) vp;
		break;

	    /*
	     * Matrix of doubles
	     */
	    case DRL_CDMATRIX_T:
		tagName = (char*) va_arg(ap, char*);
		nx      = (int*)  va_arg(ap, int*);
		ny      = (int*)  va_arg(ap, int*);
		var[0]  = (void*) va_arg(ap, double***);


		if (DrlFScanInt(fp, nx) != SUCCESS)
			ERROR("DRL_CDMATRIX_T bad format");
		if (DrlFScanInt(fp, ny) != SUCCESS)
			ERROR("DRL_CDMATRIX_T bad format");

		vp = (void*) DrlVTypeMatrAlloc(*nx, *ny, DRL_DOUBLE_T);
		if (vp == NULL) ERROR("DMATRIX_T malloc failed");

		for (i=0; i<=*nx-1; i++) 
		for (j=0; j<=*ny-1; j++) {
			if (DrlFScanDouble(fp, ((double**) vp)[i]+j)
			    != SUCCESS)
				ERROR("DMATRIX_T bad format");
		}

		*((double***) var[0]) = (double**) vp;
		break;

#ifdef	_SKIP

	    case DRL_LILARRAY_L:
		/* Format is:
		 * <# of elem> <elem1> elem2> .... <elem n>
		tagName    = (char*)    va_arg(ap, char*);
		varType[0] = (DVType) va_arg(ap, DVType);
		var[0]     = (void*)    va_arg(ap, void*);

		/* read size */
		if (FScanInt(fp, &numItems) != 


		vp = DrlVTypeVectAlloc(sz, LToCType(varType[0]));
#endif




		break;

	    /*
	     * Wrong format
	     */
	    default:
		DrlErrMsg("%s: bad arg type\n", routine);
		return(1);
	    }
	}

	/* made it through */
	status = 0;
done:
	va_end(ap);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);

#undef	ERROR
#undef	GETNEXTLINE
}



/*f---------------------------------------------------------------------
 * I/O: print structured input.
 *
 * <br><br>
 * Performs a formatted print of variables and arrays of variables
 * on a file pointer "fp". "line" is a counter that is incremented
 * at evry printed line, and "token" is the character printed
 * between arrays (usually the tabulation <i> '\\t'</i>).
 * The format and variable arguments are identical to
 * <i> VarStructScan</i>.
 * Returns 0 iff successful.
 */

int
DrlFPrintStruct(
	FILE *fp,
	int *line,
	int token,
	/* DRL_CVAR_T, char *tagName, int varType, void *var,
	 * DRL_CARRAY_T, char *tagName, int nItems, int nVar,
	 *	int varType1, void* var1,
	 *	...
	 *	int varTypeN, void* varN,
	 * DRL_CVECTOR_T, char *tagName, int nItems,
	 *	DVType varType, void* var,
	 * DRL_CDVECTOR_T, char *tagName,
	 * 	int nItems, double *v,
	 * DRL_CDMATRIX_T, char *tagName,
	 *	int nx, int ny, double **v,
	 * 0)
	 */ ...)
{
static	char		routine[] = "DrlFPrintStruct";
	int		status = FAILURE;
	va_list		ap;
	int		line1=0, 
			i, j, m,
			nItems, nx, ny,
			nVar;
	DVType	varType[MAXARG],
			strType;
	char		*tagName,
			nextLine[256],
			t0[128],
			tkn[2];
	void		*vtv[MAXARG];

#undef	ERROR
#define	ERROR(str)	{DrlErrMsg("%s: %s (name `%s', line %d)\n",\
			 routine, (str), (tagName != NULL ? tagName : "N/A"),\
			*line); status = FAILURE; goto done;}

	/* */
	if (line == NULL) line = &line1;
	tkn[0] = token; tkn[1] = '\0';

	va_start(ap, token);

	while ((strType = (DVType) va_arg(ap, DVType)) != 0) {
	    switch (strType) {
	    case DRL_CVAR_T:
		tagName    = (char*) va_arg(ap, char*);
		varType[0] = (DVType) va_arg(ap, DVType);
		vtv[0]     = (void*) va_arg(ap, void*);

		/* print the variable */
		if (tagName != NULL) {
		    DrlFPrintf(fp, "%s = %s\n", tagName,
			DrlVTypePrint(NULL, varType[0], vtv[0]));
		} else {
		    DrlFPrintf(fp, "%s\n",
			DrlVTypePrint(NULL, varType[0], vtv[0]));
		}
		break;

	    case DRL_CARRAY_T:
		/* get arguments */
		tagName = (char*) va_arg(ap, char*);
		nItems  = (int) va_arg(ap, int);
		nVar    = (int) va_arg(ap, int);
		for (i=0; i <= nVar-1; i++) {
		    varType[i] = (DVType) va_arg(ap, DVType);
		    vtv[i]     = (void*) va_arg(ap, void*);
		    if (DrlVTypeCheckValidCType(varType[i]) != 0) 
			ERROR("bad C type in DRL_CARRAY_T");
		}

		/* print the array */
		if (tagName != NULL) DrlFPrintf(fp, "%s\n", tagName);
		DrlFPrintf(fp, "\t%d\n", nItems);
		for(m=0; m<=nItems-1;m++) {
		    strcpy(nextLine,"");
		    for (i=0; i <= nVar-1; i++) {
			DrlVTypePrintSmart(t0, varType[i],
				DrlVTypeOffsetVect(vtv[i], m, varType[i]));
			strcat(nextLine, t0);
			if (i != nVar-1) {
				strcat(nextLine, " ");
			}
		    }

		    DrlFPrintf(fp, "%s\n", nextLine);
		    (*line)++;
		}
		break;

#ifdef	_OLD
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		/* print the array */
		DrlFPrintf(fp, "#\n");
		if (tagName != NULL) DrlFPrintf(fp, "%s\n", tagName);
		for(m=0; m<=nItems-1;m++) {
		    strcpy(nextLine,"\t");

		    for (i=0; i <= nVar-1; i++) {
			DrlVTypePrint(t0, varType[i],
				DrlVTypeOffsetVect(vtv[i], m, varType[i]));
			strcat(nextLine, t0);
			if (i != nVar-1) {
				strcat(nextLine, " ");
			}
		    }

		    DrlFPrintf(fp, "%s\n", nextLine);
		    (*line)++;
		}
	 	DrlFPrintf(fp, "end\n");
#endif
		break;

	    /*
	     * Vector of arbitary type
	     */
	    case DRL_CVECTOR_T:
		tagName    = (char*) va_arg(ap, char*);
		nItems     = (int) va_arg(ap, int);
		varType[0] = (DVType) va_arg(ap, DVType);
		vtv[0]     = (void*) va_arg(ap, void*);

		if (tagName != NULL) DrlFPrintf(fp, "# %s\n", tagName);
		DrlFPrintf(fp, "%d\n", nItems);

		for (i=0; i<=nItems-1; i++) {
			DrlFPrintf(fp, "\t%s\n",
			    DrlVTypePrint(NULL, varType[0],
				DrlVTypeOffsetVect(vtv[0], i, varType[0])));
		}
		break;


	    /*
	     * Print vector of doubles
	     */
	    case DRL_CDVECTOR_T:
		tagName = (char*) va_arg(ap, char*);
		nItems  = (int)   va_arg(ap, int);
		vtv[0]  = (void*) va_arg(ap, double*);

		if (tagName != NULL) DrlFPrintf(fp, "%s\n", tagName);
		DrlFPrintf(fp, "%d\n", nItems);
		for (i=0; i<=nItems-1; i++) {
			DrlFPrintf(fp, "%g\t", ((double*)vtv[0])[i]);
		}
		DrlFPrintf(fp, "\n");
		break;

	    /*
	     * Print matrix of doubles
	     */
	    case DRL_CDMATRIX_T:
		tagName = (char*) va_arg(ap, char*);
		nx      = (int) va_arg(ap, int);
		ny      = (int) va_arg(ap, int);
		vtv[0]  = (void*) va_arg(ap, double**);

		if (tagName != NULL) DrlFPrintf(fp, "%s\n", tagName);
		DrlFPrintf(fp, "%d\n%d\n", nx, ny);
		for (i=0; i<=nx-1; i++) {
		    for (j=0; j<=ny-1; j++) {
			DrlFPrintf(fp, "%8.4f ", ((double**) vtv[0])[i][j]);
		    }
		    DrlFPrintf(fp, "\n");
		}

		break;


	    /*
	     * Wrong format
	     */
	    default:
		DrlErrMsg("%s: bad var type\n", routine);
		return(1);
	    }
	}
	va_end(ap);

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f---------------------------------------------------------------------
 * I/O: print TeX output.
 *
 * <br><br>
 * Performs a LaTeX-formatted print of variables and arrays of variables
 * on a file pointer <i> fp</i>.
 * The format and variable arguments are similar to <i> VarStructPrint</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlTeXStructPrint(
	FILE *fp,
	/*
	 * DRL_CARRAY_T, int nItems, int nVar,
	 *   char *varName1, char *tabSep1, int varType1,
	 * 				char *varFmt1, void* var1,
	 *   ...
	 *   char *varNameN, char *tabSepN, int varTypeN,
	 * 				char *varFmtN, void* varN,
	 * DRL_CMATRIX_T,
	 *   int nx, int varTypeX, char *varFmtX, void* varX,
	 *   int ny, int varTypeY, char *varFmtY, void* varY,
	 *   int varTypeXY, char *varFmtXY, void* varXY,
	 * DRL_CDMATRIX_T,
	 *   int nx, int ny, double *vx, double *vy, double **v,
	 * 0)
	 */ ...)
{
static	char		routine[] = "DrlTeXStructPrint";
	int		status = FAILURE;
	va_list		ap;
	int		i, j, m,
			nItems, nx, ny,
			nVar;
	DVType		varType[MAXARG],
			strType;
	char		*varName[MAXARG],
			*varFmt[MAXARG],
			*tabSep[MAXARG],
			buf[128],
			nextLine[256];
	void		*vtv[MAXARG],
			*vx, *vy;

#undef	ERROR
#define	ERROR(str)	{DrlErrMsg("%s: %s (name `%s', line %d)\n",\
			 routine, (str), (tagName != NULL ? tagName : "N/A"),\
			*line); goto done;}

	/* scan arguments */
	va_start(ap, fp);

	while ((strType = (DVType) va_arg(ap, DVType)) != 0) {
	    switch (strType) {
	    case DRL_CARRAY_T:
		/* get arguments */
		nItems = (int) va_arg(ap, int);
		nVar = (int) va_arg(ap, int);
		for (i=0; i <= nVar-1; i++) {
		    varName[i] = (char*)    va_arg(ap, char*);
		    tabSep[i]  = (char*)    va_arg(ap, char*);
		    varType[i] = (DVType) va_arg(ap, DVType);
		    varFmt[i]  = (char*)    va_arg(ap, char*);
		    vtv[i]     = (void*)    va_arg(ap, void*);
		    if ((DrlVTypeCheckValidCType(varType[i]) != 0) &&
		        (DrlVTypeCheckValidLType(varType[i]) != 0)) {
			DrlErrMsg("%s: array argument # %d has invalid "
				"type.\n", routine, i);
			goto done;
		    }
		}

		/* print the array */

		DrlFPrintf(fp, "%% LaTex table \n");
		DrlFPrintf(fp, "\\begin{center}\\begin{tabular}{");
		for (i=0; i <= nVar-1; i++)
			DrlFPrintf(fp, tabSep[i]);
		DrlFPrintf(fp, "}\n\\hline\n");


		for (i=0; i <= nVar-1; i++) {
		    if (i != 0) DrlFPrintf(fp, " & ");
		    DrlFPrintf(fp, "%s", varName[i]);
		}
		DrlFPrintf(fp, "\\\\\n\\hline\\hline\n");



		DrlFPrintf(fp, "%\n");
		for(m=0; m<=nItems-1;m++) {
		    strcpy(nextLine, "");

		    for (i=0; i <= nVar-1; i++) {
		        if (i != 0) strcat(nextLine, " & ");
			DrlVTypePrintFmt(buf,
				varType[i],
				DrlVTypeOffsetVect(vtv[i], m, varType[i]),
				varFmt[i]);
			strcat(nextLine, buf);
		    }
		    strcat(nextLine, "\\\\\n");
		    DrlFPrintf(fp, "%s", nextLine);
		}

		DrlFPrintf(fp, "\\hline\n\\end{tabular}"
			"\\end{center}\n\\vspace{3mm}\n");


		break;

	    /*
	     * Print matrix of Variable types
	     */
	    case DRL_CMATRIX_T:
		nx         = (int)      va_arg(ap, int);
		varType[1] = (DVType) va_arg(ap, DVType);
		varFmt[1]  = (char*)    va_arg(ap, char*);
		vtv[1]     = (void*)    va_arg(ap, void*);
		ny         = (int)      va_arg(ap, int);
		varType[2] = (DVType) va_arg(ap, DVType);
		varFmt[2]  = (char*)    va_arg(ap, char*);
		vtv[2]     = (void*)    va_arg(ap, void*);
		varType[0] = (DVType) va_arg(ap, DVType);
		varFmt[0]  = (char*)    va_arg(ap, char*);
		vtv[0]     = (void*)    va_arg(ap, void*);


		DrlFPrintf(fp, "%% LaTex table \n");
		DrlFPrintf(fp, "\\begin{center}"
			"\\begin{tabular}{|l||");
		for (j=0; j<=ny-1; j++)
			DrlFPrintf(fp, "c|");
		DrlFPrintf(fp, "}\n\\hline\n");


		DrlFPrintf(fp, " ");
		for (j=0; j<=ny-1; j++) {
			DrlVTypePrintFmt(buf,
				varType[2],
				DrlVTypeOffsetVect(vtv[2], j, varType[2]),
				varFmt[2]);
			DrlFPrintf(fp, "& %s", buf);
		}
		DrlFPrintf(fp, "\\\\\n\\hline\\hline\n");

		for (i=0; i<=nx-1; i++) {
		    DrlVTypePrintFmt(buf,
				varType[1],
				DrlVTypeOffsetVect(vtv[1], i, varType[1]),
				varFmt[1]);
		    DrlFPrintf(fp, "%s ", buf);

		    for (j=0; j<=ny-1; j++) {
			DrlVTypePrintFmt(buf,
				varType[0],
				DrlVTypeOffsetMatrix(vtv[0], i, j, varType[0]),
				varFmt[0]);
			DrlFPrintf(fp, "& %s ", buf);
		    }
		    DrlFPrintf(fp, "\\\\\n");
		}

		DrlFPrintf(fp, "\\hline\n\\end{tabular}"
			"\\end{center}\n\\vspace{3mm}\n");

		break;


	    /*
	     * Print matrix of doubles
	     */
	    case DRL_CDMATRIX_T:
		nx      = (int) va_arg(ap, int);
		ny      = (int) va_arg(ap, int);
		vx      = (void*) va_arg(ap, double*);
		vy      = (void*) va_arg(ap, double*);
		vtv[0]  = (void*) va_arg(ap, double**);


		DrlFPrintf(fp, "%% LaTex table \n");
		DrlFPrintf(fp, "\\begin{center}"
			"\\begin{tabular}{|l||");
		for (j=0; j<=ny-1; j++)
			DrlFPrintf(fp, "c|");
		DrlFPrintf(fp, "}\n\\hline\n");


		DrlFPrintf(fp, " ");
		for (j=0; j<=ny-1; j++) {
			DrlFPrintf(fp, "& %8.4f ", ((double*) vy)[j]);
		}
		DrlFPrintf(fp, "\\\\\n\\hline\\hline\n");

		for (i=0; i<=nx-1; i++) {
		    DrlFPrintf(fp, "%.4f ", ((double*) vx)[i]);
		    for (j=0; j<=ny-1; j++) {
			DrlFPrintf(fp, "& %8.4f ", ((double**) vtv[0])[i][j]);
		    }
		    DrlFPrintf(fp, "\\\\\n");
		}

		DrlFPrintf(fp, "\\hline\n\\end{tabular}"
			"\\end{center}\n\\vspace{3mm}\n");

		break;


	    /*
	     * Wrong format
	     */
	    default:
		DrlErrMsg("%s: bad format.\n", routine);
		goto done;
	    }
	}
	va_end(ap);

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}



/*f-------------------------------------------------------------
 * I/O: read double vector (skip comments).
 *
 * <br><br>
 * Scans a vector of double on a file pointer <i> fp</i>.
 * The recognized format is 
 * \begin{verbatim}
 * <number of elements>
 * <element #1>
 * ...
 * <element #N>
 * \end{verbatim}
 * If the number of elements read is 0, does not allocate {tt x}
 * but returns without error.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanDoubleVect(	
	FILE *fp,		/* (I) file */
	double **x,		/* (O) allocated vector */
	int *n)			/* (O) number of elements */
{
static	char	routine[] = "DrlFScanDoubleVect";
	int	status = FAILURE;
	int	i;

	if (DrlFScanInt(fp, n) != SUCCESS) goto done;
	if (*n <= 0) {
		*x = NULL;
		status = SUCCESS;
		goto done;
	}

	if ((*x = DrlDoubleVectAlloc(0, *n-1)) == NULL) goto done;

	for (i=0; i<=*n-1; i++)
		if (DrlFScanDouble(fp, &((*x)[i])) != SUCCESS) goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f-------------------------------------------------------------
 * I/O: print double vector.
 *
 * <br><br>
 * Prints a vector of double on a file pointer <i> fp</i>
 * (prints on the error log in <i> fp</i> is NULL).
 * Same format as in <i> DrlFScanDoubleVect</i>.
 */

int
DrlFPrintDoubleVect(
	FILE *fp,		/* (O) file */
	char *fmt,		/* (I) print format (e.g. "%lf") (or NULL) */
	double *x,		/* (I) vector [0..n-1] */
	int n)			/* (I) nember of elements */
{
	int	i;

	fmt = (fmt == NULL ? "%g\t" : fmt);

	DrlFPrintf(fp, "%d\n", n);
	for (i=0; i<=n-1; i++) DrlFPrintf(fp, fmt, x[i]);
	DrlFPrintf(fp, "\n");

	return(SUCCESS);
}

/*f-------------------------------------------------------------
 * I/O: read double matrix (skip comments).
 *
 * <br><br>
 * Scans a mnatrix of double on a file pointer "fp".
 * The recognized format is 
 * \begin{verbatim}
 * <number of rows>
 * <number of column>
 * <element (1,1)>     ... <element (1,NCOLUMNS)>
 * <element (2,1)>     ...
 * ...
 * <element (NROWS,1)> ... <element (NROWS,NCOLUMNS)>
 * \end{verbatim}
 * If either the number of rows or columns read is 0,
 * does not allocate {tt x} but returns without error.
 * Returns SUCCESS/FAILURE.
 */

int
DrlFScanDoubleMatr(
	FILE *fp,		/* (I) file */
	double ***x,		/* (O) allocated matrix [0..n1-1][0..n2-1] */
	int *n1,		/* (O) number of rows */
	int *n2)		/* (O) number of columns */
{
static	char	routine[] = "DrlFScanDoubleMatr";
	int	status = FAILURE;
	int	i, j;

	if (DrlFScanInt(fp, n1) != SUCCESS) goto done;
	if (DrlFScanInt(fp, n2) != SUCCESS) goto done;
	if ((*n1 <= 0) || (*n2 <= 0)) {
		*x = NULL;
		status = SUCCESS;
		goto done;
	}


	if ((*x = DrlDoubleMatrAlloc(0, *n1-1, 0, *n2-1)) == NULL)
		goto done;

	for (i=0; i<=*n1-1; i++) {
		for (j=0; j<=*n2-1; j++)
			if (DrlFScanDouble(fp, ((*x)[i])+j) != SUCCESS)
				goto done;
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f-------------------------------------------------------------
 * I/O: print double matrix.
 *
 * <br><br>
 * Prints a matrix of double on a file pointer <i> fp</i>
 * (prints on the error log in <i> fp</i> is NULL).
 * Same format as in <i> DrlFScanDoubleMatr</i>.
 */

int
DrlFPrintDoubleMatr(
	FILE *fp,		/* (I) file */
	char *fmt,		/* (I) print format (e.g. "%lf") (or NULL) */
	double **x,		/* (I) matrix [0..n1-1][0..n2-1] */
	int n1,			/* (I) number of rows */
	int n2)			/* (I) number of columns */
{
	int	i, j;

	fmt = (fmt == NULL ? "%g\t" : fmt);
	DrlFPrintf(fp, "%d\n%d\n", n1, n2);
	for (i=0; i<=n1-1; i++) {
		for (j=0; j<=n2-1; j++) DrlFPrintf(fp, fmt, x[i][j]);
		DrlFPrintf(fp, "\n");
	}
	return(SUCCESS);
}


