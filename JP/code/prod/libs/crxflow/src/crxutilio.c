/****************************************************************************/
/*      Standard input output for yield and volatility curves.              */
/****************************************************************************/
/*      STDINPUT.c                                                          */
/****************************************************************************/


/*
$Header$
*/

#include <ctype.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include <alib/macros.h>
#include <alib/ldate.h>
#include <alib/date_sup.h>		/* C analytics */
#include <alib/convert.h>
#include <alib/strutil.h>

#include "crxutilio.h"

#undef	INSERTCHAR
#define	INSERTCHAR(s,p,c)	{char *q; for(q=s+(int)strlen(s);q>=p;q--)\
				 *(q+1)=*q; *p = c;}
#undef	DELETECHAR
#define	DELETECHAR(s,p)		{char *q; for(q=p; *(q+1)!= '\0'; q++) \
				*q=*(q+1); *q = '\0';}

/* buffers for printing */
#define	BUF_IDX	32
static	char	tmpBuf[BUF_IDX][64];
static	int	tmpBufIdx=0;
#define	BUF_NEXT(str)	str = (str != NULL ? str : tmpBuf[tmpBufIdx++]);\
			tmpBufIdx = (tmpBufIdx > BUF_IDX-1 ? 0 : tmpBufIdx)

/*f----------------------------------------------------------------------
 * Allocates a vector of type "type" of length "size".
 */
void* CrxVTypeVectAlloc(int size, TVType type)
{
static	char	routine[] = "CrxVTypeVectAlloc";
	int	i;
	void	*p = NULL;
	size_t	typeSize;

	switch (type) {
		/**
		 ** LIL counted arrays
		 **/
	case DRL_DOUBLE_L:
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		p = (void*) NEW_ARRAY(double, size+1);
		if (p) ((double*)p)[0] = (double) size;
		return(p);
	case DRL_LONG_L:
		p = (void*) NEW_ARRAY(long, size+1);
		if (p) ((long*)p)[0] = (long) size;
		return(p);
	case DRL_TDATE_L:
		p = (void*) NEW_ARRAY(TDate, size+1);
		if (p) ((TDate*)p)[0] = (TDate)size;
		return(p);

	case DRL_CHAR_BLOCK_L:
		p = (void*) NEW_ARRAY(char, 128*size+1);
		if (p) ((char*)p)[0] = (char)((unsigned char) size);
		return(p);

		/**
		 ** Other
		 **/
	case DRL_STRING_T:
		/* for DRL_STRING_T, allocate char string */
		p = (void*) NEW_ARRAY(char*, size);
		if (p == NULL) {
			DR_Error("%s: malloc DRL_STRING_T "
				"length %d failed.\n", routine, (int) (size));
			return(NULL);
		}
	    	for (i=0; i<=size-1; i++) {
			((char**)p)[i] = NEW_ARRAY(char,  WRAP_STR_BYTES);
			if (((char**)p)[i] == NULL) {
				DR_Error("%s: malloc failed.\n", routine);
				return(NULL);
			}
		}
		break;


	default:
		/* check type size */
		if (CrxVTypeCheckValidCType(type) != 0) {
			DR_Error("%s: bad C type.\n", routine);
			return(NULL);
		}

	 	/* get type size */
		if ((typeSize = CrxVTypeSizeof(type)) <= 0) {
			DR_Error("%s: bad type size.\n");
			return(NULL);
		}

		/* allocate memory */
		if ((p = (void*) MALLOC((size_t) (typeSize*size))) == NULL) {
			DR_Error("%s: malloc length %d failed.\n", routine, 
					(int) (typeSize*size));
			return(NULL);
		}
		break;
	}

	return(p);
}


/*f----------------------------------------------------------------------
 * Frees an array of type "type" of length "size".
 */

int CrxVTypeVectFree(void* p, int size, TVType type)
{
	int	i;

	/* nothing to do */
	if (p == NULL) return(0);

	switch (type) {
		/**
		 ** LIL counted arrays
		 **/
	case DRL_DOUBLE_L:
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		FREE((double*) p);
		return(SUCCESS);
	case DRL_LONG_L:
		FREE((long*) p);
		return(SUCCESS);
	case DRL_TDATE_L:
		FREE((TDate*) p);
		return(SUCCESS);
	case DRL_CHAR_BLOCK_L:
		FREE((char*) p);
		return(SUCCESS);

		/**
		 ** Other
		 **/
	case DRL_STRING_T:
	    	for (i=0; i<=size-1; i++) {
			if (((char**)p)[i] != NULL)
				FREE((void*) ((char**)p)[i]);
	    	}
		FREE((void*) p);
		break;
	default:
		if (CrxVTypeCheckValidCType(type) != 0) {
			DR_Error("CrxVTypeVectFree: bad C type\n");
			return(1);
		}

		FREE((void*) p);

		return(0);
	}
	return(SUCCESS);
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

int CrxFScanString(FILE *fp, char *s)
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

int CrxFScanVType(FILE *fp, TVType type, void *valPtr)
{
static	char	buf[1024];
	if (CrxFScanString(fp, buf) != SUCCESS)
		return(FAILURE);
	/*if (CrxFGetToken(fp, NULL, buf, sizeof(buf)) != SUCCESS)
		return(FAILURE);*/

	if (CrxVTypeScan(buf, type, valPtr) != SUCCESS)
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

int CrxFPrintVType(FILE *fp, TVType type, void *valPtr)
{
static	char	buf[256];

	strcpy(buf, CrxVTypePrint(NULL, type, valPtr));

	if (fp == NULL) {
		DR_Error(buf);
	} else {
		fputs(buf, fp);
	}

	return(SUCCESS);
}

/*f----------------------------------------------------------------------
 * Reads a file <i> fp</i> an array of <i> numVect</i> columns
 * each containing <i> numItems</i> elements,
 * of the form
 * \begin{verbatim}
 *   <elem #1        of vect #1>  ...  <elem #1        of vect #numVect> 
 *   ...
 *   <elem #numItems of vect #1>  ...  <elem #numItems of vect #numVect> 
 *   end
 * \end{verbatim}
 * and allocates LIL counted arrays conatainig the elements.
 */

int CrxLilVectArrayFpRead(
	FILE *fp,                             /* (I) file pointer */
	int numItems,                         /* (I) num elements in each column */
	int numVect,                          /* (I) num of columns */
	TVType *varType,                      /* (I) array of variable types [0..numVect-1] */
	void ***lptr)                         /* (I) array of void* pointers [0..numVect-1] */
{
static	char	routine[] = "CrxLilVectArrayFpRead";
	int	status = FAILURE;

	int	idxV, idx;
	void	*vptr;


	/* reset to NULL */
	for (idxV=0; idxV<=numVect-1; idxV++)
		*(lptr[idxV]) = (void*) NULL;

	/* malloc */
	for (idxV=0; idxV<=numVect-1; idxV++) {
		*(lptr[idxV]) = CrxVTypeVectAlloc(numItems, varType[idxV]);
		if (*(lptr[idxV]) == NULL) goto done;
	}



	/* read elements */
	for (idx=0; idx<=numItems-1; idx++) {
	    for (idxV=0; idxV<=numVect-1; idxV++) {

		/* access element #idx of array #idxV */
		vptr = CrxVTypeOffsetVect(*(lptr[idxV]), idx, varType[idxV]);
		if (vptr == NULL) {
			goto done;
		}

		/* read element */
		if (CrxFScanVType(fp, varType[idxV], vptr) != SUCCESS) {
		    DR_Error("%s: can't read element #%d/%d "
			"in column # %d/%d "
			"(expecting type %s).\n", routine,
			idx+1, numItems, idxV+1, numVect,
			CrxVTypeName(varType[idxV]));
		    goto done;
		}
	    }
	}


	/* made it through */
	status = SUCCESS;
done:
	if (status != 0) {
		/* free mem */
		for (idxV=0; idxV<=numVect-1; idxV++)
		    CrxVTypeVectFree(*(lptr[idxV]), numItems, varType[idxV]);
		DR_Error("%s: failed.\n", routine);
	}
	return(status);
}


/*f----------------------------------------------------------------------
 * A version of <i> CrxLilVectArrayFpRead</i> that allows the call
 * with variable number of arguments (one does not need to
 * set up arrays as in <i> CrxLilVectArrayFpRead</i>. \\
 * <b> Example:</b> to read 
 * \begin{verbatim}
 * # dates and notional
 *    01-Jan-1998    100.0
 *    01-Jul-1998     80.0
 *    01-Jan-1999     60.0
 *    01-Jul-1999     40.0
 *    01-Jan-1999     20.0
 * \end{verbatim}
 * call the routine as in 
 * \begin{verbatim}
 *     TDate	*dates = NULL;
 *     double	*notionals = NULL;
 *     ...
 *     status = CrxLilVectArrayFpReadV(
 *                    stdin,
 *                    5,
 *                    DRL_TDATE_L,  (void*) &dates,
 *                    DRL_DOUBLE_L, (void*) &notionals,
 *                    DRL_NULL_T);
 *     ...
 * \end{verbatim}
 */

int CrxLilVectArrayFpReadV(
	FILE *fp,                             /* (I) file pointer */
	int numItems,                         /* (I) num elements to be read */
	/* TVType varType, void *lptr,
	 * ...
	 * TVType varType, void *lptr,
	 * DRL_NULL_T (last argument MUST be DRL_NULL_T)
	 */
	...)
{
static	char	routine[] = "CrxLilVectArrayFpReadV";
	int	status = FAILURE;
	va_list		ap;
#undef	MAX_ARGS
#define	MAX_ARGS	32
	TVType		varType[MAX_ARGS];
	void		**lptr[MAX_ARGS];
	int		numVect = 0;


	/* get arguments */
	va_start(ap, numItems);
	while ((varType[numVect] = (TVType) va_arg(ap, TVType))
			!= DRL_NULL_T) {
		lptr[numVect] = (void**) va_arg(ap, void*);
		numVect++;
		if (numVect >= MAX_ARGS) {
			DR_Error("%s: too many arguments (max %d).\n",
				routine, MAX_ARGS);
			goto done;
		}
	}
	va_end(ap);

	/* call routine */
	if (CrxLilVectArrayFpRead(
		fp,
		numItems,
		numVect,
		varType,
		lptr) != SUCCESS)
			goto done;

	/* made it through */
	status = SUCCESS;
done:
	if (status != 0) {
		DR_Error("%s: failed.\n", routine);
	}
	return(status);
#undef	MAX_ARGS
}

/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */
char* CrxVTypeName(TVType type)
{
	switch (type) {
	case DRL_DOUBLE_T:		return("double");
	case DRL_FLOAT_T:		return("float");
	case DRL_INT_T:			return("int");
	case DRL_LONG_T:		return("long");
	case DRL_CHAR_T:		return("char");

	case DRL_STRING_T:		return("char*");
	case DRL_CHAR_ARRAY_T:		return("char[]");
	case DRL_TDATE_T:		return("TDate");
	case DRL_TDATEINTERVAL_T:	return("TDateInterval");
	case DRL_TDAYCOUNT_T:		return("TDayCount");

	case DRL_PERCENT_T:		return("percent");
	case DRL_CUR_T:			return("currency");
	case DRL_CURK_T:		return("currencyK");
	case DRL_CURM_T:		return("currencyM");
	case DRL_BOOLEAN_T:		return("TBoolean");

	case DRL_FLOAT_L:		return("LIL_float");
	case DRL_LONG_L:		return("LIL_long");
	case DRL_TDATE_L:		return("LIL_date");
	case DRL_TDATEINTERVAL_L:	return("LIL_date_interval");
	case DRL_CHAR_BLOCK_L:		return("LIL_char_block");
	case DRL_CHAR_L:		return("LIL_char");

	case DRL_PERCENT_L:		return("LIL_percent");
	case DRL_BPOINT_L:		return("LIL_bpoint");


	default:		return("ERROR");
	}
}



/*f----------------------------------------------------------------------
 * Returns SUCCESS if type "type" is a valid C type.
 */
int CrxVTypeCheckValidCType(TVType type)
{
	switch (type) {
	/* C types */
	case DRL_POINTER_T:
	case DRL_DOUBLE_T:
	case DRL_FLOAT_T:
	case DRL_INT_T:
	case DRL_LONG_T:
	case DRL_CHAR_T:
	case DRL_STRING_T:
	case DRL_TDATE_T:
	case DRL_TDATEINTERVAL_T:
	case DRL_TDAYCOUNT_T:

	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
	case DRL_BOOLEAN_T:
		/**
		 ** Complex data structures
		 **/

		return(0);
	default:
		return(1);
	}
}

/*f----------------------------------------------------------------------
 * Returns 0 if type <i> type</i> is a valid LIL type.
 */
int CrxVTypeCheckValidLType(TVType type)
{
	switch (type) {
	case DRL_FLOAT_L:
	case DRL_INT_L:
	case DRL_TDATE_L:
	case DRL_TDATEINTERVAL_L:
	case DRL_CHAR_BLOCK_L:
	case DRL_CHAR_L:

	case DRL_PERCENT_L:
	case DRL_BPOINT_L:

		return(0);
	default:
		return(1);
	}
}



/*f----------------------------------------------------------------------
 * Returns the size of <i> type</i>.
 */
size_t CrxVTypeSizeof(TVType type)
{
static	char	routine[] = "CrxVTypeSizeof";
	switch (type) {
	/* C types */
	case DRL_DOUBLE_T:
	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
		return sizeof(double);
	case DRL_FLOAT_T:
		return sizeof(float);
	case DRL_INT_T:
	case DRL_BOOLEAN_T:
		return sizeof(int);
	case DRL_LONG_T:
		return sizeof(long);
	case DRL_CHAR_T:
		return sizeof(char);
	case DRL_STRING_T:
		return sizeof(char*);
	case DRL_TDATE_T:
		return sizeof(TDate);
	case DRL_TDATEINTERVAL_T:
		return sizeof(TDateInterval);
	case DRL_TDAYCOUNT_T:
		return sizeof(TDayCount);

	/* LIL types */
	case DRL_FLOAT_L:
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		return sizeof(double);
	case DRL_LONG_L:
		return sizeof(long);
	case DRL_CHAR_BLOCK_L:
		return sizeof(char*);
	case DRL_CHAR_L:
		return sizeof(char);
	case DRL_TDATEINTERVAL_L:
		return sizeof(TDateInterval);
	case DRL_TDATE_L:
		return sizeof(long);
	default:
		DR_Error("%s: bad type %ld.\n", routine, (long) type);
		return((size_t) 0);
	}
}

/*f----------------------------------------------------------------------
 * Returns the element with offset "offset" in an array
 * "p" of type "type".
 */

void* CrxVTypeOffsetVect(void *p, int offset, TVType type)
{

	switch (type) {
	case DRL_POINTER_T:
		return(((void**)p) + offset);
	case DRL_DOUBLE_T:
	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
		return(((double*)p) + offset);
	case DRL_FLOAT_T:
		return(((float*)p) + offset);
	case DRL_INT_T:
	case DRL_BOOLEAN_T:
		return(((int*)p) + offset);
	case DRL_LONG_T:
		return(((long*)p) + offset);
	case DRL_CHAR_T:
		return(((char*)p) + offset);
	case DRL_STRING_T:
		return(((char**)p) + offset);
	case DRL_TDATE_T:
		return(((TDate*)p) + offset);
	case DRL_TDATEINTERVAL_T:
		return(((TDateInterval*)p) + offset);
	case DRL_TDAYCOUNT_T:
		return(((TDayCount*)p) + offset);

		/*
		 * LIL (counted) arrays
		 */

	case DRL_DOUBLE_L:
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		return(((double*)p) + offset +1);
	case DRL_LONG_L:
		return(((long*)p) + offset +1);
	case DRL_CHAR_L:
		return(((char*)p) + offset +1);
	case DRL_TDATE_L:
		return(((long*)p) + offset +1);

	case DRL_CHAR_BLOCK_L:
		/*return(((char*)p) + WRAP_STR_IDX(offset+1));*/
		return &((char*) p)[GTO_STR_POS(offset+1)];

		/**
		 ** Complex data structures
		 **/


	default:
		DR_Error("OffsetType: bad type %ld.\n", (long)type);
		return(NULL);
	}
}


/*f-------------------------------------------------------------
 * Scans a double value in a currency format
 * $\pm XXX,XXX,XXX.XX$ in the string <i> s</i>
 * and puts the result in <i> x</i>. 
 * If the string containts the percent sign, the output is 
 * divided by 100.
 * Returns 0 on success.
 */

int CrxCurScan(char *s, double *x)
{
	char	buf[256],
		*p;
	int	isPercent = FALSE;

	strcpy(buf, s);
	if (strchr(buf, '%')) isPercent = TRUE;

	for(p = buf; *(p+1) != 0; p++) {
		if (*p == ',') DELETECHAR(s, p);
	}

	if (sscanf(buf, "%lf", x) != 1) return (FAILURE);

	if (isPercent) *x *= 1e-2;
	return (SUCCESS);
}

/*f-------------------------------------------------------------
 * Prints a double value <i> x</i> in a currency format
 * $\pm XXX,XXX,XXX.XX$ in the string <i> s</i>.
 * If <i> s</i> is NULL, uses a static buffer and returns it.
 */

char* CrxCurPrint(char *s, double x, int ndigits)
{
static	char	fmt[32];
	char	*p ;
	int	i;

	/*s = (s == NULL ? buf : s);*/
	BUF_NEXT(s);

	if (ndigits <= 0) {
		sprintf(s, "%.0f", x);

		/* goto point */
		p = s + strlen(s);

	} else {
		sprintf(fmt, "%%.%df", ndigits);
		sprintf(s, fmt, x);

		/* goto point */
		p = s;
		while (*p != '.') p++;
	}


	if (p == s) goto done;

	for(;;) {
	    for (i=0; i<=2; i++) {
		p--;
		if ((p == s) || !(isdigit(*p))) goto done;
	    }
	    if (!isdigit(*(p-1))) goto done;

	    if (isdigit(*(p-1))) INSERTCHAR(s, p, ',');

		/*if (((p-1) == s) || !(isdigit(*(p-1)))) goto done;
		INSERTCHAR(s, p, ',');*/
	}


done:
	return(s);
}



/*f-------------------------------------------------------------
 * Date routines : scan date interval in char string.
 *                                                          
 * <br><br>
 * Scans for a date interval in the string "s". If scan is
 * successful, returns 0 and puts the result in "tlapse".
 * If failure, returns a non-zero value.
 * Recognized formats are of the form ``MXNY..'' where
 * M, N are number of periods (non-negative) and X, Y
 * a character defining a period (D,W,M,Q,etc\dots).
 */

int CrxTDateIntervalScan(char *s, TDateInterval *tlaps)
{
    static	char	routine[]="CrxTDateIntervalScan";
    char		*p,  buf[256];
    int		nPer;
    
    /* Special case */
    if ((!strcmp(s, "O/N")) ||
        (!strcmp(s, "ON"))  ||
        (!strcmp(s, "SN"))  ||
        (!strcmp(s, "S/N"))) {
        if (GtoMakeDateInterval(1, 'D', tlaps) != SUCCESS)
            goto done;
        return(SUCCESS);
    }
    
    /* Scan for IMMn intervals */
    if ((p = strstr(s, "IMM")) != NULL) {
        p +=3;
        if (sscanf(p, "%d", &nPer) != 1)
            goto done;
        if (GtoMakeDateInterval(nPer, 'I', tlaps) != SUCCESS)
            goto done;
        return(SUCCESS);
    }
    
    /* replace '3y' by '3A' for GtoString... */
    strcpy(buf, s);
    for (p = buf; *p != '\0'; p++) {
        if (toupper(*p) == 'Y') *p = 'A';
    }
    
    return (GtoStringToDateInterval(buf, routine, tlaps));
    
    
done:
    DR_Error("CrxTDateIntervalScan: can't scan `%s'.\n", s);
    return(FAILURE);
}


/*f-------------------------------------------------------------
 * This routine scans for strings
 * in a character string pointer by <i> p</i>.
 * The read strings can be enclosed in quotes (') or double quotes ("),
 * which are removed.
 * On exit, the result string is put in <i> arg</i> and returned (if <i> arg</i> 
 * is NULL, a static copy is returned).
 * The pointer <i> p</i> is advanced to the remaining unread portion 
 * the scanned string. It should not be changed between
 * successive calls to the routine.
 * For example,
 * \begin{verbatim}
 *          char    *p, *q, *string = "1.23 'hello world'";
 *          ...
 *          / * p points stringto be scanned * /
 *          p = &string[0];
 *          q = CrxStrScanString(&p, NULL);	/ * q is "1.23" * /
 *          / * p now points to the unread portion of string  * /
 *          ...
 *          q = CrxStrScanString(&p, NULL);	/ * q is "hello world" * /
 *          ...
 * \end{verbatim}
 */

char* CrxStrScanString(char **p, char *arg)
{
const	char	EOS = '\0';		/* end of string */
static	char	arg1[256];
	char	*s, *q, c;

	arg = (arg != NULL ? arg : arg1);
	q = arg;
	s = *p;

	while (((c = *s++) != EOS) && (isspace(c)));

	if (c == EOS) return(NULL);


	if (c == '"') {
	        while (((c = *s++) != EOS) && (c != '"'))
			*q++ = c;
		if (c == EOS) return(NULL);
		*q = '\0';
	} else if (c == '\'') {
	        while (((c = *s++) != EOS) && (c != '\''))
			*q++ = c;
		if (c == EOS) return(NULL);
		*q = '\0';
	} else {
		*q++ = c;
	        while (((c = *s++) != EOS) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	}
	*p = s;
	return(arg);
}

/*f-------------------------------------------------------------
 * Date routines : scan date from string.
 *                                                          
 * <br><br>
 * Scan the string <i> string</i> for a date. Recognized formats
 * are <i> mm/dd/yy</i>, <i> mm/dd/yyyy</i>, or <i> dd-MMM-yyyy</i>
 * or <i> yyyy.mm.dd</i>.
 * Puts the result in <i> aDate</i>. Returns 0 if success.
 */

int CrxTDateScan(char *string, TDate *aDate)
{
static	char	routine[] = "CrxTDateScan";
	char	dateCopy[32];
	int	yyyy, mm, dd;


	/* special case of yyyy.mm.dd not handled by GTO */
	if (strchr(string, '.') != NULL) {
		if (sscanf(string, "%d.%d.%d", &yyyy, &mm, &dd) != 3) {
			DR_Error("%s: can't read date `%s'.\n", 
				  routine, string);
			return(FAILURE);
		}
		sprintf(dateCopy, "%04d%02d%02d", yyyy, mm, dd);
	} else {
		strncpy(dateCopy, string, sizeof(dateCopy)-1);
	}

    /* allow zero dates - specified as 0 or nil */
    if (strcmp (string, "0") == 0 || strcmp (string, "nil") == 0)
    {
        *aDate = 0;
        return SUCCESS;
    }

	/* Use Gto routine */
	return GtoStringToDate(dateCopy, aDate);


#ifdef	_SKIP	/* Do not use */
static	char	routine[] = "CrxTDateScan";

	int	yyyy, mm, dd ;
	char	s0[255], *s1, *s2 ;
static	char	toks[] = "-/ ,";

#define	FMT_EXTENDED
#ifdef	FMT_EXTENDED

	if ((s1 = CrxStrToken(string, toks, s0, &s2)) == NULL) return(__LINE__);

	if (sscanf(s1, "%d", &mm) == 1) {
		/*
		 * scan for `dd-Mon-yyyy' or `mm-dd-yyyy'
		 */
		if ((s1 = CrxStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &dd) == 1) {
			/*
			 * format `mm-dd-yyyy'
			 */
			if ((s1 = CrxStrToken(NULL, toks, s0, &s2)) == NULL)
				return(__LINE__);
			if (sscanf(s1, "%d", &yyyy) != 1)
				return(__LINE__);
		} else {
			/*
			 * format `dd-Mon-yyyy'
			 */
			dd = mm ;
			mm = 0;
			while ((++mm <= 12)
				&& (strncmp(s1, monthName[mm], 3) != 0));
			if (mm >= 13) return(__LINE__);

			if ((s1 = CrxStrToken(NULL, toks, s0, &s2)) == NULL)
				return(__LINE__);
			if (sscanf(s1, "%d", &yyyy) != 1)
				return(__LINE__);
		}
	} else {
		/*
		 * format `Mon dd, yyy'
		 */
		mm = 0;
		while ((++mm <= 12) && (strncmp(s1, monthName[mm], 3) != 0));
		if (mm >= 13) return(__LINE__);

		if ((s1 = CrxStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &dd) != 1)
			return(__LINE__);
		if ((s1 = CrxStrToken(NULL, toks, s0, &s2)) == NULL)
			return(__LINE__);
		if (sscanf(s1, "%d", &yyyy) != 1)
			return(__LINE__);
	}

#else	/* FMT_EXTENDED */

	if (sscanf(string, "%d/%d/%d", &mm, &dd, &yyyy) == 3) {
		if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
		if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

		*aDate = CrxTDateMake(mm, dd, yyyy) ;
		return(CrxTDateCheckValid(*aDate)) ;
	}
	if (sscanf(string, "%d-%s-%d", &dd, s0, &yyyy) == 3) {
		if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
 		if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

		mm = 0;
		while ((++mm <= 12) && (strncmp(s0, monthName[mm], 3) != 3));
 		*aDate = CrxTDateMake(mm, dd, yyyy) ;
		return(CrxTDateCheckValid(*aDate)) ;
	}
#endif	/* FMT_EXTENDED */


	if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
	if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

 	*aDate = CrxTDateMake(mm, dd, yyyy) ;
	return(CrxTDateCheckValid(*aDate)) ;
#endif
}

/*f----------------------------------------------------------------------
 * Scans in a formatted string <i> str</i> a type <i> type</i> and
 * puts the result in <i> p</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> CrxVTypeOffsetVect</i>, even if it is the first element.}
 * Returns 0 if scan successful.
 */

int CrxVTypeScan(char *str, TVType type, void *p)
{
static	char	routine[] = "CrxVTypeScan";
	char	*buf2, buf[256];

	switch (type) {
	case DRL_DOUBLE_T:
		/*if (sscanf(str, "%lf", (double*)p) != 1) goto error; */
		if (CrxCurScan(str, (double*)p) != SUCCESS) goto error;
		break;
	case DRL_PERCENT_T:
		if (sscanf(str, "%lf", (double*)p) != 1) goto error;
		*((double*)p) *= 1e-2;
		break;
	case DRL_FLOAT_T:
		if (sscanf(str, "%f",  (float*)p) != 1) goto error;
		break;
	case DRL_INT_T:
		if (sscanf(str, "%d",  (int*)p) != 1) goto error;
		break;
	case DRL_LONG_T:
		if (sscanf(str, "%ld",  (long*)p) != 1) goto error;
		break;
	case DRL_CHAR_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (sscanf(buf, "%c",  (char*)p) != 1) goto error;
		break;
	case DRL_STRING_T:
		if (sscanf(str, "%s",  ((char**)p)[0]) != 1) goto error;
		break;
	case DRL_CHAR_ARRAY_T:
		strncpy(buf, str, sizeof(buf));
		buf2 = &buf[0];
		if (CrxStrScanString(&buf2, ((char*)p)) == NULL)
			goto error;
		/*if (sscanf(str, "%s",  ((char*)p)) != 1) goto error;*/
		break;
	case DRL_TDATE_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		/*if (GtoStringToDate(buf, (TDate*)p) != 0) goto error;*/
		if (CrxTDateScan(buf, (TDate*)p) != 0) goto error;
		break;
	case DRL_TDATEINTERVAL_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (CrxTDateIntervalScan(buf, (TDateInterval*)p)
				!= SUCCESS) goto error;
		break;
	case DRL_TDAYCOUNT_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (GtoStringToDayCountConv(buf, 
			(TDayCount*)p) != SUCCESS) goto error;
		break;

		/**
		 ** Derived types
		 **/
	case DRL_CUR_T:
		if (CrxCurScan(str, (double*)p) != SUCCESS) goto error;
		break;
	case DRL_CURK_T:
		if (CrxCurScan(str, (double*)p) != SUCCESS) goto error;
		*((double*)p) *= 1e3;
		break;
	case DRL_CURM_T:
		if (CrxCurScan(str, (double*)p) != SUCCESS) goto error;
		*((double*)p) *= 1e6;
		break;
	case DRL_BOOLEAN_T:
		/* T,F or integer value */
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (toupper(buf[0]) == 'T') {
			*((int*) p) = 1;
		} else if (toupper(buf[0]) == 'F') {
			*((int*) p) = 0;
		} else {
			if (sscanf(buf, "%d",  (int*)p) != 1)
				goto error;
		}
		break;

		/**
		 ** LIL counted arrays
		 **/
	case DRL_DOUBLE_L:
		if (sscanf(str, "%lf", (double*)p) != 1) goto error;
		break;
	case DRL_LONG_L:
		if (sscanf(str, "%ld",  (long*)p) != 1) goto error;
		break;
	case DRL_TDATE_L:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (GtoStringToDate(buf, (TDate*)p) != 0) goto error;
		break;
	case DRL_CHAR_BLOCK_L:
		if (sscanf(str, "%s",  ((char*)p)) != 1) goto error;
		break;

		/**
		 ** Derived LIL counted arrays
		 **/
	case DRL_PERCENT_L:
		if (sscanf(str, "%lf", (double*)p) != 1) goto error;
		*((double*)p) *= 1e-2;
		break;
	case DRL_BPOINT_L:
		if (sscanf(str, "%lf", (double*)p) != 1) goto error;
		*((double*)p) *= 1e-4;
		break;



		/**
		 ** Complex data structures
		 **/

	default:
		DR_Error("%s: bad type %ld.\n", routine, (long) type);
		return(2);
	}

	return(0);
error:
	DR_Error("%s: can't scan `%s' (type %s)\n",
		routine, str, CrxVTypeName(type));
	return(1);
}


/*f----------------------------------------------------------------------
 * Prints in a formatted string <i> str</i> a pointer <i> </i>p
 * of type <i> type</i>.  Returns <i> str</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> CrxVTypeOffsetVect</i>, even if it is the first element.}
 */

char* CrxVTypePrint(char *str, TVType type, void *p)
{
static	char	routine[] = "CrxVTypePrint";
static	char	buf[64];

	str = (str == (char*)NULL ? buf : str);

	if (p == NULL) {
		sprintf(str, "NULL");
		return(str);
	}

	switch (type) {
	case DRL_POINTER_T:
		sprintf(str, "%p", *((void**) p));
		break;
	case DRL_DOUBLE_T:
		sprintf(str, "%12.8f", *((double*) p));
		break;
	case DRL_FLOAT_T:
		sprintf(str, "%f",  *((float*) p));
		break;
	case DRL_INT_T:
		sprintf(str, "%d",  *((int*) p));
		break;
	case DRL_LONG_T:
		sprintf(str, "%ld",  *((long*) p));
		break;
	case DRL_CHAR_T:
		sprintf(str, "%c",  *((char*) p));
		break;
	case DRL_STRING_T:
		sprintf(str, "%s",  ((char**)p)[0]);
		break;
	case DRL_CHAR_ARRAY_T:
		sprintf(str, "%s",  ((char*)p));
		break;
	case DRL_TDATE_T:
		sprintf(str, "%10s", GtoFormatDate(*((TDate*) p)));
		break;
	case DRL_TDATEINTERVAL_T:
		sprintf(str, "%6s",
			GtoFormatDateInterval(((TDateInterval*) p)));
		break;
	case DRL_TDAYCOUNT_T:
		sprintf(str, "%10s",
			GtoFormatDayCountConv(*((TDayCount*) p)));
		break;

		/* derived Types */

	case DRL_PERCENT_T:
		sprintf(str, "%1.8f", *((double*) p) * 1e2);
		break;
	case DRL_CUR_T:
		CrxCurPrint(str, *((double*) p), 0);
		break;
	case DRL_CURK_T:
		CrxCurPrint(str, *((double*) p)*1e-3, 0);
		break;
	case DRL_CURM_T:
		CrxCurPrint(str, *((double*) p)*1e-6, 0);
		break;
	case DRL_BOOLEAN_T:
		sprintf(str, "%s", (*((int*) p) != 0 ? "T" : "F"));
		break;

		/**
		 ** LIL counted arrays
		 **/
	case DRL_DOUBLE_L:
		sprintf(str, "%12.8f", *((double*) p));
		break;

		/**
		 ** Derived LIL counted arrays
		 **/
	case DRL_PERCENT_L:
		sprintf(str, "%12.8f", *((double*) p) * 1e2);
		break;
	case DRL_BPOINT_L:
		sprintf(str, "%12.8f", *((double*) p) * 1e4);
		break;


		/**
		 ** Complex data structures
		 **/

	default:
		DR_Error("%s: bad type %ld.\n", routine, (long)type);
		strcpy(str, "ERR");
		break;
	}

	return(str);
}



/*f-------------------------------------------------------------
 * Date routines : scan date in YYYYMDD format.
 *
 * <br><br>
 * Scan the string "string" for a date. Recognized format
 * is YYYYMMDD. Puts the
 * result in "aDate". Returns 0 if success.
 */

int
CrxTDateScanYMD(char *string, TDate *aDate)
{
        int     yyyy, mm, dd ;
        long    ymd;

        if (strcmp (string, "0") == 0 ||
            strcmp (string, "nil") == 0)
        {
            *aDate = 0;
            return SUCCESS;
        }

        if (sscanf(string, "%ld", &ymd) != 1) return(-2);

        yyyy = (int) floor(ymd*1e-4);
        mm = (int) floor((ymd-10000*yyyy)*1e-2);
        dd = ymd - 10000*yyyy - 100*mm;

        if ((yyyy >= 50) && (yyyy <= 99)) yyyy += 1900 ;
        if ((yyyy >= 0)  && (yyyy <= 50)) yyyy += 2000 ;

        *aDate = CrxTDateMake(mm, dd, yyyy) ;
        return(SUCCESS) ;
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

static int CrxFGetToken(FILE *fp, char *sep, char *s, int sz)
{
    char    *q;
    int     c;
    static  char    buf[256];

    sep = sep; /* avoid warning on unused parameter */
    sz  = sz;  /* avoid warning on unused parameter */

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




int CrxFAdvanceToToken(FILE *fp, char *token)
{
static  char    buf[64];
        for (;;) {
            if (CrxFGetToken(fp, NULL, buf, sizeof(buf))!= SUCCESS)
                return (FAILURE);
            if (!strcmp(buf, token)) {
                return(SUCCESS);
            }
        }
        return(FAILURE);
}





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

char *
CrxFGetLine(char *s, int n, FILE *fp, int *line)
{
        int     acceptFlag = 0;
        char    *p;

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



/*----------------------------------------------------------------------
 * Internally called.
 */

int
CrxStrLongValueScanV(
        char *s,                /* (I) input string (unchanged) */
        char *errName,          /* (I) variable name for debugging */
        long *value,            /* (O) value */
        va_list ap)             /* (I) arguments list */
{
        int             status = FAILURE;
        char            *varName;
        long            varValue;
        char            errString[1024];

        sprintf(errString, "CrxStrLongValueScanV: can't read %s in `%s' "
                "(possible choices are ", errName, s);

        while ((varName = (char*) va_arg(ap, char*)) != NULL) {
            varValue = (long) va_arg(ap, long);

            /* Add choice to error string just in case */
            strcat(errString, " `");
            strcat(errString, varName);
            strcat(errString, "'");

            /* if (!strncmp(s, varName, strlen(varName))) {*/
            if (!strcmp(s, varName)) {
                *value =  varValue;
                status = SUCCESS;
                goto done;
            }
        }

        strcat(errString, ").\n");
        DR_Error(errString);
        return(FAILURE);
done:
        return(status);
}




/*f---------------------------------------------------------------------
 * Scans a long value in a string <i> s</i> by trying to match
 * the string with different values. On successful return, the
 * value is put in <i> value</i>.
 * <br> <b> Example:</b>
 * \begin{verbatim}
 * char  *buf;
 * long  value;
 * ...
 * errCode = CrxStrLongValueScan(buf, "stub type" &value,
 *              "B", GTO_BOND_STUB,
 *              "S", GTO_SIMPLE_STUB,
 *              NULL);
 * ...
 * \end{verbatim}
 */
int
CrxStrLongValueScan(
        char *s,                /* (I) input string (unchanged) */
        char *errName,          /* (I) variable name for debugging */
        long *value,            /* (O) value */
        /* char *varName1, long varValue1,
         * ...
         * char *varNameN, long varValueN,
         * NULL)                        LAST ARGUMENT MUST BE NULL
         */ ...)
{
static  char            routine[] = "CrxStrLongValueScan";
        int             status = FAILURE;
        va_list         ap;

        va_start(ap, value);
        status = CrxStrLongValueScanV(s, errName, value, ap);
        va_end(ap);

        if (status != SUCCESS) {
            DR_Error("%s: scanning `%s' failed.\n", routine, s);
        }
        return(status);
}






/*f--------------------------------------------------------------
 * Reads a TCurve "that" from a file pointer "fp".
 * The format of the file is the same as for <i> CrxTCurveFileWrite</i>.
 * Returns 0 iff OK.
 */

static int CrxTCurveFpRead(TCurve **that, FILE *fp)
{
    static  char            routine[] = "CrxTCurveFpRead";
    char            buf[255], buf2[255];
    int             i, line = 0, c,
        errCode = FAILURE;
    TDate           baseDate;
    int             nPts;
    double          basis;
    TDayCount       dayCountConv;
    
    *that = NULL;
    
#define CrxTCurveDate(object, index)    ((object)->fArray[(index)].fDate)
#define CrxTCurveRate(object, index)    ((object)->fArray[(index)].fRate)
    
    /*
     * Identify format
     */
    if ((c = fgetc(fp)) == EOF) goto done;
    ungetc(c, fp);
    
    /*
     * (1) London format:
     * # Start date
     * 19970411
     * # Money Market basis (360 or 365)
     * 360
     * # Annual or semi-annual curve ("A" or "S")
     * S
     * # Year basis for benchmark swaps ("ACT", "365" or "360")
     * ACT
     * # No of entries
     * 84
     * #zero maturity yyyymmdd rates (ACT/365F annual)
     * 19970412 5.734349908886
     * 19970511 5.913551526763
     * etc
     */
    if (CrxFGetLine(buf, sizeof(buf),fp, &line) == NULL)
        goto done;
    if (CrxTDateScanYMD(buf, &baseDate) != 0)
        goto done;
    
    if (CrxFGetLine(buf, sizeof(buf),fp, &line) == NULL)
        goto done;
    if (CrxFGetLine(buf, sizeof(buf),fp, &line) == NULL)
        goto done;
    if (CrxFGetLine(buf, sizeof(buf),fp, &line) == NULL)
        goto done;
    
    if (CrxFGetLine(buf, sizeof(buf),fp, &line) == NULL)
        goto done;
    if (sscanf(buf, "%d", &nPts) != 1)
        goto done;
    
    basis = 1L;
    dayCountConv = GTO_ACT_365F;
    
    /* create new TCurve */
    *that = GtoNewTCurve(baseDate, nPts, basis, dayCountConv);
    if (*that == NULL) goto done;
    
    if (nPts > 0) {
        for (i=0; i<=nPts-1; i++) {
            if ((CrxFGetLine(buf, sizeof(buf), fp, &line) == NULL) ||
                (sscanf(buf, "%s %lf", buf2, &CrxTCurveRate(*that,i)) != 2)
                ||
                (CrxTDateScanYMD(buf2, &CrxTCurveDate(*that, i)) != 0)) {
                goto done;
            }
            CrxTCurveRate(*that, i) *= 1e-2;
        }
    }
    
    /* made it through */
    errCode = SUCCESS;
    
    
done:
    if (errCode != SUCCESS) {
        if (*that != NULL) GtoFreeTCurve(*that);
        *that = NULL;
        DR_Error("%s: failed (line %d).\n", routine, line);
    }
    return(errCode);
}




/**---------------------------------------------------------
 */
int
CrxTCurveFileRead(TCurve **that, char *fnam)
{
        FILE            *fp;
        int             errCode ;
static  char            routine[] = "CrxTCurveFileRead";


        if ((fp = fopen(fnam, "r")) == NULL) {
                DR_Error("%s: can't open file `%s'\n",
                        routine, fnam);
                errCode = -1;
                return errCode;
        } else {
                errCode = CrxTCurveFpRead(that, fp);
                if (errCode != 0)
                        DR_Error("%s: can't read file `%s' (code %d)\n",
                                routine, fnam, errCode);
                fclose(fp);
                return(errCode);
        }
}





/*f----------------------------------------------------------------------
 * Compares two elements <i> p</i> and <i> q</i> of type <i> type</i>
 * and returns an integer less
 * than, equal to, or greater than zero, depending on
 * whether <i> p</i> is lexicographically less than, equal to, or
 * greater than <i> q</i>.
 */

int
CrxVTypeCompare(const void *p, const void *q, TVType type)
{
static	char	routine[] = "CrxVTypeCompare";
	double	val1, val2;

#undef	_COMPARE
#define	_COMPARE(type)	\
		{if (*((type*)p) == *((type*)q))  return(0); \
		else if (*((type*)p) < *((type*)q))  return(-1); \
		else return(1);}

	switch (type) {
	case DRL_DOUBLE_T:
	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
		if (IS_ALMOST_ZERO(*((double*)p) - *((double*)q))) {
			return(0);
		} else if (*((double*)p) < *((double*)q)) {
			return(-1);
		} else {
			return(1);
		}
	case DRL_FLOAT_T:
		if (IS_ALMOST_ZERO(*((float*)p) - *((float*)q))) {
			return(0);
		} else if (*((float*)p) < *((float*)q)) {
			return(-1);
		} else {
			return(1);
		}
	case DRL_INT_T:
	case DRL_BOOLEAN_T:
		_COMPARE(int);
	case DRL_LONG_T:
		_COMPARE(long);
	case DRL_CHAR_T:
		_COMPARE(char);
	case DRL_STRING_T:
		return strcmp(*((char**)p), *((char**)q));
	case DRL_TDATE_T:
		_COMPARE(TDate);
	case DRL_TDATEINTERVAL_T:
		GtoDateIntervalToYears((TDateInterval*)p, &val1);
		GtoDateIntervalToYears((TDateInterval*)q, &val2);
		return (IS_ALMOST_ZERO(val2 - val1) ? 0 :
			(val1 < val2 ? -1 : 1));
	case DRL_TDAYCOUNT_T:
		_COMPARE(TDayCount);
	default:
		DR_Error("%s: bad type %ld.\n", routine, (long)type);
		return(0);
	}
#undef	_COMPARE
}



int
CrxVTypeLet(void *p, int offset1, void *q, int offset2, TVType type)
{
static  char    routine[] = "CrxVTypeLet";
        switch (type) {
        case DRL_DOUBLE_T:
        case DRL_PERCENT_T:
        case DRL_CUR_T:
        case DRL_CURK_T:
        case DRL_CURM_T:
                ((double*)p)[offset1] = ((double*)q)[offset2];
                break;
        case DRL_FLOAT_T:
                ((float*)p)[offset1] = ((float*)q)[offset2];
                break;
        case DRL_INT_T:
        case DRL_BOOLEAN_T:
                ((int*)p)[offset1] = ((int*)q)[offset2];
                break;
        case DRL_LONG_T:
                ((long*)p)[offset1] = ((long*)q)[offset2];
                break;
        case DRL_CHAR_T:
                ((char*)p)[offset1] = ((char*)q)[offset2];
                break;
        case DRL_STRING_T:
                strncpy(((char**)p)[offset1],
                        ((char**)q)[offset2],
                        WRAP_STR_BYTES);
                break;
        case DRL_CHAR_ARRAY_T:
                strcpy((char*)p, (char*)q);
                break;
        case DRL_TDATE_T:
                ((TDate*)p)[offset1] = ((TDate*)q)[offset2];
                break;
        case DRL_TDATEINTERVAL_T:
                ((TDateInterval*)p)[offset1] = ((TDateInterval*)q)[offset2];
                break;
        case DRL_TDAYCOUNT_T:
                ((TDayCount*)p)[offset1] = ((TDayCount*)q)[offset2];
                break;
        default:
                DR_Error("%s: bad type %ld.\n", routine, (long)type);
                return(1);
        }
        return(0);
}






/* qsort -- qsort interface implemented by faster quicksort.
   J. L. Bentley and M. D. McIlroy, SPE 23 (1993) 1249-1265.
   Copyright 1993, John Wiley.
*/

    /*assume sizeof(long) is a power of 2 */
#define SWAPINIT(a, es) swaptype =         \
    (((a)-(char*)0) | (es)) % sizeof(long) ? 2 : (es) > sizeof(long);
#define swapcode(TYPE, parmi, parmj, n) {  \
    register TYPE *pi = (TYPE *) (parmi);  \
    register TYPE *pj = (TYPE *) (parmj);  \
    do {                                   \
        register TYPE t = *pi;             \
        *pi++ = *pj;                       \
        *pj++ = t;                         \
    } while ((n -= sizeof(TYPE)) > 0);     \
}
#include <stddef.h>
static void swapfunc(char *a, char *b, size_t n, int swaptype)
{   if (swaptype <= 1) swapcode(long, a, b, n)
    else swapcode(char, a, b, n)
}
#define swap(a, b)                         \
    if (swaptype == 0) {                   \
        t = *(long*)(a);                   \
        *(long*)(a) = *(long*)(b);         \
        *(long*)(b) = t;                   \
    } else                                 \
        swapfunc(a, b, es, swaptype)

#define PVINIT(pv, pm)                                \
    if (swaptype != 0) { pv = a; swap(pv, pm); }      \
    else { pv = (char*)&v; *(long*)pv = *(long*)pm; }

#define vecswap(a, b, n) if (n > 0) swapfunc(a, b, n, swaptype)

#undef min
#define min(x, y) ((x)<=(y) ? (x) : (y))

static char *med3(char *a, char *b, char *c, int (*cmp)())
{       return cmp(a, b) < 0 ?
                  (cmp(b, c) < 0 ? b : cmp(a, c) < 0 ? c : a)
                : (cmp(b, c) > 0 ? b : cmp(a, c) > 0 ? c : a);
}



void CrxQSort(
        void *base,
        size_t n,
        size_t es,
        int (*cmp)(const void *, const void *))
{
    char    *a = (char*) base;
    char *pa, *pb, *pc, *pd, *pl, *pm, *pn, *pv;
    int r, swaptype;
    long t, v;
    size_t s;
    
    SWAPINIT(a, es);
    if (n < 7) {     /* Insertion sort on smallest arrays */
        for (pm = a + es; pm < a + n*es; pm += es)
            for (pl = pm; pl > a && cmp(pl-es, pl) > 0; pl -= es)
                swap(pl, pl-es);
        return;
    }
    pm = a + (n/2)*es;    /* Small arrays, middle element */
    if (n > 7) {
        pl = a;
        pn = a + (n-1)*es;
        if (n > 40) {    /* Big arrays, pseudomedian of 9 */
            s = (n/8)*es;
            pl = med3(pl, pl+s, pl+2*s, cmp);
            pm = med3(pm-s, pm, pm+s, cmp);
            pn = med3(pn-2*s, pn-s, pn, cmp);
        }
        pm = med3(pl, pm, pn, cmp); /* Mid-size, med of 3 */
    }
    PVINIT(pv, pm);       /* pv points to partition value */
    pa = pb = a;
    pc = pd = a + (n-1)*es;
    for (;;) {
        while (pb <= pc && (r = cmp(pb, pv)) <= 0) {
            if (r == 0) { swap(pa, pb); pa += es; }
            pb += es;
        }
        while (pb <= pc && (r = cmp(pc, pv)) >= 0) {
            if (r == 0) { swap(pc, pd); pd -= es; }
            pc -= es;
        }
        if (pb > pc) break;
        swap(pb, pc);
        pb += es;
        pc -= es;
    }
    pn = a + n*es;
    s = min(pa-a,  pb-pa   ); vecswap(a,  pb-s, s);
    s = min(pd-pc, pn-pd-es); vecswap(pb, pn-s, s);
    if ((s = pb-pa) > es) CrxQSort(a,    s/es, es, cmp);
    if ((s = pd-pc) > es) CrxQSort(pn-s, s/es, es, cmp);
}



static  TVType  _CrxVTypeVectSort_type;

static  int
_CrxVTypeVectSort_compare(const void *p, const void *q)
{
        return CrxVTypeCompare(p, q, _CrxVTypeVectSort_type);

}


/*f-------------------------------------------------------------
 * Sorts a vector <i> p</i> with <i> numItem</i> elements
 * of type <i> type</i> in ascending order.
 * If the flag <i> removeMultiple</i> is TRUE, removes multiple
 * instances of a same element, and <i> numItems</i>
 * may be changed on exit.
 */


int
CrxVTypeVectSort(
	void *p,		/* (B) vector */
	int *numItems,		/* (B) num elem in vector */
	TVType type,		/* (I) vector type */
	int removeMultiple)	/* (I) TRUE=remove multiple items */
{
static	char	routine[] = "CrxVTypeVectSort";
	int	status = FAILURE;
	size_t	nel;
	size_t	size;
	int	idx, idx2;

	/* Perform sorting */
	_CrxVTypeVectSort_type = type;
	size = (size_t) CrxVTypeSizeof(type);
	nel = (size_t) *numItems;

	CrxQSort((void*) p,
		(size_t) nel,
		(size_t) size, 
		_CrxVTypeVectSort_compare);

	/* Remove double items */
	if (removeMultiple) {
	    idx = 0;
	    while (idx < *numItems-1) {
		if (CrxVTypeCompare(
			    CrxVTypeOffsetVect(p, idx,   type),
			    CrxVTypeOffsetVect(p, idx+1, type),
			    type) == 0) {
		    for (idx2=idx+1; idx2<*numItems-1; idx2++) {
		        if (CrxVTypeLet(
				p, idx2,
				p, idx2+1,
				type) != SUCCESS)
					goto done;
		    }
		    --(*numItems);
		} else {
		    idx++;
		}
	    }
	}


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DR_Error("%s: failure.\n", routine);
	}
	return(status);
}




/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */

int CrxVTypeNameScan(TVType *type, char *s)
{
static	char	routine[] = "CrxVTypeNameScan";
#undef	CHECK
#define	CHECK(str, val)		{if (!strcmp(str, s)) \
				{*type = val; return (SUCCESS);}}

	CHECK("double",	DRL_DOUBLE_T);
	CHECK("float",	DRL_FLOAT_T);
	CHECK("int",	DRL_INT_T);
	CHECK("long",	DRL_LONG_T);
	CHECK("char",	DRL_CHAR_T);

	CHECK("char*",		DRL_STRING_T);
	CHECK("string",		DRL_STRING_T);
	CHECK("char_array",	DRL_CHAR_ARRAY_T);
	CHECK("date",		DRL_TDATE_T);
	CHECK("date_interval",	DRL_TDATEINTERVAL_T);
	CHECK("daycount",	DRL_TDAYCOUNT_T);

	CHECK("percent", DRL_PERCENT_T);
	CHECK("cur", DRL_CUR_T);
	CHECK("curK", DRL_CURK_T);
	CHECK("curM", DRL_CURM_T);
	CHECK("boolean", DRL_BOOLEAN_T);

	DR_Error("%s: can't read type in `%s'.\n", routine, s);
	return(FAILURE);


#undef	CHECK
}



/*f-------------------------------------------------------------
 * Allocates a array of dimension <i> nDim</i> of type <i> type</i>
 * and sizes in each dimensions given by the values of 
 * the array <i> nSize</i> (an vector of length <i> nDim</i>).
 */

void*
CrxVTypeTensorAlloc(TVType type, int nDim, int *nSize)
{
static	char	routine[] = "CrxVTypeTensorAlloc";
	int	i;
	void	**ptr = NULL;

	if (nDim > 1) {
	    ptr = NEW_ARRAY(void*, nSize[0]);
	    if (ptr == NULL) return(NULL);
	    for (i=0; i<=nSize[0]-1; i++) {
		ptr[i] = CrxVTypeTensorAlloc(type, nDim-1, nSize+1);
		if (ptr[i] == NULL) {
		    return(NULL);
		}
	    }
	    return(ptr);
	} else if (nDim == 1) {
	    ptr = CrxVTypeVectAlloc(nSize[0], type);
	    if (ptr == NULL) {
		return(NULL);
	    }
	    return(ptr);
	} else {
	    DR_Error("%s: bad dimension %d\n", routine, nDim);
	    return(NULL);
	}
}



/*f-------------------------------------------------------------
 * Frees a array of arbitrary dimension  allocated by
 * <i> CrxVTypeTensorAlloc</i>.
 */

int
CrxVTypeTensorFree(void *ptr, TVType type, int nDim, int *nSize)
{
static	char	routine[] = "CrxVTypeTensorAlloc";
	int	i;

	if (ptr == NULL) return(SUCCESS);
	if (nDim > 1) {
	    for (i=0; i<=nSize[0]-1; i++) {
		CrxVTypeTensorFree((void*) ((void**) ptr)[i],
			type, nDim-1, nSize+1);
	    }
	    FREE(ptr);
	    return(SUCCESS);
	} else if (nDim == 1) {
	    CrxVTypeVectFree(ptr, nSize[0], type);
	    return(SUCCESS);
	} else {
	    DR_Error("%s: bad dimension %d\n", routine, nDim);
	    return(FAILURE);
	}
}


/*f-------------------------------------------------------------
 * Same a <i> CrxVTypeTensorAlloc</i>, but takes a variable
 * number of arguments instead of an array for the sizes
 * in each dimension.
 */

void*
CrxVTypeTensorAllocV(TVType type, int nDim,
		... /* int size1, ..., int sizeNDim */)
{
static	char	routine[] = "CrxVTypeTensorAllocV";
#define	NMAXDIM	32
	int	n, nSize[NMAXDIM];
	va_list	ap;

	va_start(ap, nDim);
	if (nDim > NMAXDIM) {
	    DR_Error("%s: maximum dimension %d (got %d).\n",
		routine, nDim, NMAXDIM);
	    return(NULL);
	}
	for (n=0; n<=nDim-1; n++) {
		nSize[n] = va_arg(ap, int);
	}
	va_end(ap);

	return CrxVTypeTensorAlloc(type, nDim, nSize);
#undef	NMAXDIM
}



/*f-------------------------------------------------------------
 * Frees a array of arbitrary dimension created by
 * <i> CrxVTypeTensorAllocV</i>.
 */

int
CrxVTypeTensorFreeV(void *ptr, TVType type, int nDim,
		... /* int size1, ..., sizeNDim */)
{
static	char	routine[] = "CrxVTypeTensorFreeV";
#define	NMAXDIM	32
	int	n, nSize[NMAXDIM];
	va_list	ap;

	va_start(ap, nDim);
	if (nDim > NMAXDIM) {
	    DR_Error("%s: maximum dimension %d (got %d).\n",
		routine, nDim, NMAXDIM);
	    return(-1);
	}
	for (n=0; n<=nDim-1; n++) {
		nSize[n] = va_arg(ap, int);
	}
	va_end(ap);

	return CrxVTypeTensorFree(ptr, type, nDim, nSize);
#undef	NMAXDIM
}



/*f-------------------------------------------------------------
 * In a array of dimension <i> nDim</i> of type <i> type</i>
 * and sizes in each dimensions given by the values of 
 * the array <i> nSize</i> (an vector of length <i> nDim</i>),
 * the routine accesses and return (as a <i> void</i> pointer)
 * the element having <i> idx</i> multi-index (a vector
 * of length <i> nDim</i> containing the indices in each dimension).
 */

void*
CrxVTypeTensorAccessElement(
	void *ptr, TVType type, int nDim, int *nSize,
	int *idx)		/* multi-index [0..nDim-1] */
{
static	char	routine[] = "CrxVTypeTensorAccessElement";

	if (ptr == NULL) return(NULL);
	if (nDim > 1) {
	    return CrxVTypeTensorAccessElement(
			(void*) ((void**) ptr)[idx[0]],
			type, nDim-1, nSize+1,
			idx+1);
	} else if (nDim == 1) {
	    return CrxVTypeOffsetVect(ptr, idx[0], type);
	} else {
	    DR_Error("%s: bad dimension %d.\n", routine, nDim);
	    return(NULL);
	}
}



/*--------------------------------------------------------------
 */

int*
CrxVTypeTensorMultiIndexNew(int nDim, int *nSize)
{
    int	*idx;
    int	i;
    
    nSize = nSize; /* avoid warning on unused parameter */

    if ((idx = NEW_ARRAY(int, nDim)) == NULL)
        return(NULL);
    for (i=0; i<=nDim-1; i++) {
        idx[i] = 0;
    }
    return(idx);
}

int
CrxVTypeTensorMultiIndexNext(int *idx, int nDim, int *nSize)
{
	int	i;
	for (i=nDim-1; i>=0; i--) {
		idx[i]++;
		if (idx[i] < nSize[i]) {
			return(TRUE);
		}
		idx[i] = 0;
	}
	FREE(idx);
	return(FALSE);
}


/*f-------------------------------------------------------------
 * In a array of dimension <i> nDim</i> of type <i> type</i>
 * and sizes in each dimensions given by the values of 
 * the array <i> nSize</i> (an vector of length <i> nDim</i>),
 * performs a scalar operation. Available operations
 * are:\\
 * <i> "sum"</i>:  add all elements, \\
 * <i> "avg"</i>:  average all elements, \\
 * <i> "min"</i>:  compute minimum of all elements, \\
 * <i> "max"</i>:  compute minimum of all elements, \\
 * <i> "min0"</i>: compute $l^1$-norm of the negative elements, \\
 * <i> "max0"</i>: compute $l^1$-norm of the positive elements, \\
 * <i> "l1"</i>:   compute $l^1$-norm, \\
 * <i> "l2"</i>:   compute $l^2$-norm, \\
 * <i> "="</i>:    set all elements to <i> value</i>, \\
 * <i> "+="</i>:   add <i> value</i> to all elements, \\
 * <i> "*="</i>:   multiply by <i> value</i> all elements.\\
 */

int
CrxVTypeTensorOperScalar(
	void *ptr, TVType type, int nDim, int *nSize,
	char *operName,		/* (I) see details */
	double *value)		/* (B) input or output */
{
    static	char	routine[] = "CrxVTypeTensorOperScalar";
	int	status = FAILURE;
	int	cnt, *idx = NULL;
	double	*dval;

	if (type != DRL_DOUBLE_T) {
	    DR_Error("%s: type must be double.\n", routine);
	    goto done;
	}

#define	LOOPSTART	if ((idx = CrxVTypeTensorMultiIndexNew(nDim, nSize)) == NULL)\
			goto done; do {\
			if ((dval = (double*) CrxVTypeTensorAccessElement(\
			ptr, type, nDim, nSize, idx)) == NULL) goto done;
#define	LOOPEND		} while (CrxVTypeTensorMultiIndexNext(idx, nDim, nSize) != FALSE);

#define	IF_OPER(name)	if (!strcmp(operName, name))


	if (!strcmp(operName, "sum")) {
		*value = 0e0;
		LOOPSTART
			*value += *dval;
		LOOPEND
	} else if (!strcmp(operName, "min")) {
		*value = 1e32;
		LOOPSTART
			*value = MIN(*dval, *value);
		LOOPEND
	} else if (!strcmp(operName, "max")) {
		*value = -1e32;
		LOOPSTART
			*value = MAX(*dval, *value);
		LOOPEND
	} else if (!strcmp(operName, "min0")) {
		*value = 0e0;
		LOOPSTART
			*value += MAX(-(*dval), 0e0);
		LOOPEND
	} else if (!strcmp(operName, "max0")) {
		*value = 0e0;
		LOOPSTART
			*value += MAX(*dval, 0e0);
		LOOPEND
	} else if (!strcmp(operName, "l1")) {
		*value = 0e0;
		LOOPSTART
			*value += fabs(*dval);
		LOOPEND
	} else if (!strcmp(operName, "L1")) {
		*value = 0e0;
		cnt = 0;
		LOOPSTART
			*value += fabs(*dval);
			cnt++;
		LOOPEND
		*value /= (double)cnt;
	} else if (!strcmp(operName, "l2")) {
		*value = 0e0;
		LOOPSTART
			*value += (*dval)*(*dval);
		LOOPEND
		*value = sqrt(*value);
	} else if (!strcmp(operName, "avg")) {
		*value = 0e0;
		cnt = 0;
		LOOPSTART
			*value += *dval;
			cnt++;
		LOOPEND
		*value /= (double)cnt;
	} else if (!strcmp(operName, "=")) {
		LOOPSTART
			*dval = *value;
		LOOPEND
	} else if (!strcmp(operName, "+=")) {
		LOOPSTART
			*dval = *value;
		LOOPEND
	} else if (!strcmp(operName, "*=")) {
		LOOPSTART
			*dval = *value;
		LOOPEND
	} else {
		DR_Error("%s: unknown operation %s.\n", routine, operName);
		goto done;
	}

	status = SUCCESS;
done:
	return(status);
}

/*f-------------------------------------------------------------
 * Same a <i> CrxVTypeTensorOperScalar</i>, but takes a variable
 * number of arguments instead of an array for the sizes
 * in each dimension.
 */

int
CrxVTypeTensorOperScalarV(
    void *ptr,		/* (I) tensor */
    TVType type,		/* (I) tensor elements type */
    char *operName,		/* (I) see details */
    double *value,		/* (B) input or output */
    int nDim,		/* (I) number of dimensions */
    ... /* int size1, ..., int sizeNDim */)
{
    static	char	routine[] = "CrxVTypeTensorOperScalarV";
#define	NMAXDIM	32
    int	n, nSize[NMAXDIM];
    va_list	ap;
    
    va_start(ap, nDim);
    if (nDim > NMAXDIM) {
        DR_Error("%s: maximum dimension %d (got %d).\n",
                 routine, nDim, NMAXDIM);
        return FAILURE;
    }
    for (n=0; n<=nDim-1; n++) {
        nSize[n] = va_arg(ap, int);
    }
    va_end(ap);
    
    return CrxVTypeTensorOperScalar(
        ptr, type, nDim, nSize, operName, value);
#undef	NMAXDIM
}


/*f-------------------------------------------------------------
 */

int
CrxVTypeTensorFpWrite(
	void *ptr, TVType type, int nDim, int *nSize,
	FILE *fp)		/* (I) */
{
static	char	routine[] = "CrxVTypeTensorFpWrite";
	int	status = FAILURE;
	int	i, *idx = NULL;
	void	*element;


	fprintf(fp, "%d\n", nDim);
	for (i=0; i<=nDim-1; i++) {
		fprintf(fp, "\t%d", nSize[i]);
	}
	fprintf(fp, "\n");

	if ((idx = CrxVTypeTensorMultiIndexNew(nDim, nSize)) == NULL)
		goto done;
	do {
		if ((element = CrxVTypeTensorAccessElement(ptr,
			type, nDim, nSize, idx)) == NULL) goto done;

		fprintf(fp, "\t%s",
			CrxVTypePrint(NULL, type, element));

		fprintf(fp, "\t#");
		for (i=0; i<=nDim-1; i++) {
			fprintf(fp, " %3d", idx[i]);
		}

		fprintf(fp, "\n");


	} while (CrxVTypeTensorMultiIndexNext(idx, nDim, nSize) != FALSE);


	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DR_Error("%s: failed.\n", routine);
	}
	return(status);
}


/* simple set of DR wrapper utility functions for we changed the prefix
   and compiled locally */

#define COMMENT_CHR '#'
#define MAX_TOKEN_LEN 256
#define FORMAT_STRING "%255s"

/* REQUIRE should be used to check input conditions. */
/* REQUIRE assumes existence of 'done' label for error cases. */
/* REQUIRE assumes that the function is named 'routine'. */
#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    GtoErrMsg("%s: Required condition (%s) fails!\n",routine,#cond);\
    goto done;\
}} while(0)

/*
***************************************************************************
** Reads a size type from file.
***************************************************************************
*/
int CrxReadSize(FILE *fp, size_t *val)
{
    return CrxReadGenericData (fp, "Size", 
                               (CrxTStringConverter)CrxStringToSize, val);
}

/*
***************************************************************************
** Reads an integer from file.
***************************************************************************
*/
int CrxReadInt(FILE *fp, int *val)
{
    return CrxReadGenericData (fp, "Integer", 
                               (CrxTStringConverter)CrxStringToInt, val);
}

/*
***************************************************************************
** Reads a long integer from file.
***************************************************************************
*/
int CrxReadLong(FILE *fp, long *val)
{
    return CrxReadGenericData (fp, "Long", 
                               (CrxTStringConverter)CrxStringToLong, val);
}

/*
***************************************************************************
** Reads a double from file.
***************************************************************************
*/
int CrxReadDouble(FILE *fp, double *val)
{
    return CrxReadGenericData (fp, "Double", 
                               (CrxTStringConverter)CrxStringToDouble, val);
}

/*
***************************************************************************
** Reads a date from file. The format is YYYYMMDD.
***************************************************************************
*/
int CrxReadDate(FILE *fp, TDate *val)
{
    return CrxReadGenericData (fp, "Date", 
                               (CrxTStringConverter)CrxStringToDate, val);
}

/*
***************************************************************************
** Reads a Boolean from file.
***************************************************************************
*/
int CrxReadBoolean(FILE *fp, TBoolean *val)
{
    return CrxReadGenericData (fp, "Boolean", 
                               (CrxTStringConverter)CrxStringToBoolean, val);
}

/*
***************************************************************************
** Reads a string from file.
***************************************************************************
*/
int CrxReadString(FILE *fp, char **val)
{
    return CrxReadGenericData (fp, "String", 
                               (CrxTStringConverter)CrxStringToString, val);
}

/*
***************************************************************************
** Reads a number of arrays from file where the elements of each array
** are corresponding.
**
** Before reading the arrays, the size of the array will be read from file
** via a previous call. Then the caller should allocate the arrays of
** that size and pass the empty arrays to this function.
**
** The caller should provide the type of the arrays providing a combination
** of D (date), F (double), L (long), S (string), I (integer), B (boolean),
** G (generic).
**
** For each string within the array types, either one or four parameters
** in the variable argument list.
**
** For generic data type, there should be four parameters - first should
** be the name of the type (for error handling purposes only), second should
** be the StringConverter function and second should be sizeof data, and
** third should be an array of the generic data type.
**
** For defined data types, there should be one parameter - an array of
** the defined pointer type.
**
** e.g. CrxReadArrays(fp, 24, "DF", TDate*, double*)
***************************************************************************
*/
int CrxReadArrays
(FILE  *fp,         /* (I) File pointer */
 size_t arraySize,  /* (I) Size of arrays */
 char  *arrayTypes, /* (I) Types of arrays */
 ...)               /* (O) Pointers to arrays for output */
{
    static char routine[] = "CrxReadArrays";
    int status = FAILURE;

    va_list ap; /* Variable argument list pointer */

    void              **arrays     = NULL;
    CrxTStringConverter *converters = NULL;
    size_t             *sizes      = NULL;
    char              **names      = NULL;
    size_t              numArrays  = strlen(arrayTypes);
    size_t              i;
    size_t              j;

    va_start(ap, arrayTypes);

    REQUIRE (numArrays > 0);

    if (arraySize == 0) return SUCCESS;

    arrays     = NEW_ARRAY (void*, numArrays);
    converters = NEW_ARRAY (CrxTStringConverter, numArrays);
    sizes      = NEW_ARRAY (size_t, numArrays);
    names      = NEW_ARRAY (char*, numArrays);

    for (j = 0; j < numArrays; ++j)
    {
        char           *name;
        CrxTStringConverter converter;
        size_t          size;
        switch (arrayTypes[j])
        {
        case 'D':
            name      = "Date";
            converter = (CrxTStringConverter)CrxStringToDate;
            size      = sizeof(TDate);
            break;
        case 'F':
            name      = "Double";
            converter = (CrxTStringConverter)CrxStringToDouble;
            size      = sizeof(double);
            break;
        case 'L':
            name      = "Long";
            converter = (CrxTStringConverter)CrxStringToLong;
            size      = sizeof(long);
            break;
        case 'S':
            name      = "String";
            converter = (CrxTStringConverter)CrxStringToString;
            size      = sizeof(char*);
            break;
        case 'I':
            name      = "Int";
            converter = (CrxTStringConverter)CrxStringToInt;
            size      = sizeof(int);
            break;
        case 'B':
            name      = "Boolean";
            converter = (CrxTStringConverter)CrxStringToBoolean;
            size      = sizeof(TBoolean);
            break;
        case 'G':
            name      = va_arg (ap, char*);
            converter = va_arg (ap, CrxTStringConverter);
            size      = va_arg (ap, size_t);
            break;
        default:
            GtoErrMsg ("%s: Invalid array type %c - should be D,F,L,S,I,B "
                       "or G\n", routine, arrayTypes[j]);
            goto done; /* failure */
        }
        arrays[j]     = va_arg (ap, void*);
        converters[j] = converter;
        sizes[j]      = size;
        names[j]      = name;
    }
    
    for (i = 0; i < arraySize; ++i)
    {
        for (j = 0; j < numArrays; ++j)
        {
            void *val = (void*)((char*)(arrays[j]) + (i * sizes[j]));
            if (CrxReadGenericData(fp, names[j], converters[j], 
                                   val) != SUCCESS)
                goto done; /* failure */
        }
    }

    status = SUCCESS;

done:

    va_end (ap);

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
    }

    FREE (arrays);
    FREE (converters);
    FREE (sizes);
    FREE (names);

    return status;
}

/*
***************************************************************************
** Generic read routine. Reads a token (maxlength = 255) and converts it
** to the underlying type via the types StringConverter routine.
***************************************************************************
*/
int CrxReadGenericData (FILE* fp, const char* name, CrxTStringConverter func,
                     void* val)
{
    char  buf[256];
    char *str;

    str = CrxReadToken (fp, sizeof(buf), buf);
    if (str == NULL)
    {
        GtoErrMsg ("CrxRead%s: Unexpected end of file\n", name);
        return FAILURE;
    }

    if (func(str, val) != SUCCESS)
    {
        GtoErrMsg ("CrxRead%s: Failed\n", name);
        return FAILURE;
    }
    
    return SUCCESS;
}


/**
***************************************************************************
** Reads next token from file. The token is returned as a pointer to the
** character buffer provided by the caller.
**
** Quotation marks (either single or double) are interpreted as the
** beginning of text and may be mixed freely (e.g. to include a double
** quote surround it by single quotes and vice versa). This allows white
** space to be embedded within the token.
**
** Comments are skipped and are taken to go to the end of the current line.
** The comment character is #. Comments encountered within text are not
** skipped.
**
** If the maximum token length (= sz-1) is exceeded, then the returned token
** is truncated and not included in the following token.
***************************************************************************
*/
char* CrxReadToken(FILE *fp, size_t sz, char *buf)
{
    static char routine[] = "CrxReadToken";
    int         status    = FAILURE;

    int         c;
    int         nc;
    char        commentChar = '#';

    TBoolean noToken    = TRUE;
    TBoolean endOfToken = FALSE;
    int      bufPos     = 0;
    int      maxBufPos  = sz-1;

    REQUIRE (sz > 2);

    memset (buf, 0, sz);

    while (!(endOfToken))
    {
        c = fgetc(fp);
        switch (c)
        {
        case '"':
        case '\'':
            /*
            ** We have found a quotation mark. Seek to the next quotation
            ** mark, and include everything we find in between! EOF also
            ** ends the text (this function does not report errors).
            */
            noToken = FALSE;
            for (;;)
            {
                nc = fgetc(fp);
                if (nc == c || nc == EOF)
                    break;
                if (bufPos < maxBufPos)
                {
                    buf[bufPos] = nc;
                    ++bufPos;
                }
            }
            break;

        case EOF:
            endOfToken = TRUE;
            break;

        default:
            /*
            ** Any remaining character, including white space and comment
            ** character (note the comment character is not a constant and
            ** so cannot be included in the switch).
            */
            if (c == commentChar)
            {
                /*
                ** We have found a comment character. We read to the end of
                ** the line. If we have started a token, then we are done,
                ** otherwise we continue.
                */
                do
                {                            
                    nc = fgetc(fp);
                }  while (nc != (int)'\n' && nc != EOF);
                if (!(noToken))
                {
                    endOfToken = TRUE;
                }
            }
            else if (isspace(c))
            {
                /*
                ** If we are at the beginning of the token (noToken=TRUE),
                ** then we skip white space. Otherwise white space ends the
                ** token. 
                */
                if (!(noToken))
                {
                    endOfToken = TRUE;
                }
            }
            else
            {
                /*
                ** Just a normal character - neither a comment nor space.
                */
                noToken = FALSE;
                if (bufPos < maxBufPos)
                {
                    buf[bufPos] = c;
                    ++bufPos;
                }
            }
        }
    }

    if (noToken)
        return NULL;

    if (bufPos > maxBufPos)
    {
        PROGRAM_BUG();
        bufPos = maxBufPos;
    }

    buf[bufPos] = '\0';
    status      = SUCCESS;

 done:

    if (status != SUCCESS)
        return NULL;

    return buf;
}

int CrxStringToSize (const char* str, size_t *val)
{
    long tmp;
    if (CrxStringToLong (str, &tmp) != SUCCESS) 
        return FAILURE;

    *val = (size_t)tmp;
    if ((long)(*val) != tmp)
        return FAILURE; /* out of range for size_t */

    return SUCCESS;
}
    
int CrxStringToInt (const char* str, int *val)
{
    long tmp;
    if (CrxStringToLong (str, &tmp) != SUCCESS) 
        return FAILURE;

    *val = (int)tmp;
    if ((long)(*val) != tmp)
        return FAILURE; /* out of range for int */

    return SUCCESS;
}

int CrxStringToLong (const char* str, long *val)
{
    char *ep;
    *val = strtol (str, &ep, 0);
    if (*ep != '\0')
        return FAILURE; /* not really an integer */
    
    return SUCCESS;
}

int CrxStringToDouble (const char* str, double *val)
{
    char *ep;
    *val = strtod (str, &ep);
    if (*ep != '\0')
        return FAILURE; /* not really a double */
    
    return SUCCESS;
}

int CrxStringToDate (const char* str, TDate *val)
{
    long tmp;
    if (CrxStringToLong (str, &tmp) != SUCCESS)
        return FAILURE;

    *val = GtoDate (tmp/10000, tmp/100 % 100, tmp % 100);
    if (*val < 0)
        return FAILURE; /* not really a date */

    return SUCCESS;
}

int CrxStringToBoolean (const char* str, TBoolean *val)
{
    if (strcmp(str, "1") == 0) 
        *val = TRUE;
    else if (strcmp(str, "0") == 0)
        *val = FALSE;
    else if (strcmp(str, "TRUE") == 0)
        *val = TRUE;
    else if (strcmp(str, "FALSE") == 0)
        *val = FALSE;
    else
        return FAILURE;

    return SUCCESS;
}

int CrxStringToString (const char* src, char **val)
{
    *val = GtoStringDuplicate((char*)src);
    return SUCCESS;
}
