/************************************************************************
 * Module:	DRL
 * Submodule:	VTYPE
 * File:	
 * Function:	Variable Type
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"
#include "drlerr.h"		/* DrlErrMsg */

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>		/* IEEE finite() */
#include <float.h>		/* DBL_MIN, etc. */

#include "drltime.h"		/* DrlDDateScan() */
#include "drlstr.h"
#include "drlio.h"
#include "drlsort.h"

#include "drlvtype.h"


/*f-------------------------------------------------------------
 * Returns returns 1 only when -INFINITY < x < +INFINITY.
 * Otherwise it returns 0 (i.e., when |x| = INFINITY or x is NaN).
 */

int
DrlDoubleIsFinite(double x)
{
#if defined(UNIX)
	/* UNIX */
	return finite(x);
#elif defined(_WIN32) || defined(WIN32)
	/* NT */
	return _finite(x);
#else
	/* not implemented */
	/*return (((x) == 0e0) ||
		(((x) >= DBL_MIN) && ((x) <= DBL_MAX)) ||
		((-(x) >= DBL_MIN) && (-(x) <= DBL_MAX)));*/
	return (TRUE);
#endif
}


/*f----------------------------------------------------------------------
 * Allocates a vector of type "type" of length "size".
 */

void*
DrlVTypeVectAlloc(int size, DVType type)
{
static	char	routine[] = "DrlVTypeVectAlloc";
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
		p = (void*) NEW_ARRAY(DDate, size+1);
		if (p) ((DDate*)p)[0] = (DDate)size;
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
			DrlErrMsg("%s: malloc DRL_STRING_T "
				"length %d failed.\n", routine, (int) (size));
			return(NULL);
		}
	    	for (i=0; i<=size-1; i++) {
			((char**)p)[i] = NEW_ARRAY(char,  DRL_WRAP_STR_BYTES);
			if (((char**)p)[i] == NULL) {
				DrlErrMsg("%s: malloc failed.\n", routine);
				return(NULL);
			}
		}
		break;


	default:
		/* check type size */
		if (DrlVTypeCheckValidCType(type) != 0) {
			DrlErrMsg("%s: bad C type.\n", routine);
			return(NULL);
		}

	 	/* get type size */
		if ((typeSize = DrlVTypeSizeof(type)) <= 0) {
			DrlErrMsg("%s: bad type size.\n");
			return(NULL);
		}

		/* allocate memory */
		if ((p = (void*) MALLOC((size_t) (typeSize*size))) == NULL) {
			DrlErrMsg("%s: malloc length %d failed.\n", routine, 
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

int
DrlVTypeVectFree(void* p, int size, DVType type)
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
		FREE((DDate*) p);
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
		if (DrlVTypeCheckValidCType(type) != 0) {
			DrlErrMsg("DrlVTypeVectFree: bad C type\n");
			return(1);
		}

		FREE((void*) p);

		return(0);
	}
	return(SUCCESS);
}

/*f----------------------------------------------------------------------
 * Allocates a matrix of type "type" of size "nx" by "ny".
 */

void*
DrlVTypeMatrAlloc(int nx, int ny, DVType type)
{
static	char	routine[] = "DrlVTypeMatrAlloc";
	int	i, j, errCode = 1;
	void	**p = NULL;
	size_t	sz;


	/* check type size */
	if (DrlVTypeCheckValidCType(type) != 0) {
		DrlErrMsg("%s: bad C type.\n", routine);
		return(NULL);
	}


	/* get type size */
	if ((sz = sizeof(void*)*nx) <= 0) goto done;

	/* allocate memory */
	if ((p = (void**) MALLOC((size_t) sz)) == NULL) {
		DrlErrMsg("%s: malloc failed.\n", routine);
		return(NULL);
	}

	if ((sz = DrlVTypeSizeof(type)*ny) <= 0) goto done;
	for (i=0; i<=nx-1; i++) {
	    if ((p[i] = (void*) MALLOC((size_t) sz)) == NULL) {
		DrlErrMsg("%s: malloc failed.\n", routine);
		return(NULL);
	    }
	}

	/* for DRL_STRING_T, allocate char string */
	if (type == DRL_STRING_T) {
	    for (i=0; i<=nx-1; i++) 
	    for (j=0; j<=ny-1; j++) {
		((char***)p)[i][j] = (char*) MALLOC((size_t) 
				DRL_WRAP_STR_BYTES*sizeof(char));
		if (((char***)p)[i][j] == NULL) {
			DrlErrMsg("%s: malloc failed\n", routine);
			return(NULL);
		}
	    }
	}

	errCode = 0;
done:
	if (errCode != 0) {
		DrlErrMsg("%s: failed (type %d, size %dx%d).\n",
			routine, type, nx, ny);
	}
	return((void*) p);
}



/*f----------------------------------------------------------------------
 * Frees a matrix
 * of type "type" of size "nx" by "ny".
 */

int
DrlVTypeMatrFree(void* p, int nx, int ny, DVType type)
{
	int	i, j;

	/* check type size */
	if (DrlVTypeCheckValidCType(type) != 0) {
		DrlErrMsg("DrlVTypeMatrFree: bad C type.\n");
		return(1);
	}

	if (p == NULL) return(0);

	if (type == DRL_STRING_T) {
	    for (i=0; i<=nx-1; i++)
	    for (j=0; j<=ny-1; j++) {
		FREE((void*) ((char***)p)[i][j]);
	    }
	}

	for (i=0; i<=nx-1; i++) {
		FREE(((void**) p)[i]);
	}
	FREE((void*) p);
	return(0);
}




/*f----------------------------------------------------------------------
 * Scans in a formatted string <i> str</i> a type <i> type</i> and
 * puts the result in <i> p</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> DrlVTypeOffsetVect</i>, even if it is the first element.}
 * Returns 0 if scan successful.
 */

int
DrlVTypeScan(char *str, DVType type, void *p)
{
static	char	routine[] = "DrlVTypeScan";
	char	*buf2, buf[256];

	switch (type) {
	case DRL_DOUBLE_T:
		/*if (sscanf(str, "%lf", (double*)p) != 1) goto error; */
		if (DrlCurScan(str, (double*)p) != SUCCESS) goto error;
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
		if (DrlStrScanString(&buf2, ((char*)p)) == NULL)
			goto error;
		/*if (sscanf(str, "%s",  ((char*)p)) != 1) goto error;*/
		break;
	case DRL_TDATE_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (DrlDDateScan(buf, (DDate*)p) != 0) goto error;
		break;
	case DRL_TDATEINTERVAL_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (DrlDIntervalScan(buf, (DInterval*)p)
				!= SUCCESS) goto error;
		break;
	case DRL_TDAYCOUNT_T:
		if (sscanf(str, "%s",  buf) != 1) goto error;
		if (DrlDDayCountScan(buf, 
			(DDayCount*)p) != SUCCESS) goto error;
		break;

		/**
		 ** Derived types
		 **/
	case DRL_CUR_T:
		if (DrlCurScan(str, (double*)p) != SUCCESS) goto error;
		break;
	case DRL_CURK_T:
		if (DrlCurScan(str, (double*)p) != SUCCESS) goto error;
		*((double*)p) *= 1e3;
		break;
	case DRL_CURM_T:
		if (DrlCurScan(str, (double*)p) != SUCCESS) goto error;
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
		if (DrlDDateScan(buf, (DDate*)p) != 0) goto error;
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
		DrlErrMsg("%s: bad type %ld.\n", routine, (long) type);
		return(2);
	}

	return(0);
error:
	DrlErrMsg("%s: can't scan `%s' (type %s)\n",
		routine, str, DrlVTypeName(type));
	return(1);
}


/*f----------------------------------------------------------------------
 * Prints in a formatted string <i> str</i> a pointer <i> </i>p
 * of type <i> type</i>.  Returns <i> str</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> DrlVTypeOffsetVect</i>, even if it is the first element.}
 */

char*
DrlVTypePrint(char *str, DVType type, void *p)
{
static	char	routine[] = "DrlVTypePrint";
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
		sprintf(str, "%10s", DrlDDatePrint(NULL, *((DDate*) p)));
		break;
	case DRL_TDATEINTERVAL_T:
		sprintf(str, "%6s",
			DrlDIntervalPrint(NULL, ((DInterval*) p)));
		break;
	case DRL_TDAYCOUNT_T:
		sprintf(str, "%10s",
			DrlDDayCountPrint(NULL, *((DDayCount*) p)));
		break;

		/* derived Types */

	case DRL_PERCENT_T:
		sprintf(str, "%1.8f", *((double*) p) * 1e2);
		break;
	case DRL_CUR_T:
		DrlCurPrint(str, *((double*) p), 0);
		break;
	case DRL_CURK_T:
		DrlCurPrint(str, *((double*) p)*1e-3, 0);
		break;
	case DRL_CURM_T:
		DrlCurPrint(str, *((double*) p)*1e-6, 0);
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
		DrlErrMsg("%s: bad type %ld.\n", routine, (long)type);
		strcpy(str, "ERR");
		break;
	}

	return(str);
}




/*f----------------------------------------------------------------------
 * Prints in a formatted string <i> str</i> a pointer <i> p</i>
 * of type <i> type</i>. Returns <i> str</i>.
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> DrlVTypeOffsetVect</i>, even if it is the first element.}
 */

char*
DrlVTypePrintFmt(char *str, DVType type, void *p, char *fmt)
{
static	char	routine[] = "DrlVTypePrintFmt";
static	char	buf[64];

	str = (str == (char*)NULL ? buf : str);


	switch (type) {
	case DRL_DOUBLE_T:
		sprintf(str, fmt, *((double*) p));
		break;
	case DRL_FLOAT_T:
		sprintf(str, fmt,  *((float*) p));
		break;
	case DRL_INT_T:
		sprintf(str, fmt,  *((int*) p));
		break;
	case DRL_LONG_T:
		sprintf(str, fmt,  *((long*) p));
		break;
	case DRL_CHAR_T:
		sprintf(str, fmt,  *((char*) p));
		break;
	case DRL_STRING_T:
		sprintf(str, fmt,  ((char**)p)[0]);
		break;
	case DRL_TDATE_T:
		sprintf(str, fmt, DrlDDatePrint(NULL, *((DDate*) p)));
		break;
	case DRL_TDATEINTERVAL_T:
		sprintf(str, fmt, DrlDIntervalPrint(NULL, ((DInterval*) p)));
		break;
	case DRL_TDAYCOUNT_T:
		sprintf(str, fmt,
		   DrlDDayCountPrint(NULL, *((DDayCount*) p)));
		break;

		/* derived Types */

	case DRL_PERCENT_T:
		sprintf(str, fmt, *((double*) p) * 1e2);
		break;

		/* LIL types */

	case DRL_CHAR_BLOCK_L:
		sprintf(str, fmt,  ((char*) p));
		break;

	default:
		DrlErrMsg("%s: bad type %ld.\n", routine, (long)type);
		strcpy(str, "ERR");
		break;
	}

	return(str);
}




/*f----------------------------------------------------------------------
 * Prints in a formatted string <i> str</i> a pointer <i> p</i>
 * of type <i> type</i>. Returns <i> str</i>.
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> DrlVTypeOffsetVect</i>, even if it is the first element.}
 */

char*
DrlVTypePrintSmart(char *str, DVType type, void *p)
{
static	char	routine[] = "DrlVTypePrintSmart";
static	char	buf[64];

	str = (str == (char*)NULL ? buf : str);


	switch (type) {
	case DRL_POINTER_T:
		sprintf(str, "%p", *((void**) p));
		break;
	case DRL_DOUBLE_T:
		if (fabs(*((double*) p)) < 1e2) {
			sprintf(str, "%12.6f", *((double*) p));
		} else  {
			sprintf(str, "%12.0f", *((double*) p));
		}
		break;
	case DRL_FLOAT_T:
		sprintf(str, "%f",  *((float*) p));
		break;
	case DRL_INT_T:
		sprintf(str, "%4d",  *((int*) p));
		break;
	case DRL_LONG_T:
		sprintf(str, "%4ld",  *((long*) p));
		break;
	case DRL_CHAR_T:
		sprintf(str, "%c",  *((char*) p));
		break;
	case DRL_STRING_T:
		sprintf(str, "%10s",  ((char**)p)[0]);
		break;
	case DRL_CHAR_ARRAY_T:
		sprintf(str, "%10s",  ((char*)p));
		break;
	case DRL_TDATE_T:
		sprintf(str, "%10s", DrlDDatePrint(NULL, *((DDate*) p)));
		break;
	case DRL_TDATEINTERVAL_T:
		sprintf(str, "%6s",
			DrlDIntervalPrint(NULL, ((DInterval*) p)));
		break;
	case DRL_TDAYCOUNT_T:
		sprintf(str, "%10s",
			DrlDDayCountPrint(NULL, *((DDayCount*) p)));
		break;

		/* derived TYpes */
	case DRL_PERCENT_T:
		sprintf(str, "%10.6f", *((double*) p) * 1e2);
		break;
	case DRL_CUR_T:
		sprintf(str, "%18s",
			DrlCurPrint(NULL, *((double*) p), 0));
		break;
	case DRL_CURK_T:
		sprintf(str, "%10s",
			DrlCurPrint(NULL, *((double*) p)*1e-3, 0));
		break;
	case DRL_CURM_T:
		sprintf(str, "%6s",
			DrlCurPrint(NULL, *((double*) p)*1e-6, 0));
		break;
	case DRL_BOOLEAN_T:
		sprintf(str, "%s", (*((int*) p) != 0 ? "T" : "F"));
		break;

		/**
		 ** Complex data structures
		 **/


	default:
		DrlErrMsg("%s: bad type %ld.\n", routine, (long)type);
		strcpy(str, "ERR");
		break;
	}

	return(str);
}


/*f----------------------------------------------------------------------
 * Returns the element with offset "offset" in an array
 * "p" of type "type".
 */

void*
DrlVTypeOffsetVect(void *p, int offset, DVType type)
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
		return(((DDate*)p) + offset);
	case DRL_TDATEINTERVAL_T:
		return(((DInterval*)p) + offset);
	case DRL_TDAYCOUNT_T:
		return(((DDayCount*)p) + offset);

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
		return &((char*) p)[DRL_WRAP_STR_IDX(offset+1)];

		/**
		 ** Complex data structures
		 **/


	default:
		DrlErrMsg("OffsetType: bad type %ld.\n", (long)type);
		return(NULL);
	}
}



/*f----------------------------------------------------------------------
 * Returns the element with offset "nx" , "ny" in a matrix
 * "p" of type "type".
 */

void*
DrlVTypeOffsetMatrix(void *p, int nx, int ny, DVType type)
{
#undef	ELEM
#define	ELEM(t)	((void*)&((t**)p)[nx][ny])
	switch (type) {
	case DRL_DOUBLE_T:
	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
		return ELEM(double);
	case DRL_FLOAT_T:
		return ELEM(float);
	case DRL_INT_T:
	case DRL_BOOLEAN_T:
		return ELEM(int);
	case DRL_LONG_T:
		return ELEM(long);
	case DRL_CHAR_T:
		return ELEM(char);
	case DRL_STRING_T:
		return ELEM(char*);
	case DRL_TDATE_T:
		return ELEM(DDate);
	case DRL_TDATEINTERVAL_T:
		return ELEM(DInterval);
	case DRL_TDAYCOUNT_T:
		return ELEM(DDayCount);
	default:
		DrlErrMsg("OffsetType: bad type %ld.\n", (long)type);
		return(NULL);
	}
#undef	ELEM
}



/*f----------------------------------------------------------------------
 * Given two arrays of type "type",
 * copies the element \# "offset2" of "q" to the
 * element \# "offset1" of "p".\\
 * <i> Example:</i> 
 * \begin{verbatim}
 *     double  p[10], q[10];
 *     ...
 *     DrlVTypeLet((void*)p, 2, (void*)q, 4, DRL_DOUBLE_T);
 *                         / * i.e. p[2]=q[4] * /
 * \end{verbatim}
 */

int
DrlVTypeLet(void *p, int offset1, void *q, int offset2, DVType type)
{
static	char	routine[] = "DrlVTypeLet";
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
			DRL_WRAP_STR_BYTES);
		break;
	case DRL_CHAR_ARRAY_T:
		strcpy((char*)p, (char*)q);
		break;
	case DRL_TDATE_T:
		((DDate*)p)[offset1] = ((DDate*)q)[offset2];
		break;
	case DRL_TDATEINTERVAL_T:
		((DInterval*)p)[offset1] = ((DInterval*)q)[offset2];
		break;
	case DRL_TDAYCOUNT_T:
		((DDayCount*)p)[offset1] = ((DDayCount*)q)[offset2];
		break;
	default:
		DrlErrMsg("%s: bad type %ld.\n", routine, (long)type);
		return(1);
	}
	return(0);
}

/*f----------------------------------------------------------------------
 * Compares two elements <i> p</i> and <i> q</i> of type <i> type</i>
 * and returns an integer less
 * than, equal to, or greater than zero, depending on
 * whether <i> p</i> is lexicographically less than, equal to, or
 * greater than <i> q</i>.
 */

int
DrlVTypeCompare(const void *p, const void *q, DVType type)
{
static	char	routine[] = "DrlVTypeCompare";
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
		_COMPARE(DDate);
	case DRL_TDATEINTERVAL_T:
		DrlDIntervalToYears((DInterval*)p, &val1);
		DrlDIntervalToYears((DInterval*)q, &val2);
		return (IS_ALMOST_ZERO(val2 - val1) ? 0 :
			(val1 < val2 ? -1 : 1));
	case DRL_TDAYCOUNT_T:
		_COMPARE(DDayCount);
	default:
		DrlErrMsg("%s: bad type %ld.\n", routine, (long)type);
		return(0);
	}
#undef	_COMPARE
}


/*f----------------------------------------------------------------------
 * Returns the numerical "value" of a value pointed by <i> p</i>
 * of type <i> type</i>.
 * The argument <i> arg</i> is only neede in some cases:\\
 * <br>
 * <br> If the type is <i> DRL_TDATEINTERVAL_T</i>, then <i> arg</i>
 * must  be the address of a <i> DDate</i> that is used
 * as reference date to convert an interval to a double value.
 * (the routine <i> DrlDIntervalToYears</i> is used
 * with a day count 30/360).
 * <br>
 * Not available for all types.
 */

int
DrlVTypeNumValue(void *p, DVType type, void *arg, double *value)
{
static	char	routine[] = "DrlVTypeNumValue";

	switch (type) {
	case DRL_DOUBLE_T:
	case DRL_PERCENT_T:
	case DRL_CUR_T:
	case DRL_CURK_T:
	case DRL_CURM_T:
		*value = (double) *((double*)p);
		break;
	case DRL_FLOAT_T:
		*value = (double) *((float*)p);
		break;
	case DRL_INT_T:
	case DRL_BOOLEAN_T:
		*value = (double) *((int*)p);
		break;
	case DRL_LONG_T:
		*value = (double) *((long*)p);
		break;
	case DRL_CHAR_T:
		*value = (double) *((char*)p);
		break;
	case DRL_STRING_T:
		*value = (double) (strlen(*((char**)p)) != 0);
		break;
	case DRL_TDATE_T:
		*value = (double) (*((DDate*)p) > 0L);
		break;
	case DRL_TDATEINTERVAL_T:
		IF_FAILED_DONE( DrlDIntervalToYears(
			((DInterval*)p),
			value));
		break;
	case DRL_TDAYCOUNT_T:
	default:
		DrlErrMsg("%s: can't get numerical value for type %s.\n",
			routine, DrlVTypeName(type));
		return(FAILURE);
	}
	return(SUCCESS);
done:
	DrlErrMsg("%s: failed.\n", routine);
	return(FAILURE);
}

/*f-------------------------------------------------------------
 * Prints a vector <i> p</i> with <i> numItem</i> elements
 * of type <i> type</i>.
 */


int
DrlVTypeVectPrint(
	void *p,		/* (B) vector */
	int numItems,		/* (B) num elem in vector */
	DVType type,		/* (I) vector type */
	FILE *fp)
{
static	char	routine[] = "DrlVTypeVectPrint";
	int	status = FAILURE;
	int	idx;

	for (idx=0; idx<numItems; idx++) {
		DrlFPrintf(fp, " [%3d/%3d] %s\n", idx, numItems,
			DrlVTypePrint(NULL, type,
			    DrlVTypeOffsetVect(p, idx, type)));
	}

	return(SUCCESS);
}



/*f-------------------------------------------------------------
 * Adds an element <i> q</i> to a vector <i> p</i> of length <i> numItems</i>.
 * If <i> removeDoubleItems</i> is TRUE, the element is added only
 * if its value is not already present in the vector (always added if false).
 * If <i> numItemsMax</i> is not NULL, it points to the current maximum
 * size of <i> p</i>: i.e. as long as the number of elements if the vector
 * is less than <i> numItemsMax</i>,  no reallocation is performed
 * and the element simply added at the end of the vector.
 * (<i> *numItemsMax</i> and <i> *p</i> are unchanged
 * and <i> *numItems</i> is increased by at most one).
 * If the maximum size if attained, then <i> p</i> reallocated and its
 * maximum size is increased by <i> reallocSize</i>
 * and, on exit, <i> numItemsMax</i> is increased by <i> reallocSize</i>,
 * and <i> *p</i> is changed.
 */

int
DrlVTypeVectAdd(
	void **p,		/* (B) address of vector */
	int *numItems,		/* (B) address num elem in vector */
	int *numItemsMax,	/* (B) address max num elem in vect (or NULL) */
	int reallocSize,	/* (I) */
	void *q,		/* (I) element to add */
	int removeDoubleItems,	/* (I) TRUE=remove double items within tol */
	DVType type)		/* (I) vector type */
{
static	char	routine[] = "DrlVTypeVectAdd";
	int	status = FAILURE;

	int	idx, numItemsMaxNew;
	void	*pnew = NULL;


	/* Check if element already in series */
	if (removeDoubleItems) {
	    for (idx=0; idx<=*numItems-1; idx++) {
		if (DrlVTypeCompare(
			DrlVTypeOffsetVect(*p, idx, type),
			q, type) == 0)
				return(SUCCESS);
	    }
	}

	/* Check if reallocation needed */
	if ((numItemsMax == NULL) ||
	    (*numItems >= *numItemsMax)) {

		if (reallocSize < 1) reallocSize = 1;

		if (numItemsMax != NULL)
			numItemsMaxNew = *numItemsMax + reallocSize;
		else
			numItemsMaxNew = *numItems + 1;

		if ((pnew = DrlVTypeVectAlloc(numItemsMaxNew, type)) == NULL)
			goto done;

		for (idx=0; idx<=*numItems-1; idx++) {
			if (DrlVTypeLet(pnew, idx, *p, idx, type) != SUCCESS)
				goto done;
		}


		DrlVTypeVectFree(*p, *numItems, type);

		if (numItemsMax != NULL)
			*numItemsMax = numItemsMaxNew;
		*p = pnew;
	}

	/* Add new element */
	if (DrlVTypeLet(*p, *numItems, q, 0, type) != SUCCESS)
		goto done;
	*numItems += 1;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failure.\n", routine);
	}
	return(status);
}



/*
 * Need those for the next routine.
 */

static	DVType	_DrlVTypeVectSort_type;

static	int
_DrlVTypeVectSort_compare(const void *p, const void *q)
{
	return DrlVTypeCompare(p, q, _DrlVTypeVectSort_type);

}


/*f-------------------------------------------------------------
 * Sorts a vector <i> p</i> with <i> numItem</i> elements
 * of type <i> type</i> in ascending order.
 * If the flag <i> removeMultiple</i> is TRUE, removes multiple
 * instances of a same element, and <i> numItems</i>
 * may be changed on exit.
 */


int
DrlVTypeVectSort(
	void *p,		/* (B) vector */
	int *numItems,		/* (B) num elem in vector */
	DVType type,		/* (I) vector type */
	int removeMultiple)	/* (I) TRUE=remove multiple items */
{
static	char	routine[] = "DrlVTypeSort";
	int	status = FAILURE;
	size_t	nel;
	size_t	size;
	int	idx, idx2;

	/* Perform sorting */
	_DrlVTypeVectSort_type = type;
	size = (size_t) DrlVTypeSizeof(type);
	nel = (size_t) *numItems;

	DrlQSort((void*) p,
		(size_t) nel,
		(size_t) size, 
		_DrlVTypeVectSort_compare);

	/* Remove double items */
	if (removeMultiple) {
	    idx = 0;
	    while (idx < *numItems-1) {
		if (DrlVTypeCompare(
			    DrlVTypeOffsetVect(p, idx,   type),
			    DrlVTypeOffsetVect(p, idx+1, type),
			    type) == 0) {
		    for (idx2=idx+1; idx2<*numItems-1; idx2++) {
		        if (DrlVTypeLet(
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
		DrlErrMsg("%s: failure.\n", routine);
	}
	return(status);
}




/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */

int
DrlVTypeNameScan(DVType *type, char *s)
{
static	char	routine[] = "DrlVTypeNameScan";
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

	DrlErrMsg("%s: can't read type in `%s'.\n", routine, s);
	return(FAILURE);


#undef	CHECK
}


/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */

char*
DrlVTypeName(DVType type)
{
	switch (type) {
	case DRL_DOUBLE_T:		return("double");
	case DRL_FLOAT_T:		return("float");
	case DRL_INT_T:			return("int");
	case DRL_LONG_T:		return("long");
	case DRL_CHAR_T:		return("char");

	case DRL_STRING_T:		return("char*");
	case DRL_CHAR_ARRAY_T:		return("char[]");
	case DRL_TDATE_T:		return("DDate");
	case DRL_TDATEINTERVAL_T:	return("DInterval");
	case DRL_TDAYCOUNT_T:		return("DDayCount");

	case DRL_PERCENT_T:		return("percent");
	case DRL_CUR_T:			return("currency");
	case DRL_CURK_T:		return("currencyK");
	case DRL_CURM_T:		return("currencyM");
	case DRL_BOOLEAN_T:		return("DBoolean");

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

int
DrlVTypeCheckValidCType(DVType type)
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

int
DrlVTypeCheckValidLType(DVType type)
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

size_t
DrlVTypeSizeof(DVType type)
{
static	char	routine[] = "DrlVTypeSizeof";
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
		return sizeof(DDate);
	case DRL_TDATEINTERVAL_T:
		return sizeof(DInterval);
	case DRL_TDAYCOUNT_T:
		return sizeof(DDayCount);

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
		return sizeof(DInterval);
	case DRL_TDATE_L:
		return sizeof(long);
	default:
		DrlErrMsg("%s: bad type %ld.\n", routine, (long) type);
		return((size_t) 0);
	}
}


/*f-------------------------------------------------------------
 * Allocates a array of dimension <i> nDim</i> of type <i> type</i>
 * and sizes in each dimensions given by the values of 
 * the array <i> nSize</i> (an vector of length <i> nDim</i>).
 */

void*
DrlVTypeTensorAlloc(DVType type, int nDim, int *nSize)
{
static	char	routine[] = "DrlVTypeTensorAlloc";
	int	i;
	void	**ptr = NULL;

	if (nDim > 1) {
	    ptr = NEW_ARRAY(void*, nSize[0]);
	    if (ptr == NULL) return(NULL);
	    for (i=0; i<=nSize[0]-1; i++) {
		ptr[i] = DrlVTypeTensorAlloc(type, nDim-1, nSize+1);
		if (ptr[i] == NULL) {
		    return(NULL);
		}
	    }
	    return(ptr);
	} else if (nDim == 1) {
	    ptr = DrlVTypeVectAlloc(nSize[0], type);
	    if (ptr == NULL) {
		return(NULL);
	    }
	    return(ptr);
	} else {
	    DrlErrMsg("%s: bad dimension %d\n", routine, nDim);
	    return(NULL);
	}
}



/*f-------------------------------------------------------------
 * Frees a array of arbitrary dimension  allocated by
 * <i> DrlVTypeTensorAlloc</i>.
 */

int
DrlVTypeTensorFree(void *ptr, DVType type, int nDim, int *nSize)
{
static	char	routine[] = "DrlVTypeTensorAlloc";
	int	i;

	if (ptr == NULL) return(SUCCESS);
	if (nDim > 1) {
	    for (i=0; i<=nSize[0]-1; i++) {
		DrlVTypeTensorFree((void*) ((void**) ptr)[i],
			type, nDim-1, nSize+1);
	    }
	    FREE(ptr);
	    return(SUCCESS);
	} else if (nDim == 1) {
	    DrlVTypeVectFree(ptr, nSize[0], type);
	    return(SUCCESS);
	} else {
	    DrlErrMsg("%s: bad dimension %d\n", routine, nDim);
	    return(FAILURE);
	}
}


/*f-------------------------------------------------------------
 * Same a <i> DrlVTypeTensorAlloc</i>, but takes a variable
 * number of arguments instead of an array for the sizes
 * in each dimension.
 */

void*
DrlVTypeTensorAllocV(DVType type, int nDim,
		... /* int size1, ..., int sizeNDim */)
{
static	char	routine[] = "DrlVTypeTensorAllocV";
#define	NMAXDIM	32
	int	n, nSize[NMAXDIM];
	va_list	ap;

	va_start(ap, nDim);
	if (nDim > NMAXDIM) {
	    DrlErrMsg("%s: maximum dimension %d (got %d).\n",
		routine, nDim, NMAXDIM);
	    return(NULL);
	}
	for (n=0; n<=nDim-1; n++) {
		nSize[n] = va_arg(ap, int);
	}
	va_end(ap);

	return DrlVTypeTensorAlloc(type, nDim, nSize);
#undef	NMAXDIM
}



/*f-------------------------------------------------------------
 * Frees a array of arbitrary dimension created by
 * <i> DrlVTypeTensorAllocV</i>.
 */

int
DrlVTypeTensorFreeV(void *ptr, DVType type, int nDim,
		... /* int size1, ..., sizeNDim */)
{
static	char	routine[] = "DrlVTypeTensorFreeV";
#define	NMAXDIM	32
	int	n, nSize[NMAXDIM];
	va_list	ap;

	va_start(ap, nDim);
	if (nDim > NMAXDIM) {
	    DrlErrMsg("%s: maximum dimension %d (got %d).\n",
		routine, nDim, NMAXDIM);
	    return(-1);
	}
	for (n=0; n<=nDim-1; n++) {
		nSize[n] = va_arg(ap, int);
	}
	va_end(ap);

	return DrlVTypeTensorFree(ptr, type, nDim, nSize);
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
DrlVTypeTensorAccessElement(
	void *ptr, DVType type, int nDim, int *nSize,
	int *idx)		/* multi-index [0..nDim-1] */
{
static	char	routine[] = "DrlVTypeTensorAccessElement";

	if (ptr == NULL) return(NULL);
	if (nDim > 1) {
	    return DrlVTypeTensorAccessElement(
			(void*) ((void**) ptr)[idx[0]],
			type, nDim-1, nSize+1,
			idx+1);
	} else if (nDim == 1) {
	    return DrlVTypeOffsetVect(ptr, idx[0], type);
	} else {
	    DrlErrMsg("%s: bad dimension %d.\n", routine, nDim);
	    return(NULL);
	}
}



/*--------------------------------------------------------------
 */

int*
DrlVTypeTensorMultiIndexNew(int nDim, int *nSize)
{
	int	*idx;
	int	i;

	if ((idx = NEW_ARRAY(int, nDim)) == NULL)
		return(NULL);
	for (i=0; i<=nDim-1; i++) {
		idx[i] = 0;
	}
	return(idx);
}

int
DrlVTypeTensorMultiIndexNext(int *idx, int nDim, int *nSize)
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
DrlVTypeTensorOperScalar(
	void *ptr, DVType type, int nDim, int *nSize,
	char *operName,		/* (I) see details */
	double *value)		/* (B) input or output */
{
static	char	routine[] = "DrlVTypeTensorOperScalar";
	int	status = FAILURE;
	int	cnt, *idx = NULL;
	double	*dval;

	if (type != DRL_DOUBLE_T) {
	    DrlErrMsg("%s: type must be double.\n", routine);
	    goto done;
	}

#define	LOOPSTART	if ((idx = DrlVTypeTensorMultiIndexNew(nDim, nSize)) == NULL)\
			goto done; do {\
			if ((dval = (double*) DrlVTypeTensorAccessElement(\
			ptr, type, nDim, nSize, idx)) == NULL) goto done;
#define	LOOPEND		} while (DrlVTypeTensorMultiIndexNext(idx, nDim, nSize) != FALSE);

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
		DrlErrMsg("%s: unknown operation %s.\n", routine, operName);
		goto done;
	}

	status = SUCCESS;
done:
	return(status);
}

/*f-------------------------------------------------------------
 * Same a <i> DrlVTypeTensorOperScalar</i>, but takes a variable
 * number of arguments instead of an array for the sizes
 * in each dimension.
 */

int
DrlVTypeTensorOperScalarV(
	void *ptr,		/* (I) tensor */
	DVType type,		/* (I) tensor elements type */
	char *operName,		/* (I) see details */
	double *value,		/* (B) input or output */
	int nDim,		/* (I) number of dimensions */
	... /* int size1, ..., int sizeNDim */)
{
static	char	routine[] = "DrlVTypeTensorOperScalarV";
#define	NMAXDIM	32
	int	n, nSize[NMAXDIM];
	va_list	ap;

	va_start(ap, nDim);
	if (nDim > NMAXDIM) {
	    DrlErrMsg("%s: maximum dimension %d (got %d).\n",
		routine, nDim, NMAXDIM);
	    return(FAILURE);
	}
	for (n=0; n<=nDim-1; n++) {
		nSize[n] = va_arg(ap, int);
	}
	va_end(ap);

	return DrlVTypeTensorOperScalar(
		ptr, type, nDim, nSize, operName, value);
#undef	NMAXDIM
}


/*f-------------------------------------------------------------
 */

int
DrlVTypeTensorFpWrite(
	void *ptr, DVType type, int nDim, int *nSize,
	FILE *fp)		/* (I) */
{
static	char	routine[] = "DrlVTypeTensorFpWrite";
	int	status = FAILURE;
	int	i, *idx = NULL;
	void	*element;


	DrlFPrintf(fp, "%d\n", nDim);
	for (i=0; i<=nDim-1; i++) {
		DrlFPrintf(fp, "\t%d", nSize[i]);
	}
	DrlFPrintf(fp, "\n");

	if ((idx = DrlVTypeTensorMultiIndexNew(nDim, nSize)) == NULL)
		goto done;
	do {
		if ((element = DrlVTypeTensorAccessElement(ptr,
			type, nDim, nSize, idx)) == NULL) goto done;

		DrlFPrintf(fp, "\t%s",
			DrlVTypePrint(NULL, type, element));

		DrlFPrintf(fp, "\t#");
		for (i=0; i<=nDim-1; i++) {
			DrlFPrintf(fp, " %3d", idx[i]);
		}

		DrlFPrintf(fp, "\n");


	} while (DrlVTypeTensorMultiIndexNext(idx, nDim, nSize) != FALSE);


	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


