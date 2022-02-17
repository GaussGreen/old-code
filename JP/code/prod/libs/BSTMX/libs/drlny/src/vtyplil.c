/************************************************************************
 * Module:	DRL
 * Submodule:	VTYPE
 * File:	
 * Function:	Variable Type
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <float.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "ldate.h"
#include "date_sup.h"
#include "convert.h"
#include "macros.h"
#include "date_sup.h"

#include "drlstr.h"	/* DrlCurrentDateTimePrint */
#include "drlio.h"	/* DrlFPrintf */
#include "drlvtype.h"	/* Prototype consistency */
#include "drlproc.h"		/* DrlStrerror() */

static	char*	strReplaceNonPrintChar(char *s);

#define	CASTSTR(ptr, nidx)	(&(((char*)ptr)[WRAP_STR_IDX(nidx+1)]))
#define	CASTTYPE(ptr, nidx, type)	(((type *) ptr)[nidx+1])

#define	__DEBUG__
#undef	__DEBUG__

#ifdef	_WINDLL
# undef __DEBUG__
#endif


/*f----------------------------------------------------------------------
 * Returns the C type corresponding to a LIL type "varType"
 * (or 0 if fails).
 */

DLL_EXPORT(int)
DrlLilToCType(TVType varType)		/* (I) variable type */
{
static	char	routine[] = "DrlLilToCType";

	switch (varType) {
	case DRL_FLOAT_L:		return(DRL_DOUBLE_T);
	case DRL_PERCENT_L:		return(DRL_DOUBLE_T);
	case DRL_BPOINT_L:		return(DRL_DOUBLE_T);
	case DRL_INT_L:			return(DRL_INT_T);
	case DRL_CHAR_BLOCK_L:		return(DRL_STRING_T);
	case DRL_CHAR_L:		return(DRL_CHAR_T);
	case DRL_TDATEINTERVAL_L:	return(DRL_TDATEINTERVAL_T);
	case DRL_TDATE_L:		return(DRL_TDATE_T);
	default:
		GtoErrMsg("%s: bad type %ld.\n", routine, (long)varType);
		return(DRL_NULL_T);
	}
}


/*f----------------------------------------------------------------------
 * Returns the length of a LIL vector "lptr"
 * of type "varType". "argName" is used for error messages. \\
 */

DLL_EXPORT(int)
DrlLilVectSize(TVType varType, void *lptr, char *argName)
{
static	char	routine[] = "DrlLilVectSize";

	/* Get input vector length */
	if (lptr == NULL) {
	    return(0);
	}

	switch (varType) {
	case DRL_FLOAT_L:
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		return ARGSIZE((double*) lptr);
	case DRL_INT_L:
		return ARGSIZE((long*)   lptr);
	case DRL_TDATE_L:
		return ARGSIZE((long*)   lptr);
	case DRL_TDATEINTERVAL_L:
		return ARGSIZE((double*) lptr);
	case DRL_CHAR_BLOCK_L:
	case DRL_CHAR_L:
		/*return ARGSIZE((char*)   lptr);*/
		return (int)(lptr == NULL ? 0 :
			(unsigned char)((char*) lptr)[0]);

	case DRL_DOUBLE_T:
	case DRL_FLOAT_T:
	case DRL_INT_T:
	case DRL_LONG_T:
	case DRL_CHAR_T:
	case DRL_STRING_T:
	case DRL_TDATE_T:
	case DRL_TDATEINTERVAL_T:
	case DRL_TDAYCOUNT_T:
		GtoErrMsg("%s: [%s] non LIL type (%s).\n",
			routine, argName, DrlVTypeName(varType));
		return(-1);

	default:
		GtoErrMsg("%s: [%s] bad type %ld.\n",
			routine, argName, (long) varType);
		return(-1);
	}
}

/*f----------------------------------------------------------------------
 * Returns the <i> used</i> length of a LIL vector "lptr"
 * of type "varType". "argName" is used for error messages. \\
 * <b> Remark:</b> the behaviour is different from the <i> DrlLilVectSize</i>
 * function in the sens that it returns the size occupied by
 * valid element (e.g. dates).
 */

DLL_EXPORT(int)
DrlLilVectActiveSize(TVType varType, void *lptr, char *argName)
{
static	char	routine[] = "DrlLilVectActiveSize";
	register int	i;
	int	sz=0;

	/*
	 * Get input vector length
	 */
	if (lptr == NULL) {
	    GtoErrMsg("%s: [%s] NULL pointer.\n", routine, argName);
	    return(-1);
	}


	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		sz = ARGSIZE((double*) lptr);
		break;
	case DRL_INT_L: /* long -> int */
		sz = ARGSIZE((long*) lptr);
		break;

	case DRL_TDATE_L: /* TDate: check > 1 */
		sz = ARGSIZE((long*) lptr);
		for (i=0; i<=sz-1; i++) {
			if (((long*) lptr)[i+1] <= 1L) break;
		}
		sz = i;
		break;

	case DRL_TDATEINTERVAL_L: /* TDateInterval */
		sz = ARGSIZE((double*) lptr);
		break;

	case DRL_CHAR_BLOCK_L:
	case DRL_CHAR_L:
		/*sz = ARGSIZE((char*) lptr);*/
		sz  = (int)(lptr == NULL ? 0 :
			(unsigned char)((char*) lptr)[0]);
		break;

	case DRL_DOUBLE_T:
	case DRL_FLOAT_T:
	case DRL_INT_T:
	case DRL_LONG_T:
	case DRL_CHAR_T:
	case DRL_STRING_T:
	case DRL_TDATE_T:
	case DRL_TDATEINTERVAL_T:
	case DRL_TDAYCOUNT_T:
		GtoErrMsg("%s: [%s] non LIL type (%s)\n",
			routine, argName, DrlVTypeName(varType));
		return(-1);

	default:
		GtoErrMsg("%s: [%s] bad type %ld.\n",
			routine, argName, (long) varType);
		return(-1);
	}

	/*
	 * Check length OK
	 */
	if (sz <= 0) {
	    GtoErrMsg("%s: argument `%s' bad array length %d (type %s)\n",
		routine, argName, sz, DrlVTypeName(varType));
	    return(sz);
	}

	return(sz);
}


/*f----------------------------------------------------------------------
 * Copies and converts the values of a LIL vector "lptr"
 * of type "varType" a  to C vector "cptr" of length "n"
 * (which should be of the proper C equivalent type).
 * and proper length). 
 * <b> WARNING:</b> No check is made on the arrays length.\\
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLilVectGet(
	void *lptr,		/* (I) input LIL vector */
	TVType varType,	/* (I) variable type */
	int loffset,		/* (I) offset on LIL vector */
	int lskip,		/* (I) interval bet. LIL values */
	int n,			/* (I) size of C vector */
	void *cptr)		/* (O) output C vector */
{
static	char	routine[] = "DrlLilVectGet";
	register int i;
	double	dVal;

#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL

#define	FORALL	for (i=0; i<=n-1; i++)
#define	LSTR	(&(((char*)lptr)[WRAP_STR_IDX(loffset+lskip*i+1)]))
#define	LVAL(t)	(((t) lptr)[loffset+lskip*i+1])
#define	CVAL(t)	(((t) cptr)[i])


	/*
	 * copy input in output (possibly parse)
	 */
	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		FORALL {
			CVAL(double*) = LVAL(double*);
		}
		break;

	case DRL_INT_L: /* long */
		FORALL {
			CVAL(int*) = (int) LVAL(long*);
		}
		break;

	case DRL_CHAR_BLOCK_L: /* char string */
		FORALL {
			strncpy(CVAL(char**), LSTR, 64);
		}
		break;

	case DRL_CHAR_L: /* char */
		FORALL {
			CVAL(char*) = LSTR[0];
		}
		break;

	case DRL_TDATEINTERVAL_L: /* LIL: double <-> C: TDateInterval */
		FORALL {
			dVal = LVAL(double*);
			if (!DoubleIsFinite(dVal)) {
			    GtoErrMsg("%s: bad TDateInterval value.\n",
				routine);
				return(5);
			}
			if (GtoYearsToDateInterval(dVal, 
				((TDateInterval*) cptr)+i) != SUCCESS)
					return(FAILURE);
		}
		break;

	case DRL_TDATE_L: /* TDate */
		FORALL {
			CVAL(TDate*) = LVAL(TDate*);
		}
		break;


	default:
		GtoErrMsg("%s: bad type %ld\n", routine, (long) varType);
		return(1);
	}
	return(0);
#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL
}


/*f----------------------------------------------------------------------
 * Copies and converts values to a LIL vector "lptr"
 * of type "varType" and length "n"
 * from C vector "cptr" (which should be of the proper C equivalent type
 * and proper length).
 * <b> WARNING:</b> No check is made on the arrays length.\\
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLilVectPut(
	void *lptr,		/* (O) input wrap vector */
	TVType varType,		/* (I) variable type (see header) */
	int loffset,		/* (I) offset on LIL vector */
	int lskip,		/* (I) interval bet. LIL values */
	int n,			/* (I) size */
	void *cptr)		/* (I) output C vector */
{
static	char		routine[] = "DrlLilVectPut";
	int		status = FAILURE;

	register int	i;
	double		dVal;
	int		iVal;

#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL

#define	FORALL	for (i=0; i<=n-1; i++)
#define	LSTR	(&((char*) lptr)[WRAP_STR_IDX(loffset+lskip*i+1)])
#define	LVAL(t)	(((t) lptr)[loffset+lskip*i+1])
#define	CVAL(t)	(((t) cptr)[i])

	/*
	 * copy input in output (possibly parse)
	 */
	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
	    FORALL {
		dVal = CVAL(double*);
		if (!DoubleIsFinite(dVal)) {
			GtoErrMsg("%s: double not finite.\n", routine);
			goto done;
		}
		LVAL(double*) = dVal;
	    }
	    break;

	case DRL_INT_L: /* long */
	    FORALL {
		iVal = CVAL(int*);
		LVAL(long*) = (long) iVal;
	    }
	    break;

	case DRL_CHAR_BLOCK_L: /* char string: the input should be a StringW type */
	    FORALL {
		strncpy(LSTR, CVAL(char**), 64);
		strReplaceNonPrintChar(LSTR);
	    }
	    break;

	case DRL_CHAR_L: /* char */
	    FORALL {
		LSTR[0] = CVAL(char*);
	    }
	    break;

	case DRL_TDATEINTERVAL_L: /* TDateInterval */
	    FORALL {
		if (GtoDateIntervalToYears(&(CVAL(TDateInterval*)),
		    &LVAL(double*)) != SUCCESS)
			goto done;
	    }
	    break;

	case DRL_TDATE_L: /* TDate */
	    FORALL {
		LVAL(TDate*) = CVAL(TDate*);
	    }
	    break;

	default:
		GtoErrMsg("%s: bad type %ld.\n", routine, (long) varType);
		goto done;
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL
}


/*-----------------------------------------------------------------------
 * For an array in Lotus: goes by lines first
 *	-----------> 		----------------->
 *	|     Y			0		NY-1
 *	|			----------------->
 *	| X			NY+1		2*NY-1
 *	|			  	....
 *	|
 *	V
 *
 * Element of (X=i, Y=j) is located in the
 * wrapper array at (j + NY*i).
 *
 */



/*f----------------------------------------------------------------------
 * Copies and converts the values of a LIL range (matrix) "lptr"
 * of type "varType" and size "n1" $x$ "n2"
 * to C matrix "cptr" (which should be of the proper C equivalent type
 * and proper size).
 * <b> WARNING:</b> No check is made on the arrays length.\\
 * Returns 0 iff successful.
 */


DLL_EXPORT(int)
DrlLilMatrGet(
	void *lptr,		/* (I) input wrap range */
	TVType varType,		/* (I) variable type */
	int n1,			/* (I) size */
	int n2,			/* (I) size */
	void *cptr)		/* (O) output C matrix */
{
static	char	routine[] = "DrlLilMatrGet";
	register int i, j;
	double	dVal;

#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL

#define	FORALL	for (i=0; i<=n1-1; i++) for (j=0; j<=n2-1; j++)
#define	LSTR	(&((char*)lptr)[WRAP_STR_IDX(n2*i+j+1)])
#define	LVAL(t)	(((t) lptr)[n2*i+j+1])
#define	CVAL(t)	(((t *) cptr)[i][j])


	/*
	 * copy input in output (possibly parse)
	 */
	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		FORALL {
			CVAL(double*) = LVAL(double*);
		}
		break;

	case DRL_INT_L: /* long */
		FORALL {
			CVAL(int*) = (int) LVAL(long*);
		}
		break;

	case DRL_CHAR_BLOCK_L: /* char string */
		FORALL {
			strncpy(CVAL(char**), LSTR, 64);
		}
		break;

	case DRL_CHAR_L: /* char */
		FORALL {
			CVAL(char*) = LSTR[0];
		}
		break;

	case DRL_TDATEINTERVAL_L: /* TDateInterval */
		FORALL {
			dVal = LVAL(double*);

			if (!DoubleIsFinite(dVal)) {
			    GtoErrMsg("%s: bad TDateInterval value.\n",
				routine);
				return(FAILURE);
			}
			if (GtoYearsToDateInterval(dVal, 
				&((TDateInterval**) cptr)[i][j]) != SUCCESS)
					return(FAILURE);
		}
		break;

	case DRL_TDATE_L: /* TDate */
		FORALL {
			CVAL(TDate*) = LVAL(TDate*);
		}
		break;


	default:
		GtoErrMsg("%s: bad type %ld.\n", routine, (long)varType);
		return(1);
	}
	return(0);
#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL
}


/*f----------------------------------------------------------------------
 * Copies and converts the values to a LIL range (matrix) "lptr"
 * of type "varType" and size "n1" $x$ "n2"
 * from C matrix "cptr" (which should be of the proper C equivalent type
 * and proper size).
 * <b> WARNING:</b> No check is made on the arrays length.\\
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLilMatrPut(
	void *lptr,		/* (O) input wrap vector */
	TVType varType,	/* (I) variable type (see header) */
	int n1,			/* (I) size */
	int n2,			/* (I) size */
	void *cptr)		/* (I) output C vector */
{
static	char		routine[] = "DrlLilMatrPut";
	register int	i, j;
	int		status = 1;
	double		dVal;
	int		iVal;

#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL

#define	FORALL	for (i=0; i<=n1-1; i++) for (j=0; j<=n2-1; j++)
#define	LSTR	(&((char*)lptr)[WRAP_STR_IDX(n2*i+j+1)])
#define	LVAL(t)	(((t) lptr)[n2*i+j+1])
#define	CVAL(t)	(((t *) cptr)[i][j])

	/*
	 * copy input in output (possibly parse)
	 */
	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
	    FORALL {
		dVal = CVAL(double*);
		if (!DoubleIsFinite(dVal)) {
			GtoErrMsg("%s: double not finite.\n", routine);
			goto done;
		}
		LVAL(double*) = dVal;
	    }
	    break;

	case DRL_INT_L: /* long */
	    FORALL {
		iVal = CVAL(int*);
		LVAL(long*) = (long) iVal;
	    }
	    break;

	case DRL_CHAR_BLOCK_L: /* char string: the input should be a StringW type */
	    FORALL {
		strncpy(LSTR, CVAL(char**), 64);
		strReplaceNonPrintChar(LSTR);
	    }
	    break;

	case DRL_CHAR_L: /* char */
	    FORALL {
		LSTR[0] = CVAL(char*);
	    }
	    break;

	case DRL_TDATEINTERVAL_L: /* TDateInterval */
	    FORALL {
		if (GtoDateIntervalToYears(&(CVAL(TDateInterval*)),
		    &LVAL(double*)) != SUCCESS)
			goto done;
	    }
	    break;

	case DRL_TDATE_L: /* TDate */
	    FORALL {
		LVAL(TDate*) = CVAL(TDate*);
	    }
	    break;

	default:
		GtoErrMsg("%s: bad type %ld.\n", routine, (long)varType);
		goto done;
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
#undef	FORALL
#undef	LVAL
#undef	LSTR
#undef	CVAL
}


/*f----------------------------------------------------------------------
 * Generic routine to convert LIL vectors to C variables.
 */

DLL_EXPORT(int)
DrlLilStructGet(int iFlag,
	/* DRL_LILVAR_L, 
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * DRL_LILVECT_L, int *nItems, intMinItems, int nMaxItems,
	 *	char *name,  int type,  int alloc,  void *lvar,  void *cvar, 
	 * DRL_LILMATR_L, int nx, int ny,
	 *	char *name,  int type,  int alloc,  void *lvar,  void *cvar, 
	 * DRL_LILVECTARRAY_L, int *nItems, intMinItems, intMaxItems, int nVar,
	 *	char *name1, int type1, int alloc1, void* lvar1, void *cvar1,
	 *	...
	 *	char *nameN, int typeN, int allocN, void* lvarN, void *cvarN,
	 * DRL_LILVECTRANGE_L, int *nItems, intMinItems, intMaxItems, int nVar,
	 *	int type, void* lvar,
	 *	char *name1, int alloc1, void *cvar1,
	 *	...
	 *	char *nameN, int allocN, void *cvarN,
	 * 0)
	 */
	...)
{
static	char		routine[] = "DrlLilStructGet";
#define	MAXARG		32	/* maximum # of args */
	va_list		ap;
	int		errCode = 1,
			status = FAILURE,
			i,
			*nItems, nMinItems, nMaxItems, nVar,
			nx, ny, sz, sz1,
			allocType[MAXARG];
	TVType		varType[MAXARG],
			strType;
	char		*argName[MAXARG];
	void		*vp,
			*lptr[MAXARG],
			*cptr[MAXARG];


#undef	CHECKSZ
#define	CHECKSZ(sz0)	{if (((sz0) < nMinItems) || ((sz0) > nMaxItems) ||\
		(nMinItems > nMaxItems)) {\
		GtoErrMsg("%s: [%s] expected %d <= #items <= %d, got %d\n",\
		routine, argName[0], nMinItems, nMaxItems, (sz0)); goto done;}}


	/*
	 *
	 */
	va_start(ap, iFlag);

	while ((strType = (TVType) va_arg(ap, TVType)) != 0) {
	    switch (strType) {
	    case DRL_LILVAR_L:
		/************************************************
		 * Single variable
		 ************************************************/
		argName[0]   = (char*)   va_arg(ap, char*);
		varType[0]   = (TVType)  va_arg(ap, TVType);
		lptr[0]      = (void*)   va_arg(ap, void*);
		cptr[0]      = (void*)   va_arg(ap, void*);

		/* get active size */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		if (sz <= 0) {
		    GtoErrMsg("%s: [%s] bad range size\n",
			routine, argName[0]);
		    goto done;
		}

		/* copy vector if non NULL */
		if (cptr[0] != NULL) {
		    errCode = DrlLilVectGet(lptr[0], varType[0],
			0, 1, 1, cptr[0]);
		    if (errCode != 0) goto done;
		}

#ifdef	__DEBUG__
		GtoErrMsg("%s: read `%s', type=%s\n",
			routine, argName[0], DrlVTypeName(varType[0]));
		GtoErrMsg("\t%s\n",
			    VTypePrint(NULL, DrlLilToCType(varType[0]), 
				cptr[0]));
#endif
		break;

	    case DRL_LILVECT_L:
		/************************************************
		 * Single vector
		 ************************************************/
		nItems       = (int*)     va_arg(ap, int*);
		nMinItems    = (int)      va_arg(ap, int);
		nMaxItems    = (int)      va_arg(ap, int);
		argName[0]   = (char*)    va_arg(ap, char*);
		varType[0]   = (TVType) va_arg(ap, TVType);
		allocType[0] = (int)      va_arg(ap, int);
		lptr[0]      = (void*)    va_arg(ap, void*);
		cptr[0]      = (void*)    va_arg(ap, void*);

		/* get active size */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		sz = MIN(sz, nMaxItems);
		if (sz <= 0) goto done;
		CHECKSZ(sz);
		*nItems = sz;


		/* allocate memory if necessary */
		if (allocType[0] == TRUE) {
		    vp = DrlVTypeVectAlloc(sz, DrlLilToCType(varType[0]));
		    if (vp == NULL) goto done;
		} else {
		    vp = (void*) cptr[0];
		}

		/* copy vector if non NULL */
		if (vp != NULL) {
		    if (DrlLilVectGet(lptr[0], varType[0],
			0, 1, sz, vp) != 0) goto done;
		}

#ifdef	__DEBUG__
		GtoErrMsg("%s: read `%s', type=%s, nItems=%d, allocType=%d\n",
			routine, argName[0], DrlVTypeName(varType[0]),
			*nItems, allocType[0]);
		for (i=0; i<=*nItems-1; i++) {
			GtoErrMsg("\t%d\t%s\n",
			    i, VTypePrint(NULL, DrlLilToCType(varType[0]), 
				VTypeOffsetVect(vp, i,
				DrlLilToCType(varType[0]))));
		}
#endif

		if (allocType[0] == TRUE) {
		    *((void**)cptr[0]) = vp;
		}


		break;

	    case DRL_LILMATR_L:
		/************************************************
		 * Matrix
		 ************************************************/
		nx           = (int)      va_arg(ap, int);
		ny           = (int)      va_arg(ap, int);
		argName[0]   = (char*)    va_arg(ap, char*);
		varType[0]   = (TVType) va_arg(ap, TVType);
		allocType[0] = (int)      va_arg(ap, int);
		lptr[0]      = (void*)    va_arg(ap, void*);
		cptr[0]      = (void*)    va_arg(ap, void*);


		/* get active size */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		if (sz != nx*ny) {
		    GtoErrMsg("%s: [%s] expected matrix size %dx%d=%d, "
			" got range size %d\n",
			routine, argName[0], nx, ny, nx*ny, sz);
		    goto done;
		}

		/* allocate memory if necessary */
		if (allocType[0] == TRUE) {
		    vp = DrlVTypeMatrAlloc(nx, ny, DrlLilToCType(varType[0]));
		    if (vp == NULL) goto done;
		} else {
		    vp = (void*) cptr[0];
		}

		/* copy vector if non NULL */
		if (vp != NULL) {
		    errCode = DrlLilMatrGet(lptr[0], varType[0], nx, ny, vp);
		    if (errCode != 0) goto done;
		}
#ifdef	__DEBUG__
		{
		int	j;
		GtoErrMsg("%s: read `%s', type=%s, "
			"nx=%d, ny=%d, allocType=%d\n",
			routine, argName[0], DrlVTypeName(varType[0]),
			nx, ny, allocType[0]);
		for (i=0; i<=nx-1; i++) 
		for (j=0; j<=ny-1; j++) {
			GtoErrMsg("\t%d\t%d\t%s\n",
			    i, j, VTypePrint(NULL, DrlLilToCType(varType[0]), 
				VTypeOffsetMatrix(vp, i, j,
				DrlLilToCType(varType[0]))));
		}
		}
#endif
		if (allocType[0] == TRUE) {
		    *((void**)cptr[0]) = vp;
		}

		break;

	    case DRL_LILVECTARRAY_L:
		/************************************************
		 * Array of ranges of same length and different types
		 ************************************************/
		nItems       = (int*) va_arg(ap, int*);
		nMinItems    = (int)  va_arg(ap, int);
		nMaxItems    = (int)  va_arg(ap, int);
		nVar         = (int)  va_arg(ap, int);
		for (i=0; i<=nVar-1; i++) {
		    argName[i]   = (char*)     va_arg(ap, char*);
		    varType[i]   = (TVType)  va_arg(ap, TVType);
		    allocType[i] = (int)       va_arg(ap, int);
		    lptr[i]      = (void*)     va_arg(ap, void*);
		    cptr[i]      = (void*)     va_arg(ap, void*);
		}

		/* get minimal active size */
		if (nVar <= 0) break;
		sz = 1000;
		for (i=0; i<=nVar-1; i++) {
		    sz1 = DrlLilVectActiveSize(varType[i], lptr[i], argName[i]);
		    sz = MIN(sz, sz1);
		}
		*nItems = sz;
		CHECKSZ(sz);
		if (sz <= 0) goto done;

		sz = MIN(sz, nMaxItems);
		*nItems = sz;
		if (sz <= 0) goto done;

		/* allocate memory if necessary */
		for (i=0; i<=nVar-1; i++) {
		    /* allocate memory if necessary */
		    if (allocType[i] == TRUE) {
			vp = DrlVTypeVectAlloc(sz, DrlLilToCType(varType[i]));
			if (vp == NULL) goto done;
		    } else {
			vp = (void*) cptr[i];
		    }

		    /* copy vector if non NULL */
		    /*$$$if (cptr[i] != NULL) {*/
		    if (vp != NULL) {
		        if (DrlLilVectGet(lptr[i], varType[i],
				0, 1, sz, vp) != 0) goto done;
		    }

		    if (allocType[i] == TRUE) {
			*((void**)cptr[i]) = vp;
		    }
		}
		break;


	    case DRL_LILVECTRANGE_L:
		/************************************************
		 * Single range representing a set of vertical
		 * vectors of same type (and same length!).
		 ************************************************/
		nItems       = (int*)     va_arg(ap, int*);
		nMinItems    = (int)      va_arg(ap, int);
		nMaxItems    = (int)      va_arg(ap, int);
		nVar         = (int)      va_arg(ap, int);
		varType[0]   = (TVType) va_arg(ap, TVType);
		lptr[0]      = (void*)    va_arg(ap, void*);
		for (i=0; i<=nVar-1; i++) {
		    argName[i]   = (char*) va_arg(ap, char*);
		    allocType[i] = (int)   va_arg(ap, int);
		    cptr[i]      = (void*) va_arg(ap, void*);
		}

		/* get minimal active size */
		if (nVar <= 0) break;

		/* check len consistent */
		if ((sz = DrlLilVectSize(varType[0], lptr[0], argName[0])) <= 0) {
		    GtoErrMsg("%s: [%s] LIL array of length %d\n",
			routine, sz);
		    goto done;
		}
		if ((sz - nVar*((int) (sz/nVar))) != 0) {
		    GtoErrMsg("%s: [%s] expected range of width %d, "
		    "got range of total length %d.\n",
		    routine, argName[0], nVar, sz);
		}
		*nItems = sz = (int) (sz/nVar);
		CHECKSZ(sz);
		if (sz <= 0) goto done;


		/* copy LIL data in C vectors */
		for (i=0; i<=nVar-1; i++) {
		    /* allocate memory if necessary */
		    if (allocType[i] == TRUE) {
			vp = DrlVTypeVectAlloc(sz, DrlLilToCType(varType[0]));
			if (vp == NULL) goto done;
		    } else {
			vp = (void*) cptr[i];
		    }

		    /* copy vector if non NULL */
		    if (vp != NULL) {
			if (DrlLilVectGet(lptr[0], varType[0],
				i, nVar, sz, vp) != 0) goto done;
		    }

		    if (allocType[i] == TRUE) {
			*((void**)cptr[i]) = vp;
		    }
		}

		break;
	    default:
		GtoErrMsg("%s: bad arg type (%ld).\n",
			routine, ((long) strType));
		return(1);
	    }
	}
	va_end(ap);

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f----------------------------------------------------------------------
 * Generic routine to convert LIL vectors  from C variables.
 */

DLL_EXPORT(int)
DrlLilStructPut(int iFlag,
	/* DRL_LILVAR_L, 
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * DRL_LILVECT_L, int *nItems, intMinItems, int nMaxItems,
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * DRL_LILMATR_L, int nx, int ny,
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * DRL_LILVECTARRAY_L, int *nItems, intnMinItems, intMaxItems, int nVar,
	 *	char *name1, int type1, void* lvar1, void *cvar1,
	 *	...
	 *	char *nameN, int typeN, void* lvarN, void *cvarN,
	 * 0)
	 */
	...)
{
static	char		routine[] = "DrlLilStructPut";
	int		errCode = FAILURE;

#define	MAXARG	32	/* maximum # of args */
	va_list		ap;
	int		*nItems, nMinItems, nMaxItems, nx, ny, sz;
	TVType	varType[MAXARG],
			strType;
	char		*argName[MAXARG];
	void		*lptr[MAXARG],
			*cptr[MAXARG];

	/*
	 *
	 */
	va_start(ap, iFlag);

	while ((strType = (TVType) va_arg(ap, TVType)) != 0) {
	    switch (strType) {
	    case DRL_LILVAR_L:
		argName[0]   = (char*)    va_arg(ap, char*);
		varType[0]   = (TVType) va_arg(ap, TVType);
		lptr[0]      = (void*)    va_arg(ap, void*);
		cptr[0]      = (void*)    va_arg(ap, void*);

		/* get active size */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		if (sz <= 0) {
		    GtoErrMsg("%s: [%s] bad range size\n", routine, argName[0]);
		    goto done;
		}

		/* copy vector if non NULL */
		if (DrlLilVectPut(lptr[0], varType[0],
			0, 1, 1, cptr[0]) != 0) goto done;
		break;


	    case DRL_LILVECT_L:
		nItems       = (int*)     va_arg(ap, int*);
		nMinItems    = (int)      va_arg(ap, int);
		nMaxItems    = (int)      va_arg(ap, int);
		argName[0]   = (char*)    va_arg(ap, char*);
		varType[0]   = (TVType) va_arg(ap, TVType);
		lptr[0]      = (void*)    va_arg(ap, void*);
		cptr[0]      = (void*)    va_arg(ap, void*);

		/* get active size and check enough space */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		CHECKSZ(sz);
		*nItems = sz;

		/* copy vector if non NULL */
		if (DrlLilVectPut(lptr[0], varType[0],
			0, 1, sz, cptr[0]) != 0) goto done;

		break;

	    /*
	     * Matrix:   C -> W
	     */
	    case DRL_LILMATR_L:
		nx           = (int)      va_arg(ap, int);
		ny           = (int)      va_arg(ap, int);
		argName[0]   = (char*)    va_arg(ap, char*);
		varType[0]   = (TVType) va_arg(ap, TVType);
		lptr[0]      = (void*)    va_arg(ap, void*);
		cptr[0]      = (void*)    va_arg(ap, void*);


		/* get active size */
		sz = DrlLilVectActiveSize(varType[0], lptr[0], argName[0]);
		if (sz != nx*ny) {
		    GtoErrMsg("%s: [%s] expected matrix size %dx%d=%d, "
			" got range size %d\n",
			routine, argName[0], nx, ny, nx*ny, sz);
		    goto done;
		}

		/* copy matrix */
		errCode = DrlLilMatrPut(lptr[0], varType[0], nx, ny, cptr[0]);
		if (errCode != 0) goto done;

		break;

	    case DRL_LILVECTARRAY_L:
#ifdef	_SKIP
		nItems       = (int*) va_arg(ap, int*);
		nMaxItems    = (int)  va_arg(ap, int);
		nMaxItems    = (int)  va_arg(ap, int);
		nVar         = (int)  va_arg(ap, int);
		for (i=0; i<=nVar-1; i++) {
		    argName[i]   = (char*)    va_arg(ap, char*);
		    varType[i]   = (TVType) va_arg(ap, TVType);
		    allocType[i] = (int)      va_arg(ap, int);
		    lptr[i]      = (void*)    va_arg(ap, void*);
		    cptr[i]      = (void*)    va_arg(ap, void*);
		}

		/* get minimal active size */
		if (nVar <= 0) break;
		sz = 1000;
		for (i=0; i<=nVar-1; i++) {
		    sz1 = DrlLilVectActiveSize(varType[i], lptr[i], argName[i]);
		    sz = MIN(sz, sz1);
		}
		CHECKSZ(sz);
		*nItems = sz;
		if (sz <= 0) goto done;

		sz = MIN(sz, nMaxItems);
		*nItems = sz;
		if (sz <= 0) goto done;

		/* allocate memory if necessary */
		for (i=0; i<=nVar-1; i++) {
		    /* allocate memory if necessary */
		    if (allocType[i] == TRUE) {
			cptr[i] = *((void**)cptr[i]);
			cptr[i] = DrlVTypeVectAlloc(sz, DrlLilToCType(varType[i]));
			if (cptr[i] == NULL) goto done;
		    }

		    /* copy vector if non NULL */
		    if (cptr[i] != NULL) {
		        if (DrlLilVectGet(lptr[i], varType[i],
				0, 1, sz, cptr[i]) != 0) goto done;
		    }
		}
#endif
		GtoErrMsg("%s: [XXX] not implemented\n", routine);
		goto done;
		break;

	    case DRL_LILMATR2VECT_L:
		/************************************************
		 * Output a matrix as a series of vectors.
		 * For Excel
		 ************************************************/
		GtoErrMsg("%s: MATR2VECT not implemented\n", routine);
		goto done;
		break;
	    default:
		GtoErrMsg("%s: bad arg type (%ld).\n",
			routine, ((long) strType));
		return(1);
	    }
	}
	va_end(ap);

	/* made it through */
	errCode = 0;
done:
	if (errCode != 0) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(errCode);
}




/*f----------------------------------------------------------------------
 * Copies a LIL vector "rptr" of type "varType" to "lptr".
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLilVectCopy(
	void *lptr,		/* (I) output LIL vector */
	void *rptr,		/* (I) input LIL vector */
	TVType varType,	/* (I) variable type */
	int n)			/* (I) size of vectors */
{
static	char	routine[] = "DrlLilVectTypeCopy";
	register int i;

#undef	FORALL
#undef	LSTR
#undef	RSTR
#undef	LVAL
#undef	RVAL

#define	FORALL	for (i=0; i<=n-1; i++)
#define	LSTR	(&(((char*)lptr)[WRAP_STR_IDX(i+1)]))
#define	RSTR	(&(((char*)rptr)[WRAP_STR_IDX(i+1)]))
#define	LVAL(t)	(((t) lptr)[i+1])
#define	RVAL(t)	(((t) rptr)[i+1])


	/*
	 * copy input in output (possibly parse)
	 */
	switch (varType) {
	case DRL_FLOAT_L: /* double */
	case DRL_PERCENT_L:
	case DRL_BPOINT_L:
		FORALL {
			LVAL(double*) = RVAL(double*);
		}
		i = -1;
		LVAL(double*) = RVAL(double*);
		break;

	case DRL_LONG_L: /* long */
		FORALL {
			LVAL(long*) = RVAL(long*);
		}
		i = -1;
		LVAL(long*) = RVAL(long*);
		break;

	case DRL_CHAR_BLOCK_L: /* char string */
	case DRL_CHAR_L: /* char string */
		FORALL {
			strncpy(LSTR, RSTR, 64);
		}
		((char*) lptr)[0] = ((char*) rptr)[0];
		break;

	case DRL_TDATE_L: /* TDate */
		FORALL {
			LVAL(TDate*) = RVAL(TDate*);
		}
		i = -1;
		LVAL(long*) = RVAL(long*);
		break;


	default:
		GtoErrMsg("%s: bad type %ld.\n", routine, (long)varType);
		return(1);
	}
	return(0);

#undef	LSTR
#undef	RSTR
#undef	RVAL
#undef	RVAL

}



/*f----------------------------------------------------------------------
 * Print a LIL vector <i> lptr</i> of type <i> varType</i> to 
 * a file pointer <i> fp</i> (or to the error log if <i> fp</i> is NULL).
 * Returns 0 iff successful.
 */

DLL_EXPORT(int)
DrlLilVectLogging(
	FILE *fp,		/* (I) file to print (or NULL) */
	TVType varType,	/* (I) variable type */
	void *lptr,		/* (I) pointer */
	char  *argName)		/* (I) for debuging */
{
static	char	routine[] = "DrlLilVectLogging";
	int	status = FAILURE;
	int	sz, idx, idxN = 0;

#undef	FORALL
#undef	PRINTF
#define	FORALL(statement)	for (idx=0; idx<=sz-1; idx++) {statement;\
				if (++idxN == 5){idxN=0; PRINTF("\n",0);}}
#define	PRINTF(fmt, val)	{if (fp != NULL) fprintf(fp, fmt, val);\
				else GtoErrMsg(fmt, val);}

	/* get size */
	if ((sz = DrlLilVectSize(varType, lptr, argName)) < 0)
		goto done;

	/* print */
	PRINTF("%s: ", (argName != NULL ? argName : "[NO_NAME]"));
	/*PRINTF("%d\n", (int) sz);*/

	PRINTF("\n", 0);

	/* For sz == 0 (NULL pointer), print only the name and
	 * skip the array 
	 */
	if (sz > 0) {
		/* print */
		switch (varType) {
		case DRL_FLOAT_L:
		case DRL_PERCENT_L:
		case DRL_BPOINT_L:
			FORALL(PRINTF("\t%12.8f", CASTTYPE(lptr, idx, double)));
			break;
		case DRL_LONG_L:
			FORALL(PRINTF("\t%ld", CASTTYPE(lptr, idx, long)));
			break;
		case DRL_CHAR_BLOCK_L: /* char string */
			FORALL(PRINTF("\t%s", CASTSTR(lptr, idx)));
			break;
		case DRL_TDATE_L:
			FORALL(PRINTF("\t%s",
				GtoFormatDate(CASTTYPE(lptr, idx, TDate))));
			break;

		default:
			GtoErrMsg("%s: bad type %ld.\n", routine, (long)varType);
			goto done;
		}
	}

	if (idxN != 0) PRINTF("\n", 0)
	PRINTF("\t%%ARRAY_END%% %%ARG_END%%\n", 0)


	status = SUCCESS;
done:
	return(status);

#undef	FORALL
#undef	PRINTF

}

/*f----------------------------------------------------------------------
 * A variable argument version of <i> DrlLilVectLogging</i> that
 * writes an set of LIL vectors to a file "fnam"
 * (if "fnam" is NULL, writes to the error log).
 * Returns 0 iff successful.
 * <br><b> Example:</b>
 * \begin{verbatim}
 * int WrappeRoutineL(
 *         double *dblVal,
 *         long *dateVal,
 *         long *numVal)
 * {
 * ...
 * status = DrlLilVectLoggingFile("addin.log", "w", "WRAPPER",
 *      DRL_FLOAT_L,  (void*) dblVal,  "DOUBLE_INPUT",
 *      DRL_TDATE_L,  (void*) dateVal, "DATE_INPUT",
 *      DRL_LONG_L,   (void*) numVal,  "LONG_INPUT",
 *      0);
 * ...
 * }
 * \end{verbatim}
 */

DLL_EXPORT(int)
DrlLilVectLoggingFile(
	char *fnam,			/* (I) file name */
	char *mode,			/* (I) write mode */
	char *funcName,			/* (I) function name */
	/* TVType varType, void *lptr, char  *argName,
	 * ...
	 * TVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...)
{
static	char	routine[] = "DrlLilVectLoggingFile";
	int	status = FAILURE;
	FILE	*fp = NULL;
	va_list		ap;
	TVType	varType;
	void		*lptr;
	char		*argName;

	va_start(ap, funcName);

	if (fnam != NULL) {
	    if ((fp = fopen(fnam, mode)) == NULL) {
		GtoErrMsg("%s: `%s' (%s).\n", routine, fnam,
			DrlStrerror());
		goto done;
	    }
	}

	DrlFPrintf(fp, "COM: ---------------------------------"
		"------------------------------------\n");
	DrlFPrintf(fp, "COM: %s\n",
		(funcName != NULL ? funcName : "[NO_NAME]"));
	DrlFPrintf(fp, "COM: %s\n", DrlCurrentDateTimePrint(NULL));
	DrlFPrintf(fp, "COM: ---------------------------------"
		"------------------------------------\n");


	while ((varType = (TVType) va_arg(ap, TVType)) != 0) {

	    lptr    = (void*) va_arg(ap, void*);
	    argName = (char*) va_arg(ap, char*);

	    if (DrlLilVectLogging(fp, varType, lptr, argName) != SUCCESS)
		goto done;
	}

	status = SUCCESS;
done:
	va_end(ap);
	if (fp) fclose(fp);
	return(status);
}



/*f----------------------------------------------------------------------
 * Similar to DrlLilVectLoggingFile, put write to a given file pointer.
 */

DLL_EXPORT(int)
DrlLilVectLoggingFp(
	FILE *fp,			/* (I) file to print (or NULL) */
	char *funcName,			/* (I) function name */
	/* TVType varType, void *lptr, char  *argName,
	 * ...
	 * TVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...)
{
static	char	routine[] = "DrlLilVectLoggingFp";
	int	status = FAILURE;
	va_list		ap;
	TVType	varType;
	void		*lptr;
	char		*argName;

	va_start(ap, funcName);

	DrlFPrintf(fp, "COM: ---------------------------------"
		"------------------------------------\n");
	DrlFPrintf(fp, "COM: %s\n",
		(funcName != NULL ? funcName : "[NO_NAME]"));
	DrlFPrintf(fp, "COM: %s\n", DrlCurrentDateTimePrint(NULL));
	DrlFPrintf(fp, "COM: ---------------------------------"
		"------------------------------------\n");


	while ((varType = (TVType) va_arg(ap, TVType)) != 0) {

	    lptr    = (void*) va_arg(ap, void*);
	    argName = (char*) va_arg(ap, char*);

	    if (DrlLilVectLogging(fp, varType, lptr, argName) != SUCCESS)
		goto done;
	}

	status = SUCCESS;
done:
	va_end(ap);
	return(status);
}


/*-----------------------------------------------------------------------
 * Replaces all non printable characters by ' ' in the string "s".
 * Returns "s".
 */

static	char*
strReplaceNonPrintChar(char *s)
{
	char	*q;
	for (q = s; *q != '\0'; q++) {
		if (!isprint(*q)) *q = ' ';
	}
	return (s);
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

DLL_EXPORT(int)
DrlLilVectArrayFpRead(
	FILE *fp,		/* (I) file pointer */
	int numItems,		/* (I) num elements in each column */
	int numVect,		/* (I) num of columns */
	TVType *varType,	/* (I) array of variable types [0..numVect-1] */
	void ***lptr)		/* (I) array of void* pointers [0..numVect-1] */
{
static	char	routine[] = "DrlLilVectArrayFpRead";
	int	status = FAILURE;

	int	idxV, idx;
	void	*vptr;


	/* reset to NULL */
	for (idxV=0; idxV<=numVect-1; idxV++)
		*(lptr[idxV]) = (void*) NULL;

	/* malloc */
	for (idxV=0; idxV<=numVect-1; idxV++) {
		*(lptr[idxV]) = DrlVTypeVectAlloc(numItems, varType[idxV]);
		if (*(lptr[idxV]) == NULL) goto done;
	}



	/* read elements */
	for (idx=0; idx<=numItems-1; idx++) {
	    for (idxV=0; idxV<=numVect-1; idxV++) {

		/* access element #idx of array #idxV */
		vptr = DrlVTypeOffsetVect(*(lptr[idxV]), idx, varType[idxV]);
		if (vptr == NULL) {
			goto done;
		}

		/* read element */
		if (DrlFScanVType(fp, varType[idxV], vptr) != SUCCESS) {
		    GtoErrMsg("%s: can't read element #%d/%d "
			"in column # %d/%d "
			"(expecting type %s).\n", routine,
			idx+1, numItems, idxV+1, numVect,
			DrlVTypeName(varType[idxV]));
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
		    DrlVTypeVectFree(*(lptr[idxV]), numItems, varType[idxV]);
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f----------------------------------------------------------------------
 * A version of <i> DrlLilVectArrayFpRead</i> that allows the call
 * with variable number of arguments (one does not need to
 * set up arrays as in <i> DrlLilVectArrayFpRead</i>. \\
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
 *     status = DrlLilVectArrayFpReadV(
 *                    stdin,
 *                    5,
 *                    DRL_TDATE_L,  (void*) &dates,
 *                    DRL_DOUBLE_L, (void*) &notionals,
 *                    DRL_NULL_T);
 *     ...
 * \end{verbatim}
 */

DLL_EXPORT(int)
DrlLilVectArrayFpReadV(
	FILE *fp,			/* (I) file pointer */
	int numItems,			/* (I) num elements to be read */
	/* TVType varType, void *lptr,
	 * ...
	 * TVType varType, void *lptr,
	 * DRL_NULL_T (last argument MUST be DRL_NULL_T)
	 */
	...)
{
static	char	routine[] = "DrlLilVectArrayFpReadV";
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
			GtoErrMsg("%s: too many arguments (max %d).\n",
				routine, MAX_ARGS);
			goto done;
		}
	}
	va_end(ap);

	/* call routine */
	if (DrlLilVectArrayFpRead(
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
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
#undef	MAX_ARGS
}




