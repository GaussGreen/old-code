/************************************************************************
 * Module:	Q3
 * Submodule:
 * File: aliblog.c
 * Function:	Alib-style Regtest Logging
 * Author:	Interest Rates DR
 * Revision:	$Header$
 ************************************************************************/
#include "aliblog.h"
#include "q3.h"

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#if defined(WIN32) || defined(MSDOS) || defined(_WIN32)
#define vsnprintf       _vsnprintf
#endif

#ifndef	_WIN32
# include <errno.h>
#endif

static const char*	Q3Strerror()
{
#ifdef	_WIN32
	return _strerror(NULL);
#else
	return strerror(errno);
#endif
}


/*f----------------------------------------------------------------------
 * Print a LIL vector <i> lptr</i> of type <i> varType</i> to 
 * a file pointer <i> fp</i> (or to the error log if <i> fp</i> is NULL).
 * Returns 0 iff successful.
 */

int Q3LilVectLogging(
    FILE *fp,		        /* (I) file to print (or NULL) */
    DVType varType,	        /* (I) variable type */
    void *lptr,		        /* (I) pointer */
    char  *argName)	        /* (I) for debugging */
{

    static    char     routine[] = "Q3LilVectLogging";
    int	      status = FAILURE;

    int	sz, idx, idxN = 0;

#undef	FORALL
#undef	PRINTF
#define	FORALL(statement)	for (idx=0; idx<=sz-1; idx++) {statement;\
				if (++idxN == 5){idxN=0; PRINTF("\n",0);}}
#define	PRINTF(fmt, val)	{if (fp != NULL) fprintf(fp, fmt, val);\
				else Q3ErrMsg(fmt, val);}

    /* get size */
    if ((sz = Q3LilVectSize(varType, lptr, argName)) < 0) goto RETURN;

    /* print */
    PRINTF("%s: ", (argName != NULL ? argName : "[NO_NAME]"));
    /*PRINTF("%d\n", (int) sz);*/

    PRINTF("\n", 0);

    /* For sz == 0 (NULL pointer), print only the name and skip the array */
    if (sz > 0) {
	/* print */
        switch (varType) {
        case Q3_DOUBLE_L:
        case Q3_PERCENT_L:
        case Q3_BPOINT_L:
	    FORALL(PRINTF("\t%12.20f", CASTTYPE(lptr, idx, double)));
	    break;
	case Q3_LONG_L:
	    FORALL(PRINTF("\t%ld", CASTTYPE(lptr, idx, long)));
	    break;
        case Q3_CHAR_BLOCK_L: /* char string */
	    FORALL(PRINTF("\t%s", CASTSTR(lptr, idx)));
	    break;
        case Q3_TDATE_L:
	    FORALL(PRINTF("\t%s", 
		       Q3TDatePrint(NULL, CASTTYPE(lptr, idx, TDate))));
	    break;
        default:
	    Q3ErrMsg("%s: bad type %ld.\n", routine, (long)varType);
	    goto RETURN;
        }
    }

    if (idxN != 0) PRINTF("\n", 0)
    PRINTF("\t%%ARRAY_END%% %%ARG_END%%\n", 0)

    status = SUCCESS;
  done:
  RETURN:

    return(status);

#undef	FORALL
#undef	PRINTF
} /* Q3LilVectLogging */

/*f----------------------------------------------------------------------
 * A variable argument version of <i> Q3LilVectLogging</i> that
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
 * status = Q3LilVectLoggingFile("addin.log", "w", "WRAPPER",
 *      Q3_DOUBLE_L, (void*) dblVal,  "DOUBLE_INPUT",
 *      Q3_TDATE_L,  (void*) dateVal, "DATE_INPUT",
 *      Q3_LONG_L,   (void*) numVal,  "LONG_INPUT",
 *      0);
 * ...
 * }
 * \end{verbatim}
 */

int Q3LilVectLoggingFile(
    char *fnam,			/* (I) file name */
    char *mode,			/* (I) write mode */
    char *funcName,		/* (I) function name */
    /* DVType varType, void *lptr, char  *argName,
     * ...
     * DVType varType, void *lptr, char  *argName,
     * 0)	LAST ARGUMENT MUST BE ZERO
     */
    ...)
{
    static    char     routine[] = "Q3LilVectLoggingFile";
    int	      status = FAILURE;

    FILE	*fp = NULL;
    va_list	ap;
    DVType	varType;
    void	*lptr;
    char	*argName;

    va_start(ap, funcName);

    if (fnam != NULL) {
        if ((fp = fopen(fnam, mode)) == NULL) {
	    Q3ErrMsg("%s: `%s' (%s).\n", routine, fnam, Q3Strerror());
	    goto RETURN;
        }
    }
    
    Q3FPrintf(fp, "COM: ---------------------------------"
	      "------------------------------------\n");
    Q3FPrintf(fp, "COM: %s\n",
	      (funcName != NULL ? funcName : "[NO_NAME]"));
    Q3FPrintf(fp, "COM: ---------------------------------"
	      "------------------------------------\n");
    
    
    while ((varType = (DVType) va_arg(ap, DVType)) != 0) {

        lptr    = (void*) va_arg(ap, void*);
        argName = (char*) va_arg(ap, char*);
      
        if (Q3LilVectLogging(fp, varType, lptr, argName) != SUCCESS)
	    goto RETURN;
    }
    
    status = SUCCESS;

  done:
  RETURN:

    va_end(ap);
    if (fp) fclose(fp);
    return(status);
} /* Q3LilVectLoggingFile */



/*f----------------------------------------------------------------------
 * Similar to Q3LilVectLoggingFile, put write to a given file pointer.
 */

int Q3LilVectLoggingFp(
    FILE *fp,			/* (I) file to print (or NULL) */
    char *funcName,		/* (I) function name */
    /* DVType varType, void *lptr, char  *argName,
     * ...
     * DVType varType, void *lptr, char  *argName,
     * 0)	LAST ARGUMENT MUST BE ZERO
     */
    ...)
{

    static    char     routine[] = "Q3LilVectLoggingFp";
    int	      status = FAILURE;
    va_list	ap;
    DVType	varType;
    void	*lptr;
    char	*argName;
    
    va_start(ap, funcName);
    
    Q3FPrintf(fp, "COM: ---------------------------------"
	      "------------------------------------\n");
    Q3FPrintf(fp, "COM: %s\n",
	      (funcName != NULL ? funcName : "[NO_NAME]"));
    Q3FPrintf(fp, "COM: ---------------------------------"
	      "------------------------------------\n");
    
    while ((varType = (DVType) va_arg(ap, DVType)) != 0) {
      
        lptr    = (void*) va_arg(ap, void*);
        argName = (char*) va_arg(ap, char*);
      
        if (Q3LilVectLogging(fp, varType, lptr, argName) != SUCCESS)
	    goto RETURN;
    }
    
    status = SUCCESS;

  done:
  RETURN:

    va_end(ap);
    return(status);
} /* Q3LilVectLoggingFp */


/*f----------------------------------------------------------------------
 * Returns the length of a LIL vector "lptr"
 * of type "varType". "argName" is used for error messages. 
 */

int Q3LilVectSize(
    DVType varType,        /* (I) Lil type */
    void *lptr,            /* (I) Lil vector */
    char *argName)         /* (I) for debugging */
{
    static	char	routine[] = "Q3LilVectSize";

    /* Get input vector length */
    if (lptr == NULL) return(0);
       
    switch (varType) {
    case Q3_DOUBLE_L:
    case Q3_PERCENT_L:
    case Q3_BPOINT_L:
        return ARGSIZE((double*) lptr);
    case Q3_LONG_L:
        return ARGSIZE((long*)   lptr);
    case Q3_TDATE_L:
        return ARGSIZE((long*)   lptr);
    case Q3_TDATEINTERVAL_L:
        return ARGSIZE((double*) lptr);
    case Q3_CHAR_BLOCK_L:
    case Q3_CHAR_L:
        /*return ARGSIZE((char*)   lptr);*/
        return (int)(lptr == NULL ? 0 : (unsigned char)((char*) lptr)[0]);

    default:
        Q3ErrMsg("%s: [%s] bad type %ld.\n", routine, argName, (long) varType);
        return(-1);
    }
}


int Q3FPrintf(
    FILE *fp,        /* (I) file name */
    char *fmt,       /* (I) format */
    ...)
{
    va_list	ap;
    char	buf[2048], buf2[2048];
    
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    
    if (fp != (FILE*)NULL) {
        fputs(buf, fp);
    } else {
        char	*p = buf, *q;
        int	i, imax=255;
        while (*p != '\0') {
	    q = buf2;
	    i = 0;
	    while ((i++ <= imax) && (*q++ = *p++));
	    p--; q--;
	    *q = '\0';
	    Q3ErrMsg("%s", buf2);
        }
    }

    return(0);
} /* Q3LilVectSize */


/*f----------------------------------------------------------------------
 * Prints TDates in mm/dd/yyyy format
 */

char* Q3TDatePrint(
    char *string,          /* (O) date string */
    TDate aDate)           /* (I) Alib TDate */
{

#undef	MAX_IDX
#define	MAX_IDX	8

    static	char	tmp[MAX_IDX][64] ;
    static	int	tmpIdx=0;
    char	*s ;
    TMonthDayYear mdy;
    
    if (string == NULL) {
        s = tmp[tmpIdx];
        tmpIdx++;
        if (tmpIdx > MAX_IDX-1) tmpIdx=0;
    } else {
        s = string;
    }
    
    GtoDateToMDY(aDate, &mdy);

    sprintf(s, "%02d/%02d/%04d",
        (int) mdy.month, (int) mdy.day, (int) mdy.year) ;
    
    return(s)  ;

#undef	MAX_IDX
} /* Q3TDatePrint */
