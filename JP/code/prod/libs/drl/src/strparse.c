/************************************************************************
 * Module:	DRL
 * Submodule:	STR
 * File:	
 * Function:	String Utilities
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>		/* time_t */

#include "drltime.h"		/* DrlDIntervalScan() */
#include "drlvtype.h"

#include "drlstr.h"		/* Prototype Consistency */

#define	__DEBUG__
#undef	__DEBUG__

#define	MAXVAR		32

#define	LEFT_BRAK	'('
#define	RIGHT_BRAK	')'


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




/*f---------------------------------------------------------------------
 * Parses a string of the form
 * \begin{verbatim}
 *   FUNC(varToken1=value1,....,varTokenN=valueN)
 * \end{verbatim}
 * Returns SUCCESS/FAILURE.
 */

int
DrlStrParseFuncScan(
	char *s,		/* (I) input string */
	char *funcToken,	/* (I) function token (or NULL) */
	int	nVar,		/* (I) # of variable description to follow */
	/*  char *varToken, DVType varType, void* varPtr,
			int varDefFlag, void* varDef */
	...)
{
static	char	routine[] = "DrlStrParseFuncScan";
	int	status = 1;	/* FAILURE */
	va_list	ap;
	int	i;
	char	buf[256],
		*p, *q, *r, *v0, *v1,
		tmpS0[255], *tmpS2;
	char	*varToken[MAXVAR];
	DVType	varType[MAXVAR];
	void*	*varPtr[MAXVAR],
		*varDef[MAXVAR];
	int	varDefFlag[MAXVAR],
		varSetFlag[MAXVAR];

	/*
	 *
	 */
	if (nVar > MAXVAR-1) {
	    DrlErrMsg("%s: too any arguments (max %d).\n",
		routine, MAXVAR);
	    return(FAILURE);
	}
	strcpy(buf, s);

	/*
	 *
	 */
	va_start(ap, nVar);
	for (i=0; i <= nVar-1; i++) {
		varToken[i]   = (char*) va_arg(ap, char*);
		varType[i]    = (DVType)  va_arg(ap, DVType);
		varPtr[i]     = (void*) va_arg(ap, void*);
		varDefFlag[i] = (int) va_arg(ap, int);
		varDef[i]     = (void*) va_arg(ap, void*);
		varSetFlag[i] = 0;
	}
	va_end(ap);



	/*
	 * Scan function name 
	 */
	/* Advance to 1st nonspace character */
	for (p=buf; (isspace(*p)) && (*p != '\0'); p++);

	if (funcToken != NULL) {
		for (q=p; (*q != '\0') && (*q != LEFT_BRAK); q++);
		if (q == p) goto done;
		if (*q != LEFT_BRAK) goto done;
	
		for (r=q; (*r != '\0') && (*r != RIGHT_BRAK); r++);
		if (*r != RIGHT_BRAK) goto done;

		*q = '\0'; q++;
		*(r+1) = '\0';

#ifdef	__DEBUG__
		printf("p=`%s'  q=`%s'\n", p, q);
#endif

		/* Scan function name. If get different name, return -1.
		 */
		if (strcmp(funcToken, p) != 0) return(-1);
	} else {
		q = p;
		p = NULL;
	}

	/*
	 * Now p points on function name, q on beginning of args
	 */



	/*
	 * scan arguments
	 */
	r = q;
	while ((v0 = DrlStrToken(r, ",)", tmpS0, &tmpS2)) != NULL) {
		if (r != NULL) r = NULL;
#ifdef	__DEBUG__
printf("v0=`%s'\n", v0);
#endif
		/*
		 *
		 */
		for (v1 = v0; (*v1 != '\0') && (*v1 != '='); v1++);
		if (*v1 == '\0') return(__LINE__);
		*v1 = '\0'; v1++;

#ifdef	__DEBUG__
printf("v0=`%s'  v1=`%s'\n", v0, v1);
#endif

		/*
		 * A string of the form "v0=v1" has been identified,
		 * Check if variable name "v0" is in the table of
		 * variable names and scan its value "v1"
		 */
		for (i=0; i<=nVar-1; i++) {
		    if (strcmp(v0, varToken[i]) == 0) {
			if (DrlVTypeScan(v1, varType[i], varPtr[i]) != 0) {
			    DrlErrMsg("%s: error arg `%s' (type %s)\n",
				routine, v0, DrlVTypeName(varType[i]));
			    goto done;
			}
			varSetFlag[i] = 1;
			i = (-1);
			break;
		    }
		}

		/* if no name matched, error */
		if (i > 0) {
		    DrlErrMsg("%s: can't match arg `%s'\n", routine, v0);
		    goto done;
		}
	} /* while(...*/

	/*
	 * Check all variables that do not have default values
	 * have been set.
	 */
	for (i=0; i <= nVar-1; i++) {
	    if (varSetFlag[i] == 0) {
		if (varDefFlag[i] != 0) {
			/* variable not set but has a dfault
			 * value. OK */
			DrlVTypeLet(varPtr[i], 0,
				varDef[i], 0, varType[i]);
		} else {
			/* variable not set nad has no default
			 * value. ERROR  */
			DrlErrMsg("%s: mandatory variable `%s' "
				"not set for `%s'\n",
				routine, varToken[i], funcToken);
			goto done;
		}
	    }
	}

	/* made it through OK */
	status = 0; /* SUCCESS */
done:
	/* errlog only if bad syntax */
	if ((status != 0)  && (status != -1)) {
	    DrlErrMsg("%s: scanning `%s' failed.\n", routine, s);
	}
	return(status);
}



/*f-------------------------------------------------------------------
 * Print a string of the form
 * \begin{verbatim}
 * FUNC(varToken11=value1,....,varTokenN=valueN)
 * \end{verbatim}
 */

char*
DrlStrParseFuncPrint(
	char	*s,		/* (I) input string */
	char	*funcToken,	/* (I) function token */
	int	nVar,		/* (I) # of variable description to follow */
	/*  char *varToken, DVType varType, void* varPtr,
			int varDefFlag, void* varDef */
	...
)
{
static	char	routine[] = "DrlStrParseFuncPrint";
	va_list	ap;
	int	i;
static	char	str[256];
	char	buf[64];
	char	*varToken[MAXVAR];
	DVType	varType[MAXVAR];
	void	*varPtr[MAXVAR],
		*varDef[MAXVAR];
	int	varDefFlag[MAXVAR];

	/*
	 *
	 */
	s = (s != (char*)NULL ? s : str);
	if (nVar > MAXVAR-1) {
	    DrlErrMsg("%s: too any arguments (max %d).\n",
		routine, MAXVAR);
	    strcpy(s, "<ERROR>");
	    return(s);
	}


	/*
	 *
	 */
	va_start(ap, nVar);
	for (i=0; i <= nVar-1; i++) {
		varToken[i]   = (char*)    va_arg(ap, char*);
		varType[i]    = (DVType) va_arg(ap, DVType);
		varPtr[i]     = (void*)    va_arg(ap, void*);
		varDefFlag[i] = (int)      va_arg(ap, int);
		varDef[i]     = (void*)    va_arg(ap, void*);
	}
	va_end(ap);


	/*
	 * print
	 */
	strcpy(s, "");

	if (funcToken != NULL) {
	    sprintf(buf, "%s(", funcToken);
	    strcat(s, buf);
	}

	for (i=0; i <= nVar-1; i++) {
		sprintf(buf, "%s=", varToken[i]);
		strcat(s, buf);
		DrlVTypePrint(buf, varType[i], varPtr[i]);
		strcat(s, buf);
		if (i<nVar-1) {
			sprintf(buf, ",");
			strcat(s, buf);
		}
	}

	if (funcToken != NULL) {
	    sprintf(buf, ")");
	    strcat(s, buf);
	}

	/* removes all spaces */
	DrlStrRemoveChars(s, " \t");

	return(s);
}


/*f---------------------------------------------------------------------
 * Parses a string of the form
 * \begin{verbatim}
 *   <varToken1>=<value1>:....:<varTokenN>=<valueN>
 * \end{verbatim}
 * Designed to be used to parse applications options defined in a string
 * (e.g. that can be retrieved in the environment).
 * For example, assuming that
 * one wants to scan three variables "iVal" (integer), "dVal" (double)
 * and "sVal" (char array), a call to the routine would be
 * \begin{verbatim}
 *     char   buf[256];                     // input string
 *     int    iVal;                         // variables to be set
 *     double dVal;
 *     char   sVal[256];
 *     iVal =0; dVal = 0e0; sVal[0] = '\0'; // default values.
 *     ...
 *     strcpy(buf, "iVal=1:dVal123.456:sVal=OK");   // for the example 
 *     ...
 *     errCode = DrlStrParseOptions(
 *                 DR_INT_T,        (void*) &iVal,  "iVal",
 *                 DR_DOUBLE_T,     (void*) &dVal,  "dVal",
 *                 DR_CHAR_ARRAY_T, (void*)  sVal,  "sVal",
 *                 DR_NULL_T);
 *     ...
 * \end{verbatim}
 * If some variables specified in the function call are absent from the
 * input string, their values are unchanged on exit.
 * If a variable appears in the string and is not specified
 * in the function call, an error is returned.
 * Whenever the token "nil" is encountered, parsing is stopped
 * and the routine returns with a success error code.
 */

int
DrlStrParseOptions(
	const char *s,		/* (I) input string */
	/* DVType varType1, void* varPtr1, char *varToken1,
	 *  ...
	 * DVType varTypeN, void* varPtrN, char *varTokenN,
	 * DRL_NULL_T)		LAST ARGUMENT MUST BE DRL_NULL_T
	 */ ...)
{
static	char	routine[] = "DrlStrParseOptions";
	int	status = FAILURE;

	va_list	ap;
	char	buf[256],
		*r, *v0, *v1,
		tmpS0[256], *tmpS2;

	int	numVar, idxVar;
	DVType	varType[MAXVAR];
	void*	*varPtr[MAXVAR];
	char	*varToken[MAXVAR];

	/* Get function arguments */
	numVar = 0;
	va_start(ap, s);
	while ((varType[numVar] = (DVType) va_arg(ap, DVType)) != 0L) {
		varPtr[numVar]   = (void*) va_arg(ap, void*);
		varToken[numVar] = (char*) va_arg(ap, char*);

		numVar++;
		if (numVar > MAXVAR-1) {
			DrlErrMsg("%s: too any arguments (max %d).\n",
				routine, MAXVAR-1);
			goto done;
		}
	}
	va_end(ap);


	/* Scan */
	strcpy(buf, s);
	r = buf;
	while ((v0 = DrlStrToken(r, ":; \t\n", tmpS0, &tmpS2)) != NULL) {
	    if (r != NULL) r = NULL;

#ifdef	__DEBUG__
	    printf("FIELD: v0=`%s'\n", v0);
#endif

	    /* Check for NIL
	     */
	    if ((!strcmp(v0, "nil")) || (!strcmp(v0, "NIL"))) {
		status = SUCCESS;
		goto done;
	    }


	    /* Check if string of the form "v0=v1"
	     */
	    for (v1 = v0; (*v1 != '\0') && (*v1 != '='); v1++);

	    if (*v1 == '\0') {
		/* String is not of the form "v0=v1"
		 * Set the corresponding variable to TRUE.
		 */
	        *v1 = '\0'; v1++;
#ifdef	__DEBUG__
		printf("BOOL : v0=`%s'\n");
#endif
		for (idxVar=0; idxVar<=numVar-1; idxVar++) {
		    if (strcmp(v0, varToken[idxVar]) == 0) {
			switch (varType[idxVar]) {
			case DRL_INT_T: 
				*((int*)varPtr[idxVar]) = (int) TRUE;
				break;
			case DRL_LONG_T: 
				*((long*)varPtr[idxVar]) = (long) TRUE;
				break;
			default:
				DrlErrMsg("%s: error arg `%s' (type %s).\n",
					routine, v0,
					DrlVTypeName(varType[idxVar]));
				goto done;
			}
			idxVar = (-1);
			break;
		    }
		}

		/* if no name matched, error */
		if (idxVar > 0) {
		    DrlErrMsg("%s: can't match arg `%s'\n", routine, v0);
		    goto done;
		}


	    } else {
		/* A string of the form "v0=v1" has been identified,
		 * Check if variable name "v0" is in the table of
		 * variable names and scan its value "v1"
		 */
	        *v1 = '\0'; v1++;

#ifdef	__DEBUG__
		printf("TOKEN: v0=`%s'  v1=`%s'\n", v0, v1);
#endif
		for (idxVar=0; idxVar<=numVar-1; idxVar++) {
		    if (strcmp(v0, varToken[idxVar]) == 0) {
			if (DrlVTypeScan(v1,
				varType[idxVar], varPtr[idxVar]) != 0) {
			    DrlErrMsg("%s: error arg `%s' (type %s)\n",
				routine, v0, DrlVTypeName(varType[idxVar]));
			    goto done;
			}
			idxVar = (-1);
			break;
		    }
		}

		/* if no name matched, error */
		if (idxVar > 0) {
		    DrlErrMsg("%s: can't match arg `%s'\n", routine, v0);
		    goto done;
		}
	    }
	} /* while(...*/

	/* made it through OK */
	status = SUCCESS;
done:
	/* errlog only if bad syntax */
	if (status != SUCCESS) {
	    DrlErrMsg("%s: scanning `%s' failed.\n", routine, s);
	}
	return(status);
}





/*f---------------------------------------------------------------------
 * Scans the next token of type "varType" in char string "s"
 * and puts the result in "var".
 * Behaves likes the ANSI  "strtok" called on "s" with tokens
 * SPACE, TAB, NEWLINE and ';'.
 * Returns 0 iff scan successful.\\
 * <b> Example:</b>
 * \begin{verbatim}
 *          int    iVal;
 *          double dVal;
 *          TDAte  tVal;
 *          ...
 *          / * string s is "0.123   12    5/14/1966" * /
 *          ...
 *          / * 1st call * /
 *          status = DrlStrTokScan(s, DRL_DOUBLE_T, (void*) &dVal);
 *          ...
 *          / * next calls: pass NULL * /
 *          status = DrlStrTokScan(NULL, DRL_INT_T, (void*) &iVal);
 *          ...
 *          status = DrlStrTokScan(NULL, DRL_TDATE_T, (void*) &tVal);
 *          ...
 * 
 * \end{verbatim}
 * 
 */

int
DrlStrTokScan(char *s, DVType varType, void *var)
{
static	char	routine[] = "DrlStrTokScan",
		tok[] = " \t;\n";
	char 	*p;
static	char	st0[512], *st1;


	/*if ((p = strtok(s, tok)) == NULL) {*/
	if ((p = DrlStrToken(s, tok, st0, &st1)) == NULL) {
		DrlErrMsg("%s: end of line\n", routine);
		return(1);
	}

	if (DrlVTypeScan(p, varType, var) != 0) {
		DrlErrMsg("%s: can't read type %s in `%s'\n",
			routine, DrlVTypeName(varType), p);
		return(1);
	}
	return(0);
}

/*f---------------------------------------------------------------------
 * Prints a variable <i> var</i> (cast to void*)
 * of type <i> varType</i> in char string <i> s</i>.
 * <b> Example:</b>
 * \begin{verbatim}
 *          int    iVal = 0.123;
 *          double dVal = 12;
 *          char   sVal[] = "nabuchodonosor";
 *          ...
 *          strcpy(s, "");
 *          DrlStrTokPrint(s, DRL_DOUBLE_T, (void*) &dVal);
 *          DrlStrTokPrint(s, DRL_INT_T, (void*) &iVal);
 *          DrlStrTokPrint(s, DRL_STRING_T, (void*) &sVal);
 *          ...
 *          / * string s is "0.123   12  nabuchodonosor" * /
 * 
 * \end{verbatim}
 */

char*
DrlStrTokPrint(char *s, DVType varType, void *var)
{
	strcat(s, "\t");
	strcat(s, DrlVTypePrint(NULL, varType, var));
	return(s);
}



/*f---------------------------------------------------------------------
 * Scans variables on a line of the form
 * \begin{verbatim}
 *   value1  value2 ... valueN
 * \end{verbatim}
 */

int
DrlStrLineVarScan(
	char *s,		/* (I) input string (unchanged) */
	/* DVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * DVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...)
{
static	char		routine[] = "DrlStrLineVarScan";
	int		status = FAILURE;
	va_list		ap;

	int		argNo = 0;
	char		*p, *q,
			tmpS0[256], *tmpS2;
static	char		tok[] = " \t;\n";

	char		*varName;
	DVType	varType;
	void*		*varPtr;

	/*
	 *
	 */
	q = s;
	va_start(ap, s);

	while ((varType = (DVType) va_arg(ap, DVType)) != 0L) {
	    argNo++;

	    varName = (char*) va_arg(ap, char*);
	    varPtr  = (void*) va_arg(ap, void*);

	    if ((p = DrlStrToken(q, tok, tmpS0, &tmpS2)) == NULL) {
		DrlErrMsg("%s: can't read variable # %d "
			"%s of type %s in `%s' (EOF)\n",
			routine, argNo,
			(varName != NULL ? varName : ""),
			DrlVTypeName(varType), p);
		goto done;
	    }
	    /* call strtok only the 1st with non null */
	    if (q != NULL) q = NULL;

	    if (DrlVTypeScan(p, varType, varPtr) != 0) {
		DrlErrMsg("%s: can't read variable # %d "
			"%s of type %s in `%s'\n",
			routine, argNo,
			(varName != NULL ? varName : ""),
			DrlVTypeName(varType), p);
		goto done;
	    }
	}
	va_end(ap);


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: scanning `%s' failed.\n", routine, s);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Prints variables on a line of the form
 * \begin{verbatim}
 *   value1  value2 ... valueN
 * \end{verbatim}
 * in a string <i> s</i>.
 * If <i> s</i> is NULL, uses a static buffer and returns it.
 */

char*
DrlStrLineVarPrint(
	char *s,		/* (I) input string (unchanged) */
	/* DVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * DVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...)
{
static	char		routine[] = "DrlStrLineVarPrint";
	int		status = FAILURE;
	va_list		ap;

static	char	str[256];
	char		*varName;
	DVType		varType;
	void*		*varPtr;

	/*
	 *
	 */
	va_start(ap, s);
	s = (s != (char*)NULL ? s : str);
	s[0] = '\0';

	while ((varType = (DVType) va_arg(ap, DVType)) != 0L) {
	    varName = (char*) va_arg(ap, char*);
	    varPtr  = (void*) va_arg(ap, void*);

	    strcat(s, " ");
	    strcat(s, DrlVTypePrintSmart(NULL, varType, varPtr));
	}
	va_end(ap);
	/*strcat(s, "\n");*/


	/* made it through OK */
	status = SUCCESS;
/*done:*/
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine);
	    strcpy(s, "ERR");
	}
	return(s);
}


/*----------------------------------------------------------------------
 * Internally called.
 */

int
DrlStrLongValueScanV(
	char *s,		/* (I) input string (unchanged) */
	char *errName,		/* (I) variable name for debugging */
	long *value,		/* (O) value */
	va_list ap)		/* (I) arguments list */
{
	int		status = FAILURE;
	char		*varName;
	long		varValue;
	char		errString[1024];

	sprintf(errString, "DrlStrLongValueScanV: can't read %s in `%s' "
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
	DrlErrMsg(errString);
	return(FAILURE);
done:
	return(status);
}


/*----------------------------------------------------------------------
 * Internally called.
 */

char*
DrlStrLongValuePrintV(
	char *s,		/* (O) output string */
	char *errName,		/* (I) variable name for debugging */
	long value,		/* (I) value */
	va_list ap)		/* (I) arguments list */
{
	char		*varName;
	long		varValue;

	BUF_NEXT(s);

	while ((varName = (char*) va_arg(ap, char*)) != NULL) {
	    varValue = (long) va_arg(ap, long);
	    if (value == varValue) {
		sprintf(s, "%s", varName);
		goto done;
	    }
	}
	sprintf(s, "ERROR");
done:
	va_end(ap);
	return(s);
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
 * errCode = DrlStrLongValueScan(buf, "stub type" &value,
 *              "B", DRL_BOND_STUB,
 *              "S", DRL_SIMPLE_STUB,
 *              NULL);
 * ...
 * \end{verbatim}
 */

int
DrlStrLongValueScan(
	char *s,		/* (I) input string (unchanged) */
	char *errName,		/* (I) variable name for debugging */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...)
{
static	char		routine[] = "DrlStrLongValueScan";
	int		status = FAILURE;
	va_list		ap;

	va_start(ap, value);
	status = DrlStrLongValueScanV(s, errName, value, ap);
	va_end(ap);

	if (status != SUCCESS) {
	    DrlErrMsg("%s: scanning `%s' failed.\n", routine, s);
	}
	return(status);
}

/*f---------------------------------------------------------------------
 * Prints a long value in a string <i> s</i> by trying to match
 * it with different possible values. See <i> DrlStringLongValueScan</i>
 * <br> <b> Example:</b>
 * \begin{verbatim}
 * char  *buf;
 * long  value = DRL_SIMPLE_STUB;
 * ...
 * DrlStrLongValuePrint(buf, &value,
 *              "B", DRL_BOND_STUB,
 *              "S", DRL_SIMPLE_STUB,
 *              NULL);
 * ...
 * \end{verbatim}
 */

char*
DrlStrLongValuePrint(
	char *s,		/* (O) output string */
	char *errName,		/* (I) variable name for debugging */
	long value,		/* (I) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...)
{
	va_list		ap;
	va_start(ap, value);
	DrlStrLongValuePrintV(s, errName, value, ap);
	va_end(ap);
	return(s);
}


/*f---------------------------------------------------------------------
 * Performs a substitutions of the environment variables 
 * in the strings <i> s</i> (using the <i> getenv</i> routine
 * of the C standard library).\\
 * Returns SUCCESS/FAILURE.
 */

int
DrlStrSubsEnv(char *s)
{
static	char	routine[] = "DrlStrSubsEnv";

	char	*p, *q, *v,
		*varval,
		buf[1024],
		varname[256];

#define	ISNCHAR(c)	(isalnum(c) || (c =='_'))

	for (;;) {
		p = s;
		while ((*p != '$') && (*p != '\0')) p++;
		if (*p == '\0') return(SUCCESS);
		q = p+1;

		if (*q == '{') {
			q++;
			v = varname;
			while ((*q != '}')  && (*q != '\0')) *v++ = *q++;
			if (*q++ == '\0')  goto done;
			*v = '\0';
		} else if (*q == '(') {
			q++;
			v = varname;
			while ((*q != ')')  && (*q != '\0')) *v++ = *q++;
			if (*q++ == '\0')  goto done;
			*v = '\0';
		} else {
			v = varname;
			while (ISNCHAR(*q) && (*q != '\0')) *v++ = *q++;
			*v = '\0';
		}


		strcpy(buf, q);

		if ((varval = getenv(varname)) == NULL) {
			DrlErrMsg("%s: in string `%s', variable `%s' "
				"undefined.\n", routine, s, varname);
			return(FAILURE);
		}
		strcpy(p, varval);
		while (*(++p));
		strcpy(p, buf);
	}
	PROGRAM_BUG();

done:
	DrlErrMsg("%s: syntax error in `%s'.\n",
		routine, s);
	return(FAILURE);
}




/***************************************************************
 *
 ***************************************************************/

/*f-------------------------------------------------------------
 * Prints in a char buffer <i> s</i> the current date and time.
 * Returns {s}.
 * If "s" is NULL, returns a pointer to a static char string.
 */

char*
DrlCurrentDateTimePrint(char *s)
{
static	char	buf[64];
	time_t	calTime;

	s = (s != NULL ? s : buf);

	if (time(&calTime) == -1) {
		strcpy(s, "N/A");
		return(s);
	}

	strcpy(s, ctime(&calTime));
	if (s[strlen(s)-1] == '\n') s[strlen(s)-1] = '\0';
	return(s) ;
}





/*f-------------------------------------------------------------
 * Prints a double <i> x</i> in a string <i> s</i>
 * using exactly <i> ndigits</i> (must be at least 7).
 * If <i> s</i> is NULL, uses a static buffer and returns it.
 */

char*
DrlFloatPrint(char *s, double x, int ndigits)
{
static	int	l;
	char	fmt[16];


	BUF_NEXT(s);


#if defined(UNIX)
	if (!DrlDoubleIsFinite(x)) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "NaN");
		return(s);
	}
#endif

	if (x == 0e0) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "0");
	} else if (x <= -1e20) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "--");
	} else if (x >= 1e80) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "++");
	} else if ((x < 1e-80) && (x > 0)) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "+0");
	} else if ((x > -1e-20) && (x < 0)) {
		sprintf(fmt, "%%%ds", ndigits);
		sprintf(s, fmt, "-0");
	} else if ((l = fabs(log(fabs(x))*0.4342e0)) > ndigits-4) {
		sprintf(fmt, "%% %d.%de", ndigits-4, ndigits-7);
		sprintf(s, fmt, x);
	} else {
		if (fabs(x) < 1e0) l = ndigits-3;
		else l = ndigits-3-l;

		sprintf(fmt, "%% %d.%df", ndigits-2, l);
		sprintf(s, fmt, x);
	}

	return(s);
}



/*f-------------------------------------------------------------
 * Scans in the string <i> s</i>
 * a date interval, either as a number (e.g. 1.25, etc)
 * or an number of intervals (e.g. 1A, 1W, 3S, etc).
 * The result is converted to double and put in <i> value</i>.
 * Returns SUCCESS/FAILURE.
 */

int
DrlDoubleIntervalScan(char *s, double *value)
{
	int	status = FAILURE;
	char	*buf = NULL;
	DInterval	interval;

	BUF_NEXT(buf);

	if (sscanf(s, "%s", buf) != 1)
		goto done;
	if (isalpha(s[strlen(s)-1])) {
	    if (DrlDIntervalScan(buf, &interval) != SUCCESS)
		goto done;
	    if (DrlDIntervalToYears(&interval, value) != SUCCESS)
	    	goto done;
	} else {
	    if (sscanf(buf, "%lf", value) != 1)
		goto done;
	}


	status = SUCCESS;
done:
	if (status != SUCCESS)
		DrlErrMsg("DoubleIntervalScan: failed on `%s'.\n", s);
	return(status);
}

/*f-------------------------------------------------------------
 * Prints a double variable <i> value</i> as a date interval
 * (e.g. 1A, 1W, 3S, etc) in a string <i> s</i>.
 * If <i> s</i> is NULL, uses a static buffer and returns it.
 */

char*
DrlDoubleIntervalPrint(char *s, double value)
{
	DInterval	interval;

	BUF_NEXT(s);

	if (DrlYearsToDInterval(value, &interval) != SUCCESS) {
		strcpy(s, "ERR");
		goto done;
	}
	strcpy(s, DrlDIntervalPrint(NULL, &interval));
done:
	return(s);
}


/*f-------------------------------------------------------------
 * Scans a double value in a currency format
 * $\pm XXX,XXX,XXX.XX$ in the string <i> s</i>
 * and puts the result in <i> x</i>. 
 * If the string containts the percent sign, the output is 
 * divided by 100.
 * Returns 0 on success.
 */

int
DrlCurScan(char *s, double *x)
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

char*
DrlCurPrint(char *s, double x, int ndigits)
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



