/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	STR - String Functions
 * File:	drlstr.h
 * Function:	String routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlstr_H
#define	_drlstr_H
#include "drlstd.h"
#include <stdarg.h>

/* basic string functions */
extern	char*	DrlStrAlloc(int n1);
extern	void	DrlStrFree(char *s);
extern	char**	DrlStrArrayAlloc(int n1, int n2);
extern	void	DrlStrArrayFree(char **p, int n1);
extern	int		DrlStrArraySort(char **p, int n1);
extern	char*	DrlStrInsertChar(char *s, char *p, char c);
extern	char*	DrlStrDeleteChar(char *s, char *p);
extern	char*	DrlStrRemoveChars(char *s, char *ct);
extern	char*	DrlStrReplace(char *s, char *before,
					char *after);
extern	char*	DrlStrReplaceStr(char *s, char *before,
					char *after);
extern	char*	DrlStrReplaceNonPrintChar(char *s);
extern	char*	DrlStrNCpy(char *s, char *ct, int n);
extern	char*	DrlStrRightNCpy(char *s, char *ct, int n);
extern	char*	DrlStrCut(char *s, int n);
extern	char*	DrlStrPrint(char *fmt, ...);
extern	char*	DrlStrToken(const char *s, const char *tok,
					char *s0, char **s1);
extern	char*	DrlStrScanString(char **p, char *q);
extern	char*	DrlStrToLower(char *s);
extern	char*	DrlStrToUpper(char *s);
extern	int		DrlStrCmpCaseIndif(const char *cs,
					const char *ct);
extern	int		DrlStrStartsByCaseIndif(const char *cs,
					const char *ct);


/* String operations on filenames */
extern	char*	DrlGetBaseName(char *s, char *ct);
extern	char*	DrlGetRootName(char *s, const char *ct);
extern	char*	DrlGetSuffix(char *s, char *ct);
extern	char*	DrlGetPathDir(char *s, char *ct);



extern	int 	DrlStrParseFuncScan(char *str,
					char *funcToken, int nVar,
					/*  char *varToken, char varType,
					 *  void* varPtr,
				 	 * int varDefFlag, void* varDef */
					...);
extern	char*	DrlStrParseFuncPrint(char *s,
					char *funcToken, int nVar,
					/* char *varToken, char varType,
					 * void* varPtr,
					 * int varDefFlag, void* varDef */
					...);

extern	int		DrlStrParseOptions(
		const char *s,		/* (I) input string */
		/* DVType varType1, void* varPtr1, char *varToken1,
		 *  ...
		 * DVType varTypeN, void* varPtrN, char *varTokenN,
		 * DRL_NULL_T)		LAST ARGUMENT MUST BE DRL_NULL_T
		 */ ...);


extern	int		DrlStrTokScan(char *s, DVType varType,
					void *var);
extern	char*	DrlStrTokPrint(char *s, DVType varType,
					void *var);

extern	int		DrlStrSubsEnv(char *s);

/*
 *
 */
extern	int		DrlStrLineVarScan(
	char *s,		/* (I) input string (unchanged) */
	/* DVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * DVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...);

extern	char*	 DrlStrLineVarPrint(
	char *s,		/* (I) input string (unchanged) */
	/* DVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * DVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...);



extern	int		DrlStrLongValueScanV(char *s, char *errName,
					long *value, va_list ap);
extern	char*	DrlStrLongValuePrintV(char *s, char *errName,
					long value, va_list ap);

extern	int		DrlStrLongValueScan(
	char *s,		/* (I) input string (unchanged) */
	char *errName,		/* (I) variable name for debugging */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)		LAST ARGUMENT MUST BE NULL
	 */ ...);

extern	char*	DrlStrLongValuePrint(
	char *s,		/* (O) output string */
	char *errName,		/* (I) variable name for debugging */
	long value,		/* (I) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)		LAST ARGUMENT MUST BE NULL
	 */ ...);


/*
 *
 */
extern	char*	DrlCurrentDateTimePrint(char *buf);
extern	char*	DrlFloatPrint(char *s, double x, int ndigits);
extern	char*	DrlDoubleIntervalPrint(char *s, double value);
extern	int		DrlDoubleIntervalScan(char *s, double *value);
extern	int		DrlCurScan(char *s, double *x);
extern	char*	DrlCurPrint(char *s, double x, int ndigits);



#endif	/* _drlstr_H */
