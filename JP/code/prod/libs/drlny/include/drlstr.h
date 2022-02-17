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

#if defined(CLIB) && (defined(_WINDLL) && !(defined(WIN32) || defined(WIN32)))
# include "cfileio.h"
#endif

/* sscanf and sprintf */
#if defined(CLIB) && defined(_WINDLL) && !(defined(WIN32) || defined(_WIN32))
# define	sscanf		GtoSscanf
#endif

/* basic string functions */
extern	DLL_EXPORT(char*)	DrlStrAlloc(int n1);
extern	DLL_EXPORT(void)	DrlStrFree(char *s);
extern	DLL_EXPORT(char**)	DrlStrArrayAlloc(int n1, int n2);
extern	DLL_EXPORT(void)	DrlStrArrayFree(char **p, int n1);
extern	DLL_EXPORT(int)		DrlStrArraySort(char **p, int n1);
extern	DLL_EXPORT(char*)	DrlStrInsertChar(char *s, char *p, char c);
extern	DLL_EXPORT(char*)	DrlStrDeleteChar(char *s, char *p);
extern	DLL_EXPORT(char*)	DrlStrRemoveChars(char *s, char *ct);
extern	DLL_EXPORT(char*)	DrlStrReplace(char *s, char *before,
					char *after);
extern	DLL_EXPORT(char*)	DrlStrReplaceStr(char *s, char *before,
					char *after);
extern	DLL_EXPORT(char*)	DrlStrReplaceNonPrintChar(char *s);
extern	DLL_EXPORT(char*)	DrlStrNCpy(char *s, char *ct, int n);
extern	DLL_EXPORT(char*)	DrlStrRightNCpy(char *s, char *ct, int n);
extern	DLL_EXPORT(char*)	DrlStrCut(char *s, int n);
extern	DLL_EXPORT(char*)	DrlStrPrint(char *fmt, ...);
extern	DLL_EXPORT(char*)	DrlStrToken(const char *s, const char *tok,
					char *s0, char **s1);
extern	DLL_EXPORT(char*)	DrlStrScanString(char **p, char *q);
extern	DLL_EXPORT(char*)	DrlStrToLower(char *s);
extern	DLL_EXPORT(char*)	DrlStrToUpper(char *s);
extern	DLL_EXPORT(int)		DrlStrCmpCaseIndif(const char *cs,
					const char *ct);
extern	DLL_EXPORT(int)		DrlStrStartsByCaseIndif(const char *cs,
					const char *ct);


/* String operations on filenames */
extern	DLL_EXPORT(char*)	DrlGetBaseName(char *s, char *ct);
extern	DLL_EXPORT(char*)	DrlGetRootName(char *s, const char *ct);
extern	DLL_EXPORT(char*)	DrlGetSuffix(char *s, char *ct);
extern	DLL_EXPORT(char*)	DrlGetPathDir(char *s, char *ct);



extern	DLL_EXPORT(int) 	DrlStrParseFuncScan(char *str,
					char *funcToken, int nVar,
					/*  char *varToken, char varType,
					 *  void* varPtr,
				 	 * int varDefFlag, void* varDef */
					...);
extern	DLL_EXPORT(char*)	DrlStrParseFuncPrint(char *s,
					char *funcToken, int nVar,
					/* char *varToken, char varType,
					 * void* varPtr,
					 * int varDefFlag, void* varDef */
					...);

extern	DLL_EXPORT(int)		DrlStrParseOptions(
		const char *s,		/* (I) input string */
		/* TVType varType1, void* varPtr1, char *varToken1,
		 *  ...
		 * TVType varTypeN, void* varPtrN, char *varTokenN,
		 * DRL_NULL_T)		LAST ARGUMENT MUST BE DRL_NULL_T
		 */ ...);


extern	DLL_EXPORT(int)		DrlStrTokScan(char *s, TVType varType,
					void *var);
extern	DLL_EXPORT(char*)	DrlStrTokPrint(char *s, TVType varType,
					void *var);

extern	DLL_EXPORT(int)		DrlStrSubsEnv(char *s);

/*
 *
 */
extern	DLL_EXPORT(int)		DrlStrLineVarScan(
	char *s,		/* (I) input string (unchanged) */
	/* TVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * TVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...);

extern	DLL_EXPORT(char*)	 DrlStrLineVarPrint(
	char *s,		/* (I) input string (unchanged) */
	/* TVType varType1, char *varName1, void* varPtr1,
	 * ...
	 * TVType varTypeN, char *varNameN, void* varPtrN,
	 * 0)			LAST ARGUMENT MUST BE ZERO
	 */ ...);



extern	DLL_EXPORT(int)		DrlStrLongValueScanV(char *s, char *errName,
					long *value, va_list ap);
extern	DLL_EXPORT(char*)	DrlStrLongValuePrintV(char *s, char *errName,
					long value, va_list ap);

extern	DLL_EXPORT(int)		DrlStrLongValueScan(
	char *s,		/* (I) input string (unchanged) */
	char *errName,		/* (I) variable name for debugging */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)		LAST ARGUMENT MUST BE NULL
	 */ ...);

extern	DLL_EXPORT(char*)	DrlStrLongValuePrint(
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
extern	DLL_EXPORT(char*)	DrlCurrentDateTimePrint(char *buf);
extern	DLL_EXPORT(char*)	DrlFloatPrint(char *s, double x, int ndigits);
extern	DLL_EXPORT(char*)	DrlDoubleIntervalPrint(char *s, double value);
extern	DLL_EXPORT(int)		DrlDoubleIntervalScan(char *s, double *value);
extern	DLL_EXPORT(int)		DrlCurScan(char *s, double *x);
extern	DLL_EXPORT(char*)	DrlCurPrint(char *s, double x, int ndigits);



#endif	/* _drlstr_H */
