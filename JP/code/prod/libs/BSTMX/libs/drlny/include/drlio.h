/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	IO - File I/O Routines
 * File:	drlio.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlio_H
#define	_drlio_H
#include "drlstd.h"

#include <stddef.h>
#include <stdio.h>

#if defined(CLIB) && (defined(_WINDLL) && !(defined(WIN32) || defined(WIN32)))
# include "cfileio.h"
#endif

/* Win3.1 does not have fscanf */
#if defined(CLIB) && defined(_WINDLL) && !(defined(WIN32) || defined(_WIN32))
extern	int	_drlDummyFscanf(FILE *fp, char *fmt,...);
# define	fscanf		_dummyFscanf
#endif


/*
 * Formatted files I/O utilities.
 */
extern	DLL_EXPORT(int)		DrlFPrintf(FILE *fp, char *fmt,...);

extern	DLL_EXPORT(char*)	DrlFGetLine(char *s, int n,
					FILE *fp, int *line);
extern	DLL_EXPORT(int)		DrlFGetToken(FILE *fp, char *sep,
					char *s, int sz);
extern	DLL_EXPORT(int)		DrlFAdvanceToNextChar(FILE *fp);
extern	DLL_EXPORT(int)		DrlFAdvanceToToken(FILE *fp, char *token);
extern	DLL_EXPORT(int)		DrlFScanDouble(FILE *fp, double *val);
extern	DLL_EXPORT(int)		DrlFScanLong(FILE *fp, long *val);
extern	DLL_EXPORT(int)		DrlFScanInt(FILE *fp, int *val);
extern	DLL_EXPORT(int)		DrlFScanString(FILE *fp, char *s);
extern	DLL_EXPORT(int)		DrlFScanVType(FILE *fp, TVType type,
					void *valPtr);
extern	DLL_EXPORT(int)		DrlFPrintVType(FILE *fp, TVType type,
					void *valPtr);


extern	DLL_EXPORT(int)		DrlFScanStringLongValue(
	FILE *fp,		/* (I) input file */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...);


extern	DLL_EXPORT(int)		DrlFPrintStringLongValue(
	FILE *fp,		/* (O) output file */
	long value,		/* (I) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...);

extern	DLL_EXPORT(int)		DrlFScanStruct(
	FILE *fp,
	int *line,
	char *toks,
	/* VAR_T, char *tagName, TVType varType, void *var,
	 * ARRAY_T, char *tagName, int *nItems, int nVar,
	 *	TVType varType1, TVType varAllocType1, void* var1,
	 *	...
	 *	TVType varTypeN, TVType varAllocTypeN, void* varN,
	 * DVECTOR_T, char *tagName,
	 *	int *nItems, double **v,
	 * DMATRIX_T, char *tagName,
	 *	int *nx, int *ny, double ***v,
	 * 0)
	 */ ...);

extern	DLL_EXPORT(int)		DrlFPrintStruct(
	FILE *fp,
	int *line,
	int token,
	/* VAR_T, char *tagName, TVType varType, void *var,
	 * ARRAY_T, char *tagName, int nItems, int nVar,
	 *	TVType varType1, void* var1,
	 *	...
	 *	TVType varTypeN, void* varN,
	 * DVECTOR_T, char *tagName,
	 * 	int nItems, double *v,
	 * DMATRIX_T, char *tagName,
	 *	int nx, int ny, double **v,
	 * 0)
	 */ ...);


extern	DLL_EXPORT(int)		DrlFScanDoubleVect(FILE *fp,
					double **x, int *n);
extern	DLL_EXPORT(int)		DrlFPrintDoubleVect(FILE *fp,
					char *fmt, double *x, int n);
extern	DLL_EXPORT(int)		DrlFScanDoubleMatr(FILE *fp,
					double ***x, int *n1, int *n2);
extern	DLL_EXPORT(int)		DrlFPrintDoubleMatr(FILE *fp,
					char *fmt, double **x, int n1, int n2);

/*
 * TeX output
 */

extern	DLL_EXPORT(int)		DrlTeXStructPrint(FILE *fp,
	/*
	 * ARRAY_T, int nItems, int nVar,
	 *   char *varName1, char *tabSep1, TVType varType1,
	 * 				char *varFmt1, void* var1,
	 *   ...
	 *   char *varNameN, char *tabSepN, TVType varTypeN,
	 * 				char *varFmtN, void* varN,
	 * MATRIX_T,
	 *   int nx, TVType varTypeX, char *varFmtX, void* varX,
	 *   int ny, TVType varTypeY, char *varFmtY, void* varY,
	 *   TVType varTypeXY, char *varFmtXY, void* varXY,
	 * DMATRIX_T,
	 *   int nx, int ny, double *vx, double *vy, double **v,
	 * 0)
	 */ ...);



/*
 * Provides stdout and stderr for Windows.
 */

#if defined(_WINDLL) && !(defined(WIN32) || defined(WIN32))
extern	DLL_EXPORT(FILE*)	DrlStdoutFp(void);
# define stdout	(DrlStdoutFp())
#endif
extern	DLL_EXPORT(int)		DrlStdoutFileSet(char *pathname, char *mode);
extern	DLL_EXPORT(int)		DrlStdoutClose(void);
extern	DLL_EXPORT(int)		DrlStderrFileSet(char *pathname, char *mode);
extern	DLL_EXPORT(int)		DrlStderrClose(void);


extern	DLL_EXPORT(FILE*)	DrlStdlogFp(int fpNo);
extern	DLL_EXPORT(int)		DrlStdlogFileSet(int fpNo, char *fnam);

/*
 * Other
 */

extern	DLL_EXPORT(int)		DrlFilePrintf(char *fnam, char *fmt,...);


#endif	/* _drlio_H */
