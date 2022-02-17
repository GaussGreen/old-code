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

#if defined(DRL_CLIB) && (defined(_WINDLL) && !(defined(WIN32) || defined(WIN32)))
# include "cfileio.h"
#endif

/* Win3.1 does not have fscanf */
#if defined(DRL_CLIB) && defined(_WINDLL) && !(defined(WIN32) || defined(_WIN32))
extern	int	_drlDummyFscanf(FILE *fp, char *fmt,...);
# define	fscanf		_dummyFscanf
#endif


/*
 * Formatted files I/O utilities.
 */
extern	int		DrlFPrintf(FILE *fp, char *fmt,...);

extern	char*	DrlFGetLine(char *s, int n,
					FILE *fp, int *line);
extern	int		DrlFGetToken(FILE *fp, char *sep,
					char *s, int sz);
extern	int		DrlFAdvanceToNextChar(FILE *fp);
extern	int		DrlFAdvanceToToken(FILE *fp, char *token);
extern	int		DrlFScanDouble(FILE *fp, double *val);
extern	int		DrlFScanLong(FILE *fp, long *val);
extern	int		DrlFScanInt(FILE *fp, int *val);
extern	int		DrlFScanString(FILE *fp, char *s);
extern	int		DrlFScanVType(FILE *fp, DVType type,
					void *valPtr);
extern	int		DrlFPrintVType(FILE *fp, DVType type,
					void *valPtr);


extern	int		DrlFScanStringLongValue(
	FILE *fp,		/* (I) input file */
	long *value,		/* (O) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...);


extern	int		DrlFPrintStringLongValue(
	FILE *fp,		/* (O) output file */
	long value,		/* (I) value */
	/* char *varName1, long varValue1,
	 * ...
	 * char *varNameN, long varValueN,
	 * NULL)			LAST ARGUMENT MUST BE NULL
	 */ ...);

extern	int		DrlFScanStruct(
	FILE *fp,
	int *line,
	char *toks,
	/* VAR_T, char *tagName, DVType varType, void *var,
	 * ARRAY_T, char *tagName, int *nItems, int nVar,
	 *	DVType varType1, DVType varAllocType1, void* var1,
	 *	...
	 *	DVType varTypeN, DVType varAllocTypeN, void* varN,
	 * DVECTOR_T, char *tagName,
	 *	int *nItems, double **v,
	 * DMATRIX_T, char *tagName,
	 *	int *nx, int *ny, double ***v,
	 * 0)
	 */ ...);

extern	int		DrlFPrintStruct(
	FILE *fp,
	int *line,
	int token,
	/* VAR_T, char *tagName, DVType varType, void *var,
	 * ARRAY_T, char *tagName, int nItems, int nVar,
	 *	DVType varType1, void* var1,
	 *	...
	 *	DVType varTypeN, void* varN,
	 * DVECTOR_T, char *tagName,
	 * 	int nItems, double *v,
	 * DMATRIX_T, char *tagName,
	 *	int nx, int ny, double **v,
	 * 0)
	 */ ...);


extern	int		DrlFScanDoubleVect(FILE *fp,
					double **x, int *n);
extern	int		DrlFPrintDoubleVect(FILE *fp,
					char *fmt, double *x, int n);
extern	int		DrlFScanDoubleMatr(FILE *fp,
					double ***x, int *n1, int *n2);
extern	int		DrlFPrintDoubleMatr(FILE *fp,
					char *fmt, double **x, int n1, int n2);

/*
 * TeX output
 */

extern	int		DrlTeXStructPrint(FILE *fp,
	/*
	 * ARRAY_T, int nItems, int nVar,
	 *   char *varName1, char *tabSep1, DVType varType1,
	 * 				char *varFmt1, void* var1,
	 *   ...
	 *   char *varNameN, char *tabSepN, DVType varTypeN,
	 * 				char *varFmtN, void* varN,
	 * MATRIX_T,
	 *   int nx, DVType varTypeX, char *varFmtX, void* varX,
	 *   int ny, DVType varTypeY, char *varFmtY, void* varY,
	 *   DVType varTypeXY, char *varFmtXY, void* varXY,
	 * DMATRIX_T,
	 *   int nx, int ny, double *vx, double *vy, double **v,
	 * 0)
	 */ ...);



/*
 * Provides stdout and stderr for Windows.
 */

#if defined(_WINDLL) && !(defined(WIN32) || defined(WIN32))
extern	FILE*	DrlStdoutFp(void);
# define stdout	(DrlStdoutFp())
#endif
extern	int		DrlStdoutFileSet(char *pathname, char *mode);
extern	int		DrlStdoutClose(void);
extern	int		DrlStderrFileSet(char *pathname, char *mode);
extern	int		DrlStderrClose(void);


extern	FILE*	DrlStdlogFp(int fpNo);
extern	int		DrlStdlogFileSet(int fpNo, char *fnam);

/*
 * Other
 */

extern	int		DrlFilePrintf(char *fnam, char *fmt,...);


#endif	/* _drlio_H */
