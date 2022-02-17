/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	VTYPE - Variable Type
 * File:	dcuvtype.h
 * Function:	Wrapper utility routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlvtype_H
#define	_drlvtype_H
#include "drlstd.h"	/* wrapper type def (DVType, FloatL, etc.) */
#include <stdio.h>


/*c
 * This module contains a set of routines to manage vectors
 * of different types and LIL utility functions.
 */


/* Management of arrays of variable type (cast to void*) */

extern	char*	DrlVTypeName(DVType type);
extern	int		DrlVTypeNameScan(DVType *type, char *s);
extern	int		DrlVTypeCheckValidCType(DVType type);
extern	int		DrlVTypeCheckValidLType(DVType type);
extern	size_t	DrlVTypeSizeof(DVType type);

extern	void*	DrlVTypeVectAlloc(int size, DVType type);
extern	int		DrlVTypeVectFree(void* p, int size,
					DVType type);
extern	void*	DrlVTypeMatrAlloc(int nx, int ny,
					DVType type);
extern	int		DrlVTypeMatrFree(void* p, int nx, int ny,
					DVType type);
extern	int		DrlVTypeScan(char *str, DVType type, void *p);
extern	char*	DrlVTypePrint(char *str, DVType type, void *p);
extern	char*	DrlVTypePrintFmt(char *str, DVType type,
					void *p, char *fmt);
extern	char*	DrlVTypePrintSmart(char *str, DVType type,
					void *p);
extern	void*	DrlVTypeOffsetVect(void *p, int offset,
					DVType type);
extern	void*	DrlVTypeOffsetMatrix(void *p, int nx, int ny,
					DVType type);
extern	int		DrlVTypeLet(void *p, int offset1, 
					void *q, int offset2, DVType type);
extern	int		DrlVTypeCompare(const void *p, const void *q,
					DVType type);
extern	int		DrlVTypeNumValue(void *p, DVType type,
					void *arg, double *value);


extern	int		DrlVTypeVectAdd(void **p, int *numItems,
					int *numItemsMax, int reallocSize,
					void *q, int removeDoubleItems,
					DVType type);
extern	int		DrlVTypeVectSort(void *p, int *numItems,
					DVType type, int removeMultiple);


extern	void*	DrlVTypeTensorAlloc(DVType type,
					int nDim, int *nSize);
extern	int		DrlVTypeTensorFree(void *ptr,
					DVType type, int nDim, int *nSize);
extern	void*	DrlVTypeTensorAllocV(DVType type, int nDim,
					... /* int size1, ..., sizeNDim */);
extern	int		DrlVTypeTensorFreeV(void *ptr, DVType type,
					int nDim,
					... /* int size1, ..., sizeNDim */);
extern	void*	DrlVTypeTensorAccessElement(void *ptr,
					DVType type, int nDim, int *nSize,
					int *idx);
extern	int		DrlVTypeTensorOperScalar(void *ptr,
					DVType type, int nDim, int *nSize,
					char *operName, double *value);
extern	int		DrlVTypeTensorOperScalarV(
					void *ptr, DVType type,
					char *operName, double *value,
					int nDim,
					... /* int size1,..., int sizeNDim */);




/* LIL interface */

extern	int	DrlLilToCType(DVType varType);
extern	int	DrlLilVectSize(DVType varType,
				void *lptr, char *argName);
extern	int	DrlLilVectActiveSize(DVType varType,
				void *lptr, char *argName);


extern	int	DrlLilVectGet(void *lptr, DVType varType,
				int loffset, int lskip, int n, void *cptr);
extern	int	DrlLilVectPut(void *lptr, DVType varType, 
				int loffset, int lskip, int n, void *cptr);
extern	int	DrlLilMatrGet(void *lptr, DVType varType,
				int n1, int n2, void *cptr);
extern	int	DrlLilMatrPut(void *lptr, DVType varType,
				int n1, int n2, void *cptr);

extern	int	DrlLilStructGet(int iFlag,
/* VAR_L, 
 *	char *name,  DVType type,  void *lvar,  void *cvar, 
 * VECT_L, int *nItems, intnMinItems, int nMaxItems,
 *	char *name,  DVType type,  int alloc,  void *lvar,  void *cvar, 
 * MATR_L, int nx, int ny,
 *	char *name,  DVType type,  int alloc,  void *lvar,  void *cvar, 
 * VECTARRAY_T, int *nItems, intnMinItems, intMaxItems, int nVar,
 *	char *name1, DVType type1, int alloc1, void* lvar1, void *cvar1,
 *	...
 *	char *nameN, DVType typeN, int allocN, void* lvarN, void *cvarN,
 * VECTRANGE_T, int *nItems, intMinItems, intMaxItems, int nVar,
 *	DVType type, void* lvar,
 *	char *name1, int alloc1, void *cvar1,
 *	...
 *	char *nameN, int allocN, void *cvarN,
 * 0)
 */
...);

extern	int	DrlLilStructPut(int iFlag,
	/* VAR_L, 
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * VECT_L, int *nItems, intMinItems, int nMaxItems,
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * MATR_L, int nx, int ny,
	 *	char *name,  int type,  void *lvar,  void *cvar, 
	 * VECTARRAY_T, int *nItems, intnMinItems, intMaxItems, int nVar,
	 *	char *name1, int type1, void* lvar1, void *cvar1,
	 *	...
	 *	char *nameN, int typeN, void* lvarN, void *cvarN,
	 * 0)
	 */
	...);


extern	int	DrlLilVectCopy(void *lptr, void *rptr,
				DVType varType, int n);
extern	int	DrlLilVectLogging(FILE *fp, DVType varType, void *lptr,
				char  *argName);
extern	int	DrlLilVectLoggingFile(
	char *fnam,		/* (I) file name */
	char *mode,		/* (I) write mode */
	char *funcName,		/* (I) function name */
	/* DVType varType, void *lptr, char  *argName,
	 * ...
	 * DVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...);

extern	int	DrlLilVectLoggingFp(
	FILE *fp,		/* (I) file pointer to write to (or NULL) */
	char *funcName,		/* (I) function name */
	/* DVType varType, void *lptr, char  *argName,
	 * ...
	 * DVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...);

extern	int	DrlLilVectArrayFpRead(
	FILE *fp,		/* (I) file pointer */
	int numItems,		/* (I) num elements in each column */
	int numVect,		/* (I) num of columns */
	DVType *varType,	/* (I) variable types [0..numVect-1] */
	void ***lptr);		/* (I) array of & of pointers [0..numVect-1] */

extern	int	DrlLilVectArrayFpReadV(
	FILE *fp,		/* (I) file pointer */
	int numItems,		/* (I) num elements to be read */
	/* DVType varType, void *lptr,
	 * ...
	 * DVType varType, void *lptr,
	 * DRL_NULL_T (last argument MUST be DRL_NULL_T)
	 */
	...);

extern	int	DrlDoubleIsFinite(double x);


#endif	/*_drlvtype_H*/



