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
#include "drlstd.h"	/* wrapper type def (TVType, FloatL, etc.) */


/*c
 * This module contains a set of routines to manage vectors
 * of different types and LIL utility functions.
 */


/* Management of arrays of variable type (cast to void*) */

extern	DLL_EXPORT(char*)	DrlVTypeName(TVType type);
extern	DLL_EXPORT(int)		DrlVTypeNameScan(TVType *type, char *s);
extern	DLL_EXPORT(int)		DrlVTypeCheckValidCType(TVType type);
extern	DLL_EXPORT(int)		DrlVTypeCheckValidLType(TVType type);
extern	DLL_EXPORT(size_t)	DrlVTypeSizeof(TVType type);

extern	DLL_EXPORT(void*)	DrlVTypeVectAlloc(int size, TVType type);
extern	DLL_EXPORT(int)		DrlVTypeVectFree(void* p, int size,
					TVType type);
extern	DLL_EXPORT(void*)	DrlVTypeMatrAlloc(int nx, int ny,
					TVType type);
extern	DLL_EXPORT(int)		DrlVTypeMatrFree(void* p, int nx, int ny,
					TVType type);
extern	DLL_EXPORT(int)		DrlVTypeScan(char *str, TVType type, void *p);
extern	DLL_EXPORT(char*)	DrlVTypePrint(char *str, TVType type, void *p);
extern	DLL_EXPORT(char*)	DrlVTypePrintFmt(char *str, TVType type,
					void *p, char *fmt);
extern	DLL_EXPORT(char*)	DrlVTypePrintSmart(char *str, TVType type,
					void *p);
extern	DLL_EXPORT(void*)	DrlVTypeOffsetVect(void *p, int offset,
					TVType type);
extern	DLL_EXPORT(void*)	DrlVTypeOffsetMatrix(void *p, int nx, int ny,
					TVType type);
extern	DLL_EXPORT(int)		DrlVTypeLet(void *p, int offset1, 
					void *q, int offset2, TVType type);
extern	DLL_EXPORT(int)		DrlVTypeCompare(const void *p, const void *q,
					TVType type);
extern	DLL_EXPORT(int)		DrlVTypeNumValue(void *p, TVType type,
					void *arg, double *value);


extern	DLL_EXPORT(int)		DrlVTypeVectAdd(void **p, int *numItems,
					int *numItemsMax, int reallocSize,
					void *q, int removeDoubleItems,
					TVType type);
extern	DLL_EXPORT(int)		DrlVTypeVectSort(void *p, int *numItems,
					TVType type, int removeMultiple);


extern	DLL_EXPORT(void*)	DrlVTypeTensorAlloc(TVType type,
					int nDim, int *nSize);
extern	DLL_EXPORT(int)		DrlVTypeTensorFree(void *ptr,
					TVType type, int nDim, int *nSize);
extern	DLL_EXPORT(void*)	DrlVTypeTensorAllocV(TVType type, int nDim,
					... /* int size1, ..., sizeNDim */);
extern	DLL_EXPORT(int)		DrlVTypeTensorFreeV(void *ptr, TVType type,
					int nDim,
					... /* int size1, ..., sizeNDim */);
extern	DLL_EXPORT(void*)	DrlVTypeTensorAccessElement(void *ptr,
					TVType type, int nDim, int *nSize,
					int *idx);
extern	DLL_EXPORT(int)		DrlVTypeTensorOperScalar(void *ptr,
					TVType type, int nDim, int *nSize,
					char *operName, double *value);
extern	DLL_EXPORT(int)		DrlVTypeTensorOperScalarV(
					void *ptr, TVType type,
					char *operName, double *value,
					int nDim,
					... /* int size1,..., int sizeNDim */);




/* LIL interface */

extern	DLL_EXPORT(int)	DrlLilToCType(TVType varType);
extern	DLL_EXPORT(int)	DrlLilVectSize(TVType varType,
				void *lptr, char *argName);
extern	DLL_EXPORT(int)	DrlLilVectActiveSize(TVType varType,
				void *lptr, char *argName);


extern	DLL_EXPORT(int)	DrlLilVectGet(void *lptr, TVType varType,
				int loffset, int lskip, int n, void *cptr);
extern	DLL_EXPORT(int)	DrlLilVectPut(void *lptr, TVType varType, 
				int loffset, int lskip, int n, void *cptr);
extern	DLL_EXPORT(int)	DrlLilMatrGet(void *lptr, TVType varType,
				int n1, int n2, void *cptr);
extern	DLL_EXPORT(int)	DrlLilMatrPut(void *lptr, TVType varType,
				int n1, int n2, void *cptr);

extern	DLL_EXPORT(int)	DrlLilStructGet(int iFlag,
/* VAR_L, 
 *	char *name,  TVType type,  void *lvar,  void *cvar, 
 * VECT_L, int *nItems, intnMinItems, int nMaxItems,
 *	char *name,  TVType type,  int alloc,  void *lvar,  void *cvar, 
 * MATR_L, int nx, int ny,
 *	char *name,  TVType type,  int alloc,  void *lvar,  void *cvar, 
 * VECTARRAY_T, int *nItems, intnMinItems, intMaxItems, int nVar,
 *	char *name1, TVType type1, int alloc1, void* lvar1, void *cvar1,
 *	...
 *	char *nameN, TVType typeN, int allocN, void* lvarN, void *cvarN,
 * VECTRANGE_T, int *nItems, intMinItems, intMaxItems, int nVar,
 *	TVType type, void* lvar,
 *	char *name1, int alloc1, void *cvar1,
 *	...
 *	char *nameN, int allocN, void *cvarN,
 * 0)
 */
...);

extern	DLL_EXPORT(int)	DrlLilStructPut(int iFlag,
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


extern	DLL_EXPORT(int)	DrlLilVectCopy(void *lptr, void *rptr,
				TVType varType, int n);
extern	DLL_EXPORT(int)	DrlLilVectLogging(FILE *fp, TVType varType, void *lptr,
				char  *argName);
extern	DLL_EXPORT(int)	DrlLilVectLoggingFile(
	char *fnam,		/* (I) file name */
	char *mode,		/* (I) write mode */
	char *funcName,		/* (I) function name */
	/* TVType varType, void *lptr, char  *argName,
	 * ...
	 * TVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...);

extern	DLL_EXPORT(int)	DrlLilVectLoggingFp(
	FILE *fp,		/* (I) file pointer to write to (or NULL) */
	char *funcName,		/* (I) function name */
	/* TVType varType, void *lptr, char  *argName,
	 * ...
	 * TVType varType, void *lptr, char  *argName,
	 * 0)	LAST ARGUMENT MUST BE ZERO
	 */
	...);

extern	DLL_EXPORT(int)	DrlLilVectArrayFpRead(
	FILE *fp,		/* (I) file pointer */
	int numItems,		/* (I) num elements in each column */
	int numVect,		/* (I) num of columns */
	TVType *varType,	/* (I) variable types [0..numVect-1] */
	void ***lptr);		/* (I) array of & of pointers [0..numVect-1] */

extern	DLL_EXPORT(int)	DrlLilVectArrayFpReadV(
	FILE *fp,		/* (I) file pointer */
	int numItems,		/* (I) num elements to be read */
	/* TVType varType, void *lptr,
	 * ...
	 * TVType varType, void *lptr,
	 * DRL_NULL_T (last argument MUST be DRL_NULL_T)
	 */
	...);

extern	DLL_EXPORT(int)	DoubleIsFinite(double x);


#endif	/*_drlvtype_H*/



