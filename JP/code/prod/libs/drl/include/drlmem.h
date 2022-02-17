/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	MEM - Memory Management Routines
 * File:	drlmem.h
 * Author:	Christian Daher
 * Function:	
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlmem_H
#define	_drlmem_H
#include "drlstd.h"

/*
 * NRC style allocation support routines
 */

#define	VECT_ALLOC_ROUTINE(routine_name, type_name) \
extern	type_name*	routine_name(long nl, long nh);
#define	VECT_FREE_ROUTINE(routine_name, type_name) \
extern	void	routine_name(type_name *v, long nl, long nh);

#define	MATR_ALLOC_ROUTINE(routine_name, type_name) \
extern	type_name**	routine_name(long nrl, long nrh,\
		long ncl, long nch);
#define	MATR_FREE_ROUTINE(routine_name, type_name) \
extern	void	routine_name(type_name**m, long nrl, long nrh,\
	long ncl, long nch);

#define	ARRAY3D_ALLOC_ROUTINE(routine_name, type_name) \
extern	type_name*** routine_name(long nl, long nh, \
		long nrl, long nrh, long ncl, long nch);
#define	ARRAY3D_FREE_ROUTINE(routine_name, type_name) \
extern	void	routine_name(type_name***m, long nl, long nh, \
	long nrl, long nrh, long ncl, long nch);




VECT_ALLOC_ROUTINE(DrlVoidPVectAlloc, void*)
VECT_ALLOC_ROUTINE(DrlVoidPPVectAlloc, void**)
VECT_ALLOC_ROUTINE(DrlDoubleVectAlloc, double)
VECT_ALLOC_ROUTINE(DrlDoublePVectAlloc, double*)
VECT_ALLOC_ROUTINE(DrlDoublePPVectAlloc, double**)
VECT_ALLOC_ROUTINE(DrlFloatVectAlloc, float)
VECT_ALLOC_ROUTINE(DrlFloatPVectAlloc, float*)
VECT_ALLOC_ROUTINE(DrlFloatPPVectAlloc, float**)
VECT_ALLOC_ROUTINE(DrlIntVectAlloc, int)
VECT_ALLOC_ROUTINE(DrlIntPVectAlloc, int*)
VECT_ALLOC_ROUTINE(DrlIntPPVectAlloc, int**)
VECT_ALLOC_ROUTINE(DrlLongVectAlloc, long)
VECT_ALLOC_ROUTINE(DrlLongPVectAlloc, long*)
VECT_ALLOC_ROUTINE(DrlLongPPVectAlloc, long**)
VECT_ALLOC_ROUTINE(DrlCharVectAlloc, char)
VECT_ALLOC_ROUTINE(DrlCharPVectAlloc, char*)
VECT_ALLOC_ROUTINE(DrlCharPPVectAlloc, char**)
VECT_ALLOC_ROUTINE(DrlDDateVectAlloc, DDate)
VECT_ALLOC_ROUTINE(DrlDDatePVectAlloc, DDate*)
VECT_ALLOC_ROUTINE(DrlDDatePPVectAlloc, DDate**)
VECT_ALLOC_ROUTINE(DrlDIntervalVectAlloc, DInterval)



VECT_FREE_ROUTINE(DrlVoidPVectFree, void*)
VECT_FREE_ROUTINE(DrlVoidPPVectFree, void**)
VECT_FREE_ROUTINE(DrlDoubleVectFree, double)
VECT_FREE_ROUTINE(DrlDoublePVectFree, double*)
VECT_FREE_ROUTINE(DrlDoublePPVectFree, double**)
VECT_FREE_ROUTINE(DrlFloatVectFree, float)
VECT_FREE_ROUTINE(DrlFloatPVectFree, float*)
VECT_FREE_ROUTINE(DrlFloatPPVectFree, float**)
VECT_FREE_ROUTINE(DrlIntVectFree, int)
VECT_FREE_ROUTINE(DrlIntPVectFree, int*)
VECT_FREE_ROUTINE(DrlIntPPVectFree, int**)
VECT_FREE_ROUTINE(DrlLongVectFree, long)
VECT_FREE_ROUTINE(DrlLongPVectFree, long*)
VECT_FREE_ROUTINE(DrlLongPPVectFree, long**)
VECT_FREE_ROUTINE(DrlCharVectFree, char)
VECT_FREE_ROUTINE(DrlCharPVectFree, char*)
VECT_FREE_ROUTINE(DrlCharPPVectFree, char**)
VECT_FREE_ROUTINE(DrlDDateVectFree, DDate)
VECT_FREE_ROUTINE(DrlDDatePVectFree, DDate*)
VECT_FREE_ROUTINE(DrlDDatePPVectFree, DDate**)
VECT_FREE_ROUTINE(DrlDIntervalVectFree, DInterval)





MATR_ALLOC_ROUTINE(DrlVoidPMatrAlloc, void*)
MATR_ALLOC_ROUTINE(DrlVoidPPMatrAlloc, void**)
MATR_ALLOC_ROUTINE(DrlDoubleMatrAlloc, double)
MATR_ALLOC_ROUTINE(DrlDoublePMatrAlloc, double*)
MATR_ALLOC_ROUTINE(DrlDoublePPMatrAlloc, double**)
MATR_ALLOC_ROUTINE(DrlFloatMatrAlloc, float)
MATR_ALLOC_ROUTINE(DrlFloatPMatrAlloc, float*)
MATR_ALLOC_ROUTINE(DrlFloatPPMatrAlloc, float**)
MATR_ALLOC_ROUTINE(DrlIntMatrAlloc, int)
MATR_ALLOC_ROUTINE(DrlIntPMatrAlloc, int*)
MATR_ALLOC_ROUTINE(DrlIntPPMatrAlloc, int**)
MATR_ALLOC_ROUTINE(DrlLongMatrAlloc, long)
MATR_ALLOC_ROUTINE(DrlLongPMatrAlloc, long*)
MATR_ALLOC_ROUTINE(DrlLongPPMatrAlloc, long**)
MATR_ALLOC_ROUTINE(DrlCharMatrAlloc, char)
MATR_ALLOC_ROUTINE(DrlCharPMatrAlloc, char*)
MATR_ALLOC_ROUTINE(DrlCharPPMatrAlloc, char**)
MATR_ALLOC_ROUTINE(DrlDDateMatrAlloc, DDate)
MATR_ALLOC_ROUTINE(DrlDDatePMatrAlloc, DDate*)
MATR_ALLOC_ROUTINE(DrlDDatePPMatrAlloc, DDate**)



MATR_FREE_ROUTINE(DrlVoidPMatrFree, void*)
MATR_FREE_ROUTINE(DrlVoidPPMatrFree, void**)
MATR_FREE_ROUTINE(DrlDoubleMatrFree, double)
MATR_FREE_ROUTINE(DrlDoublePMatrFree, double*)
MATR_FREE_ROUTINE(DrlDoublePPMatrFree, double**)
MATR_FREE_ROUTINE(DrlFloatMatrFree, float)
MATR_FREE_ROUTINE(DrlFloatPMatrFree, float*)
MATR_FREE_ROUTINE(DrlFloatPPMatrFree, float**)
MATR_FREE_ROUTINE(DrlIntMatrFree, int)
MATR_FREE_ROUTINE(DrlIntPMatrFree, int*)
MATR_FREE_ROUTINE(DrlIntPPMatrFree, int**)
MATR_FREE_ROUTINE(DrlLongMatrFree, long)
MATR_FREE_ROUTINE(DrlLongPMatrFree, long*)
MATR_FREE_ROUTINE(DrlLongPPMatrFree, long**)
MATR_FREE_ROUTINE(DrlCharMatrFree, char)
MATR_FREE_ROUTINE(DrlCharPMatrFree, char*)
MATR_FREE_ROUTINE(DrlCharPPMatrFree, char**)
MATR_FREE_ROUTINE(DrlDDateMatrFree, DDate)
MATR_FREE_ROUTINE(DrlDDatePMatrFree, DDate*)
MATR_FREE_ROUTINE(DrlDDatePPMatrFree, DDate**)

ARRAY3D_ALLOC_ROUTINE(DrlDoubleArray3DAlloc, double)
ARRAY3D_ALLOC_ROUTINE(DrlIntArray3DAlloc, int)
ARRAY3D_ALLOC_ROUTINE(DrlLongArray3DAlloc, long)
ARRAY3D_ALLOC_ROUTINE(DrlDDateArray3DAlloc, DDate)

ARRAY3D_FREE_ROUTINE(DrlDoubleArray3DFree, double)
ARRAY3D_FREE_ROUTINE(DrlIntArray3DFree, int)
ARRAY3D_FREE_ROUTINE(DrlLongArray3DFree, long)
ARRAY3D_FREE_ROUTINE(DrlDDateArray3DFree, DDate)


#undef	VECT_ALLOC_ROUTINE
#undef	VECT_FREE_ROUTINE
#undef	MATR_ALLOC_ROUTINE
#undef	MATR_FREE_ROUTINE
#undef	ARRAY3D_ALLOC_ROUTINE
#undef	ARRAY3D_FREE_ROUTINE


/* Useful Macros */
#if !defined(DRL_CLIB)

#ifndef	NEW
#define	NEW(t)		((t *) MALLOC(sizeof(t)))
#endif

#ifndef	DELETE
#define	DELETE(ptr)	FREE((void*) (ptr))
#endif

#ifndef NEW_ARRAY
#define NEW_ARRAY(t,n)      (t *) MALLOC(sizeof(t)*(n))
#define FREE_ARRAY(ptr)     FREE(ptr)
#endif

#endif	/*DRL_CLIB*/

/*
 * Conversion between C and NRC conventions
 */


#define	MATR_NR2C_ROUTINE(routine_name, type_name)  \
void \
routine_name(type_name*** m, long nrl, long nrh, long ncl, long nch);

#define	MATR_C2NR_ROUTINE(routine_name, type_name)  \
void \
routine_name(type_name*** m, long nrl, long nrh, long ncl, long nch);



MATR_NR2C_ROUTINE(DrlDoubleMatrNrToC, double)

MATR_C2NR_ROUTINE(DrlDoubleMatrCToNr, double)


#undef	MATR_NR2C_ROUTINE
#undef	MATR_C2NR_ROUTINE


#endif	/* _drlmem_H */
