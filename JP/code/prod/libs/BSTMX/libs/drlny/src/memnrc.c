/************************************************************************
 * Module:	DCU - MEM
 * Function:	Memory Utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include <stddef.h>
#include <stdio.h>
#ifdef	UNIX
# include <stdlib.h>		/* malloc(), free() */
#endif 

#ifdef	CLIB
# include "cerror.h"
# include "cmemory.h"
#endif

#include "drlmem.h"	/* for prototype consistency */

#ifdef	UNIX
# undef  MALLOC
# undef  FREE
# define MALLOC(sz)	malloc((size_t) (sz))
# define FREE(ptr)	free((char*) ptr)
#endif	/*UNIX*/


/* these have a different meaning in the header file */
#undef	VECT_ALLOC_ROUTINE
#undef	VECT_FREE_ROUTINE
#undef	MATR_ALLOC_ROUTINE
#undef	MATR_FREE_ROUTINE

#undef	MATR_NR2C_ROUTINE
#undef	MATR_C2NR_ROUTINE

/* offset */
#define NR_END 0



/*----------------------------------------------------------------------
 *
 */

#define	MATR_NR2C_ROUTINE(routine_name, type_name)  \
DLL_EXPORT(void) \
routine_name(type_name*** m, long nrl, long nrh, long ncl, long nch) \
{ \
	long	i; \
	for(i=nrl; i<=nrh; i++) { \
	    (*m)[i] += ncl; \
	    (*m)[i] -= NR_END; \
	} \
	(*m) += nrl; \
	(*m) -= NR_END; \
}


MATR_NR2C_ROUTINE(DrlDoubleMatrNrToC, double)


/*----------------------------------------------------------------------
 * Converts a C indexed matrix to a NRC index matrix.
 */

#define	MATR_C2NR_ROUTINE(routine_name, type_name)  \
DLL_EXPORT(void) \
routine_name(type_name*** m, long nrl, long nrh, long ncl, long nch) \
{ \
	long	i; \
	(*m) += NR_END; \
	(*m) -= nrl; \
	for(i=nrl; i<=nrh; i++) { \
	    (*m)[i] -= ncl; \
	    (*m)[i] += NR_END; \
	} \
}



MATR_C2NR_ROUTINE(DrlDoubleMatrCToNr, double)



#ifdef	_COMMENT

/*c---------------------------------------------------------------------
 * Vector and matrix memory allocation routines.
 *
 * <br><br>
 * The following table gives the name of the routine
 * corresponding to the type allocated.
 * For example, the routine name that allocates
 * a matrix of pointer to void is 
 * <i> VoidPMatrAlloc</i>, and the routine that frees a vector
 * of dates is <i> TDateVectFree</i>.
 * <br><br>
 * <table COLS=2>
 * <tr> <td>Routine name  </td> <td> C type allocated </td> </tr>
 * <tr> <td>DrlVoidP	 </td> <td> void*	</td> </tr>
 * <tr> <td>DrlVoidPP	 </td> <td> void**	</td> </tr>
 * <tr> <td>DrlDouble	 </td> <td> double	</td> </tr>
 * <tr> <td>DrlDoubleP	 </td> <td> double*	</td> </tr>
 * <tr> <td>DrlDoublePP	 </td> <td> double**	</td> </tr>
 * <tr> <td>DrlFloat	 </td> <td> float	</td> </tr>
 * <tr> <td>DrlFloatP	 </td> <td> float*	</td> </tr>
 * <tr> <td>DrlFloatPP	 </td> <td> float**	</td> </tr>
 * <tr> <td>DrlInt	 </td> <td> int		</td> </tr>
 * <tr> <td>DrlIntP	 </td> <td> int*	</td> </tr>
 * <tr> <td>DrlIntPP	 </td> <td> int**	</td> </tr>
 * <tr> <td>DrlLong	 </td> <td> long	</td> </tr>
 * <tr> <td>DrlLongP	 </td> <td> long*	</td> </tr>
 * <tr> <td>DrlLongPP	 </td> <td> long**	</td> </tr>
 * <tr> <td>DrlChar	 </td> <td> char	</td> </tr>
 * <tr> <td>DrlCharP	 </td> <td> char*	</td> </tr>
 * <tr> <td>DrlCharPP	 </td> <td> char**	</td> </tr>
 * <tr> <td>DrlTDate	 </td> <td> TDate	</td> </tr>
 * <tr> <td>DrlTDateP	 </td> <td> TDate*	</td> </tr>
 * <tr> <td>DrlTDatePP	 </td> <td> TDate**	</td> </tr>
 * </table>
 */

extern	<type>*	 Drl<Type>VectAlloc(long nl, long nh);
extern	void	 Drl<Type>VectFree(<type> *v, long nl, long nh);
extern	<type>** Drl<Type>MatrAlloc(long nrl, long nrh, long ncl, long nch);
extern	void	 Drl<Type>MatrFree(<type> **m, long nrl, long nrh,
			long ncl, long nch);

/*e*/
#endif /*_COMMENT*/

/* allocate an vector with subscript range v[nl..nh] */
#define	VECT_ALLOC_ROUTINE(routine_name, type_name) \
DLL_EXPORT(type_name*) routine_name(long nl, long nh) \
{ \
	type_name *v; \
	v = (type_name*) \
		MALLOC((size_t) ((nh-nl+1+NR_END)*sizeof(type_name))); \
	if (!v) { \
		GtoErrMsg("%s: alloc failure (size=%d)\n", \
			#routine_name, (nh-nl+1+NR_END)); \
		return (NULL); \
	} \
	return(v-nl+NR_END); \
}



/* free a vector allocated with vector() */
#define	VECT_FREE_ROUTINE(routine_name, type_name) \
DLL_EXPORT(void) routine_name(type_name *v, long nl, long nh) \
{ \
	if (v) \
		FREE((void*) (v+nl-NR_END)); \
}


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
VECT_ALLOC_ROUTINE(DrlTDateVectAlloc, TDate)
VECT_ALLOC_ROUTINE(DrlTDatePVectAlloc, TDate*)
VECT_ALLOC_ROUTINE(DrlTDatePPVectAlloc, TDate**)
VECT_ALLOC_ROUTINE(DrlTDateIntervalVectAlloc, TDateInterval)



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
VECT_FREE_ROUTINE(DrlTDateVectFree, TDate)
VECT_FREE_ROUTINE(DrlTDatePVectFree, TDate*)
VECT_FREE_ROUTINE(DrlTDatePPVectFree, TDate**)
VECT_FREE_ROUTINE(DrlTDateIntervalVectFree, TDateInterval)







/*----------------------------------------------------------------------
 * Matrix Allocation routines:
 * allocate a matrix with subscript range m[nrl..nrh][ncl..nch]
 */

#undef	_NRC_TYPE
#ifdef	_NRC_TYPE

#define	MATR_ALLOC_ROUTINE(routine_name, type_name) \
DLL_EXPORT(type_name**) routine_name(long nrl, long nrh, long ncl, long nch) \
{ \
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; \
	type_name **m; \
 \
	/* allocate pointers to rows */ \
	m=(type_name**) MALLOC((size_t)((nrow+NR_END)*sizeof(type_name*))); \
	if (!m) { \
		GtoErrMsg("%s: alloc failure 1 (nrow=%d)\n", \
			#routine_name, nrow); \
		return(NULL); \
	} \
	m += NR_END; \
	m -= nrl; \
 \
	/* allocate rows and set pointers to them */ \
	m[nrl]=(type_name*) MALLOC((size_t)((nrow*ncol+NR_END)*sizeof(type_name))); \
	if (!m[nrl]) { \
		GtoErrMsg("%s: alloc failure 2 (ncol=%d)\n", \
			#routine_name, ncol); \
		return(NULL); \
	} \
	m[nrl] += NR_END; \
	m[nrl] -= ncl; \
 \
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol; \
 \
	return(m); \
}




#define	MATR_FREE_ROUTINE(routine_name, type_name) \
DLL_EXPORT(void) routine_name(type_name**m, long nrl, long nrh, long ncl, long nch) \
{ \
	if (m) { \
		FREE((void*) (m[nrl]+ncl-NR_END)); \
		FREE((void*) (m+nrl-NR_END)); \
	} \
} 


#else
/*
 *
 */

#define	MATR_ALLOC_ROUTINE(routine_name, type_name) \
DLL_EXPORT(type_name**) routine_name(long nrl, long nrh, long ncl, long nch) \
{ \
	long	i, nrow=nrh-nrl+1,ncol=nch-ncl+1; \
	size_t	sz; \
	type_name	**m; \
 \
	sz = (size_t) (nrow+NR_END)*sizeof(type_name*); \
	if ((m = (type_name**)MALLOC(sz)) == NULL) { \
		GtoErrMsg("%s: alloc failure 1 (nrow=%d)\n", \
			#routine_name, nrow);  \
		return(NULL); \
	} \
	m += NR_END; \
	m -= nrl; \
 \
	sz = (size_t) (ncol+NR_END)*sizeof(type_name); \
	for(i=nrl; i<=nrh; i++) { \
	    if ((m[i] = (type_name*)MALLOC(sz)) == NULL) { \
		GtoErrMsg("%s: alloc failure 2 (ncol=%d)\n", \
			#routine_name, ncol);  \
		return(NULL); \
	    } \
	    m[i] += NR_END; \
	    m[i] -= ncl; \
	} \
	return(m);  \
}



#define	MATR_FREE_ROUTINE(routine_name, type_name) \
DLL_EXPORT(void) routine_name(type_name**m, long nrl, long nrh, long ncl, long nch) \
{ \
	long	i; \
	if (m) { \
	    for(i=nrl; i<=nrh; i++) { \
		FREE((void*) (m[i]+ncl-NR_END)); \
	    } \
	    FREE((void*) (m+nrl-NR_END)); \
	} \
}


#endif



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
MATR_ALLOC_ROUTINE(DrlTDateMatrAlloc, TDate)
MATR_ALLOC_ROUTINE(DrlTDatePMatrAlloc, TDate*)
MATR_ALLOC_ROUTINE(DrlTDatePPMatrAlloc, TDate**)



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
MATR_FREE_ROUTINE(DrlTDateMatrFree, TDate)
MATR_FREE_ROUTINE(DrlTDatePMatrFree, TDate*)
MATR_FREE_ROUTINE(DrlTDatePPMatrFree, TDate**)




/*----------------------------------------------------------------------
 * Allocates a 3D array of doubles (nl, nh) x (nrl, nrh) x (ncl, nch).
 *
 */


#define	ARRAY3D_ALLOC_ROUTINE(routine_name, type_name) \
DLL_EXPORT(type_name***) routine_name(long nl, long nh, \
	long nrl, long nrh, long ncl, long nch) \
{ \
	long	i, j, n=nh-nl+1,nrow=nrh-nrl+1,ncol=nch-ncl+1; \
	size_t	sz; \
	type_name	***m; \
 \
	sz = (size_t) (n+NR_END)*sizeof(type_name**); \
	if ((m = (type_name***)MALLOC(sz)) == NULL) { \
		GtoErrMsg("%s: alloc failure 1 (n=%d)\n", \
			#routine_name, n);  \
		return(NULL); \
	} \
	m += NR_END; \
	m -= nl; \
 \
	sz = (size_t) (nrow+NR_END)*sizeof(type_name*); \
	for(i=nl; i<=nh; i++) { \
	    if ((m[i] = (type_name**)MALLOC(sz)) == NULL) { \
		GtoErrMsg("%s: alloc failure 2 (nrow=%d)\n", \
			#routine_name, nrow);  \
		return(NULL); \
	    } \
	    m[i] += NR_END; \
	    m[i] -= nrl; \
 \
	    sz = (size_t) (ncol+NR_END)*sizeof(type_name); \
	    for(j=nrl; j<=nrh; j++) { \
	        if ((m[i][j] = (type_name*)MALLOC(sz)) == NULL) { \
		    GtoErrMsg("%s: alloc failure 2 (ncol=%d)\n", \
			    #routine_name, ncol);  \
		    return(NULL); \
	        } \
	        m[i][j] += NR_END; \
	        m[i][j] -= ncl; \
	    } \
	} \
	return(m);  \
}



#define	ARRAY3D_FREE_ROUTINE(routine_name, type_name) \
DLL_EXPORT(void) routine_name(type_name***m, \
	long nl, long nh, long nrl, long nrh, long ncl, long nch) \
{ \
	long	i, j; \
	if (m) { \
	    for(i=nl; i<=nh; i++) { \
	        for(j=nrl; j<=nrh; j++) { \
		    FREE((void*) (m[i][j]+ncl-NR_END)); \
		} \
		FREE((void*) (m[i]+nrl-NR_END)); \
	    } \
	    FREE((void*) (m+nl-NR_END)); \
	} \
}




ARRAY3D_ALLOC_ROUTINE(DrlDoubleArray3DAlloc, double)
ARRAY3D_ALLOC_ROUTINE(DrlIntArray3DAlloc, int)
ARRAY3D_ALLOC_ROUTINE(DrlLongArray3DAlloc, long)
ARRAY3D_ALLOC_ROUTINE(DrlTDateArray3DAlloc, TDate)

ARRAY3D_FREE_ROUTINE(DrlDoubleArray3DFree, double)
ARRAY3D_FREE_ROUTINE(DrlIntArray3DFree, int)
ARRAY3D_FREE_ROUTINE(DrlLongArray3DFree, long)
ARRAY3D_FREE_ROUTINE(DrlTDateArray3DFree, TDate)


