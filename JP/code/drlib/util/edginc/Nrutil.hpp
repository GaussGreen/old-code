//----------------------------------------------------------------------------
//
//      File           : Nrutil.hpp
//
//      Description    : Numerical Recipes utility macors.
//----------------------------------------------------------------------------

#ifndef EDR_NR_UTILS_HPP
#define EDR_NR_UTILS_HPP

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

DRLIB_BEGIN_NAMESPACE

//void nrerror(char error_text[]);
UTIL_DLL float *fvector(long nl, long nh);
UTIL_DLL int *ivector(long nl, long nh);
UTIL_DLL unsigned char *cvector(long nl, long nh);
UTIL_DLL unsigned long *lvector(long nl, long nh);
UTIL_DLL double *dvector(long nl, long nh);
UTIL_DLL float **matrix(long nrl, long nrh, long ncl, long nch);
UTIL_DLL double **dmatrix(long nrl, long nrh, long ncl, long nch);
UTIL_DLL int **imatrix(long nrl, long nrh, long ncl, long nch);
UTIL_DLL float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
UTIL_DLL float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
UTIL_DLL float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
UTIL_DLL void free_vector(float *v, long nl, long nh);
UTIL_DLL void free_ivector(int *v, long nl, long nh);
UTIL_DLL void free_cvector(unsigned char *v, long nl, long nh);
UTIL_DLL void free_lvector(unsigned long *v, long nl, long nh);
UTIL_DLL void free_dvector(double *v, long nl, long nh);
UTIL_DLL void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
UTIL_DLL void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
UTIL_DLL void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
UTIL_DLL void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
UTIL_DLL void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
UTIL_DLL void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

DRLIB_END_NAMESPACE

#endif /* _NR_UTILS_H_ */
