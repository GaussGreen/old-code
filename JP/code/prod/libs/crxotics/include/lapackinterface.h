/*********************************************************************************
 * LaPack Interface for MKL
 *
 ********************************************************************************/

#ifndef _LAPACKINTERFACE_H_
#define _LAPACKINTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

// MKL static linking: mkl_c.lib mkl_p3.lib mkl_lapack.lib
// MKL dynamic linking: mkl_c_dll.lib 
// LAPACK: clapack.lib

#ifndef MKL_NT
#if defined(_WIN32) 
#define MKL_NT 1
#else
#define MKL_NT 0
#endif
#endif

#if MKL_NT
#define BLAS_DECL(fname) fname
#define LAPACK_DECL(fname) fname

#else

#define BLAS_DECL(fname) fname##_
#define LAPACK_DECL(fname) fname##_

#endif
#define  BY_VAL(type) type *


/*
*
*   BLAS functions
*
*/
    void BLAS_DECL(dcopy) (int *, double *, int *, double *, int *);
    void BLAS_DECL(dtrmv) (char *, char *, char *, int *, double *, int *, double *, int *);

    double BLAS_DECL(ddot)(int *, double *, int *, double *, int *);
    
void BLAS_DECL(dscal) (BY_VAL(int) N, BY_VAL(double) alpha, 
    double *X, BY_VAL(int) incX);

void BLAS_DECL(daxpy) (BY_VAL(int) N, BY_VAL(double) alpha, 
    double *X, BY_VAL(int) incX, double *y, BY_VAL(int) incY);

void BLAS_DECL(dgemm) (
    BY_VAL(char) transa, BY_VAL(char) transb, 
    BY_VAL(int) m, BY_VAL(int) n, BY_VAL(int) k, BY_VAL(double) alpha, 
    double *a, BY_VAL(int) lda, double *b, BY_VAL(int) ldb, 
    BY_VAL(double) beta, double *c, BY_VAL(int) ldc);

void BLAS_DECL(dsymm) (
    BY_VAL(char) side, BY_VAL(char) uplo,
    BY_VAL(int) m, BY_VAL(int) n, BY_VAL(double) alpha, 
    double *a, BY_VAL(int) lda, double *b, BY_VAL(int) ldb, 
    BY_VAL(double) beta, double *c, BY_VAL(int) ldc);

    void BLAS_DECL(dtrmm) (
        BY_VAL(char) side, BY_VAL(char) uplo, BY_VAL(char) transa, BY_VAL(char) diag,
        BY_VAL(int) m, BY_VAL(int) n, BY_VAL(double) alpha, 
        double *a, BY_VAL(int) lda, double *b, BY_VAL(int) ldb );


/*
*
*   LAPACK functions
*
*/

    void LAPACK_DECL(dgesv)(int *,int *,double *, int *,int *, double *, int *,int *);
    
void LAPACK_DECL(dgetrs) (
	BY_VAL(char) trans, BY_VAL(int) n, BY_VAL(int) nrhs, double *a, BY_VAL(int) lda, int *ipiv, double *b, 
	BY_VAL(int) ldb, int *info);

void LAPACK_DECL(dgetrf) (
    BY_VAL(int) m, BY_VAL(int) n, double *a, BY_VAL(int) lda, int *ipiv, int *info);

void LAPACK_DECL(dgetri) (
    BY_VAL(int) n, double *a, BY_VAL(int) lda, int *ipiv, double *work, BY_VAL(int) lwork, int *info);

void LAPACK_DECL(dtrtri) (
    BY_VAL(char) uplo, BY_VAL(char) diag, BY_VAL(int) n, double *a, BY_VAL(int) lda, int *info);

void LAPACK_DECL(dsytrf) (
    BY_VAL(char) uplo, BY_VAL(int) n, double *a, BY_VAL(int) lda, int *ipiv, double *work, BY_VAL(int) lwork, int *info);

void LAPACK_DECL(dsytri) (
    BY_VAL(char) uplo, BY_VAL(int) n, double *a, BY_VAL(int) lda, int *ipiv, int *info);

void LAPACK_DECL(dgees) (
    BY_VAL(char) jobvs, BY_VAL(char) sort, void *select,
    BY_VAL(int) n, double *a, BY_VAL(int) lda, 
    int *sdim, double *wr, double *wi, double *vs,
    BY_VAL(int) ldvs, double *work, BY_VAL(int) lwork,
    void *bwork, int *info);

void LAPACK_DECL(dgeev)(
    BY_VAL(char) JOBVL, BY_VAL(char) JOBVR, BY_VAL(int) N, 
    double *A, BY_VAL(int) LDA, double *WR, double *WI, 
    double *VL, BY_VAL(int) LDVL, double *VR, BY_VAL(int) LDVR, 
    double *WORK, BY_VAL(int) LWORK, int *INFO );

void LAPACK_DECL(dsyevd)(
    BY_VAL(char) jobz, BY_VAL(char) uplo, BY_VAL(int) N, 
    double *A, BY_VAL(int) LDA, double *W,
    double *WORK, BY_VAL(int) lwork, 
    int *IWORKBY, BY_VAL(int) liwork, 
    int *INFO );

void LAPACK_DECL(dpotrf)(
    BY_VAL(char) uplo, BY_VAL(int) n, double *a,
    BY_VAL(int) lda, int *info);

#ifdef __cplusplus
}
#endif

#endif
