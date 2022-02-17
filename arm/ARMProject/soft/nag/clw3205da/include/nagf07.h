#ifndef NAGF07
#define NAGF07
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf07.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f07 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2154 (Feb 1998).
 */

#ifdef NAG_PROTO
extern void f07ad0(Integer m, Integer n, double NAG_HUGE *a, Integer lda, double NAG_HUGE piv[], 
            Integer NAG_HUGE *info);
#else
extern void f07ad0();
#endif

#ifdef NAG_PROTO
extern void f07ad1(Integer m, Integer n, double NAG_HUGE *a, Integer lda, 
            double NAG_HUGE piv[], Integer NAG_HUGE *info);
#else
extern void f07ad1();
#endif

#ifdef NAG_PROTO
extern void f07ad2(Integer n, double NAG_HUGE *a, Integer lda, Integer k1, Integer k2, 
            double NAG_HUGE piv[], Integer incx);
#else
extern void f07ad2();
#endif

#ifdef NAG_PROTO
extern void f07adg(Integer m, Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *piv,
            Integer NAG_HUGE *info);
#else
extern void f07adg();
#endif

#ifdef NAG_PROTO
extern void f07adh(Integer m, Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *piv,
            Integer NAG_HUGE *info);
#else
extern void f07adh();
#endif

#ifdef NAG_PROTO
extern void f07adj(Integer n, double NAG_HUGE *a, Integer tda, Integer k1, Integer k2,
            double NAG_HUGE *piv, Integer incx);
#else
extern void f07adj();
#endif

#ifdef NAG_PROTO
extern void f07ae0(char NAG_HUGE *trans, Integer n, Integer nrhs, double NAG_HUGE *a, Integer lda,
            double NAG_HUGE piv[], double NAG_HUGE *b, Integer ldb, Integer NAG_HUGE *info);
#else
extern void f07ae0();
#endif

#ifdef NAG_PROTO
extern void f07aeg(MatrixTranspose trans, Integer n, Integer nrhs, double NAG_HUGE *a,
            Integer tda, double NAG_HUGE *piv, double NAG_HUGE *b, Integer tdb, Integer NAG_HUGE *info);
#else
extern void f07aeg();
#endif

#ifdef NAG_PROTO
extern void f07fdc(MatrixTriangle uplo, Integer n, double NAG_HUGE a[], Integer tda,
            Integer NAG_HUGE *info);
#else
extern void f07fdc();
#endif

#ifdef NAG_PROTO
extern void f07fdz(MatrixTriangle uplo, Integer n, double NAG_HUGE a[], Integer tda,
            Integer NAG_HUGE *info);
#else
extern void f07fdz();
#endif

#ifdef NAG_PROTO
extern void f07fry(Integer n, Complex *x, Integer incx);
#else
extern void f07fry();
#endif

#ifdef NAG_PROTO
extern void f07mdx(double a, double b, double c, double *rt1,
double *rt2, double *cs1, double *sn1);
#else
extern void f07mdx();
#endif

#ifdef NAG_PROTO
extern void f07tgz(char *uplo, char *trans, char *diag, char *normin, Integer n,
double a[], Integer lda, double x[],
double *scale, double cnorm[], Integer *info);
#else
extern void f07tgz();
#endif

#ifdef NAG_PROTO
extern void f07tuz(char *uplo, char *trans, char *diag, char *normin, Integer n,
Complex a[], Integer lda, Complex x[],
double *scale, double cnorm[], Integer *info);
#else
extern void f07tuz();
#endif

#ifdef NAG_PROTO
extern Integer f07za0(char NAG_HUGE *name);
#else
extern Integer f07za0();
#endif

#ifdef NAG_PROTO
extern void f07za1(Integer ispec, char NAG_HUGE *name, Integer NAG_HUGE *ival, Integer rwflag);
#else
extern void f07za1();
#endif

#ifdef NAG_PROTO
extern Integer f07zay(char NAG_HUGE *name);
#else
extern Integer f07zay();
#endif

#ifdef NAG_PROTO
extern void f07zaz(Integer ispec, char NAG_HUGE *name, Integer NAG_HUGE *ival, Integer rwflag);
#else
extern void f07zaz();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGF07 */
