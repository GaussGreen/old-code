#ifndef NAGF08
#define NAGF08
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf08.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f08 Chapter
 *
 * Mark 5, 1997.
 */

#ifdef NAG_PROTO
extern void f08aef(Integer m, Integer n, double a[], Integer lda,
	    double tau[], double work[], Integer lwork,
	    Integer *info);
#else
extern void f08aef();
#endif

#ifdef NAG_PROTO
extern void f08aev(Integer n,  double *alpha, double x[],  Integer incx,
double *tau);
#else
extern void f08aev();
#endif

#ifdef NAG_PROTO
extern void f08aew(OperationSide Side,  Integer m, Integer n,  double v[],
Integer incv,  double tau, double c[],  Integer tdc,
double work[]);
#else
extern void f08aew();
#endif

#ifdef NAG_PROTO
extern void f08aex(char *direct, char *storev, Integer n, Integer k,
double v[], Integer ldv, double tau[],
double t[], Integer ldt);
#else
extern void f08aex();
#endif

#ifdef NAG_PROTO
extern void f08aey(char *side, char *trans, char *direct, char *storev, Integer m,
	    Integer n, Integer k, double v[], Integer ldv,
	    double t[], Integer ldt, double c[], Integer ldc,
	    double work[], Integer ldwork);
#else
extern void f08aey();
#endif

#ifdef NAG_PROTO
extern void f08aez(Integer m, Integer n, double a[], Integer lda,
double tau[], double work[], Integer *info);
#else
extern void f08aez();
#endif

#ifdef NAG_PROTO
extern void f08aff(Integer m, Integer n, Integer k, double a[],
	    Integer lda, double tau[], double work[],
	    Integer lwork, Integer *info);
#else
extern void f08aff();
#endif

#ifdef NAG_PROTO
extern void f08afz(Integer m, Integer n, Integer k, double a[],
Integer lda, double tau[], double work[],
Integer *info);
#else
extern void f08afz();
#endif

#ifdef NAG_PROTO
extern void f08agf(char *side, char *trans, Integer m, Integer n, Integer k,
double a[], Integer lda, double tau[],
double c[], Integer ldc, double work[],
Integer lwork, Integer *info);
#else
extern void f08agf();
#endif

#ifdef NAG_PROTO
extern void f08agz(char *side, char *trans, Integer m, Integer n, Integer k,
double a[], Integer lda, double tau[],
double c[], Integer ldc, double work[],
Integer *info);
#else
extern void f08agz();
#endif

#ifdef NAG_PROTO
extern void f08asf(Integer m, Integer n, Complex a[], Integer lda,
Complex tau[], Complex work[], Integer lwork,
	     Integer *info);
#else
extern void f08asf();
#endif

#ifdef NAG_PROTO
extern double f08asu(double x, double y, double z);
#else
extern double f08asu();
#endif

#ifdef NAG_PROTO
extern void f08asv(Integer n, Complex *alpha, Complex x[],
Integer incx, Complex *tau);
#else
extern void f08asv();
#endif

#ifdef NAG_PROTO
extern void f08asw(char *side, Integer m, Integer n, Complex v[],
Integer incv, Complex *tau, Complex c[],
Integer ldc, Complex work[]);
#else
extern void f08asw();
#endif

#ifdef NAG_PROTO
extern void f08asx(char *direct, char *storev, Integer n, Integer k,
Complex v[], Integer ldv, Complex tau[],
Complex t[], Integer ldt);
#else
extern void f08asx();
#endif

#ifdef NAG_PROTO
extern void f08asy(char *side, char *trans, char *direct, char *storev,
	     Integer m, Integer n, Integer k, Complex v[],
	     Integer ldv, Complex t[], Integer ldt,
	     Complex c[], Integer ldc, Complex work[],
	     Integer ldwork);
#else
extern void f08asy();
#endif

#ifdef NAG_PROTO
extern void f08asz(Integer m, Integer n, Complex a[], Integer lda,
	     Complex tau[], Complex work[], Integer *info);
#else
extern void f08asz();
#endif

#ifdef NAG_PROTO
extern void f08atf(Integer m, Integer n, Integer k, Complex a[],
	     Integer lda, Complex tau[], Complex work[],
	     Integer lwork, Integer *info);
#else
extern void f08atf();
#endif

#ifdef NAG_PROTO
extern void f08atz(Integer m, Integer n, Integer k, Complex a[],
	    Integer lda, Complex tau[], Complex work[],
	    Integer *info);
#else
extern void f08atz();
#endif

#ifdef NAG_PROTO
extern void f08auf(char *side, char *trans, Integer m, Integer n, Integer k,
	    Complex a[], Integer lda, Complex tau[],
	    Complex c[], Integer ldc, Complex work[],
	    Integer lwork, Integer *info);
#else
extern void f08auf();
#endif

#ifdef NAG_PROTO
extern void f08auz(char *side, char *trans, Integer m, Integer n, Integer k,
	     Complex a[], Integer tda, Complex tau[],
	     Complex c[], Integer tdc, Complex work[],
	     Integer *info);
#else
extern void f08auz();
#endif

#ifdef NAG_PROTO
extern void f08cfy(Integer m, Integer n, Integer k, double a[],
	     Integer lda, double tau[], double work[],
	     Integer lwork, Integer *info);
#else
extern void f08cfy();
#endif

#ifdef NAG_PROTO
extern void f08cfz(Integer m, Integer n, Integer k, double a[],
	     Integer lda, double tau[], double work[],
	     Integer *info);
#else
extern void f08cfz();
#endif

#ifdef NAG_PROTO
extern void f08fef(char *uplo, Integer n, double a[], Integer lda,
double d[], double e[], double tau[],
double work[], Integer lwork, Integer *info);
#else
extern void f08fef();
#endif

#ifdef NAG_PROTO
extern void f08fey(char *uplo, Integer n, Integer nb, double a[],
	     Integer lda, double e[], double tau[], double w[],
	     Integer ldw);
#else
extern void f08fey();
#endif

#ifdef NAG_PROTO
extern void f08fez(char *uplo, Integer n, double a[], Integer lda,
	     double d[], double e[], double tau[],
	     Integer *info);
#else
extern void f08fez();
#endif

#ifdef NAG_PROTO
extern void f08fff(char *uplo, Integer n, double a[], Integer lda,
	    double tau[], double work[], Integer lwork,
	    Integer *info);
#else
extern void f08fff();
#endif

#ifdef NAG_PROTO
extern void f08gef(MatrixTriangle UpperLower,  Integer n,  double ap[], double d[],
double e[], double tau[],  Integer *info);
#else
extern void f08gef();
#endif

#ifdef NAG_PROTO
extern void f08ggf(OperationSide Side, MatrixTriangle Uplo, MatrixTranspose Trans,  Integer m, Integer n,
double ap[], double tau[], double c[],  Integer tdc,
double work[],  Integer *info);
#else
extern void f08ggf();
#endif

#ifdef NAG_PROTO
extern void f08hew(double f, double g, double *cs, double *sn,
	    double *r);
#else
extern void f08hew();
#endif

#ifdef NAG_PROTO
extern void f08jef(char *compz, Integer n, double d[], double e[],
	    double z[], Integer ldz, double work[],
	    Integer *info);
#else
extern void f08jef();
#endif

#ifdef NAG_PROTO
extern void f08jex(double a, double b, double c, double *rt1,
double *rt2);
#else
extern void f08jex();
#endif

#ifdef NAG_PROTO
extern void f08jff(Integer n,  double d[], double e[],  Integer *info);
#else
extern void f08jff();
#endif

#ifdef NAG_PROTO
extern void f08jjf(char *range, char *order,  Integer n,  double vl,
double vu,  Integer il, Integer iu,  double abstol,
double d[], double e[],  Integer *m, Integer *nsplit,
double w[],  Integer iblock[], Integer isplit[],
double work[],  Integer iwork[], Integer *info);
#else
extern void f08jjf();
#endif

#ifdef NAG_PROTO
extern void f08jjz(Integer ijob, Integer nitmax, Integer n, Integer mmax,
		   Integer minp, Integer nbmin,  double abstol,
		   double reltol, double pivmin, double d[],
		   double e[], double e2[],  Integer nval[],
		   double ab[], double c[],  Integer *mout, Integer nab[],
		   double work[],  Integer iwork[], Integer *info);
#else
extern void f08jjz();
#endif

#ifdef NAG_PROTO
extern void f08jkf(Integer n,  double d[], double e[],  Integer m,
double w[],  Integer iblock[], Integer isplit[],
double z[],  Integer tdz,  double work[],
Integer iwork[], Integer ifail[], Integer *info);
#else
extern void f08jkf();
#endif

#ifdef NAG_PROTO
extern void f08nef(Integer n, Integer ilo, Integer ihi, double a[],
Integer lda, double tau[], double work[],
Integer lwork, Integer *info);
#else
extern void f08nef();
#endif

#ifdef NAG_PROTO
extern void f08ney(Integer n, Integer k, Integer nb, double a[],
	     Integer lda, double tau[], double t[], Integer ldt,
	     double y[], Integer ldy);
#else
extern void f08ney();
#endif

#ifdef NAG_PROTO
extern void f08nez(Integer n, Integer ilo, Integer ihi, double a[],
	    Integer lda, double tau[], double work[],
	    Integer *info);
#else
extern void f08nez();
#endif

#ifdef NAG_PROTO
extern void f08ngf(char *side, char *trans, Integer m, Integer n, Integer ilo,
	    Integer ihi, double a[], Integer lda, double tau[],
	    double c[], Integer ldc, double work[],
	    Integer lwork, Integer *info);
#else
extern void f08ngf();
#endif

#ifdef NAG_PROTO
extern void f08nhf(char *job, Integer n, double a[], Integer lda,
Integer *ilo, Integer *ihi, double scale[], Integer *info);
#else
extern void f08nhf();
#endif

#ifdef NAG_PROTO
extern void f08njf(char *job, char *side, Integer n, Integer ilo, Integer ihi,
double scale[], Integer m, double v[], Integer ldv,
Integer *info);
#else
extern void f08njf();
#endif

#ifdef NAG_PROTO
extern void f08nsf(Integer n, Integer ilo, Integer ihi, Complex a[],
	    Integer lda, Complex tau[], Complex work[],
	    Integer lwork, Integer *info);
#else
extern void f08nsf();
#endif

#ifdef NAG_PROTO
extern void f08nsy(Integer n, Integer k, Integer nb, Complex a[],
	    Integer lda, Complex tau[], Complex t[],
	    Integer ldt, Complex y[], Integer ldy);
#else
extern void f08nsy();
#endif

#ifdef NAG_PROTO
extern void f08nsz(Integer n, Integer ilo, Integer ihi, Complex a[],
	    Integer lda, Complex tau[], Complex work[],
	    Integer *info);
#else
extern void f08nsz();
#endif

#ifdef NAG_PROTO
extern void f08nuf(char *side, char *trans, Integer m, Integer n, Integer ilo,
	    Integer ihi, Complex a[], Integer lda,
	    Complex tau[], Complex c[], Integer ldc,
	    Complex work[], Integer lwork, Integer *info);
#else
extern void f08nuf();
#endif

#ifdef NAG_PROTO
extern void f08nvf(char *job, Integer n, Complex a[], Integer lda,
Integer *ilo, Integer *ihi, double scale[], Integer *info);
#else
extern void f08nvf();
#endif

#ifdef NAG_PROTO
extern void f08nwf(char *job, char *side, Integer n, Integer ilo, Integer ihi,
double scale[], Integer m, Complex v[],
Integer ldv, Integer *info);
#else
extern void f08nwf();
#endif

#ifdef NAG_PROTO
extern void f08pef(char *job, char *compz, Integer n, Integer ilo, Integer ihi,
double h[], Integer ldh, double wr[], double wi[],
double z[], Integer ldz, double work[],
Integer lwork, Integer *info);
#else
extern void f08pef();
#endif

#ifdef NAG_PROTO
extern void f08pew(char *side, Integer m, Integer n, double v[],
	    double tau, double c[], Integer ldc,
	    double work[]);
#else
extern void f08pew();
#endif

#ifdef NAG_PROTO
extern void f08pey(double *a, double *b, double *c, double *d,
	     double *rt1r, double *rt1i, double *rt2r,
	     double *rt2i, double *cs, double *sn);
#else
extern void f08pey();
#endif

#ifdef NAG_PROTO
extern void f08pez(Boolean wantt, Boolean wantz, Integer n, Integer ilo,
	     Integer ihi, double h[], Integer ldh, double wr[],
	     double wi[], Integer iloz, Integer ihiz, double z[],
	     Integer ldz, Integer *info);
#else
extern void f08pez();
#endif

#ifdef NAG_PROTO
extern void f08pkf(char *job, char *eigsrc, char *initv, Boolean select[],
	     Integer n, double h[], Integer ldh, double wr[],
	     double wi[], double vl[], Integer ldvl,
	     double vr[], Integer ldvr, Integer mm, Integer *m,
	     double work[], Integer ifaill[], Integer ifailr[],
	     Integer *info);
#else
extern void f08pkf();
#endif

#ifdef NAG_PROTO
extern void f08pkz(Boolean rightv, Boolean noinit, Integer n, double h[],
	     Integer ldh, double wr, double wi, double vr[],
	     double vi[], double b[], Integer ldb,
	     double work[], double eps3, double smlnum,
	     double bignum, Integer *info);
#else
extern void f08pkz();
#endif

#ifdef NAG_PROTO
extern void f08psf(char *job, char *compz, Integer n, Integer ilo,
Integer ihi, Complex h[], Integer ldh,
Complex w[], Complex z[], Integer ldz,
Complex work[], Integer lwork, Integer *info);
#else
extern void f08psf();
#endif

#ifdef NAG_PROTO
extern void f08psw(char *side, Integer m, Integer n, Complex v[],
Complex *tau, Complex c[], Integer ldc,
Complex work[]);
#else
extern void f08psw();
#endif

#ifdef NAG_PROTO
extern void f08psz(Boolean wantt, Boolean wantz, Integer n, Integer ilo,
	     Integer ihi, Complex h[], Integer ldh,
	     Complex w[], Integer iloz, Integer ihiz,
	     Complex z[], Integer ldz, Integer *info);
#else
extern void f08psz();
#endif

#ifdef NAG_PROTO
extern void f08pxf(char *job, char *eigsrc, char *initv, Boolean select[],
Integer n, Complex h[], Integer ldh,
Complex w[], Complex vl[], Integer ldvl,
Complex vr[], Integer ldvr, Integer mm, Integer *m,
Complex work[], double rwork[], Integer ifaill[],
Integer ifailr[], Integer *info);
#else
extern void f08pxf();
#endif

#ifdef NAG_PROTO
extern void f08pxz(Boolean rightv, Boolean noinit, Integer n,
	     Complex h[], Integer ldh, Complex *w,
	     Complex v[], Complex b[], Integer ldb,
	     double rwork[], double eps3, double smlnum,
	     Integer *info);
#else
extern void f08pxz();
#endif



#ifdef __cplusplus
}
#endif
#endif /* not NAGF08 */



