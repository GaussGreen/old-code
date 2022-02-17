#ifndef NAGF02
#define NAGF02
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf02.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f02 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2150 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02aac(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *r, 
            NagError NAG_HUGE *fail);
#else
extern void f02aac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02abc(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *r, double NAG_HUGE *v, 
            Integer tdv, NagError NAG_HUGE *fail);
#else
extern void f02abc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02adc(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, Integer tdb, 
            double NAG_HUGE *r, NagError NAG_HUGE *fail);
#else
extern void f02adc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02aec(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, Integer tdb,
            double NAG_HUGE *r, double NAG_HUGE *v, Integer tdv, NagError NAG_HUGE *fail);
#else
extern void f02aec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02afc(Integer n, double NAG_HUGE *a, Integer tda, Complex NAG_HUGE *r,
            Integer NAG_HUGE *iter, NagError NAG_HUGE *fail);
#else
extern void f02afc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02agc(Integer n, double NAG_HUGE a[],  Integer tda, Complex NAG_HUGE r[],
            Complex NAG_HUGE v[],  Integer tdv, Integer NAG_HUGE iter[], NagError NAG_HUGE *fail);
#else
extern void f02agc();
#endif

#ifdef NAG_PROTO
extern void f02amc(Integer n, double eps, double NAG_HUGE *d, double NAG_HUGE *e, double NAG_HUGE *z, 
            Integer tdz, NagError NAG_HUGE *fail);
#else
extern void f02amc();
#endif

#ifdef NAG_PROTO
extern void f02apc(Integer nn, double acc, double NAG_HUGE *h, Integer tdh, 
            Complex NAG_HUGE *eigvals, Integer NAG_HUGE *icnt, NagError NAG_HUGE *fail);
#else
extern void f02apc();
#endif

#ifdef NAG_PROTO
extern void f02aqc(Integer n, Integer low, Integer upp, 
            double machep, double NAG_HUGE *h, Integer tdh, double NAG_HUGE *vecs, 
            Integer tdvecs, Complex NAG_HUGE *eigvals, Integer NAG_HUGE *cnt, NagError NAG_HUGE *fail);
#else
extern void f02aqc();
#endif

#ifdef NAG_PROTO
extern void f02avc(Integer n, double acheps, double NAG_HUGE *d, double NAG_HUGE *e, NagError NAG_HUGE *fail);
#else
extern void f02avc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02awc(Integer n, Complex NAG_HUGE a[], Integer tda, double NAG_HUGE r[], NagError NAG_HUGE *fail);
#else
extern void f02awc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02axc(Integer n, Complex NAG_HUGE a[], Integer tda, double NAG_HUGE r[],
            Complex NAG_HUGE v[], Integer tdv, NagError NAG_HUGE *fail);
#else
extern void f02axc();
#endif

#ifdef NAG_PROTO
extern void f02ayc(Integer n, double eps, double NAG_HUGE d[], double NAG_HUGE e[],
            Complex NAG_HUGE v[], Integer tdv, NagError NAG_HUGE *fail);
#else
extern void f02ayc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02bjc(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, double tol, Complex NAG_HUGE alfa[],
            double NAG_HUGE beta[], Boolean wantv,
            double NAG_HUGE v[], Integer tdv, Integer NAG_HUGE iter[], NagError NAG_HUGE *fail);
#else
extern void f02bjc();
#endif

#ifdef NAG_PROTO
extern void f02bjw(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, Boolean matz, double NAG_HUGE z[], Integer tdz);
#else
extern void f02bjw();
#endif

#ifdef NAG_PROTO
extern void f02bjx(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, double eps1, Boolean matz,
            double NAG_HUGE z[], Integer tdz, Integer NAG_HUGE iter[], Integer NAG_HUGE *info);
#else
extern void f02bjx();
#endif

#ifdef NAG_PROTO
extern void f02bjy(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, Complex NAG_HUGE alfa[], double NAG_HUGE beta[], Boolean matz,
            double NAG_HUGE z[], Integer tdz);
#else
extern void f02bjy();
#endif

#ifdef NAG_PROTO
extern void f02bjz(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, Complex NAG_HUGE alfa[], double NAG_HUGE beta[],
            double NAG_HUGE z[], Integer tdz);
#else
extern void f02bjz();
#endif

#ifdef NAG_PROTO
extern void f02eaz(double amax, double rmin, double rmax,
double *sigma, Boolean *scale);
#else
extern void f02eaz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02ecc(Nag_Select_Eigenvalues crit, Integer n, double a[], Integer lda,
             double wl, double wu, Integer mest, Integer *m,
             Complex w[], Complex v[], Integer tdv, NagError *fail);
#else
extern void f02ecc();
#endif


#ifdef NAG_PROTO
extern void f02ecf(char *crit, Integer n, double a[], Integer lda,
	     double wl, double wu, Integer mest, Integer *m,
	     double wr[], double wi[], double vr[],
	     Integer ldvr, double vi[], Integer ldvi,
	     double work[], Integer lwork, Integer iwork[],
	     Boolean bwork[], Integer *ifail);
#else
extern void f02ecf();
#endif

#ifdef NAG_PROTO
extern void f02faf(char *job, char *uplo, Integer n, double a[],
	    Integer lda, double w[], double work[], Integer lwork,
	    Integer *ifail);
#else
extern void f02faf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02gcc(Nag_Select_Eigenvalues crit, Integer n, Complex a[], Integer lda,
	     double wl, double wu, Integer mest, Integer *m,
	     Complex w[], Complex v[], Integer ldv, NagError *fail);
#else
extern void f02gcc();
#endif



#ifdef NAG_PROTO
extern void f02gcf(char *crit, Integer n, Complex a[], Integer lda,
	     double wl, double wu, Integer mest, Integer *m,
	     Complex w[], Complex v[], Integer ldv,
	     Complex work[], Integer lwork, double rwork[],
	     Integer iwork[], Boolean bwork[], Integer *ifail);
#else
extern void f02gcf();
#endif

#ifdef NAG_PROTO
extern void f02swc(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *d, double NAG_HUGE *e, Integer ncoly,
            double NAG_HUGE *y, Integer tdy, Boolean wantq, double NAG_HUGE *q, Integer tdq, NagError NAG_HUGE *fail);
#else
extern void f02swc();
#endif

#ifdef NAG_PROTO
extern void f02sxc(Integer n, double NAG_HUGE *a, Integer tda, Integer ncoly, double NAG_HUGE *y,
            Integer tdy, NagError NAG_HUGE *fail);
#else
extern void f02sxc();
#endif

#ifdef NAG_PROTO
extern void f02syc(Integer n, double NAG_HUGE *d, double NAG_HUGE *e, Integer ncolb, double NAG_HUGE *b, Integer tdb,
            Integer nrowy, double NAG_HUGE *y, Integer tdy, Integer ncolz, double NAG_HUGE *z, Integer tdz,
            Integer NAG_HUGE *iter, Integer NAG_HUGE *failinfo, NagError NAG_HUGE *fail);
#else
extern void f02syc();
#endif

#ifdef NAG_PROTO
extern void f02szc(Integer n,  double NAG_HUGE d[], double NAG_HUGE e[], double NAG_HUGE sv[],
            Boolean wantb,  double NAG_HUGE b[],  Boolean wanty,
            double NAG_HUGE y[],  Integer tdy, Integer ly,  Boolean wantz,
            double NAG_HUGE z[],  Integer tdz, Integer ncz,  double NAG_HUGE work1[],
            NagError NAG_HUGE *fail);
#else
extern void f02szc();
#endif

#ifdef NAG_PROTO
extern void f02szz(Integer n, double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE x[]);
#else
extern void f02szz();
#endif

#ifdef NAG_PROTO
extern void f02uwc(Integer n, Complex NAG_HUGE *a, Integer tda, double NAG_HUGE *d, double NAG_HUGE *e,
            Integer ncoly, Complex NAG_HUGE *y, Integer tdy, Boolean wantq, Complex NAG_HUGE *q,
            Integer tdq, NagError NAG_HUGE *fail);
#else
extern void f02uwc();
#endif

#ifdef NAG_PROTO
extern void f02uxc(Integer n, Complex NAG_HUGE *a, Integer tda, Integer ncoly, Complex NAG_HUGE *y,
            Integer tdy, NagError NAG_HUGE *fail);
#else
extern void f02uxc();
#endif

#ifdef NAG_PROTO
extern void f02uyc(Integer n, double NAG_HUGE *d, double NAG_HUGE *e, Integer ncolb, Complex NAG_HUGE *b, Integer tdb,
            Integer nrowy, Complex NAG_HUGE *y, Integer tdy, Integer ncolz, Complex NAG_HUGE *z, Integer tdz,
            Integer NAG_HUGE *iter, Integer NAG_HUGE *failinfo, NagError NAG_HUGE *fail);
#else
extern void f02uyc();
#endif

#ifdef NAG_PROTO
extern void f02wax(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             Boolean wantb,  double NAG_HUGE b[], double NAG_HUGE sv[],
             double NAG_HUGE work[],  Integer lwork, Integer NAG_HUGE *ifail);
#else
extern void f02wax();
#endif

#ifdef NAG_PROTO
extern void f02way(Integer n,  double NAG_HUGE c[],  Integer tdc,  double NAG_HUGE pt[],
            Integer tdpt, Integer NAG_HUGE *fail);
#else
extern void f02way();
#endif

#ifdef NAG_PROTO
extern void f02waz(Integer m, Integer n, double NAG_HUGE *a, Integer tda,
            double NAG_HUGE *z, double NAG_HUGE *b, double NAG_HUGE *c, int NAG_HUGE *localerror);
#else
extern void f02waz();
#endif

#ifdef NAG_PROTO
extern void f02wbc(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             Boolean wantb,  double NAG_HUGE b[], double NAG_HUGE sv[],
             double NAG_HUGE work[],  Integer lwork, Integer NAG_HUGE *ifail);
#else
extern void f02wbc();
#endif

#ifdef NAG_PROTO
extern void f02wby(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE work[], Integer NAG_HUGE *fail);
#else
extern void f02wby();
#endif

#ifdef NAG_PROTO
extern void f02wbz(Integer m, Integer n,  double NAG_HUGE c[],  Integer tdc,
             double NAG_HUGE pt[],  Integer tdpt,  double NAG_HUGE work[],  Integer NAG_HUGE *fail);
#else
extern void f02wbz();
#endif

#ifdef NAG_PROTO
extern void f02wcy(Integer n, double NAG_HUGE c[], Integer tdc, double NAG_HUGE q[],
            Integer tdq,  double NAG_HUGE work1[], double NAG_HUGE work2[], Integer NAG_HUGE *fail);
#else
extern void f02wcy();
#endif

#ifdef NAG_PROTO
extern void f02wdc(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
            Boolean wantb,  double NAG_HUGE b[], double tol,  Boolean NAG_HUGE *svd,
            Integer NAG_HUGE *irank,  double NAG_HUGE z[], double NAG_HUGE sv[],
            Boolean wantr,  double NAG_HUGE r[],  Integer tdr,  Boolean wantpt,
            double NAG_HUGE pt[],  Integer tdpt,  double NAG_HUGE work[],
            Integer lwork, NagError NAG_HUGE *fail);
#else
extern void f02wdc();
#endif

#ifdef NAG_PROTO
extern Integer f02wdy(Integer n, double NAG_HUGE sv[], double tol);
#else
extern Integer f02wdy();
#endif

#ifdef NAG_PROTO
extern double f02wdz(Integer n,double NAG_HUGE *a,Integer tda,double NAG_HUGE *work, int NAG_HUGE *localerror);
#else
extern double f02wdz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02wec(Integer m, Integer n, double NAG_HUGE *a, Integer tda, Integer ncolb,
            double NAG_HUGE *b, Integer tdb, Boolean wantq, double NAG_HUGE *q, Integer tdq, double NAG_HUGE *sv,
            Boolean wantp, double NAG_HUGE *pt, Integer tdpt, Integer NAG_HUGE *iter, double NAG_HUGE *e,
            Integer NAG_HUGE *failinfo, NagError NAG_HUGE *fail);
#else
extern void f02wec();
#endif

#ifdef NAG_PROTO
extern void f02wex(Integer nx, Integer ny, double NAG_HUGE *alpha, double NAG_HUGE *x, Integer incx,
            double NAG_HUGE *y, Integer incy, double tol, double NAG_HUGE *zeta);
#else
extern void f02wex();
#endif

#ifdef NAG_PROTO
extern void f02wuc(Integer n, double NAG_HUGE *a, Integer tda, Integer ncolb, double NAG_HUGE *b,
            Integer tdb, Boolean wantq, double NAG_HUGE *q, Integer tdq, double NAG_HUGE *sv, Boolean wantp,
            Integer NAG_HUGE *iter, double NAG_HUGE *e, Integer NAG_HUGE *failinfo, NagError NAG_HUGE *fail);
#else
extern void f02wuc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f02xec(Integer m, Integer n, Complex NAG_HUGE *a, Integer tda, Integer ncolb,
            Complex NAG_HUGE *b, Integer tdb, Boolean wantq, Complex NAG_HUGE *q, Integer tdq,
            double NAG_HUGE *sv, Boolean wantp, Complex NAG_HUGE *ph, Integer tdph, Integer NAG_HUGE *iter,
            double NAG_HUGE *e, Integer NAG_HUGE *failinfo, NagError NAG_HUGE *fail);
#else
extern void f02xec();
#endif

#ifdef NAG_PROTO
extern void f02xex(Integer nx, Integer ny, Complex NAG_HUGE *alpha, Complex NAG_HUGE *x, Integer incx,
            Complex NAG_HUGE *y, Integer incy, double tol, Complex NAG_HUGE *theta);
#else
extern void f02xex();
#endif

#ifdef NAG_PROTO
extern void f02xuc(Integer n, Complex NAG_HUGE *a, Integer tda, Integer ncolb,
            Complex NAG_HUGE *b, Integer tdb, Boolean wantq, Complex NAG_HUGE *q, Integer tdq,
            double NAG_HUGE *sv, Boolean wantp, Integer NAG_HUGE *iter, double NAG_HUGE *e, Integer NAG_HUGE *failinfo,
            NagError NAG_HUGE *fail);
#else
extern void f02xuc();
#endif

#ifdef NAG_PROTO
extern void f02xus(Integer n, double NAG_HUGE *d, double NAG_HUGE *e, Boolean wantcs, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f02xus();
#endif

#ifdef NAG_PROTO
extern void f02xut(double test, Integer n, double NAG_HUGE *d, double NAG_HUGE *e, Boolean NAG_HUGE *force, Integer NAG_HUGE *p);
#else
extern void f02xut();
#endif

#ifdef NAG_PROTO
extern void f02xuu(Integer job, double d1, double e1, double dnm1, double dn,
            double enm2, double enm1, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f02xuu();
#endif

#ifdef NAG_PROTO
extern void f02xuv(Nag_InitRotation shift, Integer m, Integer n, double NAG_HUGE *d, double NAG_HUGE *e,
            double c, double s, Boolean wantlt, double NAG_HUGE *cl, double NAG_HUGE *sl, Boolean wantrt,
            double NAG_HUGE *cr, double NAG_HUGE *sr);
#else
extern void f02xuv();
#endif

#ifdef NAG_PROTO
extern void f02xuw(Integer m, Integer n, double NAG_HUGE *d, double NAG_HUGE *e, Boolean wantcs,
            double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f02xuw();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGF02 */
