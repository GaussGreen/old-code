#ifndef NAGF04
#define NAGF04
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf04.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f04 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2152 (Feb 1998).
 */

#ifdef NAG_PROTO
extern void f04aac(double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, Integer tdb, Integer n,
            Integer m, double NAG_HUGE *c, Integer tdc, NagError NAG_HUGE *fail);
#else
extern void f04aac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04adc(Integer n, Integer nrhs, Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *b,
            Integer tdb, Complex NAG_HUGE *x, Integer tdx, NagError NAG_HUGE *fail);
#else
extern void f04adc();
#endif

#ifdef NAG_PROTO
extern void f04aec(Integer n, Integer m, double NAG_HUGE *a, Integer tda, double NAG_HUGE *b,
            Integer tdb, double NAG_HUGE *c, Integer tdc, double NAG_HUGE *aa, Integer tdaa,
            double NAG_HUGE *bb, Integer tdbb, NagError NAG_HUGE *fail);
#else
extern void f04aec();
#endif

#ifdef NAG_PROTO
extern void f04afc(Integer n, Integer nrhs, double NAG_HUGE *a, Integer tda, double NAG_HUGE *p,
            double NAG_HUGE *b, Integer tdb,  double NAG_HUGE *x, Integer tdx, double NAG_HUGE *bb,
            Integer tdbb, Integer NAG_HUGE *k, NagError NAG_HUGE *fail);
#else
extern void f04afc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04agc(Integer n, Integer nrhs, double NAG_HUGE *a, Integer tda, double NAG_HUGE *p, 
            double NAG_HUGE *b, Integer tdb, double NAG_HUGE *x, Integer tdx, NagError NAG_HUGE *fail);
#else
extern void f04agc();
#endif

#ifdef NAG_PROTO
extern void f04ahc(Integer n, Integer nrhs, double NAG_HUGE *a, Integer tda,
            double NAG_HUGE *aa, Integer tdaa, Integer NAG_HUGE *pivot, double NAG_HUGE *b,
            Integer tdb, double NAG_HUGE *x, Integer tdx,
            double NAG_HUGE *bb, Integer tdbb, Integer NAG_HUGE *k, NagError NAG_HUGE *fail);
#else
extern void f04ahc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04ajc(Integer n, Integer nrhs, double NAG_HUGE *a, 
            Integer tda, Integer NAG_HUGE *pivot, double NAG_HUGE *b, Integer tdb, NagError NAG_HUGE *fail);
     
#else
extern void  f04ajc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04akc(Integer n, Integer nrhs, Complex NAG_HUGE *a, Integer tda, Integer NAG_HUGE *pivot,
            Complex NAG_HUGE *b, Integer tdb, NagError NAG_HUGE *fail);
#else
extern void f04akc();
#endif

#ifdef NAG_PROTO
extern void f04aqc(Integer n, double NAG_HUGE rl[], double NAG_HUGE d[],
            double NAG_HUGE b[], double NAG_HUGE x[]);
#else
extern void f04aqc();
#endif

#ifdef NAG_PROTO
extern void f04aqz(Integer n, Integer m,  double rl[], double d[],
double b[], double x[]);
#else
extern void f04aqz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04arc(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, 
            double NAG_HUGE *x, NagError NAG_HUGE *fail);
#else
extern void f04arc();
#endif

#ifdef NAG_PROTO
extern void f04asc(double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, Integer n, double NAG_HUGE *c,
            double NAG_HUGE *wk1, double NAG_HUGE *wk2, NagError NAG_HUGE *fail);
#else
extern void f04asc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04awc(Integer n, Integer nrhs, Complex NAG_HUGE *a, Integer tda, double NAG_HUGE *p,
            Complex NAG_HUGE *b, Integer tdb, Complex NAG_HUGE *x, Integer tdx, NagError NAG_HUGE *fail);
#else
extern void f04awc();
#endif

#ifdef NAG_PROTO
extern void f04jac(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             double NAG_HUGE b[], double tol, double NAG_HUGE *sigma,
             Integer NAG_HUGE *irank, NagError NAG_HUGE *fail);
#else
extern void f04jac();
#endif

#ifdef NAG_PROTO
extern void f04jay(Integer n, Integer rank, double NAG_HUGE sv[],
            double NAG_HUGE b[], double NAG_HUGE pt[], Integer tdpt, double NAG_HUGE x[],
            double NAG_HUGE work[]);
#else
extern void f04jay();
#endif

#ifdef NAG_PROTO
extern void f04jaz(Integer m, Integer n, Integer rank, double NAG_HUGE sv[],
            double NAG_HUGE b[], double NAG_HUGE pt[], Integer tdpt,
            double NAG_HUGE x[], double NAG_HUGE *sigma, double NAG_HUGE work[]);
#else
extern void f04jaz();
#endif

#ifdef NAG_PROTO
extern void f04jdc(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             double NAG_HUGE b[], double tol, double NAG_HUGE *sigma,
             Integer NAG_HUGE *irank, NagError NAG_HUGE *fail);
#else
extern void f04jdc();
#endif

#ifdef NAG_PROTO
extern void f04jgc(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
            double NAG_HUGE b[], double tol,  Boolean NAG_HUGE *svd,
            double NAG_HUGE *sigma,  Integer NAG_HUGE *irank,  double NAG_HUGE comm[],
            Integer lcomm, NagError NAG_HUGE *fail);
#else
extern void f04jgc();
#endif

#ifdef NAG_PROTO
extern double f04jgr(Integer n, double NAG_HUGE *a, Integer nra, double NAG_HUGE *work,
              int NAG_HUGE *localerror);
#else
extern double f04jgr();
#endif

#ifdef NAG_PROTO
extern double f04jgs(Integer n, double NAG_HUGE *a, Integer tda, int NAG_HUGE *localerror);
#else
extern double f04jgs();
#endif

#ifdef NAG_PROTO
extern void f04jgt(Integer n, double NAG_HUGE *x, double NAG_HUGE *scale, double NAG_HUGE *sumsq, double tiny,
            Boolean undflw);
#else
extern void f04jgt();
#endif

#ifdef NAG_PROTO
extern double f04jgu(double scale, double sumsq, double big);
#else
extern double f04jgu();
#endif

#ifdef NAG_PROTO
extern double f04jgv(Integer n, double NAG_HUGE *x, double tiny, double big);
#else
extern double f04jgv();
#endif

#ifdef NAG_PROTO
extern double f04jgw(Integer n, double NAG_HUGE x[]);
#else
extern double f04jgw();
#endif

#ifdef NAG_PROTO
extern Boolean f04jgx(double a, double b, double small1);
#else
extern Boolean f04jgx();
#endif

#ifdef NAG_PROTO
extern void f04jgy(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *b, double NAG_HUGE *x,
            Integer NAG_HUGE *failinfo, int NAG_HUGE *localerror);
#else
extern void f04jgy();
#endif

#ifdef NAG_PROTO
extern void f04jgz(Integer m, Integer n, Boolean NAG_HUGE *svd, Integer rank,
            double NAG_HUGE b[], double NAG_HUGE a[], Integer tda, double NAG_HUGE sv[],
            double NAG_HUGE x[], double NAG_HUGE *sigma, double NAG_HUGE work[], Integer NAG_HUGE *ierr);
#else
extern void f04jgz();
#endif

#ifdef NAG_PROTO
extern void f04lef(Integer job, Integer n,  double a[], double b[],
	    double c[], double d[],  Integer in[],  double y[],
	    double *tol,  Integer *ifail);
#else
     extern void f04lef();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f04mcc(Nag_SolveSystem selct, Integer n, Integer nrhs, double NAG_HUGE *al, Integer lal, 
            double NAG_HUGE *d, Integer NAG_HUGE *row, double NAG_HUGE *b, Integer tdb, double NAG_HUGE *x, 
            Integer tdx, NagError NAG_HUGE *fail);
#else
extern void f04mcc();
#endif

#ifdef NAG_PROTO
extern void f04mcv(Integer n1, Integer n2, double NAG_HUGE *u, Integer iu1, Integer iu2, 
            Integer lu, double scale, double NAG_HUGE *v, Integer iv1, Integer iv2,
            Integer lv);
#else
extern void f04mcv();
#endif

#ifdef NAG_PROTO
extern void f04mcw(Integer n,double NAG_HUGE *l,Integer ll,Integer NAG_HUGE *nrow,Integer p,
            double NAG_HUGE *b,Integer ib1,Integer ib2,Integer lb,double NAG_HUGE *x,
            Integer ix1,Integer ix2,Integer lx);
#else
extern void f04mcw();
#endif

#ifdef NAG_PROTO
extern void f04mcx(Integer n, double NAG_HUGE *l, Integer ll, Integer NAG_HUGE *nrow, Integer p,
            double NAG_HUGE *b, Integer ib1, Integer ib2, Integer lb, double NAG_HUGE *x,
            Integer ix1, Integer ix2, Integer lx);
#else
extern void f04mcx();
#endif

#ifdef NAG_PROTO
extern void f04mcy(Integer n, double NAG_HUGE *d, Integer p, double NAG_HUGE *b, Integer ib1,
            Integer ib2, Integer lb, double NAG_HUGE *x, Integer ix1, Integer ix2, 
            Integer lx);
#else
extern void f04mcy();
#endif

#ifdef NAG_PROTO
extern void f04mcz(Integer n, double NAG_HUGE *l, Integer ll, double NAG_HUGE *d, Integer NAG_HUGE *nrow,
            Integer p, double NAG_HUGE *b, Integer ib1, Integer ib2, Integer lb,
            Nag_SolveSystem iselct, double NAG_HUGE *x, Integer ix1, Integer ix2, Integer lx);
#else
extern void f04mcz();
#endif

#ifdef NAG_PROTO
extern void f04yac(Integer job, Integer p, double sigma, double NAG_HUGE a[],
            Integer tda, Boolean svd, Integer rank, double NAG_HUGE sv[],
            double NAG_HUGE cj[], NagError NAG_HUGE *fail);
#else
extern void f04yac();
#endif

#ifdef NAG_PROTO
extern void f04yay(Integer job, Integer n, double NAG_HUGE t[], Integer tdt,
            double NAG_HUGE b[], Integer NAG_HUGE *ierror);
#else
extern void f04yay();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGF04 */
