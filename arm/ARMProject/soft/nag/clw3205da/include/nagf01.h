#ifndef NAGF01
#define NAGF01
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf01.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f01 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2149 (Feb 1998).
 */

#ifdef NAG_PROTO
extern void f01aac(double NAG_HUGE a[], Integer tda, Integer n, double NAG_HUGE x[], Integer tdx,
            NagError NAG_HUGE *fail);
#else
extern void f01aac();
#endif

#ifdef NAG_PROTO
extern void f01acc(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[], Integer tdb,
            double NAG_HUGE z[], Integer NAG_HUGE *ell, NagError NAG_HUGE *fail);
#else
extern void f01acc();
#endif

#ifdef NAG_PROTO
extern void f01acz(Integer n, double eps, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, double NAG_HUGE z[], Integer NAG_HUGE *l, Integer NAG_HUGE *ifail);
#else
extern void f01acz();
#endif

#ifdef NAG_PROTO
extern void f01adc(Integer n, double NAG_HUGE a[], Integer tda, NagError NAG_HUGE *fail);
#else
extern void f01adc();
#endif

#ifdef NAG_PROTO
extern void f01aec(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE b[],
            Integer tdb, double NAG_HUGE dl[], NagError NAG_HUGE *fail);
#else
extern void f01aec();
#endif

#ifdef NAG_PROTO
extern void f01afc(Integer n, Integer im1, Integer im2, double NAG_HUGE b[], Integer tdb,
            double NAG_HUGE dl[], double NAG_HUGE z[], Integer tdz);
#else
extern void f01afc();
#endif

#ifdef NAG_PROTO
extern void f01agc(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE d[], double NAG_HUGE e[],
            double NAG_HUGE e2[]);
#else
extern void f01agc();
#endif

#ifdef NAG_PROTO
extern void f01ajc(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE d[], 
            double NAG_HUGE e[], double NAG_HUGE z[], Integer tdz);
#else
extern void f01ajc();
#endif

#ifdef NAG_PROTO
extern void f01akc(Integer n, Integer k, Integer l,  double NAG_HUGE a[],  Integer tda,
            Integer NAG_HUGE intger[]);
#else
extern void f01akc();
#endif

#ifdef NAG_PROTO
extern void f01apc(Integer n, Integer low, Integer iupp, Integer NAG_HUGE intger[],
            double NAG_HUGE h[],  Integer tdh,  double NAG_HUGE v[],  Integer tdv);
#else
extern void f01apc();
#endif

#ifdef NAG_PROTO
extern void f01atc(Integer n, Integer ib, double NAG_HUGE a[], Integer tda, Integer NAG_HUGE *low,
            Integer NAG_HUGE *lhi, double NAG_HUGE d[]);
#else
extern void f01atc();
#endif

#ifdef NAG_PROTO
extern void f01atz(Integer m, double NAG_HUGE *a, Integer tda, double NAG_HUGE *d, 
            Integer k, Integer l, Integer n, Integer j);
#else
extern void f01atz();
#endif

#ifdef NAG_PROTO
extern void f01auc(Integer n, Integer low, Integer lhi, Integer m,  double NAG_HUGE d[],
            double NAG_HUGE z[],  Integer tdz);
#else
extern void f01auc();
#endif

#ifdef NAG_PROTO
extern void f01bcc(Integer n,  Complex NAG_HUGE a[], Integer tda,
            double NAG_HUGE d[], double NAG_HUGE e[],
            double NAG_HUGE c[], double NAG_HUGE s[]);
#else
extern void f01bcc();
#endif

#ifdef NAG_PROTO
extern void f01bcy( Integer m, Integer n, Complex NAG_HUGE a[], Integer tda,
            Complex NAG_HUGE b[], Integer incb, double NAG_HUGE cr[], double NAG_HUGE ci[]);
#else
extern void f01bcy();
#endif

#ifdef NAG_PROTO
extern void f01bcz(Integer n, Complex NAG_HUGE a[], Integer tda, double NAG_HUGE br[],
            double NAG_HUGE bi[], double NAG_HUGE cr[], double NAG_HUGE ci[]);
#else
extern void f01bcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01bnc(Integer n, Complex NAG_HUGE a[], Integer tda, double NAG_HUGE p[], NagError NAG_HUGE *fail);
#else
extern void f01bnc();
#endif

#ifdef NAG_PROTO
extern void f01bqc(Integer n, double NAG_HUGE *eps, double NAG_HUGE rl[],
            double NAG_HUGE d[], Integer NAG_HUGE *row, NagError NAG_HUGE *fail);
#else
extern void f01bqc();
#endif

#ifdef NAG_PROTO
extern void f01bqz(Integer n,  double *eps, double rl[],  Integer irl,
double d[],  Integer *ifail);
#else
extern void f01bqz();
#endif

#ifdef NAG_PROTO
extern void f01brg(double NAG_HUGE mat[], 
            Integer col_num,
            Integer lead_dim,
            Integer trail_dim,
            double NAG_HUGE vec[], 
            Integer num_elements);
#else
extern void f01brg();
#endif

#ifdef NAG_PROTO
extern void f01brh(double NAG_HUGE mat[], 
            Integer col_num,
            Integer lead_dim,
            Integer trail_dim,
            double NAG_HUGE vec[], 
            Integer num_elements);
#else
extern void f01brh();
#endif

#ifdef NAG_PROTO
extern void f01bri(double NAG_HUGE smat[], 
            double lmat [],
            Integer sr,
            Integer sc,
            Integer lr,
            Integer lc);
#else
extern void f01bri();
#endif

#ifdef NAG_PROTO
extern void f01brj(double NAG_HUGE smat[], 
            double lmat [],
            Integer sr,
            Integer sc,
            Integer lr,
            Integer lc);
#else
extern void f01brj();
#endif

#ifdef NAG_PROTO
extern void f01brk(Integer NAG_HUGE mat[], 
            Integer col_num,
            Integer lead_dim,
            Integer trail_dim,
            Integer NAG_HUGE vec[], 
            Integer num_elements);
#else
extern void f01brk();
#endif

#ifdef NAG_PROTO
extern void f01brl(Integer NAG_HUGE mat[], 
            Integer col_num,
            Integer lead_dim,
            Integer trail_dim,
            Integer NAG_HUGE vec[], 
            Integer num_elements);
#else
extern void f01brl();
#endif

#ifdef NAG_PROTO
extern void f01crc(double NAG_HUGE *a, Integer m, Integer n, NagError NAG_HUGE *fail);
#else
extern void f01crc();
#endif

#ifdef NAG_PROTO
extern void f01lef(Integer n,  double a[], double lambda, double b[],
	    double c[], double tol, double d[],  Integer in[],
	    Integer *ifail);
#else
     extern void f01lef();
#endif

#ifdef NAG_PROTO
extern void f01lzc(Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE c[],
            Integer tdc,  Boolean wantb, double NAG_HUGE b[], Boolean wantq,
            Boolean wanty, double NAG_HUGE y[], Integer tdy, Integer ly,
            Boolean wantz, double NAG_HUGE z[], Integer tdz, Integer ncz,
            double NAG_HUGE d[], double NAG_HUGE e[], NagError NAG_HUGE *fail);
#else
extern void f01lzc();
#endif

#ifdef NAG_PROTO
extern void f01lzw(double t, double NAG_HUGE *c, double NAG_HUGE *s, double sqteps,
            double rsqtps, double big);
#else
extern void f01lzw();
#endif

#ifdef NAG_PROTO
extern void f01lzx(Integer n, double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE x[]);
#else
extern void f01lzx();
#endif

#ifdef NAG_PROTO
extern void f01lzy(Integer n, double c, double s, double NAG_HUGE *x, double NAG_HUGE *y);
#else
extern void f01lzy();
#endif

#ifdef NAG_PROTO
extern double f01lzz(double a, double b, double small1, double big);
#else 
extern double f01lzz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01mcc(Integer n, double NAG_HUGE *a, Integer lal, Integer NAG_HUGE *row, double NAG_HUGE *al, 
            double NAG_HUGE *d, NagError NAG_HUGE *fail);
#else
extern void f01mcc();
#endif

#ifdef NAG_PROTO
extern void f01qac(Integer m, Integer n, double NAG_HUGE a[], Integer tda, double NAG_HUGE c[],
            Integer tdc, double NAG_HUGE z[], NagError NAG_HUGE *fail);
#else
extern void f01qac();
#endif

#ifdef NAG_PROTO
extern void f01qaw(Integer n, double NAG_HUGE *x, double xmul, double NAG_HUGE *y, Boolean undflw);
#else
extern void f01qaw();
#endif

#ifdef NAG_PROTO
extern double f01qax(Integer nr, Integer n, double v, Boolean plus,
              double NAG_HUGE *x, double NAG_HUGE *y, Boolean undflw);
#else
extern double f01qax();
#endif

#ifdef NAG_PROTO
extern void f01qay(Integer n, double NAG_HUGE *x, Boolean norm, double NAG_HUGE *z1, double small1,
            double tiny, double big);
#else
extern void f01qay();
#endif

#ifdef NAG_PROTO
extern void f01qaz(Integer n, double NAG_HUGE *z, double z1, double NAG_HUGE *x);
#else
extern void f01qaz();
#endif

#ifdef NAG_PROTO
extern void f01qbc(Integer m, Integer n,  double NAG_HUGE a[],  Integer tda,
             double NAG_HUGE c[],  Integer tdc,  double NAG_HUGE work[],
             NagError NAG_HUGE *fail);
#else
extern void f01qbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01qcc(Integer m, Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *zeta,
            NagError NAG_HUGE *fail);
#else
extern void f01qcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01qdc(MatrixTranspose trans, Nag_WhereElements wheret,
            Integer m, Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *zeta,
            Integer ncolb, double NAG_HUGE *b, Integer tdb, NagError NAG_HUGE *fail);
#else
extern void f01qdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01qec(Nag_WhereElements wheret, Integer m, Integer n, Integer ncolq,
            double NAG_HUGE *a, Integer tda, double NAG_HUGE *zeta, NagError NAG_HUGE *fail);
#else
extern void f01qec();
#endif

#ifdef NAG_PROTO
extern void f01qfc(char NAG_HUGE *pivot,  Integer m, Integer n,  double NAG_HUGE a[],
            Integer lda,  double NAG_HUGE zeta[],  Integer NAG_HUGE perm[],
            double NAG_HUGE work[],  Integer NAG_HUGE *ifail);
#else
extern void f01qfc();
#endif

#ifdef NAG_PROTO
extern void f01qjc(Integer m, Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *zeta,
            NagError NAG_HUGE *fail);
#else
extern void f01qjc();
#endif

#ifdef NAG_PROTO
extern void f01qkc(Nag_WhereElements wheret, Integer m, Integer n, Integer nrowp,
            double NAG_HUGE *a, Integer tda, double NAG_HUGE *zeta, NagError NAG_HUGE *fail);
#else
extern void f01qkc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01rcc(Integer m, Integer n, Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *theta,
            NagError NAG_HUGE *fail);
#else
extern void f01rcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01rdc(MatrixTranspose trans, Nag_WhereElements wheret, Integer m, Integer n,
            Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *theta, Integer ncolb, 
            Complex NAG_HUGE *b, Integer tdb, NagError NAG_HUGE *fail);
#else
extern void f01rdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f01rec(Nag_WhereElements wheret, Integer m, Integer n, Integer ncolq, 
            Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *theta, NagError NAG_HUGE *fail);
#else
extern void f01rec();
#endif

#ifdef NAG_PROTO
extern void f01rjc(Integer m, Integer n, Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *theta,
            NagError NAG_HUGE *fail);
#else
extern void f01rjc();
#endif

#ifdef NAG_PROTO
extern void f01rkc(Nag_WhereElements wheret, Integer m, Integer n, Integer rowp,
            Complex NAG_HUGE *a, Integer tda, Complex NAG_HUGE *theta, NagError NAG_HUGE *fail);
#else
extern void f01rkc();
#endif

#ifdef NAG_PROTO
extern void f01sag(Integer mode, Integer m, Integer n,  double v[],
            double w[],  Integer lena, Integer luparm[],
            double parmlu[], double a[],  Integer indc[],
            Integer indr[], Integer ip[], Integer iq[], Integer lenc[],
            Integer lenr[], Integer locc[], Integer locr[], 
            Integer *inform);
#else
     extern void f01sag();
#endif

#ifdef NAG_PROTO
extern void f01sah(Integer m, Integer n,  double w[],  Integer lena,
	    Integer luparm[],  double parmlu[], double a[],
	    Integer indc[], Integer indr[], Integer ip[], 
	    Integer iq[], Integer lenc[], Integer lenr[], 
	    Integer locc[], Integer locr[], Integer *inform);
#else
     extern void f01sah();
#endif

#ifdef NAG_PROTO
extern void f01saj(Integer n,  Boolean reals,  Integer luparm[], Integer *ltop,
	    Integer lena,  double a[],  Integer ind[], Integer len[],
	    Integer loc[]);
#else
     extern void f01saj();
#endif

#ifdef NAG_PROTO
extern void f01sak(Integer n, Integer len[], Integer iperm[], Integer iw[],
	    Integer *nrank);
#else
     extern void f01sak();
#endif

#ifdef NAG_PROTO
extern void f01sal(Integer nzpiv, Integer *nzchng, Integer indr[], 
	    Integer lenold[], Integer lennew[], Integer iqloc[], 
	    Integer iq[], Integer iqinv[]);
#else
     extern void f01sal();
#endif

#ifdef NAG_PROTO
extern void f01sam(Integer m, Integer n, Integer len[], Integer iperm[],
             Integer loc[], Integer inv[], Integer num[]);
#else
extern void f01sam();
#endif

#ifdef NAG_PROTO
extern void f01san(Integer m, Integer melim, Integer ncold, Integer nspare,
	    Integer lpivc1, Integer lpivc2, Integer lpivr1, Integer lpivr2,
	    Integer *lrow, Integer lenc[], Integer lenr[], Integer locc[],
	    Integer locr[], Integer indc[], Integer indr[], 
	    Integer ifill[], Integer jfill[]);
#else
extern void f01san();
#endif

#ifdef NAG_PROTO
extern void f01sap(Integer m, Integer n, Integer nelem, Integer indc[],
	    Integer indr[], Integer lenc[], Integer lenr[], 
	    Integer locc[], Integer locr[]);
#else
extern void f01sap();
#endif

#ifdef NAG_PROTO
extern void f01saq(Integer m, Integer n, Integer nelem, Integer indc[],
	    Integer lenc[], Integer locc[], Integer iw[], Integer *lerr,
	    Integer *inform);
#else
extern void f01saq();
#endif

#ifdef NAG_PROTO
extern void f01sar(Integer n, Integer numa,  double a[],  Integer inum[],
	    Integer jnum[], Integer len[], Integer loc[]);
#else
extern void f01sar();
#endif

#ifdef NAG_PROTO
extern void f01sas(Integer m, Integer n, Integer nelem,  double small1,
	    double a[],  Integer indc[], Integer indr[], Integer lenc[],
	    Integer lenr[],  double *amax,  Integer *numnz,
	    Integer *lerr, Integer *inform);
#else
extern void f01sas();
#endif

#ifdef NAG_PROTO
extern void f01sat(Integer kol1, Integer kol2, Integer kol[],  double a[],
	    Integer indc[], Integer lenc[], Integer locc[]);
#else
extern void f01sat();
#endif

#ifdef NAG_PROTO
extern void f01sau(Integer m, Integer n, Integer lena, Integer maxmn,
	    double lmax,  Integer maxcol, Integer maxrow, Integer *ibest,
	    Integer *jbest, Integer *mbest,  double a[],  Integer indc[],
	    Integer indr[], Integer ip[], Integer iq[], Integer lenc[],
	    Integer lenr[], Integer locc[], Integer locr[], Integer iploc[],
	    Integer iqloc[]);
#else
extern void f01sau();
#endif

#ifdef NAG_PROTO
extern void f01sav(Integer m, Integer melim, Integer ncold, Integer nspare,
	    double small1,  Integer lpivc1, Integer lpivc2,
	    Integer *lfirst, Integer lpivr2, Integer lfree, Integer minfre,
	    Integer *lrow, Integer *lcol, Integer *lu, Integer *nfill,
	    double a[],  Integer indc[], Integer indr[], Integer lenc[],
	    Integer lenr[], Integer locc[], Integer locr[], Integer mark[],
	    double al[],  Integer markl[],  double au[],
	    Integer ifill[], Integer jfill[]);
#else
extern void f01sav();
#endif

#ifdef NAG_PROTO
extern void f01saw(Integer m, Integer n, Integer lena, Integer lend, Integer lu1,
             Integer mleft, Integer nleft, Integer nrank, Integer nrowu,
             Integer *lenl, Integer *lenu, Integer *nsing,  Boolean keeplu,
             double small1, double a[], double d[],
             Integer indc[], Integer indr[], Integer ip[], Integer iq[],
             Integer lenc[], Integer lenr[], Integer locc[], Integer ipinv[],
             Integer ipvt[]);
#else
extern void f01saw();
#endif

#ifdef NAG_PROTO
extern void f01sax(Integer m, Integer n, Integer nelem, Integer lena,
	    Integer luparm[],  double parmlu[], double a[],
	    Integer indc[], Integer indr[], Integer ip[], Integer iq[],
	    Integer lenc[], Integer lenr[], Integer locc[], Integer locr[],
	    Integer iploc[], Integer iqloc[], Integer ipinv[], Integer iqinv[],
	    Integer *inform, Integer *lenl, Integer *lenu, Integer *minlen,
	    Integer *mersum, Integer *nutri, Integer *nltri, Integer *nrank);
#else
extern void f01sax();
#endif

#ifdef NAG_PROTO
extern void f01say(Integer m, Integer n, Integer *nelem, Integer lena,
	    Integer luparm[],  double parmlu[], double a[],
	    Integer indc[], Integer indr[], Integer ip[], Integer iq[],
	    Integer lenc[], Integer lenr[], Integer locc[], Integer locr[],
	    Integer iploc[], Integer iqloc[], Integer ipinv[],
	    Integer iqinv[],  double w[],  Integer *inform);
#else
extern void f01say();
#endif

#ifdef NAG_PROTO
extern void f01saz(double a[],  Integer lda, Integer m, Integer n,
	    double small1,  Integer *nsing, Integer ipvt[], Integer iq[]);
#else
extern void f01saz();
#endif

#ifdef NAG_PROTO
extern void f01sbt(Integer mode1, Integer mode2, Integer m, Integer n,
	    Integer jrep,  double v[], double w[],  Integer lena,
	    Integer luparm[],  double parmlu[], double a[],
	    Integer indc[], Integer indr[], Integer ip[], Integer iq[],
	    Integer lenc[], Integer lenr[], Integer locc[], Integer locr[],
	    Integer *inform,  double *diag, double *vnorm);
#else
extern void f01sbt();
#endif

#ifdef NAG_PROTO
extern void f01sbu(Integer m, Integer n, Integer jzap, Integer *kzap, 
	    Integer lena, Integer *lenu, Integer *lrow, Integer nrank,  
	    double a[], Integer indr[], Integer ip[], Integer iq[], 
	    Integer lenr[], Integer locr[]);
#else
extern void f01sbu();
#endif

#ifdef NAG_PROTO
extern void f01sbv(Integer m, Integer n, Integer jsing, Integer lena,
	    Integer luparm[],  double parmlu[],  Integer *lenu,
	    Integer *lrow, Integer *nrank,  double a[],  Integer indc[],
	    Integer indr[], Integer ip[], Integer iq[], Integer lenr[],
	    Integer locc[], Integer locr[], Integer *inform,
	    double *diag);
#else
extern void f01sbv();
#endif

#ifdef NAG_PROTO
extern void f01sbw(Integer m, Integer n, Integer kfirst, Integer klast,
	    Integer lena, Integer luparm[],  double parmlu[],
	    Integer *lenl, Integer *lenu, Integer *lrow,  double a[],
	    Integer indc[], Integer indr[], Integer ip[], Integer iq[],
	    Integer lenr[], Integer locc[], Integer locr[], 
	    Integer *inform, double *diag);
#else
extern void f01sbw();
#endif

#ifdef NAG_PROTO
extern void f01sbx(Integer m, Integer n, Integer jelm,  double v[],
	    Integer lena, Integer luparm[],  double parmlu[],
	    Integer *lenl, Integer *lrow, Integer nrank,  double a[],
	    Integer indc[], Integer indr[], Integer ip[], Integer iq[],
	    Integer lenr[], Integer locc[], Integer locr[], 
	    Integer *inform, double *diag);
#else
extern void f01sbx();
#endif

#ifdef NAG_PROTO
extern void f01sby(Integer kfirst, Integer klast, Integer ip[]);
#else
extern void f01sby();
#endif

#ifdef NAG_PROTO
extern void f01sbz(Integer m, Integer jadd,  double v[],  Integer lena,
	    Integer luparm[],  double parmlu[],  Integer lenl,
	    Integer *lenu, Integer *lrow, Integer nrank,  double a[],
	    Integer indr[], Integer ip[], Integer lenr[], Integer locr[],
	    Integer *inform, Integer *klast,  double *vnorm);
#else
extern void f01sbz();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGF01 */
