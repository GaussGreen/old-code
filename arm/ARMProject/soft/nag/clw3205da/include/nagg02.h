#ifndef NAGG02
#define NAGG02
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg02.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g02 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2156 (Feb 1998).
 */

/* Message codes and extern of message list */
#include <nag_g02mesg.h>

#ifdef NAG_PROTO
extern void g02aax(MatrixTriangle uplo, Integer n, double NAG_HUGE *a, double NAG_HUGE *b);
#else
extern void g02aax();
#endif

#ifdef NAG_PROTO
extern void g02aay(MatrixTriangle uplo, MatrixUnitTriangular diag, Integer n,
            double NAG_HUGE *a);
#else
extern void g02aay();
#endif

#ifdef NAG_PROTO
extern void g02aaz(MatrixTriangle uplo, MatrixUnitTriangular diag, Integer n,
            double NAG_HUGE *a);
#else
extern void g02aaz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02brc(Integer n, Integer m, double NAG_HUGE x[], Integer tdx, 
            Integer NAG_HUGE svar[], Integer NAG_HUGE sobs[], double NAG_HUGE corr[], 
            Integer tdc, NagError NAG_HUGE *fail);
#else
extern void g02brc();
#endif

#ifdef NAG_PROTO
extern void g02buc(Nag_IncludeMean mean, Integer n, Integer m, double NAG_HUGE x[],
            Integer tdx, double NAG_HUGE wt[], double NAG_HUGE *sw, double NAG_HUGE wmean[],
            double NAG_HUGE c[], NagError NAG_HUGE *fail);
#else
extern void g02buc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02bxc(Integer n, Integer m, double NAG_HUGE x[], 
              Integer tdx, Integer NAG_HUGE sx[], double NAG_HUGE wt[], double NAG_HUGE *sw, 
              double NAG_HUGE wmean[], double NAG_HUGE std[], double NAG_HUGE r[], Integer tdr,
              double NAG_HUGE v[], Integer tdv, NagError NAG_HUGE *fail);
#else
extern void g02bxc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02cac(Nag_SumSquare mean, Integer n, 
            double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE wt[], double NAG_HUGE *a, double NAG_HUGE *b, 
            double NAG_HUGE *a_serr, double NAG_HUGE *b_serr, double NAG_HUGE *rsq, double NAG_HUGE *rss,
            double NAG_HUGE *df, NagError NAG_HUGE *fail);
#else
extern void g02cac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02cbc(Nag_SumSquare mean, Integer n, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE wt[], double clm, double clp, double NAG_HUGE yhat[],
            double NAG_HUGE yml[], double NAG_HUGE ymu[], double NAG_HUGE yl[], double NAG_HUGE yu[],
            double NAG_HUGE h[], double NAG_HUGE res[], double NAG_HUGE *rms, NagError NAG_HUGE *fail);
#else
extern void g02cbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dac(Nag_IncludeMean mean,Integer n,double NAG_HUGE x[], Integer tdx,Integer m,
            Integer NAG_HUGE sx[], Integer ip, double NAG_HUGE y[], double NAG_HUGE wt[], double NAG_HUGE *rss,
            double NAG_HUGE *df, double NAG_HUGE b[], double NAG_HUGE se[], double NAG_HUGE cov[], double NAG_HUGE res[],
            double NAG_HUGE h[], double NAG_HUGE q[], Integer tdq, Boolean NAG_HUGE *svd, Integer NAG_HUGE *rank,
            double NAG_HUGE p[], double tol, double NAG_HUGE com_ar[], NagError NAG_HUGE *fail);
#else
extern void g02dac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dcc(Nag_UpdateObserv update, Nag_IncludeMean mean, Integer m,
            Integer NAG_HUGE sx[], double NAG_HUGE q[], Integer tdq, Integer ip,
            double NAG_HUGE x[],  Integer nr, Integer tdx, Integer  ix,
            double y, double NAG_HUGE *wt, double NAG_HUGE *rss, NagError NAG_HUGE *fail);
#else
extern void g02dcc();
#endif

#ifdef NAG_PROTO
extern void g02dcz(Integer n, double alpha, double NAG_HUGE x[], Integer  incx, double NAG_HUGE a[],
            Integer tda, double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE *zeta);
#else
extern void g02dcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02ddc(Integer n,Integer ip, double NAG_HUGE q[], Integer tdq, double NAG_HUGE *rss,
            double NAG_HUGE *df, double NAG_HUGE b[], double NAG_HUGE se[], double NAG_HUGE cov[], Boolean NAG_HUGE *svd,
            Integer NAG_HUGE *rank,double NAG_HUGE p[], double tol, NagError NAG_HUGE *fail);
#else
extern void g02ddc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dec(Integer n, Integer ip, double NAG_HUGE q[], Integer tdq, double NAG_HUGE p[],
            double NAG_HUGE wt[], double NAG_HUGE x[], double NAG_HUGE *rss, double tol, NagError NAG_HUGE *fail);
#else
extern void g02dec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dfc(Integer ip, double NAG_HUGE q[], Integer tdq, Integer indx, double NAG_HUGE *rss,
            NagError NAG_HUGE *fail);
#else
extern void g02dfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dgc(Integer n,double NAG_HUGE wt[], double NAG_HUGE *rss, Integer ip, Integer rank,
            double NAG_HUGE cov[], double NAG_HUGE q[], Integer tdq, Boolean svd, double NAG_HUGE p[],
            double NAG_HUGE y[], double NAG_HUGE b[], double NAG_HUGE se[], double NAG_HUGE res[],
            double NAG_HUGE com_ar[], NagError NAG_HUGE *fail);
#else
extern void g02dgc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dkc(Integer ip, Integer iconst, double NAG_HUGE p[], double NAG_HUGE c[],
            Integer tdc, double NAG_HUGE b[], double rss, double df,
            double NAG_HUGE se[], double NAG_HUGE cov[], NagError NAG_HUGE *fail);
#else
extern void g02dkc();
#endif

#ifdef NAG_PROTO
extern void g02dkz(Integer ip, Integer iconst, double NAG_HUGE p[], Integer imp,
            double NAG_HUGE c[], Integer tdc, double NAG_HUGE b[], double s,
            double NAG_HUGE se[], double NAG_HUGE cov[], Integer NAG_HUGE *ind);
#else
extern void g02dkz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02dnc(Integer ip, Integer rank, double NAG_HUGE b[], double NAG_HUGE cov[],
            double NAG_HUGE p[], double NAG_HUGE f[], Boolean NAG_HUGE *est, double NAG_HUGE *stat,
            double NAG_HUGE *sestat, double NAG_HUGE *t, double tol, NagError NAG_HUGE *fail);
#else
extern void g02dnc();
#endif

#ifdef NAG_PROTO
extern void g02dnz(Integer ip, Integer rank, double NAG_HUGE b[], double NAG_HUGE cov[], double NAG_HUGE p[],
            Integer imp, double NAG_HUGE f[], Boolean NAG_HUGE *est, double NAG_HUGE *stat,
            double NAG_HUGE *sestat, double NAG_HUGE *t, double tol, Integer NAG_HUGE *ind);
#else
extern void g02dnz();
#endif

#ifdef NAG_PROTO
extern Integer g02eaz(Integer m, Integer n, Integer *ifail);
#else
extern Integer g02eaz();
#endif


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02fac(Integer n, Integer ip, Integer nres, double NAG_HUGE res[], double NAG_HUGE h[],
            double rms, double NAG_HUGE sres[], NagError NAG_HUGE *fail);
#else
extern void g02fac();
#endif

#ifdef NAG_PROTO

extern NAG_DLL_EXPIMP void NAG_CALL g02gac(Nag_Link link, Nag_IncludeMean mean, Integer n, double NAG_HUGE x[],
            Integer tdx, Integer m, Integer NAG_HUGE sx[], Integer ip, double NAG_HUGE y[],
            double NAG_HUGE wt[], double NAG_HUGE offset[], double NAG_HUGE *scale, double ex_power,
            double NAG_HUGE *rss, double NAG_HUGE *df, double NAG_HUGE b[], Integer NAG_HUGE *rank, double NAG_HUGE se[],
            double NAG_HUGE cov[], double NAG_HUGE v[], Integer tdv, double tol,
            Integer max_iter, Integer print_iter, char NAG_HUGE *outfile, double eps,
            NagError NAG_HUGE *fail);
#else
extern void g02gac();
#endif

#ifdef NAG_PROTO
extern double g02gav(double fv, double y, double NAG_HUGE *t);
#else
extern double g02gav();
#endif

#ifdef NAG_PROTO
extern void g02gaw(Integer n, double NAG_HUGE fv[], double NAG_HUGE t[], double NAG_HUGE var[], double NAG_HUGE wt[], 
            Integer no);
#else
extern void g02gaw();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02gbc(Nag_Link link, Nag_IncludeMean mean, Integer n, double NAG_HUGE x[],
            Integer tdx, Integer m, Integer NAG_HUGE sx[], Integer ip, double NAG_HUGE y[],
            double NAG_HUGE binom_t[], double NAG_HUGE wt[], double NAG_HUGE offset[], double NAG_HUGE *dev,
            double NAG_HUGE *df, double NAG_HUGE b[], Integer NAG_HUGE *rank, double NAG_HUGE se[],
            double NAG_HUGE cov[], double NAG_HUGE v[], Integer tdv, double tol,
            Integer max_iter, Integer print_iter, char NAG_HUGE *outfile,
            double eps, NagError NAG_HUGE *fail);
#else
extern void g02gbc();
#endif

#ifdef NAG_PROTO
extern void g02gbs(Integer ip, Integer rank, Boolean svd, double NAG_HUGE q[], Integer tdq,
            double NAG_HUGE cov[], double NAG_HUGE wk[]);
#else
extern void g02gbs();
#endif

#ifdef NAG_PROTO
extern void g02gbt(Nag_IncludeMean mean, Integer n, Integer m, double NAG_HUGE x[],
            Integer tdx, Integer NAG_HUGE isx[], Integer ip,  double NAG_HUGE q[],
            Integer tdq, Boolean svd, Integer rank, double NAG_HUGE wwt[], 
            double NAG_HUGE h[], double NAG_HUGE wk[]);
#else
extern void g02gbt();
#endif

#ifdef NAG_PROTO
extern void g02gbu(Integer n, char link, double NAG_HUGE y[], double NAG_HUGE t[], double NAG_HUGE fv[], 
            double NAG_HUGE eta[], double NAG_HUGE wt[], Integer no, Integer NAG_HUGE *ind);
#else
extern void g02gbu();
#endif

#ifdef NAG_PROTO
extern double g02gbv(double fv, double y, double NAG_HUGE *t);
#else
extern double g02gbv();
#endif

#ifdef NAG_PROTO
extern void g02gbw(Integer n,  double NAG_HUGE fv[], double NAG_HUGE t[], double NAG_HUGE var[], double NAG_HUGE wt[],
            Integer no);
#else
extern void g02gbw();
#endif

#ifdef NAG_PROTO
extern void g02gbx(Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE t[], double NAG_HUGE der[],
             double a, double NAG_HUGE wt[], Integer no);
#else
extern void g02gbx();
#endif

#ifdef NAG_PROTO
extern void g02gby(Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE fv[], double NAG_HUGE t[],
            double NAG_HUGE *a, double NAG_HUGE wt[], Integer no, Integer NAG_HUGE *ind);
#else
extern void g02gby();
#endif

#ifdef NAG_PROTO
extern void g02gbz(void (NAG_HUGE *clink) (Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE fv[],
                           double NAG_HUGE t[], double NAG_HUGE *a, double NAG_HUGE wt[], Integer no, Integer NAG_HUGE *ind),
            void (NAG_HUGE *cder) (Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE t[],
                          double NAG_HUGE der[], double a, double NAG_HUGE wt[], Integer no),
            double (NAG_HUGE *cdev) (double fv, double y, double NAG_HUGE *t),
            void (NAG_HUGE *cvar) (Integer n,  double NAG_HUGE fv[], double NAG_HUGE t[], double NAG_HUGE var[],
                          double NAG_HUGE wt[], Integer no),
            char link, Nag_IncludeMean mean, char weight, Integer n,
            double NAG_HUGE x[], Integer tdx, Integer m, Integer NAG_HUGE isx[], double NAG_HUGE y[],
            double NAG_HUGE t[], double NAG_HUGE wt[], Integer no, double NAG_HUGE *dev, Integer NAG_HUGE *rank,
            double NAG_HUGE b[], Integer ip, double NAG_HUGE fv[], double NAG_HUGE eta[], double NAG_HUGE var[],
            double NAG_HUGE wwt[], double NAG_HUGE off[], double NAG_HUGE q[], Integer tdq, double NAG_HUGE *a,
            Boolean NAG_HUGE *svd, double tol, Integer maxit, Integer NAG_HUGE *iter, 
            Integer iprint, Nag_FileSt NAG_HUGE *stream, double eps, double NAG_HUGE wk[], Integer NAG_HUGE *ierror);
#else
extern void g02gbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02gcc(Nag_Link link, Nag_IncludeMean mean, Integer n, double NAG_HUGE x[],
            Integer tdx, Integer m, Integer NAG_HUGE sx[], Integer ip, double NAG_HUGE y[], 
            double NAG_HUGE wt[], double NAG_HUGE offset[], double ex_power, double NAG_HUGE *dev,
            double NAG_HUGE *df, double NAG_HUGE b[], Integer NAG_HUGE *rank, double NAG_HUGE se[], double NAG_HUGE cov[], 
            double NAG_HUGE v[], Integer tdv, double tol, Integer max_iter, 
            Integer print_iter, char NAG_HUGE *outfile, double eps, NagError NAG_HUGE *fail);
#else
extern void g02gcc();
#endif

#ifdef NAG_PROTO
extern void g02gcu(Integer n, char link, double a, double NAG_HUGE y[], double NAG_HUGE fv[], 
            double NAG_HUGE eta[], double NAG_HUGE wt[], Integer no, double ymin, double y0);
#else
extern void g02gcu();
#endif

#ifdef NAG_PROTO
extern double g02gcv(double fv, double y, double NAG_HUGE *t);
#else
extern double g02gcv();
#endif

#ifdef NAG_PROTO
extern void g02gcw(Integer n,  double NAG_HUGE fv[], double NAG_HUGE t[], double NAG_HUGE var[],
             double NAG_HUGE wt[],  Integer no);
#else
extern void g02gcw();
#endif

#ifdef NAG_PROTO
extern void g02gcx(Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE t[],
            double NAG_HUGE der[], double a, double NAG_HUGE wt[], Integer no);
#else
extern void g02gcx();
#endif

#ifdef NAG_PROTO
extern void g02gcy(Integer n, char link, double NAG_HUGE eta[], double NAG_HUGE fv[],
            double NAG_HUGE t[], double NAG_HUGE *a, double NAG_HUGE wt[], Integer no,
            Integer NAG_HUGE *ind);
#else
extern void g02gcy();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02gdc(Nag_Link link, Nag_IncludeMean mean, Integer n, double NAG_HUGE x[],
            Integer tdx, Integer m, Integer NAG_HUGE sx[], Integer ip, double NAG_HUGE y[],
            double NAG_HUGE wt[], double NAG_HUGE offset[], double NAG_HUGE *scale, double ex_power,
            double NAG_HUGE *dev, double NAG_HUGE *df, double NAG_HUGE b[], Integer NAG_HUGE *rank, double NAG_HUGE se[],
            double NAG_HUGE cov[], double NAG_HUGE v[], Integer tdv, double tol,
            Integer max_iter, Integer print_iter, char NAG_HUGE *outfile, double eps,
            NagError NAG_HUGE *fail);
#else
extern void g02gdc();
#endif

#ifdef NAG_PROTO
extern double g02gdv(double fv, double y, double NAG_HUGE *t);
#else
extern double g02gdv();
#endif

#ifdef NAG_PROTO
extern void g02gdw(Integer n,  double NAG_HUGE fv[], double NAG_HUGE t[], double NAG_HUGE var[],
             double NAG_HUGE wt[],  Integer no);
#else
extern void g02gdw();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02gkc(Integer ip, Integer nclin, double NAG_HUGE v[], Integer tdv, double NAG_HUGE c[],
            Integer tdc, double NAG_HUGE b[], double scale, double NAG_HUGE se[], double NAG_HUGE cov[],
            NagError NAG_HUGE *fail);
#else
extern void g02gkc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02gnc(Integer ip, Integer rank, double NAG_HUGE b[], double NAG_HUGE cov[], double NAG_HUGE v[],
            Integer tdv, double NAG_HUGE f[], Boolean NAG_HUGE *est, double NAG_HUGE *stat, 
            double NAG_HUGE *sestat, double NAG_HUGE *z, double tol, NagError NAG_HUGE *fail);
#else
extern void g02gnc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02hac(Nag_RegType regtype, Nag_PsiFun psifun, 
            Nag_SigmaEst sigma_est, Nag_CovMatrixEst covmat_est, Integer n, 
            Integer m, double NAG_HUGE x[], Integer tdx, double NAG_HUGE y[], 
            double cpsi, double NAG_HUGE hpsi[], double cucv, double dchi, 
            double NAG_HUGE theta[], double NAG_HUGE *sigma, double NAG_HUGE c[], Integer tdc, 
            double NAG_HUGE rs[], double NAG_HUGE wt[], double tol, Integer max_iter,
            Integer print_iter, char NAG_HUGE *outfile, double NAG_HUGE info[], NagError NAG_HUGE *fail);
#else
extern void g02hac();
#endif

#ifdef NAG_PROTO
extern void g02hat(Integer isigma, Integer indw, Integer ipsi, Integer n, double d,
            double NAG_HUGE wgt[], double NAG_HUGE *beta, Integer maxit, double tol,
            double NAG_HUGE work[], Integer NAG_HUGE *ifail);
#else
extern void g02hat();
#endif

#ifdef NAG_PROTO
extern double g02hau(double t, Integer ips, double c, double h1, double h2, double h3);
#else
extern double g02hau();
#endif

#ifdef NAG_PROTO
extern double g02hav(double x, Integer ifun, double c);
#else
extern double g02hav();
#endif

#ifdef NAG_PROTO
extern void g02haw(double NAG_HUGE x[], double NAG_HUGE a[], double cucv, Integer n, Integer m,
            Integer mm, Integer mdx, Integer maxit, Integer nitmon,
	    Nag_FileSt NAG_HUGE *stream, double tol, Integer NAG_HUGE *nit, double NAG_HUGE sz[],
	    double NAG_HUGE sc2[], Integer ifun, Integer NAG_HUGE *ifail);
#else
extern void g02haw();
#endif

#ifdef NAG_PROTO
extern void g02hax(Integer indw, Integer ipsi, Integer isigma, Integer n,
	    Integer m,  double NAG_HUGE x[],  Integer tdx,  double NAG_HUGE y[],
	    double beta, double cpsi, double h1, double h2,
	    double h3, double cucv, double d, double NAG_HUGE wgt[],
	    double NAG_HUGE theta[],  Integer NAG_HUGE *k,  double NAG_HUGE *sigma,
	    double NAG_HUGE rs[], double tol,  Integer maxit, Integer nitmon,
	    Nag_FileSt NAG_HUGE *stream, Integer NAG_HUGE *nit,  double NAG_HUGE wrk[],
	    Integer NAG_HUGE *ifail);
#else
extern void g02hax();
#endif

#ifdef NAG_PROTO
extern void g02hay(Integer indw, Integer n, Integer m, double NAG_HUGE x[], Integer tdx,
            Integer mm, double cucv, double NAG_HUGE wgt[], Integer maxit,
            Integer nitmon, Nag_FileSt NAG_HUGE *stream, double tol, Integer NAG_HUGE *nit,
	    double NAG_HUGE sa[], double NAG_HUGE sz[], double NAG_HUGE sc2[], Integer NAG_HUGE *ifail);
#else
extern void g02hay();
#endif

#ifdef NAG_PROTO
extern void g02haz(Integer indw, Integer ipsi, Integer indc, double sigma,
             Integer n, Integer m,  double NAG_HUGE x[],  Integer tdx,
             double NAG_HUGE rs[], double NAG_HUGE wgt[], double cpsi, double h1,
             double h2, double h3, double NAG_HUGE c[],  Integer tdc,
             double NAG_HUGE d[], double NAG_HUGE e[], double NAG_HUGE wk[],
             Integer NAG_HUGE *ifail);
#else
extern void g02haz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g02hkc(Integer n, Integer m, double NAG_HUGE x[], Integer tdx, 
            double eps, double NAG_HUGE cov[], double NAG_HUGE theta[], Integer max_iter,
            Integer print_iter, const char NAG_HUGE *outfile, double tol, Integer NAG_HUGE *iter, 
            NagError NAG_HUGE *fail);
#else
extern void g02hkc();
#endif

#ifdef NAG_PROTO
extern void g02hku(double t, double NAG_HUGE b[], double NAG_HUGE *u, double NAG_HUGE *ud,
            double NAG_HUGE *w, double NAG_HUGE *wd);
#else
extern void g02hku();
#endif

#ifdef NAG_PROTO
extern double g02hkv(double tau2, double xp, double a2, double b2);
#else
extern double g02hkv();
#endif

#ifdef NAG_PROTO
extern void g02hkw(double a2, double b2,  Integer nvar, double NAG_HUGE *t2);
#else
extern void g02hkw();
#endif

#ifdef NAG_PROTO
extern void g02hkx(double eps, double NAG_HUGE *c);
#else
extern void g02hkx();
#endif

#ifdef NAG_PROTO
extern double g02hky(double cap,  Integer nvar,  double tol);
#else
extern double g02hky();
#endif

#ifdef NAG_PROTO
extern void g02hkz(double eps, Integer nvar, double NAG_HUGE *a2, double NAG_HUGE *b2);
#else
extern void g02hkz();
#endif

#ifdef NAG_PROTO
extern void g02hlc(void (NAG_HUGE *ucv) (double t, double NAG_HUGE userp[], double NAG_HUGE *u, double NAG_HUGE *ud, double NAG_HUGE *w, double NAG_HUGE *wd),  
            double NAG_HUGE userp[],  Integer indm, Integer n,
            Integer m,  double NAG_HUGE x[],  Integer tdx,  double NAG_HUGE cov[],
            double NAG_HUGE a[], double NAG_HUGE wt[], double NAG_HUGE theta[], double bl,
            double bd,  Integer maxit, Integer nitmon,  double tol,
            Integer NAG_HUGE *nit,  double NAG_HUGE wk[],  Nag_FileSt NAG_HUGE *stream, NagError NAG_HUGE *fail);
#else
extern void g02hlc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG02 */
