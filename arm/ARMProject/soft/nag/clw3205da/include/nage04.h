#ifndef NAGE04
#define NAGE04
#ifdef __cplusplus
extern "C"
{
#endif


/* <nage04.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library e04 Chapter
 *
 * Mark 4 re-issue, 1996.
 * Mark 5 revised. IER-2148 (Feb 1998).
 */

/* This header file uses the type FILE,  hence stdio.h has to be 
 * included.
 */
#include <stdio.h>

/* Message codes and extern of message list */
#include <nag_e04mesg.h>

/* Also note that
 * a) e04uc0 is old e04ucj
 * b) e04uc1 is old e04uck
*/

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04abc(NAG_E04ABC_FUN funct, double e1, double e2, double *a,
		       double *b, Integer max_fun, double *x,
		       double *f, Nag_Comm *comm, NagError *fail);
#else
extern void e04abc();
#endif

#ifdef NAG_PROTO
extern void e04abz(double eps, double t, double eta,
            double sftbnd, double xlamda, double NAG_HUGE *u,
            double NAG_HUGE *fu, double gu, double NAG_HUGE *xmin, double NAG_HUGE *fmin,
            double NAG_HUGE *xw, double NAG_HUGE *fw, double NAG_HUGE *xv, double NAG_HUGE *fv,
            double NAG_HUGE *a, double NAG_HUGE *fa, double NAG_HUGE *b, double NAG_HUGE *oldf,
            double NAG_HUGE *b1, double NAG_HUGE *scxbd, double NAG_HUGE *e, double NAG_HUGE *d,
            double NAG_HUGE *rr, double NAG_HUGE *ss, double NAG_HUGE *gtest1, double NAG_HUGE *gtest2,
            double NAG_HUGE *tol, int NAG_HUGE *iloc, int NAG_HUGE *itest);
#else
extern void e04abz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04bbc(NAG_E04BBC_FUN funct, double e1, double e2, double *a,
		       double *b, Integer max_fun, double *x,
		       double *f, double *g, Nag_Comm *comm, NagError *fail);
#else
extern void e04bbc();
#endif

#ifdef NAG_PROTO
extern void e04bbz(double eps, double t, double eta,
		   double xlamda, double NAG_HUGE *u, double NAG_HUGE *fu, 
		   double NAG_HUGE *gu, double NAG_HUGE *xmin, 
		   double NAG_HUGE *fmin, double NAG_HUGE *gmin,
		   double NAG_HUGE *xw, double NAG_HUGE *fw, 
		   double NAG_HUGE *gw, double NAG_HUGE *a,
		   double NAG_HUGE *b, double NAG_HUGE *oldf, 
		   double NAG_HUGE *b1, double NAG_HUGE *scxbd,
		   double NAG_HUGE *e, double NAG_HUGE *d, double NAG_HUGE *rr,
		   double NAG_HUGE *ss, double NAG_HUGE *gtest1, 
		   double NAG_HUGE *gtest2, double NAG_HUGE *tol,
		   int NAG_HUGE *iloc, int NAG_HUGE *itest);
#else
extern void e04bbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04ccc(Integer n, NAG_E04CCC_FUN funct, double NAG_HUGE *x, 
		       double NAG_HUGE *fmin, Nag_E04_Opt NAG_HUGE *options,
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04ccc();
#endif

#ifdef NAG_PROTO
extern void e04ccz(Integer n, double NAG_HUGE xc[], double fmin, double fmax, 
		   double NAG_HUGE simplex[], Nag_E04_Opt NAG_HUGE *opt, 
		   Nag_FileSt NAG_HUGE *stream,
		   Nag_Comm NAG_HUGE *comm);
#else
extern void e04ccz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04dgc(Integer n, NAG_E04DGC_FUN objfun, double NAG_HUGE x[], 
		       double NAG_HUGE *objf, double NAG_HUGE grad[],
		       Nag_E04_Opt NAG_HUGE *options, 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04dgc();
#endif

#ifdef NAG_PROTO
extern void e04dgt(Nag_PrintType print_level, double alfa, double objf,
            double gnorm, double NAG_HUGE objgrd[], double NAG_HUGE x[],
            double xnorm, double xkxknm, Integer n, Integer inform,
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04dgt();
#endif

#ifdef NAG_PROTO
extern void e04dgu(Boolean debug, Integer debug_iter, Integer iter,
            Boolean needfd, NAG_E04DGC_FUN objfun,
            Integer n, double NAG_HUGE *alfbnd, double alfmax,
            double NAG_HUGE *epsaf, double epsrf, double eta,
            double dxnorm, double xnorm, double NAG_HUGE dx[],
            double NAG_HUGE *alfa, double NAG_HUGE *objf, double NAG_HUGE *gdx,
            double NAG_HUGE grad[], double NAG_HUGE x[], Integer NAG_HUGE *inform,
            Integer NAG_HUGE *nfun, double NAG_HUGE ugrad[], double NAG_HUGE x1[],
            Nag_Comm NAG_HUGE *comm);
#else
extern void e04dgu();
#endif

#ifdef NAG_PROTO
extern Boolean e04dgv(double objf, double oldf, double alfa, double pnorm,
               double  rtftol, double cubert, double ftol, double xnorm,
               double gnorm, double epsaf);
#else
extern Boolean e04dgv();
#endif

#ifdef NAG_PROTO
extern double e04dgw(double objf, double fguess, double gtp, double smax);
#else
extern double e04dgw();
#endif

#ifdef NAG_PROTO
extern void e04dgx(Integer n, double gamma_, double NAG_HUGE sj[],
            double NAG_HUGE hjv[], double NAG_HUGE hjyj[], double yjsj,
            double yjhyj, double vsj, double vhyj,
            double NAG_HUGE hjp1v[], Integer iter, Boolean debug,
            Integer debug_iter);
#else
extern void e04dgx();
#endif

#ifdef NAG_PROTO
extern void e04dgy(Integer n, double NAG_HUGE hnew[], double NAG_HUGE hold[],
            double alpha, double NAG_HUGE pk[], double NAG_HUGE yk[],
            double NAG_HUGE g[], Integer iter, Boolean debug,
            Integer debug_iter);
#else
extern void e04dgy();
#endif

#ifdef NAG_PROTO
extern void e04dgz(Integer n, double NAG_HUGE x[], double NAG_HUGE *xnorm,
	    NAG_E04DGC_FUN objfun,
            Integer max_iter,
            double max_line_step, Integer NAG_HUGE *inform, double NAG_HUGE *optim_tol,
            double f_prec, double NAG_HUGE *objf, Integer debug_iter,
            Nag_PrintType print_level,
            double linesearch_tol, Boolean debug, double f_est,
            double NAG_HUGE sk[], double NAG_HUGE yk[], double NAG_HUGE diagb[], double NAG_HUGE olddb[],
            double NAG_HUGE sr[], double NAG_HUGE yr[], double NAG_HUGE oldg[], double NAG_HUGE hg[],
            double NAG_HUGE hyk[], double NAG_HUGE pk[], double NAG_HUGE hyr[],
            double NAG_HUGE wx[], double NAG_HUGE g[], double NAG_HUGE wgrad[],
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04dgz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04fcc(Integer m, Integer n, NAG_E04FCC_FUN lsqfun,
            double NAG_HUGE x[], double NAG_HUGE *fsumsq, double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdj,
            Nag_E04_Opt NAG_HUGE *options, Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04fcc();
#endif

#ifdef NAG_PROTO
extern void e04fch(Integer m, Integer n, double NAG_HUGE xc[],
            double NAG_HUGE fvecc[], double NAG_HUGE fjacc[], Integer tdj,
            double NAG_HUGE g[], double alpha, double xk_norm,
            Nag_E04_Opt opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04fch();
#endif

#ifdef NAG_PROTO
extern void e04fcv(Integer m, Integer n, Nag_Fun NAG_HUGE *fcall,
            double eps, double t,
            double eta, double sftbnd, double xlamda,
            double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
            double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE fjac[],
            Integer tdj, double NAG_HUGE fvec[], double NAG_HUGE g[], Integer NAG_HUGE *numf,
            Integer NAG_HUGE *flag, double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04fcv();
#endif

#ifdef NAG_PROTO
extern void e04fcw(Integer m, Integer n, NAG_E04FCC_FUN lsfun,
            double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdj,
            Nag_Comm NAG_HUGE *comm);
#else
extern void e04fcw();
#endif

#ifdef NAG_PROTO
extern void e04fcx(Integer m, Integer n, NAG_E04FCC_FUN lsfun,
            double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdj,
            Nag_Comm NAG_HUGE *comm);
#else
extern void e04fcx();
#endif

#ifdef NAG_PROTO
extern void e04fcy(Integer m, Integer n, Integer ns, double epsmch,
            Nag_Fun NAG_HUGE *fcall, Integer grade, double NAG_HUGE x[],
            double NAG_HUGE fvec[], double ff, Boolean p1zero,
            double NAG_HUGE v[], Integer tdv, double NAG_HUGE p[], double NAG_HUGE phesl[],
            double NAG_HUGE phesd[], double NAG_HUGE rhs[],
            double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04fcy();
#endif

#ifdef NAG_PROTO
extern void e04fcz(Integer m, Integer n, Nag_Fun NAG_HUGE *fcall,
            Nag_PrintType print_level, Integer max_iter, double eta,
            double xtol, double stepmx, double NAG_HUGE x[],
            double NAG_HUGE *fsumsq, double NAG_HUGE fvec[], double NAG_HUGE fjac[],
            Integer tdj, double NAG_HUGE s[], double NAG_HUGE vt[], Integer tdvt,
            Integer NAG_HUGE *niter, Integer NAG_HUGE *nftotl, int NAG_HUGE *error,
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04fcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04gbc(Integer m, Integer n, NAG_E04GBC_FUN lsqfun,
            double NAG_HUGE x[], double NAG_HUGE *fsumsq, double NAG_HUGE fvec[], double NAG_HUGE fjac[],
            Integer tdj, Nag_E04_Opt NAG_HUGE *options, Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04gbc();
#endif

#ifdef NAG_PROTO
extern void e04gbx(Integer m, Integer n, NAG_E04GBC_FUN lsfun,
            double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdj,
            Nag_Comm NAG_HUGE *comm);
#else
extern void e04gbx();
#endif

#ifdef NAG_PROTO
extern void e04gbz(Integer m, Integer n, double NAG_HUGE h[],
            double alpha, double NAG_HUGE p[], double NAG_HUGE gplus[],
            double NAG_HUGE g[], double NAG_HUGE fjac[], Integer tdj, double NAG_HUGE w[]);
#else
extern void e04gbz();
#endif

#ifdef NAG_PROTO
extern void e04gdq(Integer m, Integer n,
            void (NAG_HUGE *lsqlin)(Integer m, Integer n, Nag_Fun NAG_HUGE *fcall,
                           double eps, double t, double eta, double sftbnd,
                           double xlamda, double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
                           double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE fjac[], Integer tdj,
                           double NAG_HUGE fvec[], double NAG_HUGE g[], Integer NAG_HUGE *nftotl,
                           Integer NAG_HUGE *iflag, double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm),
            Nag_Fun NAG_HUGE *fcall, double eta, Integer grade, Nag_PrintType print_level,
            double stepmx, double epsmch, double NAG_HUGE *tau, double xnorm,
            Boolean NAG_HUGE *nomove, double NAG_HUGE x[], double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdj, double NAG_HUGE g[], double NAG_HUGE p[],
            double NAG_HUGE *alpha, double NAG_HUGE *pnorm,
            double NAG_HUGE *ssqnew, double NAG_HUGE *ssqold, Integer NAG_HUGE *niter,
            Integer NAG_HUGE *nftotl, Integer NAG_HUGE *nwhy, double NAG_HUGE w[],
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04gdq();
#endif

#ifdef NAG_PROTO
extern void e04gdr(Integer m, Integer n, double NAG_HUGE a[], Integer tda,
            double NAG_HUGE x[], double NAG_HUGE y[]);
#else
extern void e04gdr();
#endif

#ifdef NAG_PROTO
extern void e04gds(Integer ns, Integer grade, double epsmch, double NAG_HUGE utf[],
            double NAG_HUGE s[], double NAG_HUGE phesl[], double NAG_HUGE phesd[],
            double NAG_HUGE rhs[], Integer NAG_HUGE *nphi);
#else
extern void e04gds();
#endif

#ifdef NAG_PROTO
extern void e04gdu(double NAG_HUGE s[], Integer NAG_HUGE *rank);
#else
extern void e04gdu();
#endif

#ifdef NAG_PROTO
extern void e04gdv(Integer n, double NAG_HUGE a[], Integer tda);
#else
extern void e04gdv();
#endif

#ifdef NAG_PROTO
extern void e04gdw(Integer maxrnk, double NAG_HUGE s[], Integer NAG_HUGE *rank);
#else
extern void e04gdw();
#endif

#ifdef NAG_PROTO
extern void e04gdx(int NAG_HUGE *flag, Integer m, Integer n, double NAG_HUGE fjac[],
            Integer tdj, double NAG_HUGE vt[], Integer lvt, double NAG_HUGE fvec[],
            double NAG_HUGE s[], Boolean wantvc, double NAG_HUGE utf[],
            double NAG_HUGE w[]);
#else
extern void e04gdx();
#endif

#ifdef NAG_PROTO
extern void e04gdy(Integer m, Integer n, Integer NAG_HUGE *maxrnk, Integer NAG_HUGE *grade,
            Integer NAG_HUGE *nwhy, double tau, double epsmch,
            double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdj,
            double NAG_HUGE s[], double NAG_HUGE vt[], Integer tdv, double NAG_HUGE utf[],
            Boolean NAG_HUGE *gauss, Boolean nomove, double NAG_HUGE *w);
#else
extern void e04gdy();
#endif

#ifdef NAG_PROTO
extern void e04gdz(Integer m, Integer n, Nag_Fun NAG_HUGE *fcall,
            Integer max_iter, double eta,
            double xtol, double stepmx, double NAG_HUGE x[],
            double NAG_HUGE *ssqnew, double NAG_HUGE fvec[], double NAG_HUGE fjac[],
            Integer tdj, Integer lvt, Integer NAG_HUGE *niter, Integer NAG_HUGE *nftotl,
            Integer NAG_HUGE *nwhy, double NAG_HUGE *rteps, double NAG_HUGE *rtol,
            double NAG_HUGE *tol, Integer NAG_HUGE *grade,
            double NAG_HUGE *peps, double NAG_HUGE g[], double NAG_HUGE *ssqold,
            Boolean NAG_HUGE *gauss, double NAG_HUGE *tau, double NAG_HUGE *alpha, double NAG_HUGE *pnorm,
            double NAG_HUGE *epsmch, Boolean NAG_HUGE *nomove, Boolean deriv_check, Nag_Comm NAG_HUGE *comm);
#else
extern void e04gdz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04hbc(Integer n, NAG_E04HBC_FUN sfun, double NAG_HUGE x[], 
		       double NAG_HUGE delta[], double NAG_HUGE hesd[],
		       double NAG_HUGE *f, double NAG_HUGE g[], 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04hbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04hcc(Integer n, NAG_E04HCC_FUN sfun, double NAG_HUGE x[], 
		       double NAG_HUGE *f, double NAG_HUGE g[], 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04hcc();
#endif

#ifdef NAG_PROTO
extern void e04hcz(Integer n, double NAG_HUGE y[], double NAG_HUGE z[]);
#else
extern void e04hcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04hdc(Integer n, NAG_E04LBC_FUN sfun, NAG_E04LBC_HESS shess,
		       double x[], double g[], double hesl[], double hesd[], 
		       Nag_Comm *user_comm, NagError *fail);
#else
extern void e04hdc();
#endif

#ifdef NAG_PROTO
extern void e04hdz(Integer n, Integer lh,  double hesl[], double hesd[],
             double p[], double w[], double *y);
#else
extern void e04hdz();
#endif

#ifdef NAG_PROTO
extern void e04hev(Integer m, Integer n, Nag_Fun NAG_HUGE *fcall,
            double eps, double t, double eta, double sftbnd,
            double xlamda, double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
            double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE fjac[],
            Integer tdj, double NAG_HUGE fvec[], double NAG_HUGE g[], Integer NAG_HUGE *numf,
            Integer NAG_HUGE *flag, double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04hev();
#endif

#ifdef NAG_PROTO
extern void e04hey(Integer n, Integer ns, Integer grade,
            double NAG_HUGE v[], Integer tdv, double NAG_HUGE p[],
            double NAG_HUGE phesl[], double NAG_HUGE phesd[], double NAG_HUGE prhs[],
            double NAG_HUGE h[], double NAG_HUGE w[]);
#else
extern void e04hey();
#endif

#ifdef NAG_PROTO
extern void e04hez(Integer m, Integer n, double NAG_HUGE fvec[], double NAG_HUGE fjac[],
            Integer tdj, double NAG_HUGE g[]);
#else
extern void e04hez();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04jbc(Integer n, NAG_E04JBC_FUN sfun, Nag_BoundType bound, 
		       double NAG_HUGE bl[], double NAG_HUGE bu[],
		       double NAG_HUGE x[], double NAG_HUGE *f, 
		       double NAG_HUGE g[], Nag_E04_Opt NAG_HUGE *options, 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04jbc();
#endif

#ifdef NAG_PROTO
extern void e04jbg(Integer n, NAG_E04JBC_FUN sfun,
            Nag_PrintType print_level, Boolean locsch, int intyp,
            void (NAG_HUGE *minlin)(Integer n, NAG_E04JBC_FUN sfun,
                           double eps, double t,
                           double eta, double sftbnd, double xlamda,
                           double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
                           double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
                           Integer NAG_HUGE *nftotl, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm),
            Integer max_iter, double eta,
            double xtol, double stepmx, double fest,
            double NAG_HUGE delta[], int ibound, double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE x[], double NAG_HUGE hesl[],
            double NAG_HUGE hesd[], Integer NAG_HUGE istate[], double NAG_HUGE *f,
            double NAG_HUGE g[], Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm,
            char NAG_HUGE error_buf[], int NAG_HUGE *error);
#else
extern void e04jbg();
#endif

#ifdef NAG_PROTO
extern void e04jbh(Integer n, double NAG_HUGE xc[], double f,
            double NAG_HUGE g[], double gpjnorm, double cond,
            double alpha, double xk_norm,
            Nag_E04_Opt opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04jbh();
#endif

#ifdef NAG_PROTO
extern void e04jbm(Integer n, Integer NAG_HUGE *nfree, int ibound,
            double NAG_HUGE delta[], double rtolsq, double toleps,
            double endeps, double stepmx, double sfirst,
	    NAG_E04JBC_FUN sfun, 
            void (NAG_HUGE *minlin)(Integer n, NAG_E04JBC_FUN sfun,
                           double eps, double t,
                           double eta, double sftbnd, double xlamda,
                           double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
                           double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
                           Integer NAG_HUGE *nftotl, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm),
            double NAG_HUGE x[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE *f,
            Integer NAG_HUGE istate[], double NAG_HUGE g[], double NAG_HUGE *gtg,
            Integer NAG_HUGE *itype, double NAG_HUGE gpj[], double NAG_HUGE *alpha, double NAG_HUGE p[],
            double NAG_HUGE hesl[], Integer NAG_HUGE *lhproj,
            double NAG_HUGE hesd[], double NAG_HUGE *boundk, Integer NAG_HUGE *nftotl,
            Integer NAG_HUGE *lcount, int NAG_HUGE *error, Boolean NAG_HUGE *finish,
            double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04jbm();
#endif

#ifdef NAG_PROTO
extern void e04jbn(Integer n, Integer NAG_HUGE istate[], double NAG_HUGE *f,
            double NAG_HUGE x[], double NAG_HUGE delta[],
	    NAG_E04JBC_FUN sfun,
            Integer NAG_HUGE *itype, Integer NAG_HUGE *nftotl, double NAG_HUGE g[],
            double NAG_HUGE fplus[], double NAG_HUGE fminus[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04jbn();
#endif

#ifdef NAG_PROTO
extern void e04jbp(Integer n, double NAG_HUGE x[], double NAG_HUGE delta[],
            double NAG_HUGE *delmin, Integer NAG_HUGE *flag);
#else
extern void e04jbp();
#endif

#ifdef NAG_PROTO
extern void e04jbq(Integer n, NAG_E04JBC_FUN sfun,
            double eps, double t,
            double eta, double sftbnd, double xlamda,
            double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
            double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
            Integer NAG_HUGE *nftotl, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm);
#else
extern void e04jbq();
#endif

#ifdef NAG_PROTO
extern void e04jbr(Integer nfree, Integer jfix, double NAG_HUGE *gtp,
            double NAG_HUGE pa[], double NAG_HUGE gpjold[]);
#else
extern void e04jbr();
#endif

#ifdef NAG_PROTO
extern void e04jbs(Integer n, double tol, Integer NAG_HUGE *nfree,
            Integer NAG_HUGE istate[], double NAG_HUGE x[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE p[], Integer NAG_HUGE *jfix, double NAG_HUGE *boundk,
            Boolean NAG_HUGE *recalc);
#else
extern void e04jbs();
#endif

#ifdef NAG_PROTO
extern double e04jbt(double fnew, double fm, double gtp, double smax);
#else
extern double e04jbt();
#endif

#ifdef NAG_PROTO
extern void e04jbu(Integer nfrnew, Integer jfix, double NAG_HUGE hesl[],
            double NAG_HUGE hesd[], Integer NAG_HUGE *lhproj,
            double NAG_HUGE w[], Integer NAG_HUGE *iflag);
#else
extern void e04jbu();
#endif

#ifdef NAG_PROTO
extern void e04jbv(Integer n, Integer NAG_HUGE *nfree, double NAG_HUGE *stepmx,
            double test, double tol, Integer NAG_HUGE istate[],
            double NAG_HUGE x[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE p[],
            Integer NAG_HUGE *jfix, double NAG_HUGE *boundk);
#else
extern void e04jbv();
#endif

#ifdef NAG_PROTO
extern double e04jbw(Integer lhd);
#else
extern double e04jbw();
#endif

#ifdef NAG_PROTO
extern void e04jbx(Integer lhd, double NAG_HUGE hesd[], Integer n,
            double NAG_HUGE hesl[], double NAG_HUGE *cond);
#else
extern void e04jbx();
#endif

#ifdef NAG_PROTO
extern void e04jby(Integer n, Integer nfree, Integer NAG_HUGE istate[], double NAG_HUGE g[],
            double NAG_HUGE gproj[], double NAG_HUGE *gtg);
#else
extern void e04jby();
#endif

#ifdef NAG_PROTO
extern void e04jbz(Integer n, double eps, double NAG_HUGE x[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], Integer NAG_HUGE *nfree, Integer NAG_HUGE *nequal,
            Integer NAG_HUGE istate[], int NAG_HUGE *intype);
#else
extern void e04jbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04kbc(Integer n, NAG_E04KBC_FUN sfun, Nag_BoundType bound, 
		       double NAG_HUGE bl[], double NAG_HUGE bu[],
		       double NAG_HUGE x[], double NAG_HUGE *f, 
		       double NAG_HUGE g[], Nag_E04_Opt NAG_HUGE *options, 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04kbc();
#endif

#ifdef NAG_PROTO
extern void e04kbg(Integer n, NAG_E04KBC_FUN sfun,
            Nag_PrintType print_level, Boolean locsch, Integer intyp,
            void (NAG_HUGE *minlin)(Integer n, NAG_E04KBC_FUN sfun,
                           double eps, double t,
                           double eta, double sftbnd, double xlamda,
                           double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
                           double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
                           Integer NAG_HUGE *numf, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm),
            Integer max_iter, double eta,
            double xtol, double stepmx, double fest,
            int ibound, double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE x[], double NAG_HUGE hesl[],
            double NAG_HUGE hesd[], Integer NAG_HUGE istate[], double NAG_HUGE *f,
            double NAG_HUGE g[], Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm, int NAG_HUGE *error);
#else
extern void e04kbg();
#endif

#ifdef NAG_PROTO
extern void e04kbp(Integer n, Integer NAG_HUGE *nfree, double fm,
            double rtolsq, double toleps, double endeps,
            double stepmx, double sfirst, NAG_E04KBC_FUN sfun,
            void (NAG_HUGE *minlin)(Integer n, NAG_E04KBC_FUN sfun,
                           double eps, double t,
                           double eta, double sftbnd, double xlamda,
                           double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
                           double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
                           Integer NAG_HUGE *nftotl, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm),
            double NAG_HUGE x[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE *f, Integer NAG_HUGE istate[],
            double NAG_HUGE g[], double NAG_HUGE *gtg, double NAG_HUGE gpj[],
            double NAG_HUGE p[], double NAG_HUGE hesl[], Integer NAG_HUGE *lhproj,
            double NAG_HUGE hesd[], double NAG_HUGE *boundk, Integer NAG_HUGE *nftotl,
            Integer NAG_HUGE *lcount, int NAG_HUGE *error, Boolean NAG_HUGE *finish,
            double NAG_HUGE w[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04kbp();
#endif

#ifdef NAG_PROTO
extern void e04kbq(Integer n, double NAG_HUGE z[], double NAG_HUGE hesd[],
            double NAG_HUGE hesl[], double NAG_HUGE p[]);
#else
extern void e04kbq();
#endif

#ifdef NAG_PROTO
extern void e04kbr(Integer n, double alpha, double NAG_HUGE z[],
            double NAG_HUGE hesd[], double NAG_HUGE hesl[], Integer NAG_HUGE *flag);
#else
extern void e04kbr();
#endif

#ifdef NAG_PROTO
extern void e04kbs(Integer n, double NAG_HUGE z[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            Integer NAG_HUGE istate[], double NAG_HUGE p[], double NAG_HUGE ztemp[],
            double NAG_HUGE *spos, double NAG_HUGE *sneg);
#else
extern void e04kbs();
#endif

#ifdef NAG_PROTO
extern void e04kbt(Integer n, Integer nfree, double NAG_HUGE *step,
            double NAG_HUGE bl[], double NAG_HUGE bu[], Integer NAG_HUGE istate[],
            double NAG_HUGE y[], double NAG_HUGE p[], double NAG_HUGE *spos,
            double NAG_HUGE *sneg, Integer NAG_HUGE *flag);
#else
extern void e04kbt();
#endif

#ifdef NAG_PROTO
extern void e04kbu(Integer n, Integer nfree,
            double eps, double NAG_HUGE x[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            Integer NAG_HUGE istate[], NAG_E04JBC_FUN sfun,
            double NAG_HUGE *step, double NAG_HUGE *fnew, Integer NAG_HUGE *nfun, double NAG_HUGE p[], double NAG_HUGE y[],
            double NAG_HUGE gy[], Integer NAG_HUGE *is, Integer NAG_HUGE *ifail, Nag_Comm NAG_HUGE *comm);
#else
extern void e04kbu();
#endif

#ifdef NAG_PROTO
extern void e04kbv(Integer inew, double NAG_HUGE g[], Integer NAG_HUGE *nfree,
            Integer NAG_HUGE istate[], double NAG_HUGE gproj[], double NAG_HUGE *gtg,
            double NAG_HUGE hesd[], Integer NAG_HUGE *lhproj, double NAG_HUGE *boundk,
            double NAG_HUGE *cond);
#else
extern void e04kbv();
#endif

#ifdef NAG_PROTO
extern void e04kbw(Integer n, Integer NAG_HUGE *numneg, Integer NAG_HUGE *nftotl,
	    NAG_E04JBC_FUN sfun,
            double f, double NAG_HUGE x[], Integer NAG_HUGE istate[],
            double NAG_HUGE gtemp[], Integer NAG_HUGE *izero, Nag_Comm NAG_HUGE *comm);
#else
extern void e04kbw();
#endif

#ifdef NAG_PROTO
extern void e04kbx(Integer n, Integer nfree, double tol,
            Integer NAG_HUGE istate[], double NAG_HUGE g[], Boolean NAG_HUGE *lm1pos,
            Integer NAG_HUGE *ineglm, Integer NAG_HUGE *numneg);
#else
extern void e04kbx();
#endif

#ifdef NAG_PROTO
extern void e04kby(Integer m, double step, double goldtp,
            double NAG_HUGE gold[], double gnewtp, double NAG_HUGE gnew[],
            double NAG_HUGE y[], double NAG_HUGE z[], double NAG_HUGE hesd[], double NAG_HUGE hesl[],
            Integer NAG_HUGE *flag);
#else
extern void e04kby();
#endif

#ifdef NAG_PROTO
extern void e04kbz(Integer n, double boundk, double NAG_HUGE hesd[], double NAG_HUGE *cond);
#else
extern void e04kbz();
#endif

#ifdef NAG_PROTO
extern void e04lbb(Integer n, Integer lel, Integer is,  double el[],
double p[]);
#else
extern void e04lbb();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04lbc(Integer n, NAG_E04LBC_FUN sfun, NAG_E04LBC_HESS shess,
		       Nag_BoundType bound, double bl[], double bu[], 
		       double x[], double *f, double g[], Nag_E04_Opt *options,
		       Nag_Comm *user_comm, NagError *fail);
#else
extern void e04lbc();
#endif

#ifdef NAG_PROTO
extern void e04lbd(Integer n, Integer lel,  double el[], double b[],
double y[]);
#else
extern void e04lbd();
#endif

#ifdef NAG_PROTO
extern void e04lbe(Integer n, double xc[], double f,
            double g[], double gpjnorm, double cond, Boolean posdef,
	    double alpha, double xk_norm,
            Nag_E04_Opt opt, Nag_FileSt *stream, Nag_Comm *comm);
#else
extern void e04lbe();
#endif

#ifdef NAG_PROTO
extern void e04lbg(Integer n, Integer *nfree, Integer *numneg,
	    Integer *nftotl,  Boolean *nomove, Boolean userh,
	    NAG_E04LBC_HESS shess,
	    double eps, double delta,
	    double bl[], double bu[],
	    NAG_E04LBC_FUN sfun,
	    double *f, double x[],  Integer istate[],
	    double gfull[], double gproj[], double hdfull[],
	    double hdproj[], double hlfull[],  Integer nh,
	    double hlproj[],  Integer lh, Integer *lhproj, Integer *iq,
	    Integer *izero,  Boolean *tight, Nag_Comm *comm);
#else
extern void e04lbg();
#endif

#ifdef NAG_PROTO
extern void e04lbh(Integer n, Integer inew,  Boolean *tight,  Integer istate[],
	     Integer *nfree);
#else
extern void e04lbh();
#endif

#ifdef NAG_PROTO
extern void e04lbj(Integer n, Integer nfree,  double tol,  Boolean lsfail,
	     Boolean tight,  Integer istate[],  double g[],
	     double p[], double hesl[],  Integer lh,
	     Boolean *lm2pos,  Integer *ineglm, Integer *numneg);
#else
extern void e04lbj();
#endif

#ifdef NAG_PROTO
extern void e04lbk(Integer n, Integer klm,  double *total, double bndsgn,
	     Integer istate[],  double hesl[],  Integer lh,
	     double p[]);
#else
extern void e04lbk();
#endif

#ifdef NAG_PROTO
extern void e04lbl(Integer n, Integer *nfree,  double tol,  Integer istate[],
	     double x[], double bl[], double bu[], double p[],
	     Boolean *recalc);
#else
extern void e04lbl();
#endif

#ifdef NAG_PROTO
extern void e04lbm(Integer n,  double *stepmx, double tol, double x[],
	     double bl[], double bu[], double p[]);
#else
extern void e04lbm();
#endif

#ifdef NAG_PROTO
extern void e04lbn(Integer n, Integer nfree, Integer inew, Integer istate[],
	     double gfull[], double hdfull[], double hlfull[],
	     Integer nh,  double gproj[], double hdproj[],
	     double hlproj[],  Integer lh,  double *delta,
	     Integer *nsphi, Integer *lhproj,  double w[]);
#else
extern void e04lbn();
#endif

#ifdef NAG_PROTO
extern void e04lbp(Integer n, Integer lh,  double x[],
	    double g[],  Integer istate[],
	    NAG_E04LBC_FUN sfun,
	    double delta, double hesl[], double hesd[], Nag_Comm *comm);
#else
extern void e04lbp();
#endif

#ifdef NAG_PROTO
extern void e04lbq(Integer n, Integer lh,  double x[],
	    double g[],  Integer istate[], Integer inew,
	    NAG_E04LBC_FUN sfun,
	    double delta, double hesl[],
	    double hesd[], Nag_Comm *comm);
#else
extern void e04lbq();
#endif

#ifdef NAG_PROTO
extern void e04lbr(Integer n,
	    NAG_E04LBC_FUN sfun,
	    NAG_E04LBC_HESS shess,
	    Boolean userh,
	    Nag_PrintType print_level,


	    void (*minlin) (Integer n,
	    NAG_E04LBC_FUN sfun,
			    double eps, double t,
			    double eta, double sftbnd, double xlamda,
			    double p[], double gtp, double x[], double *f,
			    double *alpha, double g[],  Integer *nftotl,
			    Integer *iflag, Nag_Comm *comm),
	    Integer maxiter,  double eta, double xtol,
	    double delta, double stepmx,  Integer ibound,
	    double bl[], double bu[], double x[],
	    double hesl[],  Integer lh,  double hesd[],
	    Integer istate[],  double *f, double g[], Nag_E04_Opt *opt,
	    Nag_FileSt *stream, Nag_Comm *comm,
	    Integer *ifail);
#else
extern void e04lbr();
#endif

#ifdef NAG_PROTO
extern void e04lbs(Integer n, NAG_E04KBC_FUN sfun,
            double eps, double t,
            double eta, double sftbnd, double xlamda,
            double NAG_HUGE p[], double gtp, double NAG_HUGE x[],
            double NAG_HUGE *f, double NAG_HUGE *alpha, double NAG_HUGE g[],
            Integer NAG_HUGE *nftotl, Integer NAG_HUGE *iflag, Nag_Comm NAG_HUGE *comm);
#else
extern void e04lbs();
#endif

#ifdef NAG_PROTO
extern void e04lbt(Integer n, Integer nfree, double tol,
            Integer NAG_HUGE istate[], double NAG_HUGE g[], double rneglm,
            Integer NAG_HUGE *ineglm);
#else
extern void e04lbt();
#endif

#ifdef NAG_PROTO
extern void e04lbu(Integer n, Integer NAG_HUGE *numneg, Integer NAG_HUGE istate[]);
#else
extern void e04lbu();
#endif

#ifdef NAG_PROTO
extern void e04lbv(Integer n, Integer NAG_HUGE istate[], double NAG_HUGE pa[],
            double signum, double NAG_HUGE p[]);
#else
extern void e04lbv();
#endif

#ifdef NAG_PROTO
extern void e04lbw(double NAG_HUGE hesd[], Integer lhd, double boundk,
            Integer NAG_HUGE *flag, double NAG_HUGE *cond);
#else
extern void e04lbw();
#endif

#ifdef NAG_PROTO
extern void e04lbx(Integer n, Integer NAG_HUGE istate[], double NAG_HUGE gfull[],
            double NAG_HUGE hdfull[], double NAG_HUGE hlfull[], double NAG_HUGE gproj[],
            double NAG_HUGE hdproj[], double NAG_HUGE hlproj[], Integer NAG_HUGE *lhproj);
#else
extern void e04lbx();
#endif

#ifdef NAG_PROTO
extern void e04lby(Integer n, double eps, double NAG_HUGE x[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], Integer NAG_HUGE *nfree, Integer NAG_HUGE istate[],
            Integer NAG_HUGE *nequal);
#else
extern void e04lby();
#endif

#ifdef NAG_PROTO
extern void e04lbz(Integer n, int ibound, double NAG_HUGE *bndmax,
            double NAG_HUGE bl[], double NAG_HUGE bu[]);
#else
extern void e04lbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04mfc(Integer n, Integer nclin, double NAG_HUGE a[], 
		       Integer tda, double NAG_HUGE bl[], double NAG_HUGE bu[],
		       double NAG_HUGE cvec[], double NAG_HUGE x[],
		       double NAG_HUGE *obj, Nag_E04_Opt NAG_HUGE *options, 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04mfc();
#endif

#ifdef NAG_PROTO
extern void e04mfg(char NAG_HUGE *prbtyp, Integer n, Integer nclin, double NAG_HUGE a[], Integer tda,
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE cvec[], double NAG_HUGE ax[],
            double NAG_HUGE x[], double NAG_HUGE *obj, double NAG_HUGE lambda[], Integer NAG_HUGE istate[],
            Boolean vertex, Boolean cset, Nag_ae04mf NAG_HUGE *ae04mf, double NAG_HUGE w[],
            Integer NAG_HUGE iw[], Nag_ee04mf NAG_HUGE *ee04mf, Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm,
            Nag_FileSt NAG_HUGE *stream, NagError NAG_HUGE *ovflow, Nag_EndState NAG_HUGE *endstate);
#else
extern void e04mfg();
#endif

#ifdef NAG_PROTO
extern void e04mfh(Integer n, Integer nclin, Integer tda, Integer NAG_HUGE istate[],
            double bigbnd, Integer NAG_HUGE *numinf, double NAG_HUGE *suminf,
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE a[], double NAG_HUGE featol[],
            double NAG_HUGE cvec[], double NAG_HUGE x[], double NAG_HUGE wtinf[]);
#else
extern void e04mfh();
#endif

#ifdef NAG_PROTO
extern void e04mfj(Boolean NAG_HUGE *rowerr, Boolean unitq, Integer nclin,
            Integer nactiv, Integer nfree, Integer nz, Integer n,
            Integer tdq, Integer tda, Integer tdt, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], Integer NAG_HUGE *jmax,
            double NAG_HUGE *errmax, double NAG_HUGE *xnorm, double NAG_HUGE a[],
            double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE featol[],
            double NAG_HUGE t[], double NAG_HUGE x[], double NAG_HUGE q[], double NAG_HUGE p[],
            double NAG_HUGE work[]);
#else
extern void e04mfj();
#endif

#ifdef NAG_PROTO
extern void e04mfk(Nag_EndState endstate, char NAG_HUGE *prbtyp,
            double obj, Integer numinf,
            Integer nfree, Integer n,
            Integer nclin, Integer ncnln, Integer nctotl,
            double bigbnd, Integer nactiv, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE clamda[],
            double NAG_HUGE rlamda[], double NAG_HUGE x[],
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04mfk();
#endif

#ifdef NAG_PROTO
extern void e04mfl(Integer nrz, Integer nz,
            double zerolm, Integer NAG_HUGE *notopt, Integer numinf,
            double NAG_HUGE *trusml, double NAG_HUGE *smllst, Integer NAG_HUGE *jsmlst,
            double NAG_HUGE *tinyst, Integer NAG_HUGE *jtiny, double NAG_HUGE gq[], Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04mfl();
#endif

#ifdef NAG_PROTO
extern void e04mfm(Integer n, Integer tda,
            Integer tdt, Integer nactiv, Integer nfree, Integer nz,
            Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[],
            double zerolm, Integer NAG_HUGE *notopt, Integer numinf,
            double NAG_HUGE *trusml, double NAG_HUGE *smllst, Integer NAG_HUGE *jsmlst,
            Integer NAG_HUGE *ksmlst, double NAG_HUGE *tinyst, Integer NAG_HUGE *jtiny,
            Integer jinf, double NAG_HUGE *trubig, double NAG_HUGE *biggst,
            Integer NAG_HUGE *jbigst, Integer NAG_HUGE *kbigst, double NAG_HUGE a[],
            double NAG_HUGE anorms[], double NAG_HUGE gq[], double NAG_HUGE rlamda[],
            double NAG_HUGE t[], double NAG_HUGE wtinf[], Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04mfm();
#endif

#ifdef NAG_PROTO
extern void e04mfn(char NAG_HUGE *subr, char NAG_HUGE *msg, double NAG_HUGE v[], Integer lenv);
#else
extern void e04mfn();
#endif

#ifdef NAG_PROTO
extern void e04mfq(Integer n, Integer nclin, Integer NAG_HUGE istate[],
            double bigbnd, Integer NAG_HUGE *nviol, Integer NAG_HUGE *jmax,
            double NAG_HUGE *errmax, double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE featol[], double NAG_HUGE x[], Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04mfq();
#endif

#ifdef NAG_PROTO
extern void e04mfr(char NAG_HUGE *job, Integer n, Integer nclin,
            Integer NAG_HUGE *nmoved, Integer iter, Integer numinf,
            Integer NAG_HUGE istate[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE featol[],
            double NAG_HUGE featlu[], double NAG_HUGE x[], Nag_ce04mf NAG_HUGE *ce04mf);
#else
extern void e04mfr();
#endif

#ifdef NAG_PROTO
extern void e04mfs(Boolean firstv, Integer n, Integer nclin, Integer NAG_HUGE istate[],
            double bigalf, double bigbnd, double pnorm,
            Boolean NAG_HUGE *hitlow, Boolean NAG_HUGE *move, Boolean NAG_HUGE *onbnd, Boolean NAG_HUGE *unbndd,
            double NAG_HUGE *alfa, double NAG_HUGE *alfap, Integer NAG_HUGE *jhit,
            double NAG_HUGE anorm[], double NAG_HUGE ap[], double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE featol[], double NAG_HUGE featlu[],
            double NAG_HUGE p[], double NAG_HUGE x[], Nag_ce04mf NAG_HUGE *ce04mf, Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04mfs();
#endif

#ifdef NAG_PROTO
extern void e04mft(Nag_Start start, Boolean vertex, Integer nclin,
            Integer nctotl, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nartif,
            Integer NAG_HUGE *nfree, Integer n, Integer tda, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], double bigbnd,
            double tolact, double NAG_HUGE a[], double NAG_HUGE ax[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE featol[], double NAG_HUGE x[],
            double NAG_HUGE wx[], double NAG_HUGE work[], Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04mft();
#endif

#ifdef NAG_PROTO
extern void e04mfv(Boolean cset, Integer n, Integer nclin, Nag_ae04mf NAG_HUGE *ae04mf);
#else
extern void e04mfv();
#endif

#ifdef NAG_PROTO
extern void e04mfy(Integer nrz, Integer tdr, double NAG_HUGE r[], double rzz);
#else
extern void e04mfy();
#endif

#ifdef NAG_PROTO
extern void e04mfz(char NAG_HUGE *prbtyp, char NAG_HUGE *msg, Boolean cset,
            Boolean rset, Boolean NAG_HUGE *unitq,
            Integer itmax, Integer jinf, Integer NAG_HUGE *nviol, Integer n,
            Integer nclin, Integer lda, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree,
            Integer NAG_HUGE *nrz, Integer NAG_HUGE *nz, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],
            void (NAG_HUGE *print_iter)(char NAG_HUGE *prbtyp, Boolean NAG_HUGE *header, Boolean rset,
                               Nag_PrintType print_level, Integer isdel, Integer NAG_HUGE *jdel,
                               Integer jadd, Integer n, Integer nclin,
                               Integer nactiv, Integer nfree, Integer nz,
                               Integer nrz, Integer tdr, Integer tdt, Integer NAG_HUGE istate[],
                               double alfa, double condrz, double condt, double drzz,
                               double grznrm, Integer numinf,
                               double suminf, Integer notopt, double objqp,
                               double trusml, double NAG_HUGE ax[], double NAG_HUGE r[], double NAG_HUGE t[],
                               double NAG_HUGE x[], Integer NAG_HUGE kx[], Integer NAG_HUGE kactive[],
                               double NAG_HUGE rlambda[], double NAG_HUGE gq[], double NAG_HUGE work[],
                               Nag_QP_Print NAG_HUGE *print_inf, Nag_E04_Opt NAG_HUGE *opt,
                               Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm),
            double NAG_HUGE *obj, Integer NAG_HUGE *numinf, double NAG_HUGE *xnorm, double NAG_HUGE a[],
            double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE cvec[],
            double NAG_HUGE featol[], double NAG_HUGE featlu[], double NAG_HUGE x[], Nag_E04_Opt NAG_HUGE *opt,
            double NAG_HUGE w[], Nag_de04nb NAG_HUGE *de04nb, Nag_ae04mf NAG_HUGE *ae04mf,
            Nag_ce04mf NAG_HUGE *ce04mf, Nag_de04mf NAG_HUGE *de04mf, Nag_QP_Print NAG_HUGE *print_inf,
            Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm, Nag_ee04mf NAG_HUGE *ee04mf, NagError NAG_HUGE *ovflow);
#else
extern void e04mfz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04myc(double **a, Integer **ha, Integer **ka, 
		       double **bl, double **bu, double **xs);
#else
extern void e04myc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04mzc(char *mps_file, Integer *n, Integer *m, Integer *nnz, 
		       Integer *iobj, double **a, Integer **ha, Integer **ka, 
		       double **bl, double **bu, double **xs,
		       Nag_E04_Opt *options, NagError *fail);
#else
extern void e04mzc();
#endif

#ifdef NAG_PROTO
extern void e04mzg(NAG_E04NKC_HESSFUN_WRAP hx, NAG_E04NKC_HESSFUN hx1, 
		   Integer ngqp, Integer ngobj,
		   Integer ngobj0, Integer nstate,  double *objqp,
		   Integer nka,  double a[],  Integer ha[], Integer ka[],
		   double gobj[], double gobjqp[], double x0[],
		   double x[], double xdif[],  Integer iz[], Integer leniz,
		   double z[],  Integer lenz, Integer iparm[], double rparm[],
		   Nag_Comm *comm);
#else
extern void e04mzg();
#endif

#ifdef NAG_PROTO
extern void e04mzh(NAG_E04NKC_HESSFUN_WRAP hx, NAG_E04NKC_HESSFUN hx1, Boolean gotr, 
	    Boolean incres, Boolean *needpi, Boolean *posdef, Integer lenr, 
	    Integer m, Integer mbs, Integer maxs, Integer n, Integer nb, 
	    Integer ngqp0, Integer ngqp, Integer nnh, Integer *ns, 
	    Integer *jqsave, Integer newsb, Integer *nuncon, double *obj,
	    double *objqp, double featol, double *step, Integer ne, 
	    Integer nka, double a[], Integer ha[], Integer ka[], 
	    Integer hfeas[], Integer hs[], Integer kbs[], double bl[], 
	    double bu[], double blbs[], double bubs[], double gbs[], 
	    double gobjqp[], double rg[], double r[], double hy[], 
	    double xbs[], double xs[], double y[], double y2[], double yq[],
	    Integer iz[], Integer leniz, double z[], Integer lenz,
	    Integer iparm[], double rparm[], 
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04mzh();
#endif

#ifdef NAG_PROTO
extern void e04mzj(Nag_ProblemType type, NAG_E04NKC_HESSFUN_WRAP hx, 
		   NAG_E04NKC_HESSFUN hx1, Boolean *gotr, Boolean *needlu, 
		   Nag_SparseQP_LU_Type *typelu, Integer lenr, Integer m,
		   Integer maxs, Integer mbs, Integer n, Integer nb, Integer ngqp,
		   Integer ngqp0, Integer ngobj, Integer ngobj0, Integer nnh,
		   Integer *ns, Integer itmax, Integer *itqp, double objadd,
		   double *objqp, Integer iobj, double tolfp, double tolqp, 
		   double tolx, Integer ne, Integer nka, double a[], Integer ha[], 
		   Integer ka[], Integer hfeas[], Integer hs[], Integer kbs[],
		   double ascale[], double b[], Integer lenb, double bl[], 
		   double bu[], double blbs[], double bubs[], double gobj[],
		   double gobjqp[], double gbs[], double pi[], double r[], 
		   double rc[], double rg[], double xbs[], double x0[], 
		   double xs[], double xdif[], Integer iy[], Integer iy2[],
		   double y[], double y2[], double yq[], double t[], 
		   Integer iz[], Integer leniz, double z[], Integer lenz, 
		   Integer iparm[], double rparm[], 
		   Nag_SparseQP_Save_Vars *vsave,
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04mzj();
#endif

#ifdef NAG_PROTO
extern void e04mzk(NAG_E04NKC_HESSFUN_WRAP hx, NAG_E04NKC_HESSFUN hx1, Integer *inform, 
	    Integer jq, Integer jradd, Integer lenr, Integer m, Integer mbs, 
	    Integer n, Integer nb, Integer nnh, Integer ns, Integer ne, 
	    Integer nka, double a[],  Integer ha[], Integer ka[], Integer kbs[],
	    double r[], double v[], double w[], double y[], Integer iz[], 
	    Integer leniz, double z[], Integer lenz, 
	    Integer iparm[], double rparm[], Nag_Comm *comm);
#else
extern void e04mzk();
#endif

#ifdef NAG_PROTO
extern void e04mzl(Integer j1_, Integer j2,  Boolean feasbl,  Integer m, Integer n,
	    Integer nng, Integer nng0, Integer ne, Integer nka,
	    double a[],  Integer ha[], Integer ka[], Integer hs[],
	    double g[], double pi[], double rc[],
	    Integer iparm[],  double rparm[]);
#else
extern void e04mzl();
#endif

#ifdef NAG_PROTO
extern void e04mzm(Integer m, Integer nbs, Integer n, Integer ns,
	    double gbs[], double pi[], double rg[],
	    double *rgnorm,  Integer ne, Integer nka,  double a[],
	    Integer ha[], Integer ka[], Integer kbs[], Integer iparm[],
	    double rparm[]);
#else
extern void e04mzm();
#endif

#ifdef NAG_PROTO
extern void e04mzn(Boolean *posdef,  Integer *inform, Integer lenr, Integer ns,
	    double *drsq, double r[],  Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04mzn();
#endif

#ifdef NAG_PROTO
extern void e04mzp(Integer m,  double *pinorm, double rhs[],
	    double pi[],  Integer iz[], Integer leniz,  double z[],
	    Integer lenz, Integer iparm[],  double rparm[]);
#else
extern void e04mzp();
#endif

#ifdef NAG_PROTO
extern void e04mzq(Boolean reset_xbs,  Integer *inform, Integer m, Integer n, Integer nb,
	    Integer nbs,  double *rowerr,  Integer ne, Integer nka,
	    double a[],  Integer ha[], Integer ka[],  double b[],
	    Integer lenb, Integer kbs[],  double xbs[], double xs[],
	    double y2[], double y[],  Integer iz[], Integer leniz,
	    double z[],  Integer lenz, Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04mzq();
#endif

#ifdef NAG_PROTO
extern void e04mzr(Boolean *jstfes,  Integer m, Integer mbs, Integer ns,
	    double featol,  Integer hfeas[],  double blbs[],
	    double bubs[], double gbs[], double xbs[],
	    Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04mzr();
#endif

#ifdef NAG_PROTO
extern void e04nbt(Integer mode, Integer ncolt, Integer n,  double NAG_HUGE t[],
            double NAG_HUGE y[]);
#else
extern void e04nbt();
#endif

#ifdef NAG_PROTO
extern void e04nbu(Integer n, Integer nu, Integer nrank, Integer tdr,
            Integer i, Integer j,  double NAG_HUGE r[], double NAG_HUGE u[],
            double NAG_HUGE c[], double NAG_HUGE s[]);
#else
extern void e04nbu();
#endif

#ifdef NAG_PROTO
extern void e04nbv(Integer n, Integer nu, Integer nrank, Integer tdr,
            Integer lenv, Integer lenw,  double NAG_HUGE r[], double NAG_HUGE u[],
            double NAG_HUGE v[], double NAG_HUGE w[], double NAG_HUGE c[], double NAG_HUGE s[]);
#else
extern void e04nbv();
#endif

#ifdef NAG_PROTO
extern void e04nbw(Integer mode, Integer n, Integer nz, Integer nfree,
            Integer nq, Boolean unitq, Integer NAG_HUGE kx[], double NAG_HUGE v[],
            double NAG_HUGE zy[], double NAG_HUGE wrk[]);
#else
extern void e04nbw();
#endif

#ifdef NAG_PROTO
extern void e04nbx(Boolean major, char NAG_HUGE *prbtyp, Nag_PrintType print_level,
            Integer maj_iter, Integer iter, 
            Integer nfree, Integer nrowa,
            Integer n, Integer nclin, Integer nctotl, Integer inform,
	    Boolean infeas, double obj, double bigbnd,
            Integer nactiv, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],  double NAG_HUGE a[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE c[], double NAG_HUGE clamda[], double NAG_HUGE rlamda[], double NAG_HUGE featol[], 
	    double NAG_HUGE x[], Nag_E04_Opt NAG_HUGE *opt,  Nag_Comm NAG_HUGE *comm, Nag_FileSt NAG_HUGE *stream,
            Nag_Search_State NAG_HUGE *st, Integer NAG_HUGE *local_error);
#else
extern void e04nbx();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04ncc(Integer m, Integer n, Integer nclin, 
		       double NAG_HUGE a[], Integer tda, double NAG_HUGE bl[],
		       double NAG_HUGE bu[], double NAG_HUGE cvec[], 
		       double NAG_HUGE b[], double NAG_HUGE h[], Integer tdh, 
		       Integer NAG_HUGE kx[], double NAG_HUGE x[], 
		       double NAG_HUGE *obj, Nag_E04_Opt NAG_HUGE *options, 
		       Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04ncc();
#endif

#ifdef NAG_PROTO
extern void e04ncg(char NAG_HUGE *task,  Boolean unitq,  Integer nfree, Integer n,
	    Integer nrank, Integer tdq, Integer tdr, Integer NAG_HUGE kx[],
	    double NAG_HUGE r[], double NAG_HUGE q[], double NAG_HUGE v[], double NAG_HUGE w[], 
	    NagError NAG_HUGE *fail_qdc);
#else
extern void e04ncg();
#endif

#ifdef NAG_PROTO
extern void e04nch(Boolean linobj, Boolean NAG_HUGE *rowerr, Boolean unitq,
            Integer nclin, Integer nactiv, Integer nfree, Integer nrank,
            Integer nz, Integer n, Integer nctotl, Integer tdzy,
            Integer tda, Integer tdr, Integer tdt, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], Integer NAG_HUGE *jmax,
            double NAG_HUGE *errmax, double NAG_HUGE *ctx, double NAG_HUGE *xnorm,
            double NAG_HUGE a[], double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE cq[], double NAG_HUGE res[], double NAG_HUGE res0[],
            double NAG_HUGE featol[], double NAG_HUGE r[], double NAG_HUGE t[],
            double NAG_HUGE x[], double NAG_HUGE zy[], double NAG_HUGE p[],
            double NAG_HUGE work[], Nag_ce04nc NAG_HUGE *ce04nc);
#else
extern void e04nch();
#endif

#ifdef NAG_PROTO
extern void e04ncj(char NAG_HUGE *prbtyp, Integer isdel, Integer iter, Integer jadd,
            Integer jdel, Nag_PrintType print_level, Integer nactiv,
            Integer nfree,
            Integer n, Integer nclin, Integer nrank, Integer tdr,
            Integer tdt, Integer nz, Integer nrz, Integer NAG_HUGE istate[],
            double alfa, double condrz, double condt,
            double gfnorm, double gzrnrm,  Integer numinf,
            double suminf, double ctx, double ssq,
            double NAG_HUGE ax[], double NAG_HUGE r[], double NAG_HUGE t[], double NAG_HUGE x[],
            double NAG_HUGE work[], Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], 
            double NAG_HUGE gq[], double NAG_HUGE rlamda[], Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm, 
            Nag_Search_State NAG_HUGE *st, 
            Nag_QP_Print NAG_HUGE *print_inf, Nag_FileSt NAG_HUGE *stream);
#else
extern void e04ncj();
#endif

#ifdef NAG_PROTO
extern void e04nck(Integer n, Integer nactiv,
            Integer nfree, Integer tda, Integer tdt, Integer numinf,
            Integer nz, Integer nrz, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],  double dinky,  Integer NAG_HUGE *jsmlst,
            Integer NAG_HUGE *ksmlst, Integer jinf, Integer NAG_HUGE *jtiny, Integer NAG_HUGE *jbigst,
            Integer NAG_HUGE *kbigst,  double NAG_HUGE *trulam, double NAG_HUGE a[],
            double NAG_HUGE anorms[], double NAG_HUGE gq[], double NAG_HUGE rlamda[],
            double NAG_HUGE t[], double NAG_HUGE wtinf[], Nag_ce04nc NAG_HUGE *ce04nc);
#else
extern void e04nck();
#endif

#ifdef NAG_PROTO
extern void e04ncl(Boolean hitcon, Boolean hitlow, Boolean linobj,
            Boolean unitgz,  Integer nclin, Integer nrank, Integer nrz,
            Integer n, Integer tdr, Integer jadd, Integer numinf,
            double alfa, double ctp, double NAG_HUGE *ctx,
            double NAG_HUGE *xnorm, double NAG_HUGE ap[], double NAG_HUGE ax[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE gq[], double NAG_HUGE hz[],
            double NAG_HUGE p[], double NAG_HUGE res[], double NAG_HUGE r[], double NAG_HUGE x[],
            double NAG_HUGE work[]);
#else
extern void e04ncl();
#endif

#ifdef NAG_PROTO
extern void e04ncm(Nag_ProblemType prob, Integer n, Integer nclin, Integer NAG_HUGE *litotl,
		   Integer NAG_HUGE *lwtotl, Nag_ae04nc NAG_HUGE *ae04nc, Nag_be04nb NAG_HUGE *be04nb);
#else
extern void e04ncm();
#endif


#ifdef NAG_PROTO
extern void e04ncp(char NAG_HUGE *prbtyp, Boolean linobj, Boolean NAG_HUGE *singlr, Boolean NAG_HUGE *unitgz,
            Boolean unitq,  Integer n, Integer nclin, Integer nfree,
            Integer lda, Integer ldzy, Integer ldr, Integer nrank,
            Integer nz, Integer NAG_HUGE *nrz, Integer NAG_HUGE istate[], Integer NAG_HUGE kx[],
            double bigbnd, double tolrnk,  Integer NAG_HUGE *numinf,
            double NAG_HUGE *suminf, double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE a[], double NAG_HUGE res[], double NAG_HUGE featol[],
            double NAG_HUGE gq[], double NAG_HUGE cq[], double NAG_HUGE r[], double NAG_HUGE x[],
            double NAG_HUGE wtinf[], double NAG_HUGE zy[], double NAG_HUGE wrk[]);
#else
extern void e04ncp();
#endif

#ifdef NAG_PROTO
extern void e04ncq(Boolean linobj, Boolean singlr, Boolean unitgz,
            Boolean unitq,  Integer n, Integer nclin, Integer nfree,
            Integer tda, Integer tdzy, Integer tdr, Integer nrank,
            Integer numinf, Integer nrz, Integer NAG_HUGE kx[],  double NAG_HUGE *ctp,
            double NAG_HUGE *pnorm, double NAG_HUGE a[], double NAG_HUGE ap[],
            double NAG_HUGE res[], double NAG_HUGE hz[], double NAG_HUGE p[], double NAG_HUGE gq[],
            double NAG_HUGE cq[], double NAG_HUGE r[], double NAG_HUGE zy[],
            double NAG_HUGE work[], Nag_ce04nc NAG_HUGE *ce04nc);
#else
extern void e04ncq();
#endif

#ifdef NAG_PROTO
extern void e04ncr(Integer n, Integer nclin, Integer NAG_HUGE istate[],
            double bigbnd, double NAG_HUGE *cvnorm, double NAG_HUGE *errmax,
            Integer NAG_HUGE *jmax, Integer NAG_HUGE *nviol,  double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE featol[], double NAG_HUGE x[],
            double NAG_HUGE work[], Nag_ce04nc NAG_HUGE *ce04nc);
#else
extern void e04ncr();
#endif

#ifdef NAG_PROTO
extern void e04nct(Boolean unitq,  Integer n, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree,
            Integer nres, Integer ngq, Integer NAG_HUGE *nz, Integer NAG_HUGE *nrz,
            Integer tda, Integer tdzy, Integer tdr, Integer tdt,
            Integer nrank, Integer NAG_HUGE *jdel, Integer NAG_HUGE *kdel, Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],  double NAG_HUGE a[], double NAG_HUGE res[], double NAG_HUGE r[],
            double NAG_HUGE t[], double NAG_HUGE gq[], double NAG_HUGE zy[], double NAG_HUGE c[],
            double NAG_HUGE s[], Nag_ce04nc NAG_HUGE *ce04nc, Nag_de04nb NAG_HUGE *de04nb);
#else
extern void e04nct();
#endif

#ifdef NAG_PROTO
extern void e04ncu(Boolean cold, Boolean vertex,  Integer nclin, Integer nctotl,
            Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nartif, Integer NAG_HUGE *nfree, Integer n,
            Integer lda, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            double bigbnd, double tolact, double NAG_HUGE a[],
            double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE x[],
            double NAG_HUGE wx[], double NAG_HUGE work[], Nag_ce04nc NAG_HUGE *ce04nc);
#else
extern void e04ncu();
#endif

#ifdef NAG_PROTO
extern void e04ncv(Boolean NAG_HUGE *unitq,  Integer NAG_HUGE *inform, Integer ifix, Integer iadd,
            Integer jadd, Integer nactiv, Integer nz, Integer nfree,
            Integer nrank, Integer nres, Integer ngq, Integer n,
            Integer tda, Integer tdzy, Integer tdr, Integer tdt,
            Integer NAG_HUGE kx[],  double condmx, double NAG_HUGE a[], double NAG_HUGE r[],
            double NAG_HUGE t[], double NAG_HUGE res[], double NAG_HUGE gq[], double NAG_HUGE zy[],
            double NAG_HUGE w[], double NAG_HUGE c[], double NAG_HUGE s[], 
            NagError NAG_HUGE *ovflow, Nag_ce04nc NAG_HUGE *ce04nc, Nag_de04nb NAG_HUGE *de04nb);
#else
extern void e04ncv();
#endif

#ifdef NAG_PROTO
void e04ncw(char NAG_HUGE *prbtyp, Integer tdh, Integer n, Integer NAG_HUGE *nrank,  double tolrnk,
             Integer NAG_HUGE kx[],  double NAG_HUGE h[],  Nag_PrintType print_level,
             Nag_FileSt NAG_HUGE *stream, Integer NAG_HUGE *inform);
#else
extern void e04ncw();
#endif

#ifdef NAG_PROTO
extern void e04ncx(Boolean NAG_HUGE *unitq, Integer NAG_HUGE *inform, Integer NAG_HUGE *nz, Integer NAG_HUGE *nfree,
            Integer nrank, Integer nres, Integer ngq, Integer n,
            Integer tdzy, Integer tda, Integer tdr, Integer tdt,
            Integer NAG_HUGE istate[], Integer NAG_HUGE kx[], double condmx,
            double NAG_HUGE a[], double NAG_HUGE r[], double NAG_HUGE t[], double NAG_HUGE res[],
            double NAG_HUGE gq[], double NAG_HUGE zy[], double NAG_HUGE w[], double NAG_HUGE c[],
            double NAG_HUGE s[], NagError NAG_HUGE *ovflow,
            Nag_ce04nc NAG_HUGE *ce04nc, Nag_de04nb NAG_HUGE *de04nb);
#else
extern void e04ncx();
#endif

#ifdef NAG_PROTO
extern void e04ncy(Boolean NAG_HUGE *unitq, Boolean vertex, Integer NAG_HUGE *inform, Integer k1,
            Integer k2, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nartif, Integer NAG_HUGE *nz,
            Integer NAG_HUGE *nfree, Integer nrank, Integer NAG_HUGE *nrejtd, Integer nres,
            Integer ngq, Integer n, Integer tdzy, Integer tda,
            Integer tdr, Integer tdt, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],  double condmx, double NAG_HUGE a[], double NAG_HUGE r[],
            double NAG_HUGE t[], double NAG_HUGE res[], double NAG_HUGE gq[], double NAG_HUGE zy[],
            double NAG_HUGE w[], double NAG_HUGE c[], double NAG_HUGE s[], NagError NAG_HUGE *ovflow,
            Nag_ce04nc NAG_HUGE *ce04nc, Nag_de04nb NAG_HUGE *de04nb);
#else
extern void e04ncy();
#endif

#ifdef NAG_PROTO
extern void e04ncz(char NAG_HUGE *prbtyp, Boolean linobj,
            Boolean NAG_HUGE *unitq,  Integer NAG_HUGE *inform, Integer NAG_HUGE *iter, Integer jinf,
            Integer nclin, Integer nctotl, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree,
            Integer nrank, Integer NAG_HUGE *nz, Integer NAG_HUGE *nrz, Integer n,
            Integer tda, Integer tdr, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[],  double NAG_HUGE *ctx, double NAG_HUGE *ssq, double ssq1,
            double NAG_HUGE *suminf,  Integer NAG_HUGE *numinf,  double NAG_HUGE *xnorm,
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE a[],
            double NAG_HUGE clamda[], double NAG_HUGE ax[], double NAG_HUGE featol[],
            double NAG_HUGE r[], double NAG_HUGE x[], double NAG_HUGE w[], Nag_E04_Opt NAG_HUGE *opt,
            Nag_Comm NAG_HUGE *comm,  Nag_FileSt NAG_HUGE *stream, Nag_Search_State NAG_HUGE *st,
            NagError NAG_HUGE *ovflow, Nag_ae04nc NAG_HUGE *ae04nc,
            Nag_be04nb NAG_HUGE *be04nb, Nag_ce04nc NAG_HUGE *ce04nc,
            Nag_de04nb NAG_HUGE *de04nb, Nag_fe04nb NAG_HUGE *fe04nb, Integer NAG_HUGE *local_error);
#else
extern void e04ncz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP 
  void NAG_CALL e04nfc(Integer n, Integer nclin, double NAG_HUGE a[], 
		       Integer tda, double NAG_HUGE bl[], double NAG_HUGE bu[],
		       double NAG_HUGE cvec[], double NAG_HUGE h[], 
		       Integer tdh, NAG_E04NFC_FUN qphess_user,
		       double NAG_HUGE x[], double NAG_HUGE *obj, 
		       Nag_E04_Opt NAG_HUGE *options, Nag_Comm NAG_HUGE *comm,
		       NagError NAG_HUGE *fail);
#else
extern void e04nfc();
#endif

#ifdef NAG_PROTO
extern void e04nfg(char NAG_HUGE *prbtyp, Integer n, Integer nclin, double NAG_HUGE a[], Integer tda,
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE cvec[], double NAG_HUGE ax[], double NAG_HUGE h[], Integer tdh,
	    NAG_E04NFC_FUN qphess,
            double NAG_HUGE x[], double NAG_HUGE *obj, double NAG_HUGE lambda[], Integer NAG_HUGE istate[],
            Boolean vertex, Boolean cset, Nag_ae04mf NAG_HUGE *ae04mf, double NAG_HUGE w[],
            Integer NAG_HUGE iw[], Nag_ee04mf NAG_HUGE *ee04mf, Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm,
            Nag_FileSt NAG_HUGE *stream, NagError NAG_HUGE *ovflow, Nag_EndState NAG_HUGE *endstate);
#else
extern void e04nfg();
#endif

#ifdef NAG_PROTO
extern void e04nfp(Boolean unitq, Integer it, Integer n, Integer NAG_HUGE *nactiv,
            Integer NAG_HUGE *nfree, Integer ngq, Integer NAG_HUGE *nz, Integer NAG_HUGE *nrz,
            Integer tda, Integer tdq, Integer tdt, Integer jdel,
            Integer kdel, Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], double NAG_HUGE a[],
            double NAG_HUGE t[], double NAG_HUGE gqm[], double NAG_HUGE tgqm[], double NAG_HUGE q[],
            double NAG_HUGE c[], double NAG_HUGE s[], Nag_de04nb NAG_HUGE *de04nb);
#else
extern void e04nfp();
#endif

#ifdef NAG_PROTO
extern void e04nfq(Boolean NAG_HUGE *unitq, Boolean vertex, Integer k1, Integer k2,
            Integer NAG_HUGE *it, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nartif, Integer NAG_HUGE *nz,
            Integer NAG_HUGE *nfree, Integer NAG_HUGE *nrejtd, Integer ngq, Integer n,
            Integer tdq, Integer tda, Integer tdt, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[],
            double NAG_HUGE a[], double NAG_HUGE t[], double NAG_HUGE gqm[], double NAG_HUGE tgqm[], double NAG_HUGE q[],
            double NAG_HUGE w[], double NAG_HUGE c[], double NAG_HUGE s[],
            Nag_de04nb NAG_HUGE *de04nb, Nag_ee04mf NAG_HUGE *ee04mf, NagError NAG_HUGE *ovflow);
#else
extern void e04nfq();
#endif

#ifdef NAG_PROTO
extern void e04nfr(Boolean NAG_HUGE *unitq, Boolean rset, Integer NAG_HUGE *inform, Integer ifix,
            Integer iadd, Integer jadd, Integer NAG_HUGE *it, Integer nactiv,
            Integer nz, Integer nfree, Integer nrz, Integer ngq,
            Integer n, Integer tda, Integer tdq, Integer tdr,
            Integer tdt, Integer NAG_HUGE kx[], double condmx, double NAG_HUGE *drzz,
            double NAG_HUGE a[], double NAG_HUGE r[], double NAG_HUGE t[], double NAG_HUGE gqm[], double NAG_HUGE tgqm[],
            double NAG_HUGE q[], double NAG_HUGE w[], double NAG_HUGE c[], double NAG_HUGE s[],
            Nag_de04nb NAG_HUGE *de04nb, Nag_ee04mf NAG_HUGE *ee04mf, NagError NAG_HUGE *ovflow);
#else
extern void e04nfr();
#endif

#ifdef NAG_PROTO
extern void e04nfs(char NAG_HUGE *prbtyp, Boolean NAG_HUGE *header, Boolean rset,
            Nag_PrintType print_level, Integer isdel,
            Integer NAG_HUGE *jdel, Integer jadd, Integer n, Integer nclin,
            Integer nactiv, Integer nfree, Integer nz, Integer nrz,
            Integer tdr, Integer tdt, Integer NAG_HUGE istate[], double alfa,
            double condrz, double condt, double drzz,
            double grznrm, Integer numinf,
            double suminf, Integer notopt, double objp,
            double trusml, double NAG_HUGE ax[], double NAG_HUGE r[], double NAG_HUGE t[],
            double NAG_HUGE x[], Integer NAG_HUGE kx[], Integer NAG_HUGE kactive[], double NAG_HUGE rlambda[],
            double NAG_HUGE gq[], double NAG_HUGE work[], Nag_QP_Print NAG_HUGE *print_inf,
            Nag_E04_Opt NAG_HUGE *opt, Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm);
#else
extern void e04nfs();
#endif

#ifdef NAG_PROTO
extern void e04nft(Boolean cset, Integer n, Integer nclin, Nag_ae04mf NAG_HUGE *ae04mf);
#else
extern void e04nft();
#endif

#ifdef NAG_PROTO
extern void e04nfu(Integer n, Integer jthcol, double NAG_HUGE h[], Integer tdh,
            double NAG_HUGE x[], double NAG_HUGE hx[], Nag_Comm NAG_HUGE *comm);
#else
extern void e04nfu();
#endif

#ifdef NAG_PROTO
extern void e04nfv(Boolean delreg, Boolean posdef, Boolean statpt,
            Boolean unitgz, Boolean unitq, Integer n, Integer nclin,
            Integer nfree, Integer tda, Integer tdq, Integer tdr,
            Integer nrz, Integer issave, Integer jdsave,
            Integer NAG_HUGE kx[], double NAG_HUGE *dnorm, double NAG_HUGE *gzdz,
            double NAG_HUGE a[], double NAG_HUGE ad[], double NAG_HUGE d[], double NAG_HUGE gq[],
            double NAG_HUGE r[], double NAG_HUGE q[], double NAG_HUGE v[], Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04nfv();
#endif

#ifdef NAG_PROTO
extern void e04nfx(Boolean unitq, NAG_E04NFC_FUN qphess,
            Integer maxnz, Integer n, Integer ngq, Integer NAG_HUGE *nrz,
            Integer nz, Integer nfree, Integer tdq, Integer tdh,
            Integer tdr, Integer NAG_HUGE kx[], double NAG_HUGE *hsize, double tolrnk,
            double NAG_HUGE gq[], double NAG_HUGE h[], double NAG_HUGE r[], double NAG_HUGE q[],
            double NAG_HUGE hz[], double NAG_HUGE wrk[], Nag_Comm NAG_HUGE *comm, Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04nfx();
#endif

#ifdef NAG_PROTO
extern void e04nfy(Boolean NAG_HUGE *singlr, Boolean NAG_HUGE *posdef, Boolean NAG_HUGE *renewr,
            Boolean unitq, Integer n, Integer nrz, Integer nfree,
            Integer tdq, Integer tdh, Integer tdr, Integer NAG_HUGE kx[],
            double NAG_HUGE *hsize, double NAG_HUGE *drzz, double tolrnk,
	    NAG_E04NFC_FUN qphess,
            double NAG_HUGE h[], double NAG_HUGE r[], double NAG_HUGE q[],
            double NAG_HUGE hz[], double NAG_HUGE wrk[], Nag_Comm NAG_HUGE *comm, Nag_ee04mf NAG_HUGE *ee04mf);
#else
extern void e04nfy();
#endif

#ifdef NAG_PROTO
extern void e04nfz(char NAG_HUGE *prbtyp, char NAG_HUGE *msg, Boolean cset,
            Boolean NAG_HUGE *unitq, Integer itmax,
            Integer NAG_HUGE *nviol, Integer n, Integer nclin, Integer lda,
            Integer ldh, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree, Integer NAG_HUGE *nrz,
            Integer NAG_HUGE *nz, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[],
	    NAG_E04NFC_FUN qphess,
            void (NAG_HUGE *print_iter)(char NAG_HUGE *prbtyp, Boolean NAG_HUGE *header, Boolean rset,
                               Nag_PrintType print_level, Integer isdel, Integer NAG_HUGE *jdel,
                               Integer jadd, Integer n, Integer nclin, Integer nactiv,
                               Integer nfree, Integer nz, Integer nrz, Integer tdr,
                               Integer tdt, Integer NAG_HUGE istate[], double alfa, double condrz,
                               double condt, double drzz, double grznrm,
                               Integer numinf, double suminf, Integer notopt,
                               double objqp, double trusml, double NAG_HUGE ax[], double NAG_HUGE r[],
                               double NAG_HUGE t[], double NAG_HUGE x[], Integer NAG_HUGE kx[], Integer NAG_HUGE kactive[],
                               double NAG_HUGE rlambda[], double NAG_HUGE gq[], double NAG_HUGE work[],
                               Nag_QP_Print NAG_HUGE *print_inf, Nag_E04_Opt NAG_HUGE *opt,
                               Nag_FileSt NAG_HUGE *stream, Nag_Comm NAG_HUGE *comm),
            double NAG_HUGE *objqp, double NAG_HUGE *xnorm, double NAG_HUGE *hsize, double NAG_HUGE a[],
            double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE cvec[],
            double NAG_HUGE featol[], double NAG_HUGE featlu[], double NAG_HUGE h[],
            double NAG_HUGE x[], Nag_E04_Opt NAG_HUGE *opt, double NAG_HUGE w[],
            Nag_de04nb NAG_HUGE *de04nb, Nag_ae04mf NAG_HUGE *ae04mf, Nag_ce04mf NAG_HUGE *ce04mf,
            Nag_de04mf NAG_HUGE *de04mf, Nag_QP_Print NAG_HUGE *print_inf, Nag_FileSt NAG_HUGE *stream,
            Nag_Comm NAG_HUGE *comm, Nag_ee04mf NAG_HUGE *ee04mf, NagError NAG_HUGE *ovflow);
#else
extern void e04nfz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL 
  e04nkc(Integer n, Integer m, Integer nnz, Integer iobj, Integer ncolh, 
	 NAG_E04NKC_HESSFUN qphx, double a[], Integer ha[], Integer ka[], 
	 double bl[], double bu[], double xs[], 
	 Integer *ninf, double *suminf, double *obj, 
	 Nag_E04_Opt *options, Nag_Comm *user_comm, NagError *fail);
#else
extern void e04nkc();
#endif

#ifdef NAG_PROTO
extern void e04nkg(Integer nbs, Integer ngobj, Integer ngobj0, Integer kbs[],
	    double gobj[], double gbs[], Integer iparm[], double rparm[]);
#else
extern void e04nkg();
#endif

#ifdef NAG_PROTO
extern void e04nkh(Integer nr, Integer lenr,  double r[], double v[],
	    double w[],  Integer lastv);
#else
extern void e04nkh();
#endif

#ifdef NAG_PROTO
extern void e04nkj(Integer m, Integer maxs, Integer mbs, Integer n, Integer nb,
	    Integer *ns, Integer *nssave, Integer hs[], Integer kbs[],
	    double bl[], double bu[], double xs[],
	    Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nkj();
#endif

#ifdef NAG_PROTO
extern void e04nkk(Integer j, Integer iobj, Integer n, char **names, char *name);
#else
extern void e04nkk();
#endif

#ifdef NAG_PROTO
extern void e04nkl(Boolean ondisk, Integer m, Integer n, Integer nb,
	    Integer nnobj, Integer ne, Integer nka,
	    double a[], Integer ha[], Integer ka[], Integer hs[],
	    double ascale[], double bl[], double bu[], double gobj[], 
	    double pi[], double rc[], double xs[], double y[], 
	    char **names, char *istate, Integer iz[], Integer leniz,
	    double z[], Integer lenz, Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nkl();
#endif

#ifdef NAG_PROTO
extern void e04nkm(/* Boolean ondisk, */double bplus, double tolfea,
	    double tolopt,  Integer *js,  double d1, double d2,
	    double *djtest,  Integer j, char *name,
	    double xj, double cj, double b1, double b2,
	    double dj,  Integer k, Integer number,  double clamda,
	    Integer n, Integer m, Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nkm();
#endif

#ifdef NAG_PROTO
extern void e04nkn(Nag_ProblemType type, Integer k, char *istate);
#else
extern void e04nkn();
#endif

#ifdef NAG_PROTO
extern void e04nkr(Boolean initialize, Integer iparm[], double rparm[], 
		   Nag_E04_Opt *opt);
#else
extern void e04nkr();
#endif

#ifdef NAG_PROTO
extern void e04nks(NAG_E04NKC_HESSFUN_WRAP hx, NAG_E04NKC_HESSFUN hx1, 
	    Integer snmodx, Nag_Start start, Integer lenb, 
	    Integer mxx, Integer nxx, Integer nbxx, Integer nexx, 
	    Integer nkax, Integer iobjxx, double objadd, double *obj, 
	    double a[], Integer ha[], Integer ka[], double b[], double bl[],
	    double bu[], double gobj[], char **names, Integer hs[],  
	    double x0[], double xs[], double pi[], double rc[],
	    Integer *inform, Integer *ns, Integer iz[], Integer leniz,
	    double z[], Integer lenz, Nag_ProblemType *type, Integer *itqp, 
	    Integer *mxitqp, Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nks();
#endif

#ifdef NAG_PROTO
extern void e04nkt(Boolean print,  Integer m, Integer n, Integer nb, Integer ngobj,
		   Integer nscl, Integer ne, Integer nka, double a[],  
		   Integer ha[], Integer ka[], Integer hs[], double ascale[], 
		   double bl[], double bu[], double gobj[], char **names,
		   double pi[], double rc[], double xs[], double y[],  Integer iz[], 
		   Integer leniz,  double z[], Integer lenz,  Nag_ProblemType type,  
		   double *pnorm1, double *pnorm2,
		   Integer iparm[], double rparm[],
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nkt();
#endif

#ifdef NAG_PROTO
extern void e04nkv(NAG_E04NKC_HESSFUN_WRAP sqhx, NAG_E04NKC_HESSFUN qphx, 
	    Nag_ProblemType type, Boolean needb, Integer lenb, Integer lenb0, 
	    Integer lenr, Integer m, Integer maxs, Integer mbs, Integer n, 
	    Integer nb, Integer ngqp0, Integer ngobj0, Integer ngobj, 
	    Integer nnh, Integer *ns, Integer nscl, Integer numleq, 
	    Integer itmax, Integer *itqp, double objadd, double *obj,
	    Integer iobj, double tolfp, double tolqp, double tolx, Integer ne,
	    Integer nka, double a[], Integer ha[], Integer ka[], Integer hfeas[], 
	    Integer hs[], Integer kbs[], double ascale[], double b[], 
	    double bl[], double bu[], double blbs[], double bubs[], 
	    double blslk[], double buslk[], double gobj[], double gobjqp[], 
	    double gbs[], double pi[], double r[], double rc[], double rg[], 
	    double rhs[], double xbs[], double x0[], double xs[], double xdif[],
	    Integer iy[], Integer iy2[], double y[], double y2[], double yq[], 
	    double t[], Integer iz[], Integer leniz, double z[], Integer lenz, 
	    Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nkv();
#endif

#ifdef NAG_PROTO
extern void e04nkw(NAG_E04NKC_HESSFUN qphx,  Integer nnh, Integer ne, Integer nka,
	    double a[],  Integer ha[], Integer ka[],  double x[],
	    double hx[],  Integer nstate, Integer iz[], Integer leniz,
	    double z[],  Integer lenz, Integer iparm[],
	    double rparm[], Nag_Comm *comm);
#else
extern void e04nkw();
#endif

#ifdef NAG_PROTO
extern void e04nkx(Integer nrealloc, Integer minfac, Integer *miniz, Integer *minz, 
	    Integer iparm[], double rparm[]);
#else
extern void e04nkx();
#endif

#ifdef NAG_PROTO
extern void e04nlg(Boolean pivot,  Integer *inform, Integer n,  double hdmin,
             double *d_max,  Integer *irank, Integer perm[], Integer lenh,
             double h[]);
#else
extern void e04nlg();
#endif

#ifdef NAG_PROTO
extern void e04nlh(Integer nr, Integer lenr,  double r[], double *drmin,
	    double *rmax,  Integer iparm[],  double rparm[]);
#else
extern void e04nlh();
#endif

#ifdef NAG_PROTO
extern void e04nlj(Integer jq, Integer jr, Integer lenr, Integer nr,
	    double r[], double w[]);
#else
extern void e04nlj();
#endif

#ifdef NAG_PROTO
extern void e04nlk(Integer n, Integer lenr,  double r[], double v[],
	    double w[],  Integer lastv,  double vnorm,
	    double tolz);
#else
extern void e04nlk();
#endif

#ifdef NAG_PROTO
extern void e04nll(MatrixTranspose transpose,  Integer lenr, Integer nr,  double r[],
	    double y[]);
#else
extern void e04nll();
#endif

#ifdef NAG_PROTO
extern void e04nlm(Integer *inform, Integer nbs, Integer jqsave, Integer kbs[],
	    double *gtd, double d[],  Integer iparm[], double rparm[], 
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nlm();
#endif

#ifdef NAG_PROTO
extern void e04nln(Integer m, Integer mbs, Integer n, Integer nb, Integer ns,
	    Integer *kbsq,  double *pivot,  Integer ne, Integer nka,
	    double a[],  Integer ha[], Integer ka[], Integer kbs[],
	    double bl[], double bu[], double xbs[],
	    double y[],  Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nln();
#endif

#ifdef NAG_PROTO
extern void e04nlp(Integer nbs,  double featol, double stepmx,
	    Integer hfeas[],  double blbs[], double bubs[],
	    double xbs[], double y[],  Boolean *hitlow, Boolean *move,
	    Boolean *onbnd, Boolean *unbndd,  double *bound,
	    double *exact, double *alpha, double *alphap,
	    Integer iparm[],  double rparm[]);
#else
extern void e04nlp();
#endif

#ifdef NAG_PROTO
extern void e04nlq(Nag_DegenJob job, Integer *inform, Integer nb, double *featol,
	    double tolx, Integer hs[], double bl[], double bu[], double xs[],
	    Integer iparm[], double rparm[], Nag_SparseQP_Save_Vars *vsave,
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nlq();
#endif

#ifdef NAG_PROTO
extern void e04nlr(Boolean fix_sbasics,  Integer m, Integer maxs, Integer mbs, Integer n,
	    Integer nb, Integer *ns, Integer hs[], Integer kbs[],
	    double bl[], double bu[], double blbs[],
	    double bubs[], double xs[], double xbs[],
	    Integer iparm[],  double rparm[]);
#else
extern void e04nlr();
#endif

#ifdef NAG_PROTO
extern void e04nls(Boolean qpstep, Boolean posdef,  Integer lenr, Integer n,
	    double r[], double g[], double d[], double *gd,
	    double *dhd);
#else
extern void e04nls();
#endif

#ifdef NAG_PROTO
extern void e04nlt(Boolean pivot,  Integer lenr, Integer m, Integer maxr,
	    Integer mbs, Integer nb, Integer *ns,  double hdmax,
	    Integer hs[], Integer kbs[], Integer perm[],  double bl[],
	    double bu[], double blbs[], double bubs[],
	    double xs[], double xbs[], double r[],
	    Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nlt();
#endif

#ifdef NAG_PROTO
extern void e04nlu(Nag_SetStateMode state_mode,  Integer nb,  double bl[], double bu[],
	    Integer hs[],  double xs[]);
#else
extern void e04nlu();
#endif

#ifdef NAG_PROTO
extern void e04nlv(NAG_E04NKC_HESSFUN_WRAP hx, NAG_E04NKC_HESSFUN hx1, Integer lenr, Integer m,
	    Integer mbs, Integer n, Integer nb, Integer nnh, Integer ns,
	    Integer ne, Integer nka,  double a[],  Integer ha[],
	    Integer ka[],  double *hdmax,  Integer kbs[],
	    double r[], double v[], double w[], double y[],
	    Integer iz[], Integer leniz,  double z[],  Integer lenz,
	    Integer iparm[],  double rparm[], Nag_Comm *comm);
#else
extern void e04nlv();
#endif

#ifdef NAG_PROTO
extern void e04nlw(Nag_ProblemType type, Boolean prtfes, Boolean jstfes, Integer mbs,
	    Integer ns, Integer itn, double step, Integer ninf, 
	    double suminf, double obj, Integer kbs[], double xbs[], 
	    Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nlw();
#endif

#ifdef NAG_PROTO
extern void e04nlx(Boolean incres, Boolean *needpi,  Integer m1, Integer m,
	    Integer nb,  double featol, double *step,
	    Integer hfeas[], Integer hs[], Integer kbs[],  double bl[],
	    double bu[], double blbs[], double bubs[],
	    double xbs[], double xs[], double y[],
	    double yq[],  Integer iz[], Integer leniz,  double z[],
	    Integer lenz, Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nlx();
#endif

#ifdef NAG_PROTO
extern void e04nly(Nag_ProblemType type,  Boolean *needlu,  Integer m, Integer n,
		   Integer nb, Integer itmax, Integer *itlp,  double objadd,
		   Integer iobj,  double tolfp, double tollp,
		   double tolx,  Integer ne, Integer nka,  double a[],
		   Integer ha[], Integer ka[], Integer hfeas[], Integer hs[],
		   Integer kbs[],  double ascale[], double b[],
		   Integer lenb,  double bl[], double bu[],
		   double blbs[], double bubs[], double gbs[],
		   double pi[], double rc[], double xbs[],
		   double xs[],  Integer iy[], Integer iy2[],  double y[],
		   double y2[], double yq[],  Integer iz[], Integer leniz,
		   double z[],  Integer lenz, Integer iparm[], double rparm[],
		   Nag_SparseQP_Save_Vars *vsave,
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nly();
#endif

#ifdef NAG_PROTO
extern void e04nlz(Integer m, Integer n, Integer nb, Integer nng, Integer nng0,
	    Integer nnh, Integer ns,  Boolean *incres,  Integer ne,
	    Integer nka,  double a[],  Integer ha[], Integer ka[],
	    Integer *newsb, Integer *nonopt, Integer hs[],  double bl[],
	    double bu[], double g[], double pi[], double rc[],
	    Integer iparm[],  double rparm[]);
#else
extern void e04nlz();
#endif

#ifdef NAG_PROTO
extern void e04nmg(Boolean stats_requd,  Integer m, Integer n, Integer nb, Integer ne,
	    Integer nka,  double a[],  Integer ha[], Integer ka[],
	    double bl[], double bu[],  Integer hrtype[],
	    Integer iparm[],  double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nmg();
#endif

#ifdef NAG_PROTO
extern void e04nmh(MatrixTranspose transpose,  Integer ne, Integer lenka,  double a[],
	    Integer ha[], Integer ka[],  double alpha, double x[],
	    Integer lenx,  double beta, double y[],  Integer leny,
	    Integer iparm[],  double rparm[]);
#else
extern void e04nmh();
#endif

#ifdef NAG_PROTO
extern void e04nmj(Boolean *needlu,  Nag_SparseQP_LU_Type typelu,  Integer m, Integer mbs,
		   Integer n, Integer nb, Integer nnl, Integer ns, Integer *nswap,
		   Integer ne, Integer nka,  double a[],  Integer ha[],
		   Integer ka[], Integer kbs[], Integer hs[],  double b[],
		   Integer lenb,  double bl[], double bu[],
		   double blbs[], double bubs[], double xbs[],
		   double xs[],  Integer iy[], Integer iy2[],  double y[],
		   double y2[],  Integer iz[], Integer leniz,  double z[],
		   Integer lenz, Integer iparm[],  double rparm[],
		   Nag_SparseQP_Save_Vars *vsave,
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nmj();
#endif

#ifdef NAG_PROTO
extern void e04nmk(Integer nb,  double bl[], double bu[], double xs[],
	    double *binf,  Integer *jbinf);
#else
extern void e04nmk();
#endif

#ifdef NAG_PROTO
extern void e04nml(Nag_SparseQP_BasisFactType basis_fact,  Integer *inform, 
		   Integer m, Integer n, Integer nbs,
		   Integer ne, Integer nka,  double a[],  Integer ha[],
		   Integer ka[], Integer kbs[], Integer ip[],  double alu[],
		   Integer indc[], Integer indr[], Integer lena, Integer iy[],
		   Integer iy2[],  double y[],  Integer iz[], Integer leniz,
		   Integer iparm[], double rparm[], 
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nml();
#endif

#ifdef NAG_PROTO
extern void e04nmm(Integer mode, Integer nrealloc, Integer minfac, Integer m, 
	    Integer n, Integer ne, Integer miniz, Integer minz, 
	    Integer *nizmin, Integer *nzmin, Integer iparm[],  double rparm[]);
#else
extern void e04nmm();
#endif

#ifdef NAG_PROTO
extern void e04nmn(Integer *inform, Integer jrep, Integer m,  double w[],
	    Integer iz[], Integer leniz,  double z[],  Integer lenz,
	    Integer iparm[],  double rparm[]);
#else
extern void e04nmn();
#endif

#ifdef NAG_PROTO
extern void e04nmp(MatrixTranspose transpose, Integer m, Integer n, Integer lenkbs, 
	    Integer kbs[], Integer ne, Integer nka, double a[], Integer ha[],
	    Integer ka[], double alpha, double x[], Integer lenx, double beta, 
	    double y[], Integer leny, Integer iparm[], double rparm[]);
#else
extern void e04nmp();
#endif

#ifdef NAG_PROTO
extern void e04nmq(Nag_SparseQP_LU_SolveType solve_type,  Integer *inform, Integer m,  double w[],
	    double y[],  Integer iz[], Integer leniz,  double z[],
	    Integer lenz, Integer iparm[],  double rparm[]);
#else
extern void e04nmq();
#endif

#ifdef NAG_PROTO
extern void e04nmr(Integer lcrash, Integer m, Integer n, Integer nb, Integer ne,
	    Integer nka,  double a[],  Integer ha[], Integer ka[],
	    Integer hpiv[], Integer hs[], Integer hrtype[],  double bl[],
	    double bu[], double xs[], Integer iparm[], double rparm[],
	    Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nmr();
#endif

#ifdef NAG_PROTO
extern void e04nms(Integer nb, Integer jobj,  double bl[], double bu[],
	    double rc[], double xs[], double *dinf,
	    Integer *jdinf);
#else
extern void e04nms();
#endif

#ifdef NAG_PROTO
extern void e04nmt(Integer nbs, Integer m, Integer nb, Integer hs[], Integer ip[],
	    Integer kbs[], Integer kbsold[], Integer locr[], Integer *nswap);
#else
extern void e04nmt();
#endif

#ifdef NAG_PROTO
extern void e04nmu(Boolean feasbl,  double featol,  Integer minimz, Integer m,
	    Integer n, Integer nb, Integer nnobj, Integer ne, Integer nka,
	    double a[],  Integer ha[], Integer ka[], Integer hs[],
	    double bl[], double bu[], double gobj[],
	    double pi[], double rc[], double xs[],
	    Integer iparm[],  double rparm[]);
#else
extern void e04nmu();
#endif

#ifdef NAG_PROTO
extern void e04nmv(Integer m, Integer n, Integer nb, Integer nnl, Integer nncon,
		   Integer nnjac, Integer hrtype[], Integer ne, Integer nka,
		   double a[],  Integer ha[], Integer ka[], double ascale[], 
		   double bl[], double bu[], double rmin[], double rmax[],  
		   Integer iparm[], double rparm[], 
		   Nag_E04_Opt *opt, Nag_Comm *comm, Nag_FileSt *stream);
#else
extern void e04nmv();
#endif

#ifdef NAG_PROTO
extern void e04nmw(Boolean scale,  Integer m, Integer n, Integer nb, Integer ne,
	    Integer nka,  double a[],  Integer ha[], Integer ka[],
	    double ascale[], double bl[], double bu[],
	    double pi[], double xs[],  Integer iparm[],
	    double rparm[]);
#else
extern void e04nmw();
#endif

#ifdef NAG_PROTO
extern void e04nmx(Integer mbs, Integer m, Integer n, Integer nb, double w[],
	    Integer ip[], Integer iq[], double bl[], double bu[],
	    Integer hs[], Integer kbs[], double xs[], Integer iparm[], 
	    double rparm[], Nag_E04_Opt *opt, Nag_Comm *comm,
	    Nag_FileSt *stream);
#else
extern void e04nmx();
#endif

#ifdef NAG_PROTO
extern void e04nmy(Integer jq, Integer m, Integer n, Integer ne, Integer nka,
             double a[],  Integer ha[], Integer ka[],  double y[]);
#else
extern void e04nmy();
#endif

#ifdef NAG_PROTO
extern void e04nmz(Boolean xbs_to_xs,  Integer ms, Integer nb, Integer kbs[],
	    double xbs[], double xs[]);
#else
extern void e04nmz();
#endif

#ifdef NAG_PROTO
extern void e04paz(Nag_Search_State NAG_HUGE *st);
#else
extern void e04paz();
#endif

#ifdef NAG_PROTO
extern void e04uc0(Boolean debug, Boolean NAG_HUGE *done, Boolean NAG_HUGE *first, Boolean NAG_HUGE *imprvd,
            Integer NAG_HUGE *inform, double alfmax, double alfsml, double epsaf,
            double eta, double NAG_HUGE *xtry, double ftry, double oldf, double oldg,
            double NAG_HUGE *tolabs, double tolrel, double toltny,
            double NAG_HUGE *alfa, double NAG_HUGE *alfbst, double NAG_HUGE *fbest, Nag_Uc01St NAG_HUGE *st);
#else
extern void e04uc0();
#endif

#ifdef NAG_PROTO
extern void e04uc1(Boolean debug, Boolean NAG_HUGE *done, Boolean NAG_HUGE *first, Boolean NAG_HUGE *imprvd,
            Integer NAG_HUGE *inform, double alfmax, double epsaf, double eta,
            double NAG_HUGE *xtry, double ftry, double gtry, double oldf,
            double oldg, double NAG_HUGE *tolabs, double tolrel, double toltny,
            double NAG_HUGE *alfa, double NAG_HUGE *alfbst, double NAG_HUGE *fbest, double NAG_HUGE *gbest,
            Nag_Uc01St NAG_HUGE *st);
#else
extern void e04uc1();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL 
  e04ucc(Integer n, Integer nclin, Integer ncnlin, double NAG_HUGE a[], 
	 Integer tda, double NAG_HUGE bl[], double NAG_HUGE bu[],
	 NAG_E04UCC_FUN objfun, NAG_E04UCC_CONFUN confun,
	 double NAG_HUGE x[], double NAG_HUGE *objf, double NAG_HUGE objgrad[],
	 Nag_E04_Opt NAG_HUGE *options, Nag_Comm NAG_HUGE *user_comm, 
	 NagError NAG_HUGE *fail);
#else
extern void e04ucc();
#endif

#ifdef NAG_PROTO
extern void e04ucg(Boolean firstv, Boolean NAG_HUGE *hitlow,  Integer NAG_HUGE istate[],
            Integer NAG_HUGE *inform, Integer NAG_HUGE *jadd, Integer n, Integer nctotl,
            Integer numinf,  double NAG_HUGE *alfa, double NAG_HUGE *palfa,
            double NAG_HUGE *atphit, double bigalf, double bigbnd,
            double pnorm, double NAG_HUGE anorm[], double NAG_HUGE ap[],
            double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE featol[], double NAG_HUGE p[], double NAG_HUGE x[], Nag_fe04nb NAG_HUGE *fe04nb);
#else
extern void e04ucg();
#endif

#ifdef NAG_PROTO
extern void e04uch(Boolean firstv, Boolean negstp,  double bigalf,
            double bigbnd, double pnorm,  Integer NAG_HUGE *jadd1,
            Integer NAG_HUGE *jadd2,  double NAG_HUGE *palfa1, double NAG_HUGE *palfa2,
            Integer NAG_HUGE istate[], Integer n, Integer nctotl,
            double NAG_HUGE anorm[], double NAG_HUGE ap[], double NAG_HUGE ax[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE featol[],
            double NAG_HUGE p[], double NAG_HUGE x[], Nag_fe04nb NAG_HUGE *fe04nb);
#else
extern void e04uch();
#endif

#ifdef NAG_PROTO
extern void e04ucj(Boolean NAG_HUGE *first, Boolean debug, Boolean NAG_HUGE *done, Boolean NAG_HUGE *imprvd,
            Integer NAG_HUGE *inform, Integer maxf, Integer NAG_HUGE *numf,
            double alfmax, double alfsml, double epsaf,
            double g0, double targtg, double ftry,
            double NAG_HUGE *tolabs, double tolrel, double toltny,
            double NAG_HUGE *alfa, double NAG_HUGE *alfbst, double NAG_HUGE *fbest, Nag_UcjkSt NAG_HUGE *st);
#else
extern void e04ucj();
#endif

#ifdef NAG_PROTO
extern void e04uck(Boolean NAG_HUGE *first, Boolean debug, Boolean NAG_HUGE *done, Boolean NAG_HUGE *imprvd,
            Integer NAG_HUGE *inform, Integer maxf, Integer NAG_HUGE *numf,
            double alfmax, double epsaf, double g0,
            double targtg, double ftry, double gtry,
            double NAG_HUGE *tolabs, double tolrel, double toltny,
            double NAG_HUGE *alfa, double NAG_HUGE *alfbst, double NAG_HUGE *fbest,
            double NAG_HUGE *gbest, Nag_UcjkSt NAG_HUGE *st);
#else
extern void e04uck();
#endif

#ifdef NAG_PROTO
extern void e04ucl(char NAG_HUGE *lsumry, Boolean unitq, Integer n, Integer ncnln,
            Integer nfree, Integer nz, Integer tdcj1, Integer tdcj2,
            Integer tdzy, Integer tdr, Integer NAG_HUGE kx[], double alfa,
            double glf1, double glf2, double qpcurv,
            double NAG_HUGE cjac1[], double NAG_HUGE cjac2[], double NAG_HUGE cjdx1[],
            double NAG_HUGE cjdx2[], double NAG_HUGE cs1[], double NAG_HUGE cs2[],
            double NAG_HUGE gq1[], double NAG_HUGE gq2[], double NAG_HUGE hpq[],
            double NAG_HUGE rpq[], double NAG_HUGE qpmul[], double NAG_HUGE r[],
            double NAG_HUGE omega[], double NAG_HUGE zy[], double NAG_HUGE wrk1[],
            double NAG_HUGE wrk2[], Nag_de04uc NAG_HUGE *de04uc, Nag_ee04nb NAG_HUGE *ee04nb,
            Nag_fe04uc NAG_HUGE *fe04uc);
#else
extern void e04ucl();
#endif

#ifdef NAG_PROTO
extern void e04ucm(Boolean unitq, Integer ncqp, Integer nactiv, Integer nfree,
            Integer nz, Integer n, Integer nlnx, Integer nctotl,
            Integer tdzy, Integer tdaqp, Integer tdr, Integer tdt,
            Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[],
            double NAG_HUGE *dxnorm, double NAG_HUGE *gdx, double NAG_HUGE aqp[],
            double NAG_HUGE adx[], double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE rpq[], double NAG_HUGE rpq0[], double NAG_HUGE dx[],
            double NAG_HUGE gq[], double NAG_HUGE r[], double NAG_HUGE t[], double NAG_HUGE zy[],
            double NAG_HUGE work[], Nag_fe04uc NAG_HUGE *fe04uc);
#else
extern void e04ucm();
#endif

#ifdef NAG_PROTO
extern void e04ucn(Boolean feasqp, Integer n, Integer nclin, Integer ncnln,
            double NAG_HUGE *objalf, double NAG_HUGE *grdalf, double qpcurv,
            Integer NAG_HUGE istate[],  double NAG_HUGE cjdx[], double NAG_HUGE cmul[],
            double NAG_HUGE cs[], double NAG_HUGE dlam[], double NAG_HUGE rho[],
            double NAG_HUGE violn[], double NAG_HUGE work1[], double NAG_HUGE work2[], 
            Nag_de04uc NAG_HUGE *de04uc, Nag_fe04uc NAG_HUGE *fe04uc);
#else
extern void e04ucn();
#endif

#ifdef NAG_PROTO
extern void e04ucp(Integer n, Integer nclin, Integer ncnln, Integer nctotl,
            Integer NAG_HUGE *litotl, Integer NAG_HUGE *lwtotl, Nag_ae04nc NAG_HUGE *ae04nc,
            Nag_ae04uc NAG_HUGE *ae04uc, Nag_be04nb NAG_HUGE *be04nb);
#else
extern void e04ucp();
#endif

#ifdef NAG_PROTO
extern void e04ucr(Boolean needfd,  Integer NAG_HUGE *inform, Integer n, Integer ncnln,
            Integer tdcj, Integer tdcju, Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad,
            Integer NAG_HUGE needc[],
            NAG_E04UCC_CONFUN confun, NAG_E04UCC_FUN objfun,
            double NAG_HUGE *alfa, double NAG_HUGE *alfbnd, double NAG_HUGE *alfmax,
            double alfsml, double dxnorm, double epsrf,
            double eta, double NAG_HUGE *gdx, double grdalf,
            double glf1, double NAG_HUGE *glf, double NAG_HUGE *objf,
            double NAG_HUGE *objalf, double qpcurv, double xnorm,
            double NAG_HUGE c[], double NAG_HUGE c2[], double NAG_HUGE cjac[],
            double NAG_HUGE cjacu[], double NAG_HUGE cjdx[], double NAG_HUGE cjdx2[],
            double NAG_HUGE cmul1[], double NAG_HUGE cmul[], double NAG_HUGE cs1[],
            double NAG_HUGE cs[], double NAG_HUGE dx[], double NAG_HUGE dlam[],
            double NAG_HUGE dslk[], double NAG_HUGE grad[], double NAG_HUGE gradu[],
            double NAG_HUGE qpmul[], double NAG_HUGE rho[], double NAG_HUGE slk1[],
            double NAG_HUGE slk[], double NAG_HUGE x1[], double NAG_HUGE x[],
            double NAG_HUGE work[], double NAG_HUGE w[],
            Nag_de04uc NAG_HUGE *de04uc, Nag_fe04uc NAG_HUGE *fe04uc, Nag_Comm NAG_HUGE *comm);
#else
extern void e04ucr();
#endif

#ifdef NAG_PROTO
extern void e04ucs(Boolean cold,  Integer n, Integer nclin, Integer ncnln,
            Integer nctotl, Integer NAG_HUGE *nactiv, Integer nfree, Integer NAG_HUGE *nz,
            Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[], double bigbnd,
            double tolact, double NAG_HUGE bl[], double NAG_HUGE bu[],
            double NAG_HUGE c[], Nag_fe04uc NAG_HUGE *fe04uc);
#else
extern void e04ucs();
#endif

#ifdef NAG_PROTO
extern void e04uct(Boolean NAG_HUGE ktcond[], Boolean convrg, char NAG_HUGE *lsumry,
            Integer tdr, Integer tdt,
            Integer n, Integer nclin, Integer ncnln, Integer nctotl,
            Integer nactiv, Integer linact, Integer nlnact, Integer nz,
            Integer nfree, Integer majit0, Integer majits, Integer minits,
            Integer NAG_HUGE istate[], double alfa, Integer nfun,
            double condhz, double condh, double condt,
            double objalf, double objf, double gfnorm,
            double gznorm, double cvnorm, double NAG_HUGE ax[],
            double NAG_HUGE c[], double NAG_HUGE r[], double NAG_HUGE t[], double NAG_HUGE violn[],
            double NAG_HUGE x[], double NAG_HUGE work[], Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm, 
            Nag_Search_State NAG_HUGE *st, Nag_de04uc NAG_HUGE *de04uc,
            Nag_FileSt NAG_HUGE *stream);
#else
extern void e04uct();
#endif

#ifdef NAG_PROTO
extern void e04ucu(Boolean NAG_HUGE *feasqp, Boolean NAG_HUGE *unitq, Integer NAG_HUGE *nqperr,
            Integer NAG_HUGE *majits, Integer NAG_HUGE *minits, Integer n, Integer nclin,
            Integer ncnln, Integer tdcj, Integer tdaqp, Integer tdr,
            Integer NAG_HUGE *linact, Integer NAG_HUGE *nlnact, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree,
            Integer NAG_HUGE *nz, Integer NAG_HUGE *numinf, Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[],
            Integer NAG_HUGE kx[], double NAG_HUGE *dxnorm, double NAG_HUGE *gdx,
            double NAG_HUGE *qpcurv, double NAG_HUGE aqp[], double NAG_HUGE adx[],
            double NAG_HUGE anorm[], double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE cjac[],
            double NAG_HUGE clamda[], double NAG_HUGE cmul[], double NAG_HUGE cs[],
            double NAG_HUGE dlam[], double NAG_HUGE dslk[], double NAG_HUGE dx[],
            double NAG_HUGE qpbl[], double NAG_HUGE qpbu[], double NAG_HUGE qptol[],
            double NAG_HUGE r[], double NAG_HUGE rho[], double NAG_HUGE slk[],
            double NAG_HUGE violn[], double NAG_HUGE x[], double NAG_HUGE wtinf[],
            double NAG_HUGE w[], Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm, 
            Nag_FileSt NAG_HUGE *stream, Nag_Search_State NAG_HUGE *st, NagError NAG_HUGE *ovflow,
            Nag_ae04nc NAG_HUGE *ae04nc, Nag_be04nb NAG_HUGE *be04nb, Nag_ce04nc NAG_HUGE *ce04nc,
            Nag_de04nb NAG_HUGE *de04nb, Nag_de04uc NAG_HUGE *de04uc, Nag_fe04nb NAG_HUGE *fe04nb,
            Nag_fe04uc NAG_HUGE *fe04uc, Integer NAG_HUGE *local_error);
#else
extern void e04ucu();
#endif

#ifdef NAG_PROTO
extern void e04ucw(Integer n, Integer nclin, Integer ncnln, Integer NAG_HUGE istate[],
            double bigbnd, double NAG_HUGE *cvnorm, double NAG_HUGE *errmax,
            Integer NAG_HUGE *jmax, Integer NAG_HUGE *nviol,  double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE featol[],
            double NAG_HUGE x[], double NAG_HUGE work[], Nag_fe04uc NAG_HUGE *fe04uc);
#else
extern void e04ucw();
#endif

#ifdef NAG_PROTO
extern void e04ucy(Nag_DerivSet NAG_HUGE *deriv_level,
            Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad, Integer tdcj, Integer tdcju,
            Integer n, Integer ncnln,
            NAG_E04UCC_CONFUN confun, NAG_E04UCC_FUN objfun,
            Integer NAG_HUGE needc[], double bigbnd,
            double epsrf, double cdint, double fdint,
            double fdchk, double NAG_HUGE *fdnorm, double NAG_HUGE *objf,
            double xnorm, double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[],
            double NAG_HUGE c1[], double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE cjdx[],
            double NAG_HUGE dx[], double NAG_HUGE grad[], double NAG_HUGE gradu[], double NAG_HUGE hforwd[],
            double NAG_HUGE hcntrl[], double NAG_HUGE x[], double NAG_HUGE wrk1[], double NAG_HUGE wrk2[],
            Nag_Comm NAG_HUGE *comm, Nag_Deriv_Inf NAG_HUGE *diff,
            const Nag_DebugSt NAG_HUGE *debug, Nag_Grad_Chk_St g_chk, int NAG_HUGE *error,
            Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04ucy();
#endif

#ifdef NAG_PROTO
extern void e04ucz(Boolean NAG_HUGE *unitq, Integer NAG_HUGE *inform, Integer NAG_HUGE *majits,
            Integer n, Integer nclin, Integer ncnln,
            Integer nctotl, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree,
            Integer NAG_HUGE *nz, Integer tdcj, Integer tdcju, Integer tdaqp,
            Integer tdr, Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad, Integer NAG_HUGE istate[],
            Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[], double NAG_HUGE *objf,
            double NAG_HUGE *fdnorm, double NAG_HUGE *xnorm,
            NAG_E04UCC_CONFUN confun, NAG_E04UCC_FUN objfun,
            double NAG_HUGE aqp[], double NAG_HUGE ax[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[],
            double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE clamda[],
            double NAG_HUGE featol[], double NAG_HUGE grad[], double NAG_HUGE gradu[],
            double NAG_HUGE r[], double NAG_HUGE x[],  Integer NAG_HUGE iw[], double NAG_HUGE w[],
            Nag_Deriv_Inf NAG_HUGE *diff, Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm,
            Nag_FileSt NAG_HUGE *stream, Nag_Search_State NAG_HUGE *st, NagError NAG_HUGE *ovflow,
            Nag_ae04nc NAG_HUGE *ae04nc, Nag_ae04uc NAG_HUGE *ae04uc,
            Nag_be04nb NAG_HUGE *be04nb, Nag_ce04nc NAG_HUGE *ce04nc,
            Nag_de04nb NAG_HUGE *de04nb, Nag_de04uc NAG_HUGE *de04uc,
            Nag_ee04nb NAG_HUGE *ee04nb, Nag_fe04nb NAG_HUGE *fe04nb,
            Nag_fe04uc NAG_HUGE *fe04uc, Integer NAG_HUGE *local_error);
#else
extern void e04ucz();
#endif

#ifdef NAG_PROTO
extern void e04udr(Boolean unitq, Integer n, Integer nfree, Integer nz,
            Integer nq, Integer nrowr, Integer NAG_HUGE iperm[], Integer NAG_HUGE kx[],
            double NAG_HUGE gq[], double NAG_HUGE r[], double NAG_HUGE zy[],
            double NAG_HUGE work[], double NAG_HUGE qrwork[], Nag_ee04nb NAG_HUGE *ee04nb);
#else
extern void e04udr();
#endif

#ifdef NAG_PROTO
extern void e04uds(Boolean central, Integer ncolj,
            Integer ncoluj, Integer n, Integer ncnln,
            double bigbnd, double cdint, double fdint,
            double NAG_HUGE *fdnorm, double objf,
            NAG_E04UCC_CONFUN confun, NAG_E04UCC_FUN objfun,
            Integer NAG_HUGE needc[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE c1[],
            double NAG_HUGE c2[], double NAG_HUGE cjac[], double NAG_HUGE ujac[], double NAG_HUGE grad[],
            double NAG_HUGE ugrad[], double NAG_HUGE hforwd[], double NAG_HUGE hcntrl[], double NAG_HUGE x[],
            Nag_Comm NAG_HUGE *comm, Nag_Deriv_Inf diff);
#else
extern void e04uds();
#endif

#ifdef NAG_PROTO
extern void e04udt(Integer NAG_HUGE *inform, Integer n, Integer nclin, Integer ncnln,
            double NAG_HUGE *alfa, double alfmin, double alfmax,
            double bigbnd, double dxnorm, double NAG_HUGE anorm[],
            double NAG_HUGE adx[], double NAG_HUGE ax[], double NAG_HUGE bl[],
            double NAG_HUGE bu[], double NAG_HUGE dslk[], double NAG_HUGE dx[],
            double NAG_HUGE slk[], double NAG_HUGE x[], Nag_fe04nb NAG_HUGE *fe04nb);
#else
extern void e04udt();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL 
  e04unc(Integer m, Integer n, Integer nclin, Integer ncnlin, 
	 double NAG_HUGE a[], Integer tda, double NAG_HUGE bl[], 
	 double NAG_HUGE bu[], double NAG_HUGE y[], NAG_E04UNC_OBJFUN objfun, 
	 NAG_E04UNC_CONFUN confun, double NAG_HUGE x[], double NAG_HUGE *objf,
	 double NAG_HUGE f[], double NAG_HUGE fjac[], Integer tdfju, 
	 Nag_E04_Opt NAG_HUGE *options, Nag_Comm NAG_HUGE *user_comm,
	 NagError NAG_HUGE *fail);
#else
extern void e04unc();
#endif

#ifdef NAG_PROTO
void e04unp(Integer mode, Integer n, Integer m,  double NAG_HUGE y[],
	    double NAG_HUGE f[], double NAG_HUGE fjac[],  Integer tdfj,
	    double NAG_HUGE *objf, double NAG_HUGE grad[], double NAG_HUGE yf[]);
#else
extern void e04unp();
#endif

#ifdef NAG_PROTO
void e04unt(Boolean needfd,  Integer NAG_HUGE *inform, Integer m, Integer n,
	    Integer ncnln, Integer tdcj, Integer tdcju, Integer tdfj,
	    Integer tdfju, Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad, Integer NAG_HUGE needc[],
	    NAG_E04UNC_CONFUN confun, NAG_E04UNC_OBJFUN objfun,  double NAG_HUGE *alfa,
	    double NAG_HUGE *alfbnd, double NAG_HUGE *alfmax, double alfsml,
	    double dxnorm, double epsrf, double eta,
	    double NAG_HUGE *gdx, double grdalf, double gl1, double NAG_HUGE *gl,
	    double NAG_HUGE *objf, double NAG_HUGE *objalf, double curvqp,
	    double xnorm, double NAG_HUGE c[], double NAG_HUGE cjac[],
	    double NAG_HUGE cjacu[], double NAG_HUGE cjdx[], double NAG_HUGE cmul1[],
	    double NAG_HUGE cmul[], double NAG_HUGE cs1[], double NAG_HUGE cs[],
	    double NAG_HUGE dx[], double NAG_HUGE dlam[], double NAG_HUGE dslk[],
	    double NAG_HUGE y[], double NAG_HUGE f[], double NAG_HUGE fjac[],
	    double NAG_HUGE fjacu[], double NAG_HUGE grad[], double NAG_HUGE gradu[],
	    double NAG_HUGE qpmul[], double NAG_HUGE rho[], double NAG_HUGE slk1[],
	    double NAG_HUGE slk[], double NAG_HUGE x1[], double NAG_HUGE x[], double NAG_HUGE w[],
	    Nag_ae04nc NAG_HUGE *ae04nc, Nag_ae04uc NAG_HUGE *ae04uc, Nag_ae04up NAG_HUGE *ae04up,
	    Nag_de04uc NAG_HUGE *de04uc, Nag_Comm NAG_HUGE *comm);
#else
extern void e04unt();
#endif

#ifdef NAG_PROTO
void e04unv(Integer m, Integer n, Integer NAG_HUGE *litotl, Integer NAG_HUGE *lwtotl, 
	    Nag_ae04up NAG_HUGE *ae04up);
#else
extern void e04unv();
#endif

#ifdef NAG_PROTO
void e04unz(Boolean NAG_HUGE *unitq,  Integer NAG_HUGE *inform,
	    Integer NAG_HUGE *majits, Integer m, Integer n, Integer nclin,
	    Integer ncnln, Integer NAG_HUGE *nactiv, Integer NAG_HUGE *nfree, Integer NAG_HUGE *nz,
	    Integer tdcj, Integer tdcju, Integer tdfj, Integer tdfju,
	    Integer tdaqp, Integer tdr, Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad,
	    Integer NAG_HUGE istate[], Integer NAG_HUGE kactiv[], Integer NAG_HUGE kx[],
	    double NAG_HUGE *objf, double NAG_HUGE *fdnorm, double NAG_HUGE *xnorm,
	    NAG_E04UNC_CONFUN confun, NAG_E04UNC_OBJFUN objfun,  double NAG_HUGE aqp[],
	    double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[],
	    double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE clamda[],
	    double NAG_HUGE y[], double NAG_HUGE f[], double NAG_HUGE fjac[],
	    double NAG_HUGE fjacu[], double NAG_HUGE featol[], double NAG_HUGE grad[],
	    double NAG_HUGE r[], double NAG_HUGE x[],  Integer NAG_HUGE iw[],  double NAG_HUGE w[],
	    Nag_Deriv_Inf NAG_HUGE *diff, Nag_E04_Opt NAG_HUGE *opt, Nag_Comm NAG_HUGE *comm,
	    Nag_FileSt NAG_HUGE *stream, Nag_Search_State NAG_HUGE *st, NagError NAG_HUGE *ovflow,
	    Nag_ae04nc NAG_HUGE *ae04nc, Nag_ae04uc NAG_HUGE *ae04uc, Nag_ae04up NAG_HUGE *ae04up, 
	    Nag_be04nb NAG_HUGE *be04nb, Nag_ce04nc NAG_HUGE *ce04nc, Nag_de04nb NAG_HUGE *de04nb, 
	    Nag_de04uc NAG_HUGE *de04uc, Nag_ee04nb NAG_HUGE *ee04nb, Nag_fe04nb NAG_HUGE *fe04nb, 
	    Nag_fe04uc NAG_HUGE *fe04uc, int NAG_HUGE *fail_code);
#else
extern void e04unz();
#endif

#ifdef NAG_PROTO
void e04upg(Nag_Grad_Chk_St g_chk, Integer m, Integer n, double NAG_HUGE x[], 
	    double NAG_HUGE fjac[], Integer tdfj, Integer ncnln, double NAG_HUGE cjac[], 
	    Integer tdcj, const Nag_GPrintSt NAG_HUGE *gprint, Nag_FileSt NAG_HUGE *stream);
#else
extern void e04upg();
#endif

#ifdef NAG_PROTO
void e04upq(Nag_DerivSet deriv_level,
	    Integer m, Integer n, Integer tdfjac, Integer tdfjacu,
	    double bigbnd, double epsrf, double oktol,
	    double fdchk, double xnorm, NAG_E04UNC_OBJFUN objfun,
	    double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE f[], double NAG_HUGE f1[],
	    double NAG_HUGE fjac[], double NAG_HUGE fjacu[], double NAG_HUGE fjdx[],
	    double NAG_HUGE dx[], double NAG_HUGE err[], double NAG_HUGE x[], double NAG_HUGE y[],
	    Nag_Grad_Chk_St g_chk, Nag_Comm NAG_HUGE *comm, int NAG_HUGE *fj_error,
	    Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04upq();
#endif

#ifdef NAG_PROTO
void e04upr(Nag_DerivSet NAG_HUGE *deriv_level, 
	    Integer m, Integer n, Integer ncnln, Integer tdcj, Integer tdcju,
	    Integer tdfj, Integer tdfju,  double bigbnd, double epsrf,
	    double NAG_HUGE *fdnorm, NAG_E04UNC_CONFUN confun, NAG_E04UNC_OBJFUN objfun,
	    Integer NAG_HUGE needc[],  double NAG_HUGE bl[], double NAG_HUGE bu[],
	    double NAG_HUGE c[], double NAG_HUGE c1[], double NAG_HUGE c2[],
	    double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE f[],
	    double NAG_HUGE f1[], double NAG_HUGE f2[], double NAG_HUGE fjac[],
	    double NAG_HUGE fjacu[], double NAG_HUGE hforwd[], double NAG_HUGE hcntrl[],
	    double NAG_HUGE x[], double NAG_HUGE y[], Nag_Deriv_Inf NAG_HUGE *diff, 
	    Nag_Comm NAG_HUGE *comm, Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04upr();
#endif

#ifdef NAG_PROTO
void e04ups(Integer n, Integer tdzy, Integer nfree, Boolean unitq,
	    Integer NAG_HUGE kx[], Integer m, double NAG_HUGE fjac[], Integer tdfj,
	    double NAG_HUGE r[], Integer tdr, double NAG_HUGE zy[], double NAG_HUGE work[],
	    NagError NAG_HUGE *fail_qdc);
#else
extern void e04ups();
#endif

#ifdef NAG_PROTO
void e04upw(Boolean centrl, Integer tdcj, Integer tdcju,
	    Integer tdfj, Integer tdfju, Integer m, Integer n, Integer ncnln,
	    double bigbnd, double cdint, double fdint,
	    double NAG_HUGE *fdnorm, NAG_E04UNC_CONFUN confun, NAG_E04UNC_OBJFUN objfun,
	    Integer NAG_HUGE needc[], double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE c1[], 
	    double NAG_HUGE c2[], double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE f[], double NAG_HUGE f1[],
	    double NAG_HUGE f2[], double NAG_HUGE fjac[], double NAG_HUGE fjacu[], double NAG_HUGE hforwd[], 
	    double NAG_HUGE hcntrl[], double NAG_HUGE x[], Nag_Comm NAG_HUGE *comm, Nag_Deriv_Inf NAG_HUGE *diff);
#else
void e04upw();
#endif

#ifdef NAG_PROTO
void e04upy(Nag_DerivSet NAG_HUGE *deriv_level,
	    Integer NAG_HUGE *nfun, Integer NAG_HUGE *ngrad, Integer tdcjac, Integer tdcjacu,
	    Integer tdfjac, Integer tdfjacu, Integer m, Integer n, Integer ncnln,
	    NAG_E04UNC_CONFUN confun, NAG_E04UNC_OBJFUN objfun,
            Integer NAG_HUGE needc[], double bigbnd, 
	    double epsrf, double cdint, double fdint, 
	    double fdchk, double NAG_HUGE *fdnorm,     
	    double xnorm, double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[],
	    double NAG_HUGE c1[], double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE cjdx[], 
	    double NAG_HUGE f[], double NAG_HUGE f1[], double NAG_HUGE fjac[], double NAG_HUGE fjacu[], 
	    double NAG_HUGE fjdx[], double NAG_HUGE dx[], double NAG_HUGE hforwd[], 
	    double NAG_HUGE hcntrl[], double NAG_HUGE x[], double NAG_HUGE wrk1[], double NAG_HUGE wrk2[],
            double NAG_HUGE wrk4[], Nag_Comm NAG_HUGE *comm, Nag_Deriv_Inf NAG_HUGE *diff,
	    Nag_Grad_Chk_St g_chk, int NAG_HUGE *error, Nag_GPrintSt NAG_HUGE *gprint);
#else
void e04upy();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL 
  e04xac(Integer n, double x[], NAG_E04UCC_FUN objfun, double *objf, 
	 double g[], double h_forward[], double h_central[],
	 double h[], Integer tdh, Nag_DerivInfo deriv_info[],
	 Nag_E04_Opt *options, Nag_Comm *comm, NagError *fail);
#else
extern void e04xac();
#endif

#ifdef NAG_PROTO
extern void e04xap(Nag_Grad_Chk_St g_chk, Integer n, double NAG_HUGE x[], double NAG_HUGE grad[],
            Integer ncnln, double NAG_HUGE cjac[], Integer tdcj,
            const Nag_GPrintSt NAG_HUGE *gprint, Nag_FileSt NAG_HUGE *stream);
#else
extern void e04xap();
#endif

#ifdef NAG_PROTO
extern void e04xaq(Nag_Grad_Chk_St g_chk, Integer n, double NAG_HUGE x[], double NAG_HUGE grad[],
            Integer ncnln, double NAG_HUGE cjac[], Integer tdcj,
            const Nag_GPrintSt NAG_HUGE *gprint, Nag_FileSt NAG_HUGE *stream);
#else
extern void e04xaq();
#endif

#ifdef NAG_PROTO
extern void e04xaw(Nag_DerivSet deriv_level,
            Integer n, Integer ncnln, Integer tdcj,
            Integer tdcju, double bigbnd, double epsrf,
            double oktol, double fdchk, double xnorm,
            NAG_E04UCC_CONFUN confun,
            Integer NAG_HUGE needc[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE c1[],
            double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE cjdx[],
            double NAG_HUGE dx[], double NAG_HUGE err[], double NAG_HUGE x[], double NAG_HUGE y[],
            Nag_Grad_Chk_St g_chk, const Nag_DebugSt NAG_HUGE *debug,
            Nag_Comm NAG_HUGE *comm, int NAG_HUGE *j_error, Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04xaw();
#endif

#ifdef NAG_PROTO
extern void e04xax(Integer n,
            double bigbnd, double epsrf, double oktol,
            double fdchk, double NAG_HUGE *objf, double xnorm,
	    NAG_E04UCC_FUN objfun,
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE grad[], double NAG_HUGE gradu[],
            double NAG_HUGE dx[], double NAG_HUGE x[], double NAG_HUGE y[],
            Nag_Grad_Chk_St g_chk,
            const Nag_DebugSt NAG_HUGE *debug, Nag_Comm NAG_HUGE *comm, int NAG_HUGE *g_error,
            Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04xax();
#endif

#ifdef NAG_PROTO
extern void e04xay(Nag_DerivSet NAG_HUGE *deriv_level,
            Integer n, Integer ncnln, Integer tdcj, Integer tdcju,
            double bigbnd, double epsrf, double NAG_HUGE *fdnorm, double objf,
            NAG_E04UCC_CONFUN confun, NAG_E04UCC_FUN objfun,
            Integer NAG_HUGE needc[],
            double NAG_HUGE bl[], double NAG_HUGE bu[], double NAG_HUGE c[], double NAG_HUGE c1[], double NAG_HUGE c2[],
            double NAG_HUGE cjac[], double NAG_HUGE cjacu[], double NAG_HUGE grad[], double NAG_HUGE gradu[],
            double NAG_HUGE hforwd[], double NAG_HUGE hcntrl[], double NAG_HUGE x[], double NAG_HUGE y[],
            Nag_Deriv_Inf NAG_HUGE *diff, const Nag_DebugSt NAG_HUGE *debug, Nag_Comm NAG_HUGE *comm,
            Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04xay();
#endif

#ifdef NAG_PROTO
extern void e04xaz(Boolean debug, Boolean NAG_HUGE *done, Boolean NAG_HUGE *first, double epsa,
            double epsr, double fx, Integer NAG_HUGE *inform,
            Integer NAG_HUGE *iter, Integer itmax, double NAG_HUGE *cdest,
            double NAG_HUGE *fdest, double NAG_HUGE *sdest, double NAG_HUGE *errbnd,
            double f1, double f2, double NAG_HUGE *h, double NAG_HUGE *hopt,
            double NAG_HUGE *hphi, Nag_XazSt NAG_HUGE *st);
#else
extern void e04xaz();
#endif

#ifdef NAG_PROTO
extern Boolean e04xbp(Integer n, Integer ncnln, Nag_GPrintSt NAG_HUGE *gprint);
#else
extern Boolean e04xbp();
#endif

#ifdef NAG_PROTO
Boolean e04xbq(Integer m, Integer n, Integer ncnln, Nag_GPrintSt NAG_HUGE *gprint);
#else
extern Boolean e04xbq();
#endif

#ifdef NAG_PROTO
extern void e04xcp(Nag_GPrintSt NAG_HUGE *gprint);
#else
extern void e04xcp();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04xxc(Nag_E04_Opt NAG_HUGE *opt);
#else
extern void e04xxc();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxj(char str[], int nc, int *field_code);
#else
extern Boolean e04xxj();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxk(char str[], int nc, int *field_code);
#else
extern Boolean e04xxk();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxl(char str[], int nc, int *field_code);
#else
extern Boolean e04xxl();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxm(char str[], int nc, int *field_code);
#else
extern Boolean e04xxm();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxn(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxn();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxp(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxp();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxq(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxq();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxr(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxr();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxs(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxs();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxt(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxt();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxu(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxu();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxv(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxv();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxw(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxw();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxx(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxx();
#endif

#ifdef NAG_PROTO
extern Boolean e04xxy(char NAG_HUGE str[], int nc, int NAG_HUGE *field_code);
#else
extern Boolean e04xxy();
#endif

#ifdef NAG_PROTO
extern void e04xya(int field_code, Nag_E04_Opt NAG_HUGE *options, char NAG_HUGE buf[]);
#else
extern void e04xya();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04xyc(const char NAG_HUGE *name, 
	const char NAG_HUGE *opt_file, Nag_E04_Opt NAG_HUGE *opt,
	Boolean print, const char NAG_HUGE *outfile, NagError NAG_HUGE *fail);
#else
extern void e04xyc();
#endif

#ifdef NAG_PROTO
extern void e04xyj(int field_code, Nag_E04_Opt *options,
            Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);
#else
extern void e04xyj();
#endif

#ifdef NAG_PROTO
extern void e04xyk(int field_code, Nag_E04_Opt *options,
            Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);
#else
extern void e04xyk();
#endif

#ifdef NAG_PROTO
extern void e04xyl(int field_code, Nag_E04_Opt *options,
            Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);
#else
extern void e04xyl();
#endif

#ifdef NAG_PROTO
extern void e04xym(int field_code, Nag_E04_Opt *options,
            Nag_FileSt *stream, Nag_Mesg *mesg, NagError *fail);
#else
extern void e04xym();
#endif

#ifdef NAG_PROTO
extern void e04xyn(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyn();
#endif

#ifdef NAG_PROTO
extern void e04xyp(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyp();
#endif

#ifdef NAG_PROTO
extern void e04xyq(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyq();
#endif

#ifdef NAG_PROTO
extern void e04xyr(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyr();
#endif

#ifdef NAG_PROTO
extern void e04xys(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xys();
#endif

#ifdef NAG_PROTO
extern void e04xyt(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyt();
#endif

#ifdef NAG_PROTO
extern void e04xyu(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyu();
#endif

#ifdef NAG_PROTO
extern void e04xyv(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyv();
#endif

#ifdef NAG_PROTO
extern void e04xyw(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyw();
#endif

#ifdef NAG_PROTO
extern void e04xyx(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyx();
#endif

#ifdef NAG_PROTO
extern void e04xyy(int field_code, Nag_E04_Opt NAG_HUGE *options,
            Nag_FileSt NAG_HUGE *stream, Nag_Mesg NAG_HUGE *mesg, NagError NAG_HUGE *fail);
#else
extern void e04xyy();
#endif

#ifdef NAG_PROTO
extern void e04xyz(Boolean (NAG_HUGE *valid_field)(char NAG_HUGE *str, int nc, int NAG_HUGE *field_code),
            int NAG_HUGE *state, int NAG_HUGE *i, char NAG_HUGE line[], Integer NAG_HUGE *linenum,
            int NAG_HUGE *field_code, FILE NAG_HUGE *fp, Nag_E04_Opt NAG_HUGE *opt, char NAG_HUGE str[],
            Nag_Opt_Found NAG_HUGE *found);
#else
extern void e04xyz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04xzc(Nag_E04_Opt NAG_HUGE *opt, char NAG_HUGE *name, NagError NAG_HUGE *fail);
#else
extern void e04xzc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04yac(Integer m, Integer n, NAG_E04YAC_FUN lsqfun,
            double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdj,
            Nag_Comm NAG_HUGE *user_comm, NagError NAG_HUGE *fail);
#else
extern void e04yac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e04ycc(Integer job, Integer m, Integer n, double fsumsq,
            double NAG_HUGE cj[], Nag_E04_Opt NAG_HUGE *options, NagError NAG_HUGE *fail);
#else
extern void e04ycc();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGE04 */
