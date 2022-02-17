#ifndef NAGF11
#define NAGF11
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf11.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f11 Chapter
 *
 * Mark 5, 1997.
 */

#ifdef NAG_PROTO
extern void f11baz(Integer action, Integer idata[], double rdata[], Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11baz();
#endif

#ifdef NAG_PROTO
extern double f11bbu(Integer norm, Integer n,  double x[], double wgt[], Nag_Sparse_Comm *comm);
#else
extern double f11bbu();
#endif

#ifdef NAG_PROTO
extern void f11bby(Boolean *next,  Integer *irevcm, Integer norm,
             double *anorm,  Integer n,  double work[],
             double u[], double v[], Nag_Sparse_Comm *comm);
#else
extern void f11bby();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void f11gac(Nag_SparseSym_Method method, Nag_SparseSym_PrecType precon, 
             Nag_SparseSym_Bisection sigcmp, Nag_SparseSym_Norm norm,
             Nag_SparseSym_Weight weight,  Nag_SparseSym_Term iterm, Integer n,  double tol,
             Integer maxitn,  double anorm, double sigmax,
             double sigtol,  Integer maxits, Integer monit,
             Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11gac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11gbc(Integer *irevcm,  double u[], double v[], double work[],
             Integer lwork, Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11gbc();
#endif

#ifdef NAG_PROTO
extern void f11gbp(Integer resid,  Boolean usexc,  Integer iterm, Integer n,
             double zetab, double bpigam, double beta1,
             double beta2, double rho1, double rho2, double x[],
             double v1[], double v2[], double w1[], double u[],
             double v[]);
#else
extern void f11gbp();
#endif


#ifdef NAG_PROTO
extern void f11gbq(Boolean *next,  Integer *irevcm, Integer itn, Integer *resid,
             Nag_SparseSym_PrecType precon,  Integer norm, Integer iterm, Integer n,
             double tol, double sigmax, double anorm,
             double xnorm0, double bnorm, double bnorm2,
             double rnorm2, double betain, double beta1,
             double beta2, double pi, double gammab,
             double rho1, double rho2, double zeta,
             double zetab, double b[], double x[], double v1[],
             double v2[], double w1[], double wgt[], double u[],
             double v[], double *stplhs, double *stprhs,
             Boolean *usexc, Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern void f11gbq();
#endif

#ifdef NAG_PROTO
extern void f11gbr(Boolean *next, Boolean *floop,  Integer *irevcm, Integer *itn,
             Nag_SparseSym_PrecType precon,  Integer n,  double *alpha,
             double *beta1, double *beta2, double *pi,
             double *gammab, double *rho1, double *rho2,
             double *zeta, double *zetab, double x[],
             double v1[], double v2[], double w1[], double w2[],
             double u[], double v[], Nag_Sparse_Comm *comm,  Integer *info);
#else
extern void f11gbr();
#endif

#ifdef NAG_PROTO
extern void f11gbs(Boolean *next,  Integer *irevcm,  Nag_SparseSym_PrecType precon,  Integer n,
             double *alpha, double *beta1, double *bnorm2,
             double *rnorm2, double x[], double v1[],
             double v2[], double w1[], double w2[], double u[],
             double v[], double *beta2, double *zetab, Nag_Sparse_Comm *comm,
             Integer *info);
#else
extern void f11gbs();
#endif

#ifdef NAG_PROTO
extern Integer f11gbt(Integer n,  double eps, double d[], double e2[],
                double x);
#else
extern Integer f11gbt();
#endif

#ifdef NAG_PROTO
extern void f11gbu(Boolean *next,  Integer *irevcm, Integer itn,  Nag_SparseSym_PrecType precon,
             Integer norm, Integer iterm,  double tol, double anorm,
             double bnorm, double xnorm, double sigmax,
             Integer n,  double b[], double x[], double wgt[],
             double u[], double v[], double *stplhs,
             double *stprhs, Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern void f11gbu();
#endif

#ifdef NAG_PROTO
extern void f11gbv(Boolean *next, Boolean *floop,  Integer *irevcm, Integer *itn,
             Integer maxitn, Integer iterm,  Nag_SparseSym_PrecType precon,  Integer norm,
             Integer n,  double x[], double r[], double p[],
             double w[], double wgt[], double u[], double v[],
             double *talpha, double *tbeta, double *xnorm,
             double *stplhs, Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11gbv();
#endif

#ifdef NAG_PROTO
extern void f11gbw(Integer *sigcmx, Integer *its,  double talpha,
             double tbeta, double d[], double e2[],
             double *sigmax, Nag_Sparse_Comm *comm, double *sigerr);
#else
extern void f11gbw();
#endif

#ifdef NAG_PROTO
extern void f11gbx(Boolean *next,  Integer *irevcm,  Nag_SparseSym_PrecType precon,
             Integer norm, Integer iterm,  double *bnorm,
             double *xnorm,  Integer n,  double b[], double x[],
             double r[], double p[], double w[], double wgt[],
             double u[], double v[], double *talpha,
             double *tbeta, double *stplhs, Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern void f11gbx();
#endif

#ifdef NAG_PROTO
extern void f11gby(Integer *irevcm, Integer idata[],  double rdata[],
             double u[], double v[], double work[],  Integer lwork,
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11gby();
#endif

#ifdef NAG_PROTO
extern void f11gbz(Integer *irevcm, Integer idata[],  double rdata[],
             double u[], double v[], double work[],  Integer lwork,
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11gbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP  void NAG_CALL f11gcc(Integer *itn,  double *stplhs, double *stprhs,
             double *anorm, double *sigmax,  Integer *its,
             double *sigerr, Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11gcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11jac(Integer n, Integer nnz,  double **a,  Integer *la,
             Integer **irow, Integer **icol, Integer lfill,  double dtol,
             Nag_SparseSym_Fact mic,  double dscale,  Nag_SparseSym_Piv pstrat,
             Integer ipiv[], Integer istr[], Integer *nnzc, Integer *npivm,
             Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11jac();
#endif

#ifdef NAG_PROTO
extern void f11jaw(Integer n, Integer *nnz,  double a[],  Integer la,
             Integer irow[], Integer icol[],  double dscale,
             Integer iwork[], Integer *ierror);
#else
extern void f11jaw();
#endif

#ifdef NAG_PROTO
extern void f11jax(Integer m, Integer n, Integer ind[], Integer ist[],
             Integer iwork[]);
#else
extern void f11jax();
#endif

#ifdef NAG_PROTO
extern void f11jay(Integer *la, Integer offset, Integer n, Integer nnz, Integer *nnzc, Integer maxf,
             double **a,  Integer **irow, Integer **icol, Integer lfill,
             double dtol, double dscale, double alpha,
             Nag_SparseSym_Piv pstrat,  Integer ipiv[], Integer *npivm, Integer nnzr[],
             Integer **levf, Integer idlevf, Integer ir[], Integer **ic,
             Integer istr[], Integer istc[], Integer istll[], Integer llnnzf[],
             Integer llnnzb[], Integer iwork[], Integer *ierror);
#else
extern void f11jay();
#endif

#ifdef NAG_PROTO
extern void f11jaz(Integer n, Integer nnz, Integer irow[], Integer icol[],
             Boolean sym,  Integer *ibad, Integer *ierror);
#else
extern void f11jaz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11jbc(Integer n,  double a[],  Integer la, Integer irow[],
             Integer icol[], Integer ipiv[], Integer istr[],  Nag_SparseSym_CheckData check,
             double y[], double x[],  NagError *fail);
#else
extern void f11jbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11jcc(Nag_SparseSym_Method method,  Integer n, Integer nnz,  double a[],
             Integer la, Integer irow[], Integer icol[], Integer ipiv[],
             Integer istr[],  double b[], double tol,
             Integer maxitn,  double x[], double *rnorm,
             Integer *itn, Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11jcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11jdc(Integer n, Integer nnz,  double a[],  Integer irow[],
             Integer icol[],  double rdiag[], double omega,
             Nag_SparseSym_CheckData check,  double y[], double x[],  Integer iwork[],
             NagError *fail);
#else
extern void f11jdc();
#endif

#ifdef NAG_PROTO
extern void f11jdz(Integer n, Integer nnz,  double a[],  Integer icol[],
             Integer istr[],  double rdiag[], double omega,
             double y[], double x[]);
#else
extern void f11jdz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11jec(Nag_SparseSym_Method method, Nag_SparseSym_PrecType precon, 
             Integer n, Integer nnz,
             double a[],  Integer irow[], Integer icol[],
             double omega, double b[], double tol,
             Integer maxitn,  double x[], double *rnorm,
             Integer *itn, Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11jec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11xec(Integer n, Integer nnz,  double a[],  Integer irow[],
             Integer icol[],  Nag_SparseSym_CheckData check,  double x[], double y[],
             NagError *fail);
#else
extern void f11xec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11zbc(Integer n, Integer *nnz,  double a[],  Integer irow[],
             Integer icol[],  Nag_SparseSym_Dups dup, Nag_SparseSym_Zeros zero, Integer istr[],
             NagError *fail);
#else
extern void f11zbc();
#endif

#ifdef NAG_PROTO
extern void f11zby(Integer n, Integer nnz,  double a[],  Integer irow[],
             Integer icol[], Integer istc[], Integer iwork[], Integer *ierr);
#else
extern void f11zby();
#endif

#ifdef NAG_PROTO
extern void f11zbz(Integer n, Integer nnz,  double a[],  Integer irow[],
             Integer icol[], Integer istr[], Integer iwork[], Integer *ierr);
#else
extern void f11zbz();
#endif


#ifdef NAG_PROTO
extern void f11ddz(Integer itrans, Integer n, Integer nnz,  double a[],
             Integer irow[], Integer icol[], Integer istr[], Integer idiag[],
             double rdiag[], double omega, double y[],
             double x[]);
#else
extern     void f11ddz();
#endif




#ifdef NAG_PROTO
extern void f11ddf(char *trans,  Integer n, Integer nnz,  double a[],
             Integer irow[], Integer icol[],  double rdiag[],
             double omega,  char *check,  double y[], double x[],
             Integer iwork[], Integer *ifail);
#else
extern     void f11ddf();
#endif

#ifdef NAG_PROTO
extern void f11dbf(char *trans,  Integer n,  double a[],  Integer la,
             Integer irow[], Integer icol[], Integer ipivp[], Integer ipivq[],
             Integer istr[], Integer idiag[],  char *check,  double y[],
             double x[],  Integer *ifail);
#else
extern     void f11dbf();
#endif


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11zac(Integer n, Integer *nnz,  double a[],  Integer irow[],
             Integer icol[],  Nag_SparseNsym_Dups dup,
             Nag_SparseNsym_Zeros zero,  Integer istr[], NagError *fail);
#else
extern     void f11zac();
#endif

#ifdef NAG_PROTO
extern void f11bcf(Integer *itn,  double *stplhs, double *stprhs,
             double *anorm, double *sigmax, Nag_Sparse_Comm *comm, Integer *ifail);
#else
extern     void f11bcf();
#endif


#ifdef NAG_PROTO
extern void f11baf(char *method, char *precon, char *norm, char *weight,
             Integer iterm, Integer n, Integer m,  double tol,
             Integer maxitn,  double anorm, double sigmax,
             Integer monit, Integer *lwreq, Nag_Sparse_Comm *comm, Integer *ifail);
#else
extern     void f11baf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11dec(Nag_SparseNsym_Method method, Nag_SparseNsym_PrecType precon,  
             Integer n, Integer nnz,
             double a[],  Integer irow[], Integer icol[],
             double omega, double b[],  Integer m,  double tol,
             Integer maxitn,  double x[], double *rnorm,
             Integer *itn, Nag_Sparse_Comm *comm, NagError *fail);
#else
extern void f11dec();
#endif

#ifdef NAG_PROTO
extern void f11daz(Integer m, Integer ind[], Integer indmin, Integer indmax,
             Integer istll[], Integer ll[], Integer iwork[]);
#else
extern     void f11daz();
#endif


#ifdef NAG_PROTO
extern void f11dax(Integer n, Integer istr[], Integer ind, Integer *ir);
#else
extern void f11dax();
#endif

#ifdef NAG_PROTO
extern void f11day(Integer n, Integer nnz, Integer *nnzc,  double **a,
             Integer *la, Integer **irow, Integer **icol, Integer lfill,
             double dtol,  Integer kindp,  double alpha,
             Integer ipivp[], Integer ipivq[], Integer *npivm, Integer istr[],
             Integer istra[], Integer istca[], Integer nnzr[], Integer *istll,
             Integer llnnzf[], Integer llnnzb[], Integer idiag[],
             Integer iwork[], Integer ibad, Integer *ierror);
#else
extern      void f11day();
#endif

#ifdef NAG_PROTO
  extern NAG_DLL_EXPIMP void NAG_CALL f11dac(Integer n, Integer nnz,  double **a,  Integer *la,
             Integer **irow, Integer **icol, Integer lfill,  double dtol,
             Nag_SparseNsym_Piv pstrat, Nag_SparseNsym_Fact milu, 
             Integer ipivp[], Integer ipivq[],
             Integer istr[], Integer idiag[], Integer *nnzc, Integer *npivm,
             NagError *fail);
#else
   extern void f11dac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f11dcc(Nag_SparseNsym_Method method,  Integer n, Integer nnz,  double a[],
             Integer la, Integer irow[], Integer icol[], Integer ipivp[],
             Integer ipivq[], Integer istr[], Integer idiag[],
             double b[],  Integer m,  double tol,  Integer maxitn,
             double x[], double *rnorm,  Integer *itn,
             Nag_Sparse_Comm *comm, NagError *fail);
#else
extern  void f11dcc();
#endif



#ifdef NAG_PROTO
 extern void f11bbx(Boolean *next,  Integer *irevcm,  Boolean precon,  Integer n,
             double xnorm, double *rnorm2, double b[],
             double x[], double r[], double q[], double t[],
             double u[], double v[], Nag_Sparse_Comm *comm);
#else
extern     void f11bbx();
#endif



#ifdef NAG_PROTO
extern void f11bbz(Integer *irevcm, Integer idata[],  double rdata[],
             double u[], double v[], double work[],  Integer lwork,
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern     void f11bbz();
#endif

#ifdef NAG_PROTO
extern void f11bbw(Boolean *next, Boolean *done,  Integer *irevcm, Integer kill,
             Boolean precon,  Integer iterm, Integer *itn, Integer maxitn,
             double tol,  Integer m, Integer n,  double *sigmax,
             double rnrm20, double *xnorm2, double *rnorm2,
             double x[], double q[], double h[], double c[],
             double s[], double t[], double u[], double v[],
             Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern      void f11bbw();
#endif

#ifdef NAG_PROTO
extern void f11bbv(Boolean *next,  Integer *irevcm,  Boolean *fcall,
             Boolean precon,  Integer iterm, Integer norm,  double tol,
             Integer itn, Integer maxitn, Integer n,  double anorm,
             double sigmax, double bnorm, double rnrm20,
             double *xnorm, double xnorm2, double rnorm2,
             double *stplhs, double *stprhs, double b[],
             double x[], double wgt[], double q[], double t[],
             double r[], double u[], double v[], Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern     void f11bbv();
#endif

#ifdef NAG_PROTO
extern void f11bbn(Boolean *next,  Integer *irevcm, Integer kill,  Boolean *restrt,
             Boolean *finish, Boolean precon,  Integer iterm, Integer *itn,
             Integer maxitn, Integer m, Integer n,  double x[],
             double r[], double rbar[], double q[], double u[],
             double v[], double gamma1[], double gammap[],
             double tau[], double dx[], Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11bbn(); 
#endif

#ifdef NAG_PROTO
extern void f11bbt(Integer *irevcm, Integer idata[],  double rdata[],
             double u[], double v[], double work[],  Integer lwork,
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern     void f11bbt();
#endif

#ifdef NAG_PROTO
extern void f11bbs(Boolean *next,  Integer *irevcm,  Boolean restrt,
             Boolean precon,  Integer iterm, Integer n,  double bnorm,
             double xnorm, double b[], double x[], double r[],
             double y[], double u[], double v[],  Integer infoch);
#else
extern     void f11bbs();
#endif

#ifdef NAG_PROTO
extern void f11bbr(Boolean *next,  Integer *irevcm,  Boolean *restrt,
             Boolean precon,  Integer iterm, Integer *itn, Integer n,
             double x[], double r[], double rbar[], double p[],
             double t[], double y[], double u[], double v[],
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern     void f11bbr();
#endif



#ifdef NAG_PROTO
extern void f11bbq(Boolean *next,  Integer *irevcm,  Boolean *restrt,
             Integer kill,  Boolean precon,  Integer iterm, Integer norm,
             double tol,  Integer itn,  Boolean *fcall,  Integer n,
             double anorm, double sigmax, double bnorm,
             double *xnorm, double *stplhs, double *stprhs,
             double b[], double x[], double r[], double y[],
             double wgt[], double u[], double v[],
             Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern     void f11bbq();
#endif

#ifdef NAG_PROTO
extern void f11bbp(Integer *irevcm, Integer idata[],  double rdata[],
             double u[], double v[], double work[],  Integer lwork,
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern     void f11bbp();
#endif


#ifdef NAG_PROTO
extern void f11baz(Integer action, Integer idata[],  double rdata[],
             Nag_Sparse_Comm *comm, Integer *info);
#else
extern void f11baz();
#endif

#ifdef NAG_PROTO
extern double f11bbu(Integer norm, Integer n,  double x[], double wgt[],
                      Nag_Sparse_Comm *comm);
#else
extern double f11bbu();
#endif

#ifdef NAG_PROTO
extern void f11bbm(Boolean *next,  Integer *irevcm,  Boolean *restrt,
             Boolean *fcall, Boolean *doxn,  Integer kill,  Boolean precon,
             Integer iterm, Integer norm,  double tol,  Integer itn,
             Integer m,  double anorm, double sigmax, double bnorm,
             double *xnorm,  Integer n,  double b[], double x[],
             double wgt[], double r[], double dx[], double px[],
             double *stplhs, double *stprhs, double u[],
             double v[], Nag_Sparse_Comm *comm, Integer *infoch);
#else
extern     void f11bbm();
#endif

#ifdef NAG_PROTO
extern void f11bbf(Integer *irevcm,  double u[], double v[],
             double work[],  Integer lwork, Nag_Sparse_Comm *comm, Integer *ifail);
#else
extern     void f11bbf();
#endif

#ifdef NAG_PROTO
extern void f11bby(Boolean *next,  Integer *irevcm, Integer norm,
             double *anorm,  Integer n,  double work[],
             double u[], double v[], Nag_Sparse_Comm *comm);
#else
extern     void f11bby();
#endif

#ifdef NAG_PROTO
extern void f11xaf(char *trans,  Integer n, Integer nnz,  double a[],
             Integer irow[], Integer icol[],  char *check,  double x[],
             double y[],  Integer *ifail);
#else
extern void f11xaf();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGF11 */








