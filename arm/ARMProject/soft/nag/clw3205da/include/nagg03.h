#ifndef NAGG03
#define NAGG03
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg03.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g03 Chapter
 *
 * Mark 5, 1997.
 */

#include <nag_stddef.h>
#include <nag_g03mesg.h>

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03aac(Nag_PrinCompMat pcmatrix, Nag_PrinCompScores scores, Integer n,
	      Integer m, double x[],  Integer tdx, Integer isx[],  double s[],
	      double *wtptr,  Integer nvar,  double e[],  Integer tde,
	      double p[],  Integer tdp,  double v[],  Integer tdv,
	      NagError *fail);
#else
extern void g03aac();
#endif

#ifdef NAG_PROTO
extern void g03aaf(char *matrix, char *std, char *weight,  Integer n, Integer m,
	    double x[],  Integer tdx, Integer isx[],  double s[],
	    double wt[],  Integer nvar,  double e[],  Integer tde,
	    double p[],  Integer tdp,  double v[],  Integer tdv,
	    double wk[],  Integer *ifail);
#else
extern void g03aaf();
#endif

#ifdef NAG_PROTO
extern void g03aaz(char *matrix, char *weight,  Integer n,  double x[],
	    Integer tdx, Integer m, Integer isx[], Integer ivar,
	    double wt[], double t, double v[],  Integer tdv,
	    double s[], double e[]);
#else
extern void g03aaz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03acc(Nag_Weightstype weights, Integer n, Integer m,  double x[],
	      Integer tdx, Integer isx[], Integer nx, Integer ing[],
	      Integer ng,  double *wtptr,  Integer nig[],
	      double cvm[],  Integer tdcvm,  double e[],
	      Integer tde, Integer *ncv,  double cvx[],  Integer tdcvx,
	      double tol,  Integer *irankx, NagError *fail);
#else
extern void g03acc();
#endif

#ifdef NAG_PROTO
extern void g03acf(char *weight,  Integer n, Integer m,  double x[],
	    Integer tdx, Integer isx[], Integer nx, Integer ing[],
	    Integer ng,  double wt[],  Integer nig[],
	    double cvm[],  Integer tdcvm,  double e[],
	    Integer tde, Integer *ncv,  double cvx[],  Integer tdcvx,
	    double tol,  Integer *irankx,  double wk[],
	    Integer iwk, Integer *ifail);
#else
extern void g03acf();
#endif

#ifdef NAG_PROTO
extern void g03acz(double e[],  Integer tde,  double wsum,  Integer ncv,
	    Integer nx, Integer ny, Integer *ierror);
#else
extern void g03acz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03adc(Integer n, Integer m,  double z[],
	    Integer tdz, Integer isz[], Integer nx, Integer ny,
	    double *wtptr, double e[],  Integer tde, Integer *ncv,
	    double cvx[],  Integer tdcvx,
	    double cvy[],  Integer tdcvy,  double tol, NagError *fail);
#else
extern void g03adc();
#endif

#ifdef NAG_PROTO
extern void g03adf(char *weight,  Integer n, Integer m,  double z[],
	    Integer tdz, Integer isz[], Integer *nx, Integer *ny,
	    double wt[], double e[],  Integer tde, Integer *ncv,
	    double cvx[],  Integer tdcvx, Integer mcv,
	    double cvy[],  Integer tdcvy,  double tol,
	    double wk[],  Integer iwk, Integer *ifail);
#else
extern void g03adf();
#endif

#ifdef NAG_PROTO
extern void g03adz(Integer n, Integer nx, Integer ny,  double cvx[],
	    Integer tdcvx,  double cvy[],  Integer tdcvy,
	    double qx[], double qy[],  Integer lqy,  double rdf,
	    double tol,  Integer *irankx, Integer *iranky, Integer *ncv,
	    double e[], double wk[],  Integer *ierror);
#else
extern void g03adz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03bac(Nag_RotationLoading stand, double g,  Integer nvar, Integer k,
	    double fl[],  Integer tdf,  double flr[],
	    double r[],  Integer tdr,  double acc,  Integer maxit,
	    Integer *iter,  NagError *fail);
#else
extern void g03bac();
#endif

#ifdef NAG_PROTO
extern void g03baf(char *stand,  double g,  Integer nvar, Integer k,
	    double fl[],  Integer tdf,  double flr[],
	    double r[],  Integer tdr,  double acc,  Integer maxit,
	    Integer *iter,  double wk[],  Integer *ifail);
#else
extern void g03baf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03bcc(Nag_TransNorm stand, Nag_RotationScale pscale, Integer n,
	    Integer m, double x[],  Integer tdx,  double y[],
	    Integer tdy, double yhat[], double r[],  Integer tdr,
	    double *alpha, double *rss, double res[], NagError *fail);
#else
extern void g03bcc();
#endif

#ifdef NAG_PROTO
extern void g03bcf(char *stand, char *pscale,  Integer n, Integer m,
	    double x[],  Integer tdx,  double y[],  Integer tdy,
	    double yhat[], double r[],  Integer tdr,
	    double *alpha, double *rss, double res[],
	    double wk[],  Integer *ifail);
#else
extern void g03bcf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03cac(Nag_FacMat matrix, Integer n, Integer m,
	    double x[],  Integer tdx, Integer nvar, Integer isx[],
	    Integer nfac,  double *wtptr, double e[], double stat[],
	    double com[], double psi[], double res[],
	    double fl[],  Integer tdfl, Nag_E04_Opt *options, double eps,
	    /* Note eps is an argument to g03cac and not to e04lbc */
	    NagError *fail);
#else
extern void g03cac();
#endif

#ifdef NAG_PROTO
extern void g03caf(char *matrix, char *weight,  Integer n, Integer m,
	    double x[],  Integer tdx, Integer nvar, Integer isx[],
	    Integer nfac,  double wt[], double e[], double stat[],
	    double com[], double psi[], double res[],
	    double fl[],  Integer tdfl, Nag_E04_Opt *user_opt, double eps,
	    Integer iwk[],
	    double wk[],  Integer lwk, NagError *fail, Integer *ifail);
#else
extern void g03caf();
#endif

#ifdef NAG_PROTO
extern void g03caw(char *weight,  Integer n,  double x[],  Integer tdx,
	    Integer m, Integer isx[], Integer ivar,  double wt[],
	    double t, double v[],  Integer tdv,  double s[],
	    double e[]);
#else
extern void g03caw();
#endif

#ifdef NAG_PROTO
extern void g03cax(Integer k,  double x[], double *f,
	    double g[],  Nag_Comm *comm);
#else
extern void g03cax();
#endif

#ifdef NAG_PROTO
extern void g03cay(Integer k,  double x[], double h[],
	    double hd[], Nag_Comm *comm);
#else
extern void g03cay();
#endif

#ifdef NAG_PROTO
extern void g03caz(const Nag_Search_State *st, Nag_Comm *comm);
#else
extern void g03caz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03ccc(Nag_FacScoreMethod method, Nag_FacRotation rotate,  Integer nvar,
	    Integer nfac,
	    double fl[],  Integer tdfl,  double psi[],
	    double e[], double r[],  Integer tdr,  double fs[],
	    Integer tdfs, NagError *fail);
#else
extern void g03ccc();
#endif

#ifdef NAG_PROTO
extern void g03ccf(char *method, char *rotate,  Integer nvar, Integer nfac,
	    double fl[],  Integer tdfl,  double psi[],
	    double e[], double r[],  Integer tdr,  double fs[],
	    Integer tdfs,  double wk[],  Integer *ifail);
#else
extern void g03ccf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03dac(Integer n, Integer m,  double x[],
	    Integer tdx, Integer isx[], Integer nvar, Integer ing[],
	    Integer ng,  double *wtptr,  Integer nig[],
	    double gmean[],  Integer tdg,  double det[],
	    double gc[], double *stat, double *df, double *sig,
	    NagError *fail);
#else
extern void g03dac();
#endif

#ifdef NAG_PROTO
extern void g03daf(char *weight,  Integer n, Integer m,  double x[],
	    Integer tdx, Integer isx[], Integer nvar, Integer ing[],
	    Integer ng,  double wt[],  Integer nig[],
	    double gmean[],  Integer tdg,  double det[],
	    double gc[], double *stat, double *df, double *sig,
	    double wk[],  Integer iwk[], Integer *ifail);
#else
extern void g03daf();
#endif

#ifdef NAG_PROTO
extern void g03daz(Integer n, Integer block[], Integer m,  double a[],
	    Integer tda);
#else
extern void g03daz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03dbc(Nag_GroupCovars equal, Nag_MahalDist mode,  Integer nvar, Integer ng,
	    double gmean[],  Integer tdg,  double gc[],
	    Integer nobs, Integer m, Integer isx[],  double x[],
	    Integer tdx,  double d[],  Integer tdd,
	    NagError *fail);
#else
extern void g03dbc();
#endif

#ifdef NAG_PROTO
extern void g03dbf(char *equal, char *mode,  Integer nvar, Integer ng,
	    double gmean[],  Integer tdg,  double gc[],
	    Integer nobs, Integer m, Integer isx[],  double x[],
	    Integer tdx,  double d[],  Integer tdd,  double wk[],
	    Integer *ifail);
#else
extern void g03dbf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03dcc(Nag_DiscrimMethod type, Nag_GroupCovars equal, Nag_PriorProbability priors, Integer nvar,
	    Integer ng, Integer nig[],  double gmean[],  Integer tdg,
	    double gc[], double det[],  Integer nobs, Integer m,
	    Integer isx[],  double x[],  Integer tdx,
	    double prior[], double p[],  Integer tdp, Integer iag[],
	    Boolean atiq,  double ati[], NagError *fail);
#else
extern void g03dcc();
#endif

#ifdef NAG_PROTO
extern void g03dcf(char *type, char *equal, char *priors,  Integer nvar,
	    Integer ng, Integer nig[],  double gmean[],  Integer tdg,
	    double gc[], double det[],  Integer nobs, Integer m,
	    Integer isx[],  double x[],  Integer tdx,
	    double prior[], double p[],  Integer tdp, Integer iag[],
	    Boolean atiq,  double ati[], double wk[],
	    Integer *ifail);
#else
extern void g03dcf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03eac(Nag_MatUpdate update, Nag_DistanceType dist, Nag_VarScaleType scale,
	    Integer n, Integer m, double x[],  Integer tdx, Integer isx[],
	    double s[], double d[],  NagError *fail);
#else
extern void g03eac();
#endif

#ifdef NAG_PROTO
extern void g03eaf(char *update, char *dist, char *scale,  Integer n, Integer m,
	    double x[],  Integer tdx, Integer isx[],  double s[],
	    double d[],  Integer *ifail);
#else
extern void g03eaf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03ecc(Nag_ClusterMethod method, Integer n,  double d[],  Integer ilc[],
	    Integer iuc[],  double cd[],  Integer iord[],
	    double dord[],  NagError *fail);
#else
extern void g03ecc();
#endif

#ifdef NAG_PROTO
extern void g03ecf(Integer method, Integer n,  double d[],  Integer ilc[],
	    Integer iuc[],  double cd[],  Integer iord[],
	    double dord[],  Integer iwk[], Integer *ifail);
#else
extern void g03ecf();
#endif

#ifdef NAG_PROTO
extern double g03ecw(Integer itype,  double dki, double dkj, double dij,
	      Integer ni, Integer nj, Integer nk);
#else
extern double g03ecw();
#endif

#ifdef NAG_PROTO
extern void g03ecx(Integer n,  double d[],  Integer inc[], Integer *ith,
	    Integer *jth,  double *dmin_);
#else
extern void g03ecx();
#endif

#ifdef NAG_PROTO
extern void g03ecy(Integer itype,
		   double (*g03ecw_) (Integer itype,  double dki, double dkj, double dij,
				     Integer ni, Integer nj, Integer nk),
		   Integer n,
		   double d[],  Integer inc[], Integer nc[], Integer join1,
		   Integer join2,  double dmin_);
#else
extern void g03ecy();
#endif

#ifdef NAG_PROTO
extern void g03ecz(Integer n, Integer ilc[], Integer iuc[],  double cd[],
	    Integer iord[],  double dord[],  Integer ind[],
	    Integer *ierror);
#else
extern void g03ecz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03efc(Integer n, Integer m,  double x[],
	    Integer tdx, Integer isx[], Integer nvar, Integer k,
	    double cmeans[],  Integer tdc,  double *wtptr,
	    Integer inc[], Integer nic[],  double css[],
	    double csw[],  Integer maxit,  NagError *fail);
#else
extern void g03efc();
#endif

#ifdef NAG_PROTO
extern void g03eff(char *weight,  Integer n, Integer m,  double x[],
	    Integer tdx, Integer isx[], Integer nvar, Integer k,
	    double cmeans[],  Integer tdc,  double wt[],
	    Integer inc[], Integer nic[],  double css[],
	    double csw[],  Integer maxit, Integer iwk[],
	    double wk[],  Integer *ifail);
#else
extern void g03eff();
#endif

#ifdef NAG_PROTO
extern void g03efu(double a[],  Integer tda, Integer n, Integer m,
	    double c[],  Integer tdc, Integer k, Integer isx[],
	    double wt[],  Integer ic1[], Integer ic2[], Integer ncp[],
	    double d[],  Integer itran[], Integer *index,
	    double csw[]);
#else
extern void g03efu();
#endif

#ifdef NAG_PROTO
extern void g03efv(double a[],  Integer tda, Integer n, Integer m,
	    double c[],  Integer tdc, Integer k, Integer isx[],
	    double wt[],  Integer ic1[], Integer ic2[], Integer ncp[],
	    double d[],  Integer itran[], Integer live[],
	    Integer *index,  double csw[]);
#else
extern void g03efv();
#endif

#ifdef NAG_PROTO
extern void g03efw(Integer n, Integer m,  double a[],  Integer tda,
	    Integer isx[], Integer nvar, Integer k,  double c[],
	    Integer tdc,  double wt[],  Integer ic1[],
	    double css[], double csw[],  Integer maxit,
	    Integer ic2[], Integer ncp[],  double d[],  Integer itran[],
	    Integer live[], Integer *ifail);
#else
extern void g03efw();
#endif

#ifdef NAG_PROTO
extern void g03efx(double a[],  Integer tda, Integer n, Integer m,
	    double c[],  Integer tdc, Integer k, Integer isx[],
	    Integer ic1[], Integer ic2[],  double an1[],
	    double an2[],  Integer ncp[],  double d[],
	    Integer itran[], Integer *index, Integer nic[]);
#else
extern void g03efx();
#endif

#ifdef NAG_PROTO
extern void g03efy(double a[],  Integer tda, Integer n, Integer m,
	    double c[],  Integer tdc, Integer k, Integer isx[],
	    Integer ic1[], Integer ic2[],  double an1[],
	    double an2[],  Integer ncp[],  double d[],
	    Integer itran[], Integer live[], Integer *index, Integer nic[]);
#else
extern void g03efy();
#endif

#ifdef NAG_PROTO
extern void g03efz(Integer n, Integer m,  double a[],  Integer tda,
	    Integer isx[], Integer nvar, Integer k,  double c[],
	    Integer tdc, Integer ic1[], Integer nic[],  double css[],
	    Integer maxit, Integer ic2[],  double an1[],
	    double an2[],  Integer ncp[],  double d[],
	    Integer itran[], Integer live[], Integer *ifail);
#else
extern void g03efz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03ehc(Nag_DendOrient orient,  Integer n,  double dord[], double dmin_,
	    double dstep,  Integer nsym,  char ***c, NagError *fail);
#else
extern void g03ehc();
#endif

#ifdef NAG_PROTO
extern void g03ehf(char *orient,  Integer n,  double dord[], double dmin_,
	    double dstep,  Integer nsym,  char ***c,
	    Integer *ifail);
#else
extern void g03ehf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03ejc(Integer n,  double cd[],  Integer iord[],
	    double dord[],  Integer *k,  double *dlevel,
	    Integer ic[], NagError *fail);
#else
extern void g03ejc();
#endif

#ifdef NAG_PROTO
extern void g03ejf(Integer n,  double cd[],  Integer iord[],
	    double dord[],  Integer *k,  double *dlevel,
	    Integer ic[], Integer *ifail);
#else
extern void g03ejf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03fac(Nag_Eigenvalues roots,  Integer n,  double d[],  Integer ndim,
	    double x[],  Integer tdx,  double eval[],
	    NagError *fail);
#else
extern void g03fac();
#endif

#ifdef NAG_PROTO
extern void g03faf(char *roots,  Integer n,  double d[],  Integer ndim,
	    double x[],  Integer tdx,  double eval[],
	    double wk[],  Integer iwk[], Integer *ifail);
#else
extern void g03faf();
#endif

#ifdef NAG_PROTO
extern void g03fay(char *roots,  double a[],  Integer n, Integer ndim,
	    double eval[],  Integer *m,  double d[], double e[],
	    double tau[],  Integer iblock[], Integer isplit[],
	    double x[],  Integer tdx,  double wk[],  Integer iwk[],
	    Integer *ierror);
#else
extern void g03fay();
#endif

#ifdef NAG_PROTO
extern void g03faz(Integer n,  double x[], double row[], double *total);
#else
extern void g03faz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03fcc(Nag_ScaleCriterion type,  Integer n, Integer ndim,  double d[],
	    double x[],  Integer tdx,  double *stress,
	    double dfit[],  
	    Nag_E04_Opt *options, NagError *fail);
#else
extern void g03fcc();
#endif

#ifdef NAG_PROTO
extern void g03fcf(char *type,  Integer n, Integer ndim,  double d[],
	    double x[],  Integer tdx,  double *stress,
	    double dfit[],  
	    double wk[],  Integer iwk[], Nag_E04_Opt *user_opt, NagError *fail, Integer *ifail);
#else
extern void g03fcf();
#endif

#ifdef NAG_PROTO
extern void g03fcw(Integer n, double x[], double *objf, double g[],
	    Nag_Comm *comm);
#else
extern void g03fcw();
#endif

#ifdef NAG_PROTO
extern void g03fcx_(Integer n, Integer m,  double x[], double d[],
	     Integer nn,  Boolean srd);
#else
extern void g03fcx_();
#endif

#ifdef NAG_PROTO
extern void g03fcy_(Integer k,  double x[], double xhat[], double tol,
	     double w[]);
#else
extern void g03fcy_();
#endif

#ifdef NAG_PROTO
extern void g03fcz_(Integer mode, Integer nm,  double x[], double *stress,
	     double der[],  Integer nstate, Integer iwk[],
	     double dfit[], size_t m01_rank[]);
#else
extern void g03fcz_();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03xzc(char ***c);
#else
extern void g03xzc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g03zac(Integer n, Integer m,  double x[],  Integer tdx,
	    Integer nvar, Integer isx[],  double s[], double e[],
	    double z[],  Integer tdz, NagError *fail);
#else
extern void g03zac();
#endif

#ifdef NAG_PROTO
extern void g03zaf(Integer n, Integer m,  double x[],  Integer tdx,
	    Integer nvar, Integer isx[],  double s[], double e[],
	    double z[],  Integer tdz, Integer *ifail);
#else
extern void g03zaf();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGG03 */
