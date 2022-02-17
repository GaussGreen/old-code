#ifndef NAGD01
#define NAGD01
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagd01.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library d01 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2144 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01ajc(NAG_D01AJC_FUN f , double a, double b, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            NagError NAG_HUGE *fail);
#else
extern void d01ajc();
#endif

#ifdef NAG_PROTO
extern void d01ajv(NAG_D01_FUN f, double a, double b,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01ajv();
#endif

#ifdef NAG_PROTO
extern void d01ajx(Integer limit, Integer last, Integer NAG_HUGE *maxerr,
            double NAG_HUGE *ermax, double NAG_HUGE elist[], Integer NAG_HUGE iord[],
            Integer NAG_HUGE *nrmax);
#else
extern void d01ajx();
#endif

#ifdef NAG_PROTO
extern void d01ajy(Integer n, double NAG_HUGE epstab[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, double NAG_HUGE res3la[], Integer NAG_HUGE *nres);
#else
extern void d01ajy();
#endif

#ifdef NAG_PROTO
extern void d01ajz(NAG_D01_FUN f, double a, double b,
            double NAG_HUGE *result, double NAG_HUGE *abserr, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *funeval);
#else
extern void d01ajz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01akc(NAG_D01AKC_FUN f, double a, double b, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            NagError NAG_HUGE *fail);
#else
extern void d01akc();
#endif

#ifdef NAG_PROTO
extern void d01akv(NAG_D01_FUN f, double a, double b,
            Integer key, double epsabs, double epsrel,
            double NAG_HUGE alist[], double NAG_HUGE blist[], double NAG_HUGE elist[],
            double NAG_HUGE rlist[], Integer limit, Integer NAG_HUGE iord[],
            double NAG_HUGE *result, double NAG_HUGE *abserr, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *last, Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01akv();
#endif

#ifdef NAG_PROTO
extern void d01akz(NAG_D01_FUN f, double a, double b,
            double NAG_HUGE *result, double NAG_HUGE *abserr, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *funeval);
#else
extern void d01akz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01alc(NAG_D01ALC_FUN f, double a, double b, 
            Integer nbrkpts, double NAG_HUGE brkpts[], double epsabs,
            double epsrel, Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp, NagError NAG_HUGE *fail);
#else
extern void d01alc();
#endif

#ifdef NAG_PROTO
extern void d01alv(NAG_D01_FUN f, double a, double b,
            Integer npts2, double NAG_HUGE points[], double NAG_HUGE *pts,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], Integer NAG_HUGE *level,
            Integer NAG_HUGE *ndin, double NAG_HUGE *result, double NAG_HUGE *abserr,
            Integer NAG_HUGE *neval, Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01alv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01amc(NAG_D01AMC_FUN f, Nag_BoundInterval boundinf, double bound, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            NagError NAG_HUGE *fail);
#else
extern void d01amc();
#endif

#ifdef NAG_PROTO
extern void d01amv(NAG_D01_FUN f, double bound, Integer boundinf,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval,
            char NAG_HUGE buf[]);
#else
extern void d01amv();
#endif

#ifdef NAG_PROTO
extern void d01amz(NAG_D01_FUN f, double boun, Integer inf,
            double a, double b, double NAG_HUGE *result,
            double NAG_HUGE *abserr, double NAG_HUGE *resabs, double NAG_HUGE *resasc, Integer NAG_HUGE *funeval);
#else
extern void d01amz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01anc(NAG_D01ANC_FUN g, double a, double b, 
            double omega, Nag_TrigTransform wt_func, double epsabs,
            double epsrel, Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,  NagError NAG_HUGE *fail);
#else
extern void d01anc();
#endif

#ifdef NAG_PROTO
extern void d01anv(NAG_D01_FUN g, double a, double b,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            double NAG_HUGE chebmo[], Integer  maxp1, Integer  limit,
            Integer NAG_HUGE iord[], Integer NAG_HUGE nnlog[], double NAG_HUGE *result,
            double NAG_HUGE *abserr,double omega,
            Integer integr, Integer icall, 
            Integer NAG_HUGE *momcom, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01anv();
#endif

#ifdef NAG_PROTO
extern void d01anw(NAG_D01_FUN g, double a, double b,
            double omega, Integer integr, Integer nrmom,
            Integer maxp1, Integer ksave, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *momcom, double NAG_HUGE chebmo[],
            Integer NAG_HUGE *funeval);
#else
extern void d01anw();
#endif

#ifdef NAG_PROTO
extern void d01anx(double NAG_HUGE *x, double NAG_HUGE fval[], double NAG_HUGE cheb12[],
            double NAG_HUGE cheb24[]);
#else
extern void d01anx();
#endif

#ifdef NAG_PROTO
extern double d01any(double x, double omega, double p2, double p3,
              double p4, Integer integr);
#else
extern double d01any();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01apc(NAG_D01APC_FUN g, double a, double b, 
            double alfa, double beta, Nag_QuadWeight  wt_func,
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            NagError NAG_HUGE *fail);
#else
extern void d01apc();
#endif

#ifdef NAG_PROTO
extern void d01apv(NAG_D01_FUN g, double a, double b,
            double alfa, double beta, double epsabs,
            double epsrel, double NAG_HUGE alist[], double NAG_HUGE blist[],
            double NAG_HUGE elist[], double NAG_HUGE rlist[], Integer  limit,
            Integer NAG_HUGE iord[], Integer integr,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01apv();
#endif

#ifdef NAG_PROTO
extern void d01apw(double alfa, double beta, double NAG_HUGE ri[],
            double NAG_HUGE rj[], double NAG_HUGE rg[], double NAG_HUGE rh[], Integer integr);
#else
extern void d01apw();
#endif

#ifdef NAG_PROTO
extern void d01apx(NAG_D01_FUN g, double a, double b,
            double bl, double br, double alfa,
            double beta, double NAG_HUGE ri[], double NAG_HUGE rj[], double NAG_HUGE rg[],
            double NAG_HUGE rh[], double NAG_HUGE *result, double NAG_HUGE *abserr,
            double NAG_HUGE *resasc, Integer  integr, Integer NAG_HUGE *nev,
            Integer NAG_HUGE *funeval);
#else
extern void d01apx();
#endif

#ifdef NAG_PROTO
extern double d01apy(double x, double a, double b, double alfa,
              double beta, Integer integr);
#else
extern double d01apy();
#endif

#ifdef NAG_PROTO
extern void d01apz(NAG_D01_FUN g,
            double (NAG_HUGE *w) (double x, double omega, double p2, double p3,
                         double p4, Integer integr),
            double p1,
            double p2, double p3, double p4, Integer kp,
            double a, double b, double NAG_HUGE *result,
            double NAG_HUGE *abserr, double NAG_HUGE *resabs, double NAG_HUGE *resasc,
            Integer NAG_HUGE *funeval);
#else
extern void d01apz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01aqc(NAG_D01AQC_FUN g, double a, double b, 
            double c, double epsabs, double epsrel,
            Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp, NagError NAG_HUGE *fail);
#else
extern void d01aqc();
#endif

#ifdef NAG_PROTO
extern void d01aqv(NAG_D01_FUN g, double a, double b, double c,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer  limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *last, Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[]);
#else
extern void d01aqv();
#endif

#ifdef NAG_PROTO
extern void d01aqy(NAG_D01_FUN g, double a, double b, double c,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Integer NAG_HUGE *krul,
            Integer NAG_HUGE *neval, Integer NAG_HUGE *funeval);
#else
extern void d01aqy();
#endif

#ifdef NAG_PROTO
extern double d01aqz(double x, double c, double p2, double p3,
              double p4, Integer  kp);
#else
extern double d01aqz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01asc(NAG_D01ASC_FUN g, double a, 
            double omega, Nag_TrigTransform wt_func,
            Integer maxintervals,
            Integer maxsubint_per_int,
            double epsabs,
            double NAG_HUGE *result, double NAG_HUGE *abserr,
            Nag_QuadSubProgress NAG_HUGE *qpsub,
            NagError NAG_HUGE *fail);
#else
extern void d01asc();
#endif

#ifdef NAG_PROTO
extern void d01asv(NAG_D01_FUN g, double a, double epsabs,
            double NAG_HUGE alist[], double NAG_HUGE blist[], double NAG_HUGE elist[],
            double NAG_HUGE rlist[], double NAG_HUGE chebmo[], Integer maxp1,
            double NAG_HUGE erlst[], double NAG_HUGE rslst[], Integer NAG_HUGE ierlst[],
            Integer limlst, Integer limit, Integer NAG_HUGE iord[], 
            Integer NAG_HUGE nnlog[], double NAG_HUGE *result, double NAG_HUGE *abserr,
            double omega, Integer integr, Integer NAG_HUGE *lst, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[], NagError NAG_HUGE *fail);
#else
extern void d01asv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL d01bac(Nag_GaussFormulae quadrule, NAG_D01BAC_FUN fun,
              double a, double b, Integer npts, NagError NAG_HUGE *fail);
#else
extern double d01bac();
#endif

#ifdef NAG_PROTO
extern void d01baw(double a, double b, Integer itype, Integer npts,
            double NAG_HUGE weight[], double NAG_HUGE abscis[], Integer NAG_HUGE *fail);
#else
extern void d01baw();
#endif

#ifdef NAG_PROTO
extern void d01bax(double a, double b, Integer itype, Integer npts,
            double NAG_HUGE weight[], double NAG_HUGE abscis[], Integer NAG_HUGE *fail);
#else
extern void d01bax();
#endif

#ifdef NAG_PROTO
extern void d01bay(double a, double b, Integer itype, Integer npts,
            double NAG_HUGE weight[], double NAG_HUGE abscis[], Integer NAG_HUGE *fail);
#else
extern void d01bay();
#endif

#ifdef NAG_PROTO
extern void d01baz(double a, double b, Integer itype, Integer npts,
            double NAG_HUGE weight[], double NAG_HUGE abscis[], Integer NAG_HUGE *fail);
#else
extern void d01baz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01bbc(Nag_GaussFormulae quadrule, double a, double b,
            Integer itype, Integer npts, double NAG_HUGE weight[],
            double NAG_HUGE abscis[], NagError NAG_HUGE *fail);
#else
extern void d01bbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01fcc(Integer ndim, NAG_D01FCC_FUN functn, double NAG_HUGE a[], 
            double NAG_HUGE b[], Integer NAG_HUGE *minpts, Integer maxpts,  double eps,
            double NAG_HUGE *finval, double NAG_HUGE *acc, NagError NAG_HUGE *fail);
#else
extern void d01fcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01gac(Integer n, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE *ans,
            double NAG_HUGE *er, NagError NAG_HUGE *fail);
#else
extern void d01gac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01gbc(Integer numvar, NAG_D01GBC_FUN functn, 
            Nag_MCMethod method, Nag_Start cont, double NAG_HUGE a[], double NAG_HUGE b[],
            Integer NAG_HUGE *minpts, Integer maxpts, double releps, double NAG_HUGE *finval,
            double NAG_HUGE *relerr, double NAG_HUGE *NAG_HUGE *wrkstr, NagError NAG_HUGE *fail);
#else
extern void d01gbc();
#endif

#ifdef NAG_PROTO
extern void d01gbz(Integer numvar, double NAG_HUGE a[], double NAG_HUGE b[],
            Integer NAG_HUGE *minpts, Integer maxpts,
            NAG_D01GBC_FUN functn,
            double releps, double NAG_HUGE *relerr, double NAG_HUGE cntlim[],
            double NAG_HUGE countr[], double NAG_HUGE lower[], double NAG_HUGE pointr[],
            double NAG_HUGE bmnusa[], double NAG_HUGE width[], double NAG_HUGE x[],
            Integer lenwrk, double NAG_HUGE wrkstr[], double NAG_HUGE *finval,
            Integer NAG_HUGE *ierror);
#else
extern void d01gbz();
#endif

/* Multi-threading versions of d01s above */
#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01sjc(NAG_D01SJC_FUN f , double a, double b, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01sjc();
#endif

#ifdef NAG_PROTO
extern void d01sjv(NAG_D01_TS_FUN f, double a, double b,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01sjv();
#endif

#ifdef NAG_PROTO
extern void d01sjz(NAG_D01_TS_FUN f, double a, double b,
            double NAG_HUGE *result, double NAG_HUGE *abserr, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01sjz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01skc(NAG_D01SKC_FUN f, double a, double b, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01skc();
#endif

#ifdef NAG_PROTO
extern void d01skv(NAG_D01_TS_FUN f, double a, double b,
            Integer key, double epsabs, double epsrel,
            double NAG_HUGE alist[], double NAG_HUGE blist[], double NAG_HUGE elist[],
            double NAG_HUGE rlist[], Integer limit, Integer NAG_HUGE iord[],
            double NAG_HUGE *result, double NAG_HUGE *abserr, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *last, Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval,
            char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01skv();
#endif

#ifdef NAG_PROTO
extern void d01skz(NAG_D01_TS_FUN f, double a, double b,
            double NAG_HUGE *result, double NAG_HUGE *abserr, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01skz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01slc(NAG_D01SLC_FUN f, double a, double b, 
            Integer nbrkpts, double NAG_HUGE brkpts[], double epsabs,
            double epsrel, Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp, 
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01slc();
#endif

#ifdef NAG_PROTO
extern void d01slv(NAG_D01_TS_FUN f, double a, double b,
            Integer npts2, double NAG_HUGE points[], double NAG_HUGE *pts,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], Integer NAG_HUGE *level,
            Integer NAG_HUGE *ndin, double NAG_HUGE *result, double NAG_HUGE *abserr,
            Integer NAG_HUGE *neval, Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01slv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01smc(NAG_D01SMC_FUN f, Nag_BoundInterval boundinf, double bound, 
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01smc();
#endif

#ifdef NAG_PROTO
extern void d01smv(NAG_D01_TS_FUN f, double bound, Integer boundinf,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *ier, 
            Integer NAG_HUGE *funeval,
            char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01smv();
#endif

#ifdef NAG_PROTO
extern void d01smz(NAG_D01_TS_FUN f, double boun, Integer inf,
            double a, double b, double NAG_HUGE *result,
            double NAG_HUGE *abserr, double NAG_HUGE *resabs, double NAG_HUGE *resasc, 
            Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01smz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01snc(NAG_D01SNC_FUN g, double a, double b, 
            double omega, Nag_TrigTransform wt_func, double epsabs,
            double epsrel, Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01snc();
#endif

#ifdef NAG_PROTO
extern void d01snv(NAG_D01_TS_FUN g, double a, double b,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            double NAG_HUGE chebmo[], Integer  maxp1, Integer  limit,
            Integer NAG_HUGE iord[], Integer NAG_HUGE nnlog[], double NAG_HUGE *result,
            double NAG_HUGE *abserr,double omega,
            Integer integr, Integer icall, 
            Integer NAG_HUGE *momcom, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, 
            char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01snv();
#endif

#ifdef NAG_PROTO
extern void d01snw(NAG_D01_TS_FUN g, double a, double b,
            double omega, Integer integr, Integer nrmom,
            Integer maxp1, Integer ksave, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, double NAG_HUGE *resabs,
            double NAG_HUGE *resasc, Integer NAG_HUGE *momcom, double NAG_HUGE chebmo[],
            Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01snw();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01spc(NAG_D01SPC_FUN g, double a, double b, 
            double alfa, double beta, Nag_QuadWeight  wt_func,
            double epsabs, double epsrel, Integer max_num_subint,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp,
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01spc();
#endif

#ifdef NAG_PROTO
extern void d01spv(NAG_D01_TS_FUN g, double a, double b,
            double alfa, double beta, double epsabs,
            double epsrel, double NAG_HUGE alist[], double NAG_HUGE blist[],
            double NAG_HUGE elist[], double NAG_HUGE rlist[], Integer  limit,
            Integer NAG_HUGE iord[], Integer integr, double NAG_HUGE *result, 
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01spv();
#endif

#ifdef NAG_PROTO
extern void d01spx(NAG_D01_TS_FUN g, double a, double b,
            double bl, double br, double alfa,
            double beta, double NAG_HUGE ri[], double NAG_HUGE rj[], double NAG_HUGE rg[],
            double NAG_HUGE rh[], double NAG_HUGE *result, double NAG_HUGE *abserr,
            double NAG_HUGE *resasc, Integer  integr, Integer NAG_HUGE *nev,
            Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01spx();
#endif

#ifdef NAG_PROTO
extern void d01spz(NAG_D01_TS_FUN g,
            double (NAG_HUGE *w) (double x, double omega, double p2, double p3,
                         double p4, Integer integr),
            double p1,
            double p2, double p3, double p4, Integer kp,
            double a, double b, double NAG_HUGE *result,
            double NAG_HUGE *abserr, double NAG_HUGE *resabs, double NAG_HUGE *resasc,
            Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01spz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01sqc(NAG_D01SQC_FUN g, double a, double b, 
            double c, double epsabs, double epsrel,
            Integer max_num_subint, double NAG_HUGE *result,
            double NAG_HUGE *abserr, Nag_QuadProgress NAG_HUGE *qp, 
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01sqc();
#endif

#ifdef NAG_PROTO
extern void d01sqv(NAG_D01_TS_FUN g, double a, double b, double c,
            double epsabs, double epsrel, double NAG_HUGE alist[],
            double NAG_HUGE blist[], double NAG_HUGE elist[], double NAG_HUGE rlist[],
            Integer  limit, Integer NAG_HUGE iord[], double NAG_HUGE *result,
            double NAG_HUGE *abserr, Integer NAG_HUGE *neval, Integer NAG_HUGE *last,Integer NAG_HUGE *ier,
            Integer NAG_HUGE *funeval, char NAG_HUGE buf[], Nag_User NAG_HUGE *comm);
#else
extern void d01sqv();
#endif

#ifdef NAG_PROTO
extern void d01sqy(NAG_D01_TS_FUN g, double a, double b, double c,
            double NAG_HUGE *result, double NAG_HUGE *abserr, Integer NAG_HUGE *krul,
            Integer NAG_HUGE *neval, Integer NAG_HUGE *funeval, Nag_User NAG_HUGE *comm);
#else
extern void d01sqy();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01ssc(NAG_D01SSC_FUN g, double a, 
            double omega, Nag_TrigTransform wt_func,
            Integer maxintervals,
            Integer maxsubint_per_int,
            double epsabs,
            double NAG_HUGE *result, double NAG_HUGE *abserr,
            Nag_QuadSubProgress NAG_HUGE *qpsub, Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01ssc();
#endif

#ifdef NAG_PROTO
extern void d01ssv(NAG_D01_TS_FUN g, double a, double epsabs,
            double NAG_HUGE alist[], double NAG_HUGE blist[], double NAG_HUGE elist[],
            double NAG_HUGE rlist[], double NAG_HUGE chebmo[], Integer maxp1,
            double NAG_HUGE erlst[], double NAG_HUGE rslst[], Integer NAG_HUGE ierlst[],
            Integer limlst, Integer limit, Integer NAG_HUGE iord[], 
            Integer NAG_HUGE nnlog[], double NAG_HUGE *result, double NAG_HUGE *abserr,
            double omega, Integer integr, Integer NAG_HUGE *lst, Integer NAG_HUGE *neval,
            Integer NAG_HUGE *ier, Integer NAG_HUGE *funeval, char NAG_HUGE buf[], 
            NagError NAG_HUGE *fail, Nag_User *comm);
#else
extern void d01ssv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL d01tac(Nag_GaussFormulae quadrule, NAG_D01TAC_FUN fun,
              double a, double b, Integer npts, Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern double d01tac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01wcc(Integer ndim, NAG_D01WCC_FUN functn, double NAG_HUGE a[], 
            double NAG_HUGE b[], Integer NAG_HUGE *minpts, Integer maxpts,  double eps,
            double NAG_HUGE *finval, double NAG_HUGE *acc, Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01wcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL d01xbc(Integer numvar, NAG_D01XBC_FUN functn, 
            Nag_MCMethod method, Nag_Start cont, double NAG_HUGE a[], double NAG_HUGE b[],
            Integer NAG_HUGE *minpts, Integer maxpts, double releps, double NAG_HUGE *finval,
            double NAG_HUGE *relerr, double NAG_HUGE *NAG_HUGE *wrkstr, 
            Nag_User NAG_HUGE *comm, NagError NAG_HUGE *fail);
#else
extern void d01xbc();
#endif

#ifdef NAG_PROTO
extern void d01xbz(Integer numvar, double NAG_HUGE a[], double NAG_HUGE b[],
            Integer NAG_HUGE *minpts, Integer maxpts,
            NAG_D01XBC_FUN functn,
            double releps, double NAG_HUGE *relerr, double NAG_HUGE cntlim[],
            double NAG_HUGE countr[], double NAG_HUGE lower[], double NAG_HUGE pointr[],
            double NAG_HUGE bmnusa[], double NAG_HUGE width[], double NAG_HUGE x[],
            Integer lenwrk, double NAG_HUGE wrkstr[], double NAG_HUGE *finval,
            Integer NAG_HUGE *ierror, Nag_User NAG_HUGE *comm);
#else
extern void d01xbz();
#endif




#ifdef __cplusplus
}
#endif
#endif /* not NAGD01 */
