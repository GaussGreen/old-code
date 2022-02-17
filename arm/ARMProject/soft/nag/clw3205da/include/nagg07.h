#ifndef NAGG07
#define NAGG07
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg07.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g07 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2158 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g07cac(Nag_TailProbability tail, Nag_PopVar equal,  Integer nx, Integer ny,
             double xmean, double ymean, double xstd,
             double ystd, double clevel, double NAG_HUGE *t, double NAG_HUGE *df,
             double NAG_HUGE *prob, double NAG_HUGE *dl, double NAG_HUGE *du,
             NagError NAG_HUGE *fail);
#else
extern void g07cac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g07dac(Integer n, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE *xme,
            double NAG_HUGE *xmd, double NAG_HUGE *xsd, NagError NAG_HUGE *fail);
#else
extern void g07dac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g07dbc(Nag_SigmaSimulEst sigma_est, Integer n,  double NAG_HUGE x[],  Nag_PsiFun psifun,
             double c, double h1, double h2, double h3,
             double dchi, double NAG_HUGE *theta, double NAG_HUGE *sigma,
             Integer maxit,  double tol, double NAG_HUGE rs[],  Integer NAG_HUGE *nit,
             double NAG_HUGE wrk[],  NagError NAG_HUGE *fail);
#else
extern void g07dbc();
#endif

#ifdef NAG_PROTO
extern double g07dbx(double t, Integer ipsi, double c, double h1, double h2, 
              double h3);
#else
extern double g07dbx();
#endif

#ifdef NAG_PROTO
extern double g07dby(double t,  Integer ipsi,  double dchi);
#else
extern double g07dby();
#endif

#ifdef NAG_PROTO
extern void g07dbz(double c, double NAG_HUGE *beta);
#else
extern void g07dbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g07ddc(Integer n,  double NAG_HUGE x[], double alpha, double NAG_HUGE *tmean,
             double NAG_HUGE *wmean, double NAG_HUGE *tvar, double NAG_HUGE *wvar,
             Integer NAG_HUGE *k,  double NAG_HUGE sx[],  NagError NAG_HUGE *fail);
#else
extern void g07ddc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG07 */
