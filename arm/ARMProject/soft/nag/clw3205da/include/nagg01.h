#ifndef NAGG01
#define NAGG01
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg01.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g01 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2155 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01aac(Integer n, double NAG_HUGE x[], double NAG_HUGE wt[], Integer NAG_HUGE *nvalid, 
             double NAG_HUGE *xmean, double NAG_HUGE *xsd, double NAG_HUGE *xskew, double NAG_HUGE *xkurt,
             double NAG_HUGE *xmin, double NAG_HUGE *xmax, double NAG_HUGE *wsum, NagError NAG_HUGE *fail);
#else
extern void  g01aac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01alc(Integer n,  double NAG_HUGE x[], double NAG_HUGE res[],  NagError NAG_HUGE *fail);
#else
extern void g01alc();
#endif

#ifdef NAG_PROTO
extern double g01bcc(double x, Integer n, NagError NAG_HUGE *fail);
#else
extern double g01bcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01bjc(Integer n,  double p,  Integer k,  double NAG_HUGE *plek,
             double NAG_HUGE *pgtk, double NAG_HUGE *peqk, NagError NAG_HUGE *fail);
#else
extern void g01bjc();
#endif

#ifdef NAG_PROTO
extern double g01bju(double x, Integer NAG_HUGE *local_error);
#else
extern double g01bju();
#endif

#ifdef NAG_PROTO
extern double g01bjv(double x, Integer NAG_HUGE *local_error);
#else
extern double g01bjv();
#endif

#ifdef NAG_PROTO
extern void g01bjw(double x, double a, double b, double NAG_HUGE *lpabx,
            Integer NAG_HUGE *local_error);
#else
extern void g01bjw();
#endif

#ifdef NAG_PROTO
extern void g01bjx(double x, double a, double b, double NAG_HUGE *jabx);
#else
extern void g01bjx();
#endif

#ifdef NAG_PROTO
extern void g01bjy(double x, double a, double b, double NAG_HUGE *iabx);
#else
extern void g01bjy();
#endif

#ifdef NAG_PROTO
extern void g01bjz(double x, double a, double b, double NAG_HUGE *pltx,
             double NAG_HUGE *pgtx, double NAG_HUGE *pabx, Integer NAG_HUGE *local_error);
#else
extern void g01bjz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01bkc(double rlamda,  Integer k,  double NAG_HUGE *plek,
             double NAG_HUGE *pgtk, double NAG_HUGE *peqk,  NagError NAG_HUGE *fail);
#else
extern void g01bkc();
#endif

#ifdef NAG_PROTO
extern void g01bkw(double x, double a, double NAG_HUGE *lpax, Integer NAG_HUGE *local_error);
#else
extern void g01bkw();
#endif

#ifdef NAG_PROTO
extern void g01bkx(double x, double a, double NAG_HUGE *jax);
#else
extern void g01bkx();
#endif

#ifdef NAG_PROTO
extern void g01bky(double x, double a, double NAG_HUGE *iax);
#else
extern void g01bky();
#endif

#ifdef NAG_PROTO
extern void g01bkz(double x, double a, double NAG_HUGE *pltx, double NAG_HUGE *pgtx,
             double NAG_HUGE *pax, Integer NAG_HUGE *local_error);
#else
extern void g01bkz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01blc(Integer n, Integer l, Integer m, Integer k,  double NAG_HUGE *plek,
             double NAG_HUGE *pgtk, double NAG_HUGE *peqk,  NagError NAG_HUGE *fail);
#else
extern void g01blc();
#endif

#ifdef NAG_PROTO
extern void g01blw(Integer k, Integer n, Integer l, Integer m,  double NAG_HUGE *lpk,
	     Integer NAG_HUGE *local_error);
#else
extern void g01blw();
#endif

#ifdef NAG_PROTO
extern void g01blx(Integer k, Integer n, Integer l, Integer m,  double NAG_HUGE *jk);
#else
extern void g01blx();
#endif

#ifdef NAG_PROTO
extern void g01bly(Integer k, Integer n, Integer l, Integer m,  double NAG_HUGE *ik);
#else
extern void g01bly();
#endif

#ifdef NAG_PROTO
extern void g01blz(Integer k, Integer n, Integer l, Integer m,  double NAG_HUGE *plek,
             double NAG_HUGE *pgtk, double NAG_HUGE *pk, Integer NAG_HUGE *local_error);
#else
extern void g01blz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01cec(double p, NagError NAG_HUGE *fail);
#else
extern double g01cec();
#endif

#ifdef NAG_PROTO
extern void g01dac(Integer n,  double NAG_HUGE pp[], double etol, double NAG_HUGE *errest,
             NagError NAG_HUGE *fail);
#else
extern void g01dac();
#endif

#ifdef NAG_PROTO
extern void g01dax(Integer i, Integer n,  double crln, double rinc,
             double NAG_HUGE *score, double NAG_HUGE *err);
#else
extern void g01dax();
#endif

#ifdef NAG_PROTO
extern double g01day(Integer is, Integer n, Integer NAG_HUGE *ierror);
#else
extern double g01day();
#endif

#ifdef NAG_PROTO
extern void g01daz(Integer n, Integer ndiv2,  double NAG_HUGE crln[], double NAG_HUGE prob1[],
             double NAG_HUGE pbest[], double NAG_HUGE pp[], double etol,
             double NAG_HUGE *errest,  Integer NAG_HUGE *ifail);
#else
extern void g01daz();
#endif

#ifdef NAG_PROTO
extern void g01dbc(Integer n,  double NAG_HUGE pp[],  NagError NAG_HUGE *fail);
#else
extern void g01dbc();
#endif

#ifdef NAG_PROTO
extern double g01dbz(Integer i, Integer n);
#else
extern double g01dbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01ddc(Integer n,  double NAG_HUGE x[],  Boolean calc_wts,  double NAG_HUGE a[],
             double NAG_HUGE *w, double NAG_HUGE *pw,  NagError NAG_HUGE *fail);
#else
extern void g01ddc();
#endif

#ifdef NAG_PROTO
extern void g01ddy(Integer n,  double NAG_HUGE a[], double NAG_HUGE *eps, Integer NAG_HUGE *local_error);
#else
extern void g01ddy();
#endif

#ifdef NAG_PROTO
extern double g01ddz(double NAG_HUGE *c,  Integer n,  double x);
#else
extern double g01ddz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01dhc(Nag_Scores scores, Nag_Ties ties,  Integer n,  double NAG_HUGE x[],
             double NAG_HUGE r[], NagError NAG_HUGE *fail);
#else
extern void g01dhc();
#endif

#ifdef NAG_PROTO
extern double g01dht(char NAG_HUGE *scores,  Integer i, Integer n, Integer NAG_HUGE *local_error);
#else
extern double g01dht();
#endif

#ifdef NAG_PROTO
extern void g01dhu(char NAG_HUGE *ties, char NAG_HUGE *scores,  Integer n,  double NAG_HUGE x[],
             double NAG_HUGE r[],  Integer NAG_HUGE iwrk[], Integer NAG_HUGE *local_error);
#else
extern void g01dhu();
#endif

#ifdef NAG_PROTO
extern void g01dhv(Integer n,  double NAG_HUGE x[], double NAG_HUGE r[],  Integer NAG_HUGE iwrk[]);
#else
extern void g01dhv();
#endif

#ifdef NAG_PROTO
extern void g01dhw(char NAG_HUGE *ties,  Integer n,  double NAG_HUGE x[], double NAG_HUGE r[],
             Integer NAG_HUGE iwrk[]);
#else
extern void g01dhw();
#endif

#ifdef NAG_PROTO
extern void g01dhx(Integer n,  double NAG_HUGE x[],  Integer NAG_HUGE iwrk[], Integer NAG_HUGE *local_error);
#else
extern void g01dhx();
#endif

#ifdef NAG_PROTO
extern void g01dhy(Integer n,  double NAG_HUGE x[], double NAG_HUGE r[],  Integer NAG_HUGE iwrk[]);
#else
extern void g01dhy();
#endif

#ifdef NAG_PROTO
extern void g01dhz(char NAG_HUGE *scores,  Integer n,  double NAG_HUGE r[], Integer NAG_HUGE *local_error);
#else
extern void g01dhz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01eac(Nag_TailProbability tail,  double x,  NagError NAG_HUGE *fail);
#else
extern double g01eac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01ebc(Nag_TailProbability tail, double t, double df, NagError NAG_HUGE *fail);
#else
extern double g01ebc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01ecc(Nag_TailProbability tail, double x, double df, NagError NAG_HUGE *fail);
#else
extern double g01ecc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01edc(Nag_TailProbability tail, double f, double df1, double df2,
              NagError NAG_HUGE *fail);
#else
extern double g01edc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g01eec(double x, double a, double b, double tol, double NAG_HUGE *p, 
            double NAG_HUGE *q, double NAG_HUGE *pdf, NagError NAG_HUGE *fail);
#else
extern void g01eec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01efc(Nag_TailProbability tail, double g, double a, double b,
              NagError NAG_HUGE *fail);
#else
extern double g01efc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01fac(Nag_TailProbability tail,  double p,  NagError NAG_HUGE *fail);
#else
extern double g01fac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01fbc(Nag_TailProbability tail, double p, double df, NagError NAG_HUGE *fail);
#else
extern double g01fbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01fcc(double p, double df, NagError NAG_HUGE *fail);
#else
extern double g01fcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01fdc(double p, double df1, double df2, NagError NAG_HUGE *fail);
#else
extern double g01fdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01fec(double p, double a, double b, double tol, NagError NAG_HUGE *fail);
#else
extern double g01fec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01ffc(double p, double a, double b, double tol, NagError NAG_HUGE *fail);
#else
extern double g01ffc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g01hac(double x, double y, double rho, NagError NAG_HUGE *fail);
#else
extern double g01hac();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG01 */
