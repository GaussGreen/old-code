#ifndef NAGG04
#define NAGG04
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg04.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g04 Chapter
 *
 * Mark 5, 1997.
 */


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g04bbc(Integer n, double y[], Nag_Blocks blocks, Integer iblock, Integer nt,
		   Integer it[], double *gmean, double bmean[],
		   double tmean[], double table[], double c[], Integer tdc,
		   Integer irep[], double r[], double ef[], double tol,
		   Integer irdf, NagError *fail);
#else
extern void g04bbc();
#endif



#ifdef NAG_PROTO
extern void g04bbf(Integer n, double y[], Integer iblock, Integer nt,
	     Integer it[], double *gmean, double bmean[],
	     double tmean[], double table[], Integer tdt,
	     double c[], Integer tdc, Integer irep[], double r[],
	     double ef[], double tol, Integer irdf,
	     double wk[], Integer *ifail);
#else
extern void g04bbf();
#endif

#ifdef NAG_PROTO
extern void g04bbx(Integer n, Integer iblock, Integer kblock, Integer nt,
Integer it[], double c[], Integer tdc, Integer irep[],
double ef[], Integer *ntdf, double tmean[],
double acc, double wk[], Integer *iwarn);
#else
extern void g04bbx();
#endif

#ifdef NAG_PROTO
extern void g04bby(Integer n, Integer iblock, Integer kblock, double r[],
	    double bmean[], double *rss, double *ess);
#else
extern void g04bby();
#endif

#ifdef NAG_PROTO
extern void g04bbz(Integer n, double y[], Integer ifac[], Integer levels,
	     double b[], double *rss);
#else
extern void g04bbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g04cac(Integer n, double y[], Integer nfac, Integer lfac[],
		   Integer nblock, Integer inter, Integer irdf, Integer *mterm,
		   double **table, double **tmean,
		   Integer *maxt, double **e, Integer **imean,
		   double **semean, double **bmean, double r[],
		   NagError *fail);
#else
extern void g04cac();
#endif

#ifdef NAG_PROTO
extern void g04caf_(Integer n, double y[], Integer nfac, Integer lfac[],
	     Integer nblock, Integer inter, Integer irdf, Integer mterm,
	     double table[], Integer *itotal, double tmean[],
	     Integer maxt, double e[], Integer imean[],
	     double semean[], double bmean[], double r[],
	     Integer iwk[], Integer *ifail);
#else
extern void g04caf_();
#endif

#ifdef NAG_PROTO
extern void g04cat(Integer itab[], Integer idim);
#else
extern void g04cat();
#endif

#ifdef NAG_PROTO
extern void g04cau(Boolean qind, Boolean qfor, Integer iprod[], Integer kdim,
	    Integer *jsub, Integer ivec[], Integer *ifault);
#else
extern void g04cau();
#endif

#ifdef NAG_PROTO
extern void g04cav(Integer n, Integer m, Integer ind[], Integer *ierror);
#else
extern void g04cav();
#endif

#ifdef NAG_PROTO
extern void g04caw(Integer n, Integer nlevel, double b[], double t[],
	     double *sem, double *ess, double r[],
	     Integer ifac[]);
#else
extern void g04caw();
#endif

#ifdef NAG_PROTO
extern void g04cax(Integer n, Integer nblock, Integer kblock, double r[],
	    double y[], double t[], double *ess);
#else
extern void g04cax();
#endif

#ifdef NAG_PROTO
extern void g04cay(Integer nfac, Integer lfac[], Integer nrepl, Integer n,
	    double y[], double r[], double e[],
	    double tmean[], Integer *ie, double semean[],
	    Integer imean[], Integer *ir, Integer ifac[],
	    double table[], Integer tdt, Integer *iterm,
	    double rss, Integer *nrdf);
#else
extern void g04cay();
#endif

#ifdef NAG_PROTO
extern void g04caz(Integer inter, Integer nfac, Integer lfac[], Integer mpf,
	     Integer nblock, Integer mrep, Integer itc, Integer nrepl,
	     Integer n, double y[], double r[], double e[],
	     double tmean[], Integer *ie, double semean[],
	     Integer imean[], Integer *ir, Integer ifac[],
	     double table[], Integer tdt, Integer *iterm,
	     double rss, Integer *nrdf, Integer iwk[]);
#else
extern void g04caz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g04czc(double **table, double **tmean, double **e,
	    Integer **imean, double **semean, double **bmean);
#else
     void g04czc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG04 */
