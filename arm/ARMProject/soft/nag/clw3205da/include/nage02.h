#ifndef NAGE02
#define NAGE02
#ifdef __cplusplus
extern "C"
{
#endif


/* <nage02.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library e02 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2147 (Feb 1998).
 */

#ifdef NAG_PROTO 
extern NAG_DLL_EXPIMP void NAG_CALL e02adc(Integer m,Integer kplus1,
            Integer nrows,double x[],double y[],
            double w[],double a[],double s[],NagError NAG_HUGE *fail);
#else
extern void e02adc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02afc(Integer nplus1,double f[],
                                            double a[],NagError NAG_HUGE *fail);
#else
extern void e02afc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02aec(Integer nplus1,double a[],double xcap,double *p,NagError NAG_HUGE *fail);
#else
extern void e02aec();
#endif


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02bac(Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE weights[], double NAG_HUGE *ss, Nag_Spline NAG_HUGE *spline,
            NagError NAG_HUGE *fail);
#else
extern void e02bac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02bbc(double x, double NAG_HUGE *s, Nag_Spline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02bbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02bcc(Nag_DerivType derivs, double x, double s[4], 
            Nag_Spline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02bcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02bdc(Nag_Spline NAG_HUGE *spline, double NAG_HUGE *integral, 
            NagError NAG_HUGE *fail);
#else
extern void e02bdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02bec(Nag_Start start, Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE weights[], double s, Integer nest,
            double NAG_HUGE *fp, Nag_Comm NAG_HUGE *warmstartinf,
            Nag_Spline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02bec();
#endif

#ifdef NAG_PROTO
extern void e02bev(double NAG_HUGE t[], Integer k, double x,
            Integer l, double NAG_HUGE h[]);
#else
extern void e02bev();
#endif

#ifdef NAG_PROTO
extern void e02bew(double NAG_HUGE t[], Integer n, Integer k2, double NAG_HUGE b[]);
#else
extern void e02bew();
#endif

#ifdef NAG_PROTO
extern double e02bex(double NAG_HUGE *p1, double NAG_HUGE *f1, double p2, double f2,
              double NAG_HUGE *p3, double NAG_HUGE *f3);
#else
extern double e02bex();
#endif

#ifdef NAG_PROTO
extern void e02bey(double NAG_HUGE *x,double NAG_HUGE t[], Integer NAG_HUGE *n,
            double NAG_HUGE *fpint, Integer NAG_HUGE nrdata[], Integer NAG_HUGE *nrint);
#else                                     
extern void e02bey();
#endif

#ifdef NAG_PROTO
extern void e02bez(Integer NAG_HUGE *iopt, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE w[],
            Integer m, double NAG_HUGE *xb, double NAG_HUGE *xe, Integer NAG_HUGE *k,
            double s, Integer nest, double NAG_HUGE *tol,
            Integer NAG_HUGE *maxit, Integer k1, Integer k2, Integer NAG_HUGE *n,
            double NAG_HUGE t[], double NAG_HUGE c[], double NAG_HUGE *fp,
            double NAG_HUGE fpint[], double NAG_HUGE z[], double NAG_HUGE a[], double NAG_HUGE b[],
            double NAG_HUGE g[], double NAG_HUGE q[], Integer NAG_HUGE nrdata[], Integer NAG_HUGE *info);
#else
extern void e02bez();
#endif

#ifdef NAG_PROTO
extern void e02dbc(Integer m, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE ff[], 
            Integer NAG_HUGE point[], Nag_2dSpline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02dbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02dcc(Nag_Start start, Integer mx, double NAG_HUGE x[], Integer my,
            double NAG_HUGE y[], double NAG_HUGE f[], double s, Integer nxest,
            Integer nyest, double NAG_HUGE *fp, Nag_Comm NAG_HUGE *warmstartinf,
            Nag_2dSpline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02dcc();
#endif

#ifdef NAG_PROTO
extern void e02dcy(Integer NAG_HUGE *ifsx, Integer NAG_HUGE *ifsy, Integer NAG_HUGE *ifbx, Integer NAG_HUGE *ifby,
            double NAG_HUGE x[], Integer mx, double NAG_HUGE y[], Integer my,
            double NAG_HUGE z[], Integer kx, Integer ky, double NAG_HUGE tx[],
            Integer nx, double NAG_HUGE ty[], Integer ny, double p,
            double NAG_HUGE c[], double NAG_HUGE *fp, double NAG_HUGE fpx[],
            double NAG_HUGE fpy[], Integer kx1,
            Integer kx2, Integer ky1, Integer ky2, double NAG_HUGE spx[],
            double NAG_HUGE spy[], double NAG_HUGE right[], double NAG_HUGE q[], double NAG_HUGE ax[],
            double NAG_HUGE ay[], double NAG_HUGE bx[], double NAG_HUGE by[], Integer NAG_HUGE nrx[],
            Integer NAG_HUGE nry[]);
#else
extern void e02dcy();
#endif

#ifdef NAG_PROTO
extern void e02dcz(Integer NAG_HUGE *iopt, double NAG_HUGE x[], Integer mx, double NAG_HUGE y[],
            Integer my, double NAG_HUGE z[], double NAG_HUGE *xb, double NAG_HUGE *xe,
            double NAG_HUGE *yb, double NAG_HUGE *ye, Integer kx, Integer ky,
            double s, Integer nxest, Integer nyest,
            double tol, Integer maxit, Integer NAG_HUGE *nx,
            double NAG_HUGE tx[], Integer NAG_HUGE *ny, double NAG_HUGE ty[], double NAG_HUGE c[],
            double NAG_HUGE *fp, double NAG_HUGE *fp0, double NAG_HUGE *fpold,
            double NAG_HUGE *reducx, double NAG_HUGE *reducy, double NAG_HUGE fpintx[],
            double NAG_HUGE fpinty[], Integer NAG_HUGE *lastdi, Integer NAG_HUGE *nplusx,
            Integer NAG_HUGE *nplusy, Integer NAG_HUGE nrx[], Integer NAG_HUGE nry[], Integer NAG_HUGE nrdatx[],
            Integer NAG_HUGE nrdaty[],double NAG_HUGE work2[],
            Integer NAG_HUGE *info);
#else
extern void e02dcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02ddc(Nag_Start start, Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE f[], double NAG_HUGE weights[], double s, 
            Integer nxest, Integer nyest, double NAG_HUGE *fp,
            Integer NAG_HUGE *rank, double NAG_HUGE *warmstartinf,
            Nag_2dSpline NAG_HUGE *spline,
            NagError NAG_HUGE *fail);
#else
extern void e02ddc();
#endif

#ifdef NAG_PROTO
extern void e02ddy(double NAG_HUGE a[], double NAG_HUGE f[], Integer n, Integer m,
            Integer na, double tol, double NAG_HUGE c[], double NAG_HUGE *sq,
            Integer NAG_HUGE *rank, double NAG_HUGE aa[], double NAG_HUGE ff[], double NAG_HUGE h[]);
#else
extern void e02ddy();
#endif

#ifdef NAG_PROTO
extern void e02ddz(Integer iopt, Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE z[], double NAG_HUGE w[], double xb, double xe,
            double yb, double ye, Integer kxx, Integer kyy,
            double s, Integer nxest, Integer  nyest,
            double eta, double tol, Integer maxit,
            Integer km1, Integer km2, Integer ib1,
            Integer ib3, Integer nc, 
            Integer NAG_HUGE *nx0, double NAG_HUGE tx[], Integer NAG_HUGE *ny0, double NAG_HUGE ty[],
            double NAG_HUGE c[], double NAG_HUGE *fp, double NAG_HUGE *fp0, Integer NAG_HUGE *rank,
            double NAG_HUGE fpint[], double NAG_HUGE coord[], double NAG_HUGE f[],
            double NAG_HUGE ff[], double NAG_HUGE a[], double NAG_HUGE q[], double NAG_HUGE bx[],
            double NAG_HUGE by[], double NAG_HUGE spx[], double NAG_HUGE spy[], double NAG_HUGE h[],
            Integer NAG_HUGE index[], Integer NAG_HUGE adres[], Integer NAG_HUGE *ier);
#else
extern void e02ddz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02dec(Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE ff[], Nag_2dSpline NAG_HUGE *spline,
            NagError NAG_HUGE *fail);
#else
extern void e02dec();
#endif

#ifdef NAG_PROTO
extern void e02dez(Integer m, Integer nxknts,
            double xmax, double NAG_HUGE xlam[],
            Integer nyknts, double ymax,
            double NAG_HUGE ylam[], double NAG_HUGE c[],
            Integer ic1, Integer ic2, double NAG_HUGE x[],
            double NAG_HUGE y[], double NAG_HUGE s[], double NAG_HUGE ydep[], Integer NAG_HUGE iybsi[]);
#else
extern void e02dez();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e02dfc(Integer mx, Integer my, double NAG_HUGE x[], double NAG_HUGE y[], 
            double NAG_HUGE ff[], Nag_2dSpline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e02dfc();
#endif

#ifdef NAG_PROTO
extern void e02dfv(Integer nknots, double NAG_HUGE lambda[],
            double xmax, double x, Integer NAG_HUGE *intvl);
#else
extern void e02dfv();
#endif

#ifdef NAG_PROTO
extern void e02dfw(Integer nknots, double xmin, double xmax,
            double NAG_HUGE lambda[], Integer dlmbda, Integer NAG_HUGE *error);
#else
extern void e02dfw();
#endif

#ifdef NAG_PROTO
extern void e02dfx(Integer norder, Integer q, double NAG_HUGE lambda[],
            double xmax, Integer m,
            double NAG_HUGE x[], Integer NAG_HUGE intvl[], Integer NAG_HUGE *nindex,
            Integer NAG_HUGE index[]);
#else
extern void e02dfx();
#endif

#ifdef NAG_PROTO
extern void e02dfy(Integer nxordr, Integer nxknts,
            double xmax, double NAG_HUGE xlam[],
            Integer nyordr, Integer nyknts,
            double ymax, double NAG_HUGE ylam[], 
            double NAG_HUGE c[], Integer ic1, Integer ic2,
            Integer mx, double NAG_HUGE x[], Integer my, double NAG_HUGE y[],
            double NAG_HUGE s[], Integer is1, Integer is2,
            double NAG_HUGE ybasis[], double NAG_HUGE xbasis[],
            double NAG_HUGE ydep[], Integer NAG_HUGE iyint[], Integer NAG_HUGE iybsi[]);
#else
extern void e02dfy();
#endif

#ifdef NAG_PROTO
extern void e02dfz(Integer nxordr, Integer nxknts,
            double xmax, double NAG_HUGE xlam[], 
            Integer nyordr, Integer nyknts,
            double ymax, double NAG_HUGE ylam[], 
            double NAG_HUGE c[], Integer ic1, Integer ic2, 
            Integer mx, double NAG_HUGE x[], Integer my, double NAG_HUGE y[],
            double NAG_HUGE s[], Integer is1, Integer is2, 
            double NAG_HUGE work1[], Integer NAG_HUGE work2[]);
#else
extern void e02dfz();
#endif

#ifdef NAG_PROTO
extern void e02zac(Integer nx, Integer ny, double NAG_HUGE lamda[],
            double NAG_HUGE mu[], Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            Integer NAG_HUGE point[], Integer NAG_HUGE adres[], 
            NagError NAG_HUGE *fail);
#else
extern void e02zac();
#endif

#ifdef NAG_PROTO
extern Integer e02zaz(Integer n, double NAG_HUGE x[], double t);
#else
extern Integer e02zaz();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGE02 */
