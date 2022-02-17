#ifndef NAGE01
#define NAGE01
#ifdef __cplusplus
extern "C"
{
#endif


/* <nage01.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library e01 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2146 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01bac(Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            Nag_Spline NAG_HUGE *spline, NagError NAG_HUGE *fail);
#else
extern void e01bac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01bec(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[], NagError NAG_HUGE *fail);
#else
extern void e01bec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01bfc(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[], Integer m, double NAG_HUGE px[],
            double NAG_HUGE pf[], NagError NAG_HUGE *fail);
#else
extern void e01bfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01bgc(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[],
            Integer m, double NAG_HUGE px[], double NAG_HUGE pf[], double NAG_HUGE pd[],
            NagError NAG_HUGE *fail);
#else
extern void e01bgc();
#endif

#ifdef NAG_PROTO
extern void e01bgy(double x1, double x2, double f1, double f2,
            double d1, double d2, Integer ne, double NAG_HUGE *xe,
            double NAG_HUGE *fe, double NAG_HUGE *de, Integer NAG_HUGE * next, Integer NAG_HUGE *info);
#else
extern void e01bgy();
#endif

#ifdef NAG_PROTO
extern void e01bgz( Integer n, double NAG_HUGE x[],  double NAG_HUGE f[],   double NAG_HUGE d[],
            Integer incfd,  Integer m, double NAG_HUGE px[],  double NAG_HUGE pf[],
            double NAG_HUGE pd[], Integer NAG_HUGE *info);
#else
extern void e01bgz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01bhc(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[],
            double a, double b, double NAG_HUGE *integral, NagError NAG_HUGE *fail);
#else
extern void e01bhc();
#endif

#ifdef NAG_PROTO
extern double e01bhx(double x1, double x2, double f1, double f2,
              double d1, double d2, double a, double b,
              Integer NAG_HUGE *info);
#else
extern double e01bhx();
#endif

#ifdef NAG_PROTO
extern double e01bhy(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[],
              Integer incfd, Boolean NAG_HUGE *skip, Integer ia, Integer ib,
              Integer NAG_HUGE *info);
#else
extern double e01bhy();
#endif

#ifdef NAG_PROTO
extern double e01bhz(Integer n, double NAG_HUGE x[], double NAG_HUGE f[], double NAG_HUGE d[],
              Integer incfd, double a, double b, Integer NAG_HUGE *info);
#else
extern double e01bhz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01dac(Integer mx, Integer my, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE f[],Nag_2dSpline NAG_HUGE *spline, NagError NAG_HUGE *fail);
     
#else
extern void e01dac();
#endif

#ifdef NAG_PROTO
extern void e01dal(Integer n, Integer ibandw, double NAG_HUGE ufctr[],
            Integer nsets, Integer ntzero,
            double NAG_HUGE theta[], Integer it1, Integer it2);
#else
extern void e01dal();
#endif

#ifdef NAG_PROTO
extern void e01dam(Integer norder, Integer nknots, double xmin,
            double xmax, double NAG_HUGE lambda[],
            Integer jintvl, double NAG_HUGE knot[]);
#else
extern void e01dam();
#endif

#ifdef NAG_PROTO
extern void e01dan(Integer nknots, double NAG_HUGE lambda[],
            double x, Integer NAG_HUGE *jintvl);
#else
extern void e01dan();
#endif

#ifdef NAG_PROTO
extern void e01dap(Integer nxdata, double NAG_HUGE xdata[], double x,
            Integer NAG_HUGE *jintvl);
#else
extern void e01dap();
#endif

#ifdef NAG_PROTO
extern void e01daq(Integer n, Integer nsets, Integer nufctr, double NAG_HUGE ufctr[],
            double NAG_HUGE theta[], Integer it1,
            Integer it2, double NAG_HUGE sfrss[],
            Integer isf1, double NAG_HUGE srss[], Integer isr1);
#else
extern void e01daq();
#endif

#ifdef NAG_PROTO
extern void e01dar(Integer n, Integer ibandw, Integer ixnzst,
            double NAG_HUGE xrow[], Integer nsets, double NAG_HUGE yrow[],
            Integer iy1,  Integer NAG_HUGE *lastcl,
            double NAG_HUGE ufctr[],  double NAG_HUGE theta[],
            Integer it1, Integer it2, double NAG_HUGE wrk[]);
#else
extern void e01dar();
#endif

#ifdef NAG_PROTO
extern void e01das(Integer n, Integer ibandw, double NAG_HUGE ufctr[], Integer NAG_HUGE *info);
#else
extern void e01das();
#endif

#ifdef NAG_PROTO
extern void e01dat(Integer norder, double NAG_HUGE *knot,
            double x, double NAG_HUGE *basis);
#else
extern void e01dat();
#endif

#ifdef NAG_PROTO
extern void e01dau(Integer nxdata, double NAG_HUGE xdata[], double x,
            Integer NAG_HUGE *j);
#else
extern void e01dau();
#endif

#ifdef NAG_PROTO
extern void e01dav(Integer norder, double NAG_HUGE knot[],
            double x, double NAG_HUGE basis[]);
#else
extern void e01dav();
#endif

#ifdef NAG_PROTO
extern void e01daw(Integer norder, double xmin, double xmax,
            Integer m, double NAG_HUGE x[], Integer nsets, double NAG_HUGE f[],
            Integer if1, Integer if2, Integer nknots,
            double NAG_HUGE lambda[], double NAG_HUGE c[], Integer ic1,
            Integer ic2, double NAG_HUGE ufctr[],  double NAG_HUGE knot[],
            double NAG_HUGE *xrow, double NAG_HUGE *wrk,   Integer NAG_HUGE *info);
#else
extern void e01daw();
#endif

#ifdef NAG_PROTO              
extern void e01dax(double xmin, double xmax, Integer m, double NAG_HUGE x[],
            double NAG_HUGE w[], Integer iw1, Integer NAG_HUGE *ifnzwt,
            Integer NAG_HUGE *mnzwt, Integer NAG_HUGE *mdnzwt, Integer NAG_HUGE *info);
#else
extern void e01dax();
#endif

#ifdef NAG_PROTO
extern void e01day(Integer norder, Integer nknots, Integer m, double NAG_HUGE x[],
            double NAG_HUGE w[], Integer iw1, Integer ifnzwt,
            Integer mdnzwt, double NAG_HUGE lambda[]);
#else
extern void e01day();
#endif

#ifdef NAG_PROTO
extern void e01daz(Integer  nxordr, Integer nxknts, double xmin,
            double xmax, Integer nyordr, Integer nyknts,
            double ymin, double ymax, Integer mx,
            double NAG_HUGE x[], Integer my, double NAG_HUGE y[], double NAG_HUGE f[],
            Integer if1, Integer if2, double NAG_HUGE xlam[],
            double NAG_HUGE ylam[],double NAG_HUGE c[], Integer ic1, Integer ic2, 
            double NAG_HUGE xufctr[],  double NAG_HUGE yufctr[],
            double NAG_HUGE e[], Integer ie1, Integer ie2,
            double NAG_HUGE xknot[], double NAG_HUGE yknot[], double NAG_HUGE xrow[],
            double NAG_HUGE yrow[], double NAG_HUGE wrk[], Integer NAG_HUGE *info);
#else
extern void e01daz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01sac(Nag_2d_Scat_Method method, Integer m, double NAG_HUGE x[], double NAG_HUGE y[],
            double NAG_HUGE f[], Nag_Scat_Struct NAG_HUGE *comm, Nag_E01_Opt NAG_HUGE *optional,
            NagError NAG_HUGE *fail);
#else
extern void e01sac();
#endif

#ifdef NAG_PROTO
extern void e01saq(Integer nfrst, Integer nlast, Integer kk, Integer NAG_HUGE iarr[]);
#else
extern void e01saq();
#endif

#ifdef NAG_PROTO
extern Integer e01sar(Integer nvertx, Integer nabor, Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[]);
#else
extern Integer e01sar();
#endif

#ifdef NAG_PROTO
extern void e01sas(Integer nin1, Integer nin2, Integer nout1, Integer nout2,
            Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[]);
#else
extern void e01sas();
#endif

#ifdef NAG_PROTO
extern Boolean e01sat(Integer in1, Integer in2, Integer io1, Integer io2,
               double NAG_HUGE x[], double NAG_HUGE y[]);
#else
extern Boolean e01sat();
#endif

#ifdef NAG_PROTO
extern void e01sau(Integer kk, Integer i1, Integer i2, Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[]);
#else
extern void e01sau();
#endif

#ifdef NAG_PROTO
extern void e01sav(Integer kk, Integer i1, Integer i2, Integer i3,
            Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[]);
#else
extern void e01sav();
#endif

#ifdef NAG_PROTO
extern void e01saw(Integer nst, double px, double py, double NAG_HUGE x[], double NAG_HUGE y[],
            Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[], Integer NAG_HUGE *i1, Integer NAG_HUGE *i2,
            Integer NAG_HUGE *i3);
#else
extern void e01saw();
#endif

#ifdef NAG_PROTO
extern void e01sax(Integer kk, double NAG_HUGE x[], double NAG_HUGE y[], Integer NAG_HUGE iadj[],
            Integer NAG_HUGE iend[], Integer NAG_HUGE *ier);
#else
extern void e01sax();
#endif

#ifdef NAG_PROTO
extern void e01say(Integer n, double NAG_HUGE x[], double NAG_HUGE y[], Integer NAG_HUGE iadj[],
            Integer NAG_HUGE iend[], Integer NAG_HUGE *ier);
#else
extern void e01say();
#endif

#ifdef NAG_PROTO
extern void e01saz(Integer n, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE z[], Integer NAG_HUGE iadj[],
            Integer NAG_HUGE iend[], double eps, Integer NAG_HUGE *nit, double NAG_HUGE zxzy[],
            Integer NAG_HUGE *coinc1, Integer NAG_HUGE *coinc2, Integer NAG_HUGE *ier);
#else
extern void e01saz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01sbc(Nag_Scat_Struct NAG_HUGE *comm, Integer n, double NAG_HUGE px[], double NAG_HUGE py[],
            double NAG_HUGE pf[], NagError NAG_HUGE *fail);
#else
extern void e01sbc();
#endif

#ifdef NAG_PROTO
extern void e01sby(double x, double y, double x1, double x2,
            double x3, double y1, double y2, double y3,
            double z1, double z2, double z3, double zx1,
            double zx2, double zx3, double zy1, double zy2,
            double zy3, Integer iflag, double NAG_HUGE *w, double NAG_HUGE *wx,
            double NAG_HUGE *wy, Integer NAG_HUGE *ier);
#else
extern void e01sby();
#endif

#ifdef NAG_PROTO
extern void e01sbz(Integer n, double px, double py, double NAG_HUGE x[],
            double NAG_HUGE y[], double NAG_HUGE z[], Integer NAG_HUGE iadj[], Integer NAG_HUGE iend[],
            double NAG_HUGE zxzy[], Integer NAG_HUGE *ist, Integer iflag,
            double NAG_HUGE *pz, double NAG_HUGE *dzx, double NAG_HUGE *dzy, Integer NAG_HUGE *ier);
#else
extern void e01sbz();
#endif

#ifdef NAG_PROTO
extern void e01sez(Integer n, double NAG_HUGE x[], double NAG_HUGE y[], double NAG_HUGE z[], double rnq,
            double NAG_HUGE fnodes[], Integer NAG_HUGE *minnq, double NAG_HUGE c[], double NAG_HUGE b[]);
#else
extern void e01sez();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL e01szc(Nag_Scat_Struct NAG_HUGE *comm);
#else
extern void e01szc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGE01 */
