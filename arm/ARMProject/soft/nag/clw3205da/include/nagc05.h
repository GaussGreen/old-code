#ifndef NAGC05
#define NAGC05
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagc05.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library c05 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2142 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05adc(double a, double b, double NAG_HUGE *x, 
            NAG_C05ADC_FUN f, double xtol, 
            double ftol, NagError NAG_HUGE *fail);
#else
extern void c05adc();
#endif

#ifdef NAG_PROTO
extern void c05azc(double NAG_HUGE *x, double NAG_HUGE *y, double NAG_HUGE *fx, double xtol,
            Integer ir, double NAG_HUGE c[], Integer NAG_HUGE * ind, NagError NAG_HUGE *fail);
#else
extern void c05azc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05nbc(Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[], 
            NAG_C05NBC_FUN f,
            double xtol, NagError NAG_HUGE *fail);
#else
extern void c05nbc();
#endif

#ifdef NAG_PROTO
extern void c05ncs(NAG_C05NBC_FUN f,
              Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
              double xtol, Integer maxfev, Integer ml, Integer mu,
              double epsfcn, double NAG_HUGE diag[], Integer mode,
              double factor, Integer nprint, Integer NAG_HUGE *info,
              Integer NAG_HUGE *nfev, double NAG_HUGE fjac[], Integer tdfjac,
              double NAG_HUGE r[],  double NAG_HUGE qtf[], double NAG_HUGE wa1[],
              double NAG_HUGE wa2[], double NAG_HUGE wa3[], double NAG_HUGE wa4[]);
#else
extern void c05ncs();
#endif

#ifdef NAG_PROTO
extern double c05nct(Integer n, double NAG_HUGE x[]);
#else
extern double c05nct();
#endif

#ifdef NAG_PROTO
extern void c05ncu(Integer n, double NAG_HUGE r[], double NAG_HUGE diag[],
            double NAG_HUGE qtb[], double delta, double NAG_HUGE x[],
            double NAG_HUGE wa1[], double NAG_HUGE wa2[]);
#else
extern void c05ncu();
#endif

#ifdef NAG_PROTO
extern void c05ncv(NAG_C05NBC_FUN f,
            Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdfjac, Integer NAG_HUGE *flag,
            Integer ml, Integer mu, double epsfcn, double NAG_HUGE wa1[],
            double NAG_HUGE wa2[]);
#else
extern void c05ncv();
#endif

#ifdef NAG_PROTO
extern void c05ncw(Integer n, double NAG_HUGE q[], Integer tdq, double NAG_HUGE wa[]);
#else
extern void c05ncw();
#endif

#ifdef NAG_PROTO
extern void c05ncx(Integer n, double NAG_HUGE a[], Integer tda,
            double NAG_HUGE rdiag[], double NAG_HUGE acnorm[]);
#else
extern void c05ncx();
#endif

#ifdef NAG_PROTO
extern void c05ncy(Integer m, Integer n, double NAG_HUGE a[], Integer tda,
            double NAG_HUGE v[], double NAG_HUGE w[]);
#else
extern void c05ncy();
#endif

#ifdef NAG_PROTO
extern void c05ncz(Integer m, Integer n, double NAG_HUGE s[],
            double NAG_HUGE u[], double NAG_HUGE v[], double NAG_HUGE w[], Boolean NAG_HUGE * sing);
#else
extern void c05ncz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05pbc( Integer n, double NAG_HUGE x[],
            double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdfjac,
            NAG_C05PBC_FUN f,
            double xtol, NagError NAG_HUGE *fail);
#else
extern void c05pbc();
#endif

#ifdef NAG_PROTO
extern void c05pcz(NAG_C05PBC_FUN f,
            Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdfjac, double xtol,
            Integer maxfev, double NAG_HUGE diag[], Integer mode,
            double factor, Integer nprint, Integer NAG_HUGE *info,
            Integer NAG_HUGE *nfev, Integer NAG_HUGE *njev, double NAG_HUGE r[],
            double NAG_HUGE qtf[], double NAG_HUGE wa1[], double NAG_HUGE wa2[],
            double NAG_HUGE wa3[], double NAG_HUGE wa4[]);
#else
extern void c05pcz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05zbc( Integer n, double NAG_HUGE x[],
            double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdfjac,
            NAG_C05ZBC_FUN lsqfun, NagError NAG_HUGE *fail);
#else
extern void c05zbc();
#endif

/* Multi-threading equivalents of above */
#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05sdc(double a, double b, double NAG_HUGE *x, 
            NAG_C05SDC_FUN f, double xtol, 
            double ftol, Nag_User *comm, NagError NAG_HUGE *fail);
#else
extern void c05sdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05tbc(Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[], 
            NAG_C05TBC_FUN f,
            double xtol, Nag_User *comm, NagError NAG_HUGE *fail);
#else
extern void c05tbc();
#endif

#ifdef NAG_PROTO
extern void c05tcs(NAG_C05TBC_FUN f,
              Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
              double xtol, Integer maxfev, Integer ml, Integer mu,
              double epsfcn, double NAG_HUGE diag[], Integer mode,
              double factor, Integer nprint, Integer NAG_HUGE *info,
              Integer NAG_HUGE *nfev, double NAG_HUGE fjac[], Integer tdfjac,
              double NAG_HUGE r[],  double NAG_HUGE qtf[], double NAG_HUGE wa1[],
              double NAG_HUGE wa2[], double NAG_HUGE wa3[], double NAG_HUGE wa4[], Nag_User *comm);
#else
extern void c05tcs();
#endif

#ifdef NAG_PROTO
extern void c05tcv(NAG_C05TBC_FUN f,
            Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdfjac, Integer NAG_HUGE *flag,
            Integer ml, Integer mu, double epsfcn, double NAG_HUGE wa1[],
            double NAG_HUGE wa2[], Nag_User *comm);
#else
extern void c05tcv();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05ubc( Integer n, double NAG_HUGE x[],
            double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdfjac,
            NAG_C05UBC_FUN f,
            double xtol, Nag_User *comm, NagError NAG_HUGE *fail);
#else
extern void c05ubc();
#endif

#ifdef NAG_PROTO
extern void c05ucz(NAG_C05UBC_FUN f,
            Integer n, double NAG_HUGE x[], double NAG_HUGE fvec[],
            double NAG_HUGE fjac[], Integer tdfjac, double xtol,
            Integer maxfev, double NAG_HUGE diag[], Integer mode,
            double factor, Integer nprint, Integer NAG_HUGE *info,
            Integer NAG_HUGE *nfev, Integer NAG_HUGE *njev, double NAG_HUGE r[],
            double NAG_HUGE qtf[], double NAG_HUGE wa1[], double NAG_HUGE wa2[],
            double NAG_HUGE wa3[], double NAG_HUGE wa4[], Nag_User *comm);
#else
extern void c05ucz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c05zcc( Integer n, double NAG_HUGE x[],
            double NAG_HUGE fvec[], double NAG_HUGE fjac[], Integer tdfjac,
            NAG_C05ZCC_FUN lsqfun, Nag_User *comm, NagError NAG_HUGE *fail);
#else
extern void c05zcc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGC05 */
