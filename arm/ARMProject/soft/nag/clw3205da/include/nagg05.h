#ifndef NAGG05
#define NAGG05
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg05.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g05 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2157 (Feb 1998).
 */

/* Note:
 * g05dgc is called by g05c{b,c,f} stringents only
 */

#include <naglg05.h>    /* Needed for Nag_ag05ca_type etc. */

#ifdef NAG_THREAD_SAFE
#include <pthread.h>
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g05cac(void);
#else
extern double g05cac();
#endif

#ifdef NAG_PROTO
extern void g05cay(Boolean reinit);
#else
extern void g05cay();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05caz(Boolean NAG_HUGE *init,
                                           Nag_ag05ca_type NAG_HUGE **ag05ca,
                                           Nag_bg05ca_type NAG_HUGE **bg05ca,
                                           Nag_cg05ca_type NAG_HUGE **cg05ca,
                                           Nag_dg05ca_type NAG_HUGE **dg05ca);
#else
extern void g05caz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05cbc(Integer seed);
#else
extern void g05cbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ccc(void);
#else
extern void g05ccc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05cfc(Integer NAG_HUGE istate[], double NAG_HUGE xstate[]);
#else
extern void g05cfc();
#endif

#ifdef NAG_PROTO
extern void g05cfz(Integer NAG_HUGE *ia, double NAG_HUGE *xb);
#else
extern void g05cfz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05cgc(Integer NAG_HUGE istate[], double NAG_HUGE xstate[], NagError NAG_HUGE *fail);
#else
extern void g05cgc();
#endif

#ifdef NAG_PROTO
extern void g05cgz(Integer NAG_HUGE *ia, double NAG_HUGE *xb);
#else
extern void g05cgz();
#endif

#ifdef NAG_PROTO
extern void g05cxc(); 
#else
extern void g05cxc();
#endif


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g05dac(double mu, double sd);
#else
extern double g05dac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g05dbc(double mu);
#else
extern double g05dbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL g05ddc(double mu, double sd);
#else
extern double g05ddc();
#endif

#ifdef NAG_PROTO
extern double g05dfc(double a, double b);
#else
extern double g05dfc();
#endif

#ifdef NAG_PROTO
extern double g05dgc(double a, double b, NagError NAG_HUGE *fail);
#else
extern double g05dgc();
#endif

#ifdef NAG_PROTO
extern Integer g05drc(double alamda,  Integer NAG_HUGE *ifail);
#else
extern Integer g05drc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP Integer NAG_CALL g05dyc(Integer m, Integer n);
#else
extern Integer g05dyc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05eac(double NAG_HUGE a[], Integer n, double NAG_HUGE c[], Integer tdc, 
            double eps, double NAG_HUGE *NAG_HUGE *r,
            NagError NAG_HUGE *fail);
#else
extern void g05eac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ecc(double t, double NAG_HUGE *NAG_HUGE *r, NagError NAG_HUGE *fail);
#else
extern void g05ecc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05edc(Integer n, double p, double NAG_HUGE *NAG_HUGE *r, 
            NagError NAG_HUGE *fail);
#else
extern void g05edc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ehc(Integer NAG_HUGE index[], Integer n, NagError NAG_HUGE *fail);
#else
extern void g05ehc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ejc(Integer NAG_HUGE ia[], Integer n, Integer NAG_HUGE iz[], Integer m, 
            NagError NAG_HUGE *fail);
#else
extern void g05ejc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05exc(double NAG_HUGE p[], Integer np, Integer sizep, Nag_DiscreteDistrib df, 
            double NAG_HUGE *NAG_HUGE *r, NagError NAG_HUGE *fail);
#else
extern void g05exc();
#endif

#ifdef NAG_PROTO
extern void g05exz(Integer m, Integer n, double NAG_HUGE r[], Integer nr);
#else
extern void g05exz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP Integer NAG_CALL g05eyc(double NAG_HUGE r[]);
#else
extern Integer g05eyc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ezc(double NAG_HUGE z[], double NAG_HUGE r[]);
#else
extern void g05ezc();
#endif

#ifdef NAG_PROTO
extern void g05fac(double a, double b,  Integer n,  double x[]);
#else
     extern void g05fac();
#endif

#ifdef NAG_PROTO
extern void g05fdc(double a, double b, Integer n, double NAG_HUGE *x);
#else
extern void g05fdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05fec(double a, double b, Integer n, double NAG_HUGE x[],
            NagError NAG_HUGE *fail);
#else
extern void g05fec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05ffc(double a, double b, Integer n, double NAG_HUGE x[],
            NagError NAG_HUGE *fail);
#else
extern void g05ffc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g05hac(Boolean start, Integer p, Integer q, double NAG_HUGE phi[],
            double NAG_HUGE theta[], double mean, double vara, Integer n,
            double NAG_HUGE w[], double NAG_HUGE ref[], NagError NAG_HUGE *fail);
#else
extern void g05hac();
#endif

/* These functions are defined in g05czzt.c, they
 * are for use by g05c* stringents.
 */
#ifdef NAG_PROTO
extern  NAG_DLL_EXPIMP Nag_bg05ca_type NAG_CALL get_bg05ca(void);
extern  NAG_DLL_EXPIMP Nag_dg05ca_type NAG_CALL get_dg05ca(void);
#else
extern  Nag_bg05ca_type NAG_CALL get_bg05ca();
extern  Nag_dg05ca_type NAG_CALL get_dg05ca();
#endif



#ifdef __cplusplus
}
#endif
#endif /* not NAGG05 */
