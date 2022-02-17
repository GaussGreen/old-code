#ifndef NAGC02
#define NAGC02
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagc02.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library c02 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2141 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c02afc(Integer n, Complex NAG_HUGE a[], Boolean scale, Complex NAG_HUGE z[], NagError NAG_HUGE *fail);
#else
extern void c02afc();
#endif

#ifdef NAG_PROTO
extern void c02afw(Complex NAG_HUGE *a, Complex NAG_HUGE *b, Complex NAG_HUGE *c, Boolean NAG_HUGE *fail);
#else
extern void c02afw();
#endif

#ifdef NAG_PROTO
extern void c02afx(Complex afx, Complex bfx, Complex cfx, Complex NAG_HUGE *zsm, Complex NAG_HUGE *zlg,
            Ac02af NAG_HUGE *c02afz_1, Bc02af NAG_HUGE *c02afz_2, Cc02ag  NAG_HUGE *c02agx_1);
#else
extern void c02afx();
#endif

#ifdef NAG_PROTO
extern void c02afy(Complex dfy, Integer ndeg, Complex NAG_HUGE *a, Complex NAG_HUGE *p, Complex NAG_HUGE *pprime,
            Complex NAG_HUGE *pdprim, double NAG_HUGE *error, Complex NAG_HUGE *deflat, Ac02af NAG_HUGE *c02afz_1,
            Bc02af c02afz_2);
#else
extern void c02afy();
#endif

#ifdef NAG_PROTO
extern void c02afz(Complex NAG_HUGE *a, Integer ndeg, Boolean NAG_HUGE *scale, Complex NAG_HUGE *z,
            int NAG_HUGE *localerror);
#else
extern void c02afz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL c02agc(Integer n, double NAG_HUGE a[], Boolean scale, Complex NAG_HUGE z[], NagError NAG_HUGE *fail);
#else
extern void c02agc();
#endif

#ifdef NAG_PROTO
extern double c02agr(double x, Integer exp_agr, Cc02ag NAG_HUGE *c02agx_1);
#else
extern double c02agr();
#endif

#ifdef NAG_PROTO
extern Boolean c02ags(double x);
#else
extern Boolean c02ags();
#endif

#ifdef NAG_PROTO
extern void c02agt(double a0, double b0, double c0, Complex NAG_HUGE *zsm, Complex NAG_HUGE *zlg,
            Ac02ag NAG_HUGE *c02agz_1, Bc02ag NAG_HUGE *c02agz_2, Cc02ag NAG_HUGE *c02agx_1);
#else
extern void c02agt();
#endif

#ifdef NAG_PROTO
extern void c02agu(Complex a0, Complex b0, Complex c0, Complex NAG_HUGE *zsm, Complex NAG_HUGE *zlg,
            Ac02ag NAG_HUGE *c02agz_1, Bc02ag NAG_HUGE *c02agz_2, Cc02ag NAG_HUGE *c02agx_1);
#else
extern void c02agu();
#endif

#ifdef NAG_PROTO
extern void c02agv(double dx, Integer ndeg, double NAG_HUGE *a, double NAG_HUGE *p, double NAG_HUGE *pprime,
            double NAG_HUGE *pdprim, double NAG_HUGE *error, double NAG_HUGE *deflat, Ac02ag NAG_HUGE *c02agz_1,
            Bc02ag NAG_HUGE *c02agz_2);
#else
extern void c02agv();
#endif

#ifdef NAG_PROTO
extern void c02agw(Complex xgw, Integer ndeg, double NAG_HUGE *a, Complex NAG_HUGE *p, Complex NAG_HUGE *pprime,
            Complex NAG_HUGE *pdprim, double NAG_HUGE *error, double NAG_HUGE *deflat, Ac02ag NAG_HUGE *c02agz_1,
            Bc02ag NAG_HUGE *c02agz_2);
#else
extern void c02agw();
#endif

#ifdef NAG_PROTO
extern Integer c02agx(double dx, Cc02ag NAG_HUGE *c02agx_1);
#else
extern Integer c02agx();
#endif

#ifdef NAG_PROTO
extern double c02agy(double dx, Integer exp_agy, Cc02ag NAG_HUGE *c02agx_1);
#else
extern double c02agy();
#endif

#ifdef NAG_PROTO
extern void c02agz(double NAG_HUGE *a, Integer ndeg, Boolean NAG_HUGE *scale, Complex NAG_HUGE *z,
            int NAG_HUGE *localerror);
#else
extern void c02agz();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGC02 */
