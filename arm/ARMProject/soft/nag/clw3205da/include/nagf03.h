#ifndef NAGF03
#define NAGF03
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf03.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f03 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2151 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f03abc(double NAG_HUGE *a, Integer tda, Integer n, double NAG_HUGE *det, double NAG_HUGE *p,
            NagError NAG_HUGE *fail);
#else
extern void f03abc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f03aec(Integer n, double NAG_HUGE *a, Integer tda, double NAG_HUGE *p, double NAG_HUGE *detf, 
            Integer NAG_HUGE *dete, NagError NAG_HUGE *fail);
#else
extern void f03aec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f03afc(Integer n, double NAG_HUGE *a, Integer tda,  Integer NAG_HUGE *pivot,
            double NAG_HUGE *detf, Integer NAG_HUGE *dete, NagError NAG_HUGE *fail);
#else
extern void f03afc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f03ahc(Integer n, Complex NAG_HUGE *a, Integer tda, Integer NAG_HUGE *pivot,
            Complex NAG_HUGE *det, Integer NAG_HUGE *dete, NagError NAG_HUGE *fail);
#else
extern void f03ahc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGF03 */
