#ifndef NAGG10
#define NAGG10
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg10.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g10 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2159 (Feb 1998).
 */


#include <nag_stddef.h>

#ifdef NAG_PROTO
extern void g10bbc(Nag_Smooth_Type smoother, Integer n, double NAG_HUGE y[], double NAG_HUGE smooth[],
            double NAG_HUGE rough[], NagError NAG_HUGE *fail);
#else
extern void g10bbc();
#endif

#ifdef NAG_PROTO
extern void g10bbo(double NAG_HUGE y[], Integer n, Integer m, double NAG_HUGE run[], Integer NAG_HUGE ind[]);
#else
extern void g10bbo();
#endif

#ifdef NAG_PROTO
extern void g10bbp(double x1, double x2, double x3, double NAG_HUGE *xmed,
            Boolean NAG_HUGE * change);
#else
extern void g10bbp();
#endif

#ifdef NAG_PROTO
extern void g10bbq(double NAG_HUGE y[], Integer n, Boolean NAG_HUGE *change);
#else
extern void g10bbq();
#endif

#ifdef NAG_PROTO
extern void g10bbr(double NAG_HUGE y[], Integer n);
#else
extern void g10bbr();
#endif

#ifdef NAG_PROTO
extern void g10bbs(double NAG_HUGE y[], Integer n);
#else
extern void g10bbs();
#endif

#ifdef NAG_PROTO
extern void g10bbt(double NAG_HUGE y[], Integer n);
#else
extern void g10bbt();
#endif

#ifdef NAG_PROTO
extern void g10bbu(double NAG_HUGE y[], Integer n, double NAG_HUGE work[]);
#else
extern void g10bbu();
#endif

#ifdef NAG_PROTO
extern void g10bbv(double NAG_HUGE y[], Integer n, double NAG_HUGE *endsav, double NAG_HUGE work[]);
#else
extern void g10bbv();
#endif

#ifdef NAG_PROTO
extern void g10bbw(double NAG_HUGE y[], Integer  n, Boolean NAG_HUGE *change);
#else
extern void g10bbw();
#endif

#ifdef NAG_PROTO
extern void g10bbx(double NAG_HUGE y[], Integer n, double endsav);
#else
extern void g10bbx();
#endif

#ifdef NAG_PROTO
extern void g10bby(double NAG_HUGE y[], Integer n);
#else
extern void g10bby();
#endif

#ifdef NAG_PROTO
extern void g10bbz(double NAG_HUGE y[], Integer n);
#else
extern void g10bbz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g10cac(Nag_Smooth_Type smoother, Integer n, double NAG_HUGE y[], double NAG_HUGE smooth[],
            double NAG_HUGE rough[], NagError NAG_HUGE *fail);
#else
extern void g10cac();
#endif

#ifdef NAG_PROTO
extern void g10cao(double NAG_HUGE y[], Integer n, Integer m, double NAG_HUGE run[], size_t NAG_HUGE ind[]);
#else
extern void g10cao();
#endif

#ifdef NAG_PROTO
extern void g10cap(double x1, double x2, double x3, double NAG_HUGE *xmed,
            Boolean NAG_HUGE * change);
#else
extern void g10cap();
#endif

#ifdef NAG_PROTO
extern void g10caq(double NAG_HUGE y[], Integer n, Boolean NAG_HUGE *change);
#else
extern void g10caq();
#endif

#ifdef NAG_PROTO
extern void g10car(double NAG_HUGE y[], Integer n);
#else
extern void g10car();
#endif

#ifdef NAG_PROTO
extern void g10cas(double NAG_HUGE y[], Integer n);
#else
extern void g10cas();
#endif

#ifdef NAG_PROTO
extern void g10cat(double NAG_HUGE y[], Integer n);
#else
extern void g10cat();
#endif

#ifdef NAG_PROTO
extern void g10cau(double NAG_HUGE y[], Integer n, double NAG_HUGE work[]);
#else
extern void g10cau();
#endif

#ifdef NAG_PROTO
extern void g10cav(double NAG_HUGE y[], Integer n, double NAG_HUGE *endsav, double NAG_HUGE work[]);
#else
extern void g10cav();
#endif

#ifdef NAG_PROTO
extern void g10caw(double NAG_HUGE y[], Integer  n, Boolean NAG_HUGE *change);
#else
extern void g10caw();
#endif

#ifdef NAG_PROTO
extern void g10cax(double NAG_HUGE y[], Integer n, double endsav);
#else
extern void g10cax();
#endif

#ifdef NAG_PROTO
extern void g10cay(double NAG_HUGE y[], Integer n);
#else
extern void g10cay();
#endif

#ifdef NAG_PROTO
extern void g10caz(double NAG_HUGE y[], Integer n);
#else
extern void g10caz();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG10 */
