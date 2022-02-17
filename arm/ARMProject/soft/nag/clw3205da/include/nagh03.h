#ifndef NAGH03
#define NAGH03
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagh03.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library h03 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2163 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL h03abc(double NAG_HUGE cost[], Integer tdcost, double NAG_HUGE avail[],
              Integer navail, double NAG_HUGE req[], Integer nreq, Integer maxit,
              Integer NAG_HUGE *numit, double NAG_HUGE optq[], Integer NAG_HUGE source[],
              Integer NAG_HUGE dest[], double NAG_HUGE *optcost, double NAG_HUGE unitcost[],
              NagError NAG_HUGE *fail);
#else
extern void h03abc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGH03 */
