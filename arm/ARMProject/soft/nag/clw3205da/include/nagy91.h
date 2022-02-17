#ifndef NAGY91
#define NAGY91
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagy91.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library y91 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2172 (Feb 1998).
 */

#include <stdio.h>



/* Structure template cmrand used in e04nfcs and y91nfz */
struct cmrand
{
    Integer one, two, three;
    Integer ix, iy, iz;
};

#ifdef NAG_PROTO
extern void y91nfz(const char *outfile, Integer NAG_HUGE *lprob, Integer NAG_HUGE *inform,
					  Integer NAG_HUGE *msglvl, Integer lda, Integer ldh, Integer n, Integer nclin,
					  Integer nfixed, Integer nactiv, Integer negevh, Integer nulevh,
					  struct cmrand NAG_HUGE *seed, double NAG_HUGE *condhy, double NAG_HUGE *condhz,
					  double NAG_HUGE *condhx, double NAG_HUGE *conda1, double NAG_HUGE *conda2,
					  double NAG_HUGE *conda3, double NAG_HUGE *conda4, double NAG_HUGE *bndlow,
					  double NAG_HUGE *bndupp, double NAG_HUGE *spread, double NAG_HUGE a[],
					  double NAG_HUGE ax[], double NAG_HUGE bl[], double NAG_HUGE bu[],
					  double NAG_HUGE c[], double NAG_HUGE lambda[], double NAG_HUGE h[], double NAG_HUGE x[],
					  double NAG_HUGE w[], Integer lw);
#else
extern void y91nfz();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGY91 */
