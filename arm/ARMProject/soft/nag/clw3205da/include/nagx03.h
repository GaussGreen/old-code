#ifndef NAGX03
#define NAGX03
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagx03.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library x03 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2168 (Feb 1998).
 */

#ifdef NAG_PROTO
extern void x03aac(double NAG_HUGE a[], Integer sizea, double NAG_HUGE b[], Integer sizeb,
	    Integer n, Integer stepa, Integer stepb, double c1,
	    double c2, double NAG_HUGE *d1, double NAG_HUGE *d2, Boolean sw, NagError NAG_HUGE *fail);
#else
extern void x03aac();
#endif

#ifdef NAG_PROTO
extern void x03aay(double NAG_HUGE a[], double NAG_HUGE b[], Integer n, Integer stepa, Integer stepb,
            double c1, double c2, double NAG_HUGE *d1, double NAG_HUGE *d2);
#else
extern void x03aay();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGX03 */
