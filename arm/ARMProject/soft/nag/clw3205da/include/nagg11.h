#ifndef NAGG11
#define NAGG11
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg11.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g11 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2160 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g11aac(Integer nrow, Integer ncol, Integer NAG_HUGE nobst[], Integer tdt,
             double NAG_HUGE expt[], double NAG_HUGE chist[], double NAG_HUGE *prob,
             double NAG_HUGE *chi, double NAG_HUGE *g, double NAG_HUGE *df,  NagError NAG_HUGE *fail);
#else
extern void g11aac();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG11 */
