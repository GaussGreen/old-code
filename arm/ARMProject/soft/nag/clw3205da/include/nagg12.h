#ifndef NAGG12
#define NAGG12
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagg12.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library g12 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2161 (Feb 1998).
 */

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL g12aac(Integer n,  double NAG_HUGE t[],  Integer NAG_HUGE ic[],  Integer NAG_HUGE ifreq[],
	    Integer NAG_HUGE *nd,  double NAG_HUGE tp[], double NAG_HUGE p[], double NAG_HUGE psig[],
	    NagError NAG_HUGE *fail);
#else
extern void g12aac();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGG12 */
