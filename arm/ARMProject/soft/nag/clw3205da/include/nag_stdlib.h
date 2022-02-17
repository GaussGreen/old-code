#ifndef NAG_STDLIB
#define NAG_STDLIB  

/* <nag_stdlib.h>
 *
 * Copyright 1989 Numerical Algorithms Group
 *
 * NAG routine interface to ANSI stdlib routines.
 *
 * Malcolm Cohen, 1989, NAG Ltd, Oxford, U.K.
 *
 * Mark 1, 1990.
 * Mark 1a revised. IER-864 (Oct 1990).
 * Mark 5 revised. IER-2137 (Feb 1998).
 */


#ifdef NAG_YES_STDLIB
#include <stdlib.h>
#include <nag_stddef.h>
#include <nagx04.h>
#else

#include <nag_stddef.h>
#ifdef NAG_PROTO
extern void abort(void);
extern void exit(int status);
extern int atoi(const char *s);
extern long atol(const char *s);
extern void *calloc(size_t nobj, size_t size);
extern void *malloc(size_t size);
extern void far *farmalloc(unsigned long size);
extern void *realloc(void *p, size_t size);
extern void far *farrealloc(void far *p, unsigned long size);
extern void free(void *p);
extern void farfree(void far *p);
extern char *alloca(int size);
#else
extern int abort(),atoi(),free(),farfree(),exit();
extern long atol();
extern char *calloc(),*malloc(),*realloc(),*alloca();
#endif

#endif

#ifdef NAG_YES_MALLOC
/*#define NAG_ALLOC(n,type) (type *)malloc((size_t)(n)*sizeof(type))
#define NAG_REALLOC(pointer,n,type) (type *)realloc((Pointer)(pointer), \
(size_t)(n)*sizeof(type))
#define NAG_FREE(x) free((Pointer)x)
*/

#define NAG_ALLOC(n,type) (type *)x04bbc((size_t)(n)*sizeof(type))
#define NAG_REALLOC(pointer,n,type) (type *)x04bcc((Pointer)(pointer), \
(size_t)(n)*sizeof(type))
#define NAG_FREE(x) x04bdc((Pointer *)&(x))


#else
#ifdef IBPA_NAG_IMP
#include <malloc.h>
#define NAG_ALLOC(n,type) (type NAG_HUGE *)halloc((long)(n), sizeof(type))
#define NAG_REALLOC(pointer,n,type) (type NAG_HUGE *)realloc((Pointer)(pointer), \
(size_t)(n)*sizeof(type))
#define NAG_FREE(x) hfree((Pointer)x)

#else
#ifdef IBPE_NAG_IMP
#include <alloc.h>
#include <dos.h>
#define NAG_CALLOC(n,type) (type huge *)farcalloc((long)(n), sizeof(type))
#define NAG_ALLOC(n,type) (type far *)farmalloc((long)(n)*sizeof(type))
#define NAG_HALLOC(n,type) (type huge *)farmalloc((long)(n)*sizeof(type))
#define NAG_REALLOC(pointer,n,type) (type far *)farrealloc((Pointer)(pointer),\
(long)(n)*sizeof(type))
#define NAG_FREE(x) farfree((Pointer)x)
#endif
#endif
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif


#endif  /* not NAG_STDLIB */
