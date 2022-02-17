#ifndef NAG_STDDEF
#define NAG_STDDEF  

/* <nag_stddef.h>
 * 
 * Copyright 1990 Numerical Algorithms Group Ltd.
 *
 * Include file for NAG C Library 
 *
 * Mark 1, 1990
 *
 * Mark 2, revised.
 * Revised handling of implementation specifics.  AJS jan92.
 *
 * Purpose : Ensures that typedefs for size_t and ptrdiff_t are included.
 */

#ifdef NAG_YES_STDDEF
#include <stddef.h>
#else

#ifndef SIZE_T
#define SIZE_T
typedef unsigned size_t;
#endif /* not SIZE_T */
#ifndef PTRDIFF_T
#define PTRDIFF_T
typedef int		ptrdiff_t;	/* result of subtracting two pointers */
#endif  /* not PTRDIFF_T */

#endif

#endif /* not NAG_STDDEF */
