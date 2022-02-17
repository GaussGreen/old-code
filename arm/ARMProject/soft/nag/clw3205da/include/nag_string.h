#ifndef NAG_STRING
#define NAG_STRING  

/* <nag_string.h>
 *
 * Copyright 1989 Numerical Algorithms Group
 *
 * NAG interface to ANSI string.h
 *
 * Malcolm Cohen, 1989, NAG Ltd., Oxford, U.K.
 *
 * Mark 1, 1990
 * Modified handling of implementation specifics Feb 1992, ajs
 */

#ifdef NAG_YES_STRING
#include <string.h>
#else
#include <strings.h>
#endif

#ifndef NAG_NO_INC_MEMORY_DOT_H
#include <memory.h>
#endif

#ifndef NAG_NO_BC_MEM_FUN
#define memcpy(s,ct,n) bcopy(ct,s,(int)(n))
#define memcmp bcmp
#ifdef NAG_PROTO
extern int bcopy (char *ct, char *s, int n);
extern int bcmp (char *b1, char *b2, int length);
#else
extern int bcopy(),bcmp();
#endif /* NAG_PROTO */
#endif

#ifndef NAG_IND_STR
#define strchr index
#define strrchr rindex
#endif

#endif  /* not NAG_STRING */
