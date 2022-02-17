/*---------------------------------------------------------------------
HEADER FILE:    macros.h

CREATED BY:     2/25/93 D. Gallager

PURPOSE:        Miscellaneous macros

$Header: /home/alibcvs/alibcvs/.cvsroot-alib/platdep/include/macros.h,v 1.85 2004/06/28 19:32:33 ptaylor Exp $
----------------------------------------------------------------------
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management.

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */
#ifndef IRX_MACROS_H
#define IRX_MACROS_H

#include "error.h"
#include "memutils.h"                   /* FREE, MALLOC */
#include <string.h>                     /* memcpy */
#include <float.h>                      /* DBL_EPSILON, etc */


#ifdef __cplusplus
extern "C"
{
#endif

#ifndef BIT_SET
#define BIT_SET(flags,bit) ((flags) |= (bit))
#endif

#ifndef BIT_CLR
#define BIT_CLR(flags,bit) ((flags) &= ~(bit))
#endif

#ifndef IS_BIT_SET
#define IS_BIT_SET(flags,bit) ((flags) & (bit))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) ((x) > 0 ? (x) : -(x))
#endif

#ifndef IFA
#define IFA(a,b,c) ((a) ? (b) : (c))
#endif

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

#ifndef CUBE
#define CUBE(x) ((x)*(x)*(x))
#endif

#ifndef IRX_IS_ODD
#define IRX_IS_ODD(x) ( (int) ( ( (unsigned) x ) & 1 ) )
#endif

#ifndef IRX_IS_EVEN
#define IRX_IS_EVEN(x) ( !IRX_IS_ODD(x) )
#endif

/* This should be correct on every platform that supports the IEEE
 * floating point standard, i.e., everything except DEC VAXes.
 *
 * Note that it is only defined here if it was _not_ defined in float.h
 *
 */
#ifndef DBL_EPSILON
#define DBL_EPSILON     2.2204460492503131E-16
#endif

#define FCLOSE(fp) do { if (fp) fclose(fp); fp = NULL; } while(0)


/* IS_ALMOST_ZERO returns TRUE if and only if 1.0 + x = 1.0, assuming x is a
 * double. In other words, it compares x with 1.0
 */
#define IS_ALMOST_ZERO(x) ((x) < DBL_EPSILON && (x) > -DBL_EPSILON)
#define ARE_ALMOST_EQUAL(x,y) IS_ALMOST_ZERO((x)-(y))
#define IS_BETWEEN(x,a,b) ((a) < (b) ? \
                           ((x) >= (a) && (x) <= (b)) : \
                           ((x) >= (b) && (x) <= (a)))

#define IS_ZERO_TO_TOL(x,tol) ((x) < (tol) && (x) > -1.*(tol) )
#define ARE_EQUAL_TO_TOL(x,y,tol) IS_ZERO_TO_TOL((x)-(y),(tol))

#ifndef PROGRAM_BUG
#define PROGRAM_BUG() irxError("Program bug:%s line %d\n",__FILE__,__LINE__)
#endif

#ifndef NEW
#define NEW            IRX_NEW
#define NEW_ARRAY      IRX_NEW_ARRAY
#define NEW_ARRAY_2D   IRX_NEW_ARRAY_2D
#define FREE           IRX_FREE
#define FREE_ARRAY     IRX_FREE
#define FREE_PTR_ARRAY IRX_FREE_PTR_ARRAY

#endif

#ifndef COPY_ARRAY
#define COPY_ARRAY(dst,src,t,n) memcpy((char*)(dst),(char*)(src),(n)*sizeof(t))
#endif

/* STRDUP does a string duplication where memory is allocated using NEW_ARRAY*/
#undef STRDUP
#define STRDUP(s) ((s)==NULL ? NULL : strcpy(NEW_ARRAY(char,strlen(s)+1),(s)))

/* REQUIRE should be used to check input conditions. */
/* REQUIRE assumes existence of 'RETURN' label for error cases. */
/* REQUIRE assumes that the function is named 'routine'. */
#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    irxError("%s: Required condition (%s) fails!\n",routine,#cond);\
    goto RETURN;\
}} while(0)

/* ASSERT should be used to check things should be true from prior code      */
/* ASSERT should not be used to check input conditions - use REQUIRE instead */
/* ASSERT assumes existence of 'RETURN' label for error cases.                 */
/* ASSERT assumes that the function is named 'routine'.                      */
/* ASSERT has same syntax as the system macro assert, but failures are not   */
/*        quite so drastic as assert().                                      */
/* ASSERT is checked in release code, unlike assert() which is only checked  */
/*        in debug code.                                                     */
#undef ASSERT
#define ASSERT(cond) do { if (!(cond))\
{\
    irxError("%s: Assertion (%s) fails: %s line %d\n",routine,#cond,__FILE__,__LINE__);\
    goto RETURN;\
}} while(0)

/*
** With gcc we use -Wfloat-equal.
**
** In some cases though we really want to allow exact float comparisons.
** In these cases use these macros to avoid the warning messages.
*/

#undef IS_NOT_EQUAL
#undef IS_EQUAL
#undef IS_NOT_ZERO
#undef IS_ZERO

#ifdef __GNUC__

#define IS_NOT_EQUAL(x,y) ((x)<(y) || (x)>(y))
#define IS_EQUAL(x,y) !IS_NOT_EQUAL((x),(y))
#define IS_NOT_ZERO(x) IS_NOT_EQUAL((x),0.0)
#define IS_ZERO(x) IS_EQUAL((x),0.0)

#else

#define IS_NOT_EQUAL(x,y) ((x)!=(y))
#define IS_EQUAL(x,y) ((x)==(y))
#define IS_NOT_ZERO(x) ((x)!=0.0)
#define IS_ZERO(x) (x)==0.0

#endif


#ifdef __cplusplus
}
#endif

#endif    /* IRX_MACROS_H */
