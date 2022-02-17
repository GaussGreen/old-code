/*
***************************************************************************
** HEADER FILE: cxmacros.h
**
** Various useful macros. These are for use within the implementation files
** and not within the CX namespace.
**
** $Header$
***************************************************************************
*/

#ifndef CX_MACROS_H
#define CX_MACROS_H

/* Various macros */
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <alib/cerror.h>
#include <alib/cgeneral.h>
#include <alib/macros.h>

#include "memutils.h"

/* NEW and NEW_ARRAY allocate memory for data types and initialises to 0 */
#undef NEW
#undef NEW_ARRAY
#undef NEW_ARRAY_2D
#define NEW CX_NEW
#define NEW_ARRAY CX_NEW_ARRAY
#define NEW_ARRAY_2D CX_NEW_ARRAY_2D

/* FREE and FREE_ARRAY safely frees memory as allocated by NEW and NEW_ARRAY */
#undef FREE
#undef FREE_ARRAY
#define FREE CX_FREE
#define FREE_ARRAY CX_FREE_ARRAY
#undef FREE_PTR_ARRAY
#define FREE_PTR_ARRAY CX_FREE_PTR_ARRAY

/* REQUIRE should be used to check input conditions. */
/* REQUIRE assumes existence of 'done' label for error cases. */
/* REQUIRE assumes that the function is named 'routine'. */
#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    GtoErrMsg("%s: Required condition (%s) fails!\n",routine,#cond);\
    goto done;\
}} while(0)

/* ASSERT should be used to check things should be true from prior code      */
/* ASSERT should not be used to check input conditions - use REQUIRE instead */
/* ASSERT assumes existence of 'done' label for error cases.                 */
/* ASSERT assumes that the function is named 'routine'.                      */
/* ASSERT has same syntax as the system macro assert, but failures are not   */
/*        quite so drastic as assert().                                      */
/* ASSERT is checked in release code, unlike assert() which is only checked  */
/*        in debug code.                                                     */
#undef ASSERT
#define ASSERT(cond) do { if (!(cond))\
{\
    GtoErrMsg("%s: Assertion (%s) fails: %s line %d\n",routine,#cond,__FILE__,__LINE__);\
    goto done;\
}} while(0)

#undef PROGRAM_BUG
#define PROGRAM_BUG() GtoErrMsg("Program bug:%s line %d\n",__FILE__,__LINE__)

/* COPY_ARRAY copies an array when type conversion is not required. */
/* COPY_ARRAY uses memcpy. */
#undef COPY_ARRAY
#define COPY_ARRAY(dst,src,T,n) memcpy((char*)(dst),\
                                       (char*)(src),\
                                       (size_t)(n)*sizeof(T))

/* CLEAR_ARRAY sets each element in an array to zero */
#undef CLEAR_ARRAY
#define CLEAR_ARRAY(array,T,n) memset((char*)array,(size_t)0,(size_t)(n)*sizeof(T))

/* STRDUP does a string duplication where memory is allocated using NEW_ARRAY*/
#undef STRDUP
#define STRDUP(s) ((s)==NULL ? NULL : strcpy(NEW_ARRAY(char,strlen(s)+1),(s)))

#undef MAX
#undef MIN
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#undef ABS
#define ABS(a) ((a)>0 ? (a) : -(a))

#ifndef DBL_EPSILON
#define DBL_EPSILON     2.2204460492503131E-16
#endif

#undef IS_ALMOST_ZERO
#undef ARE_ALMOST_EQUAL
#define IS_ALMOST_ZERO(x) ((x) < DBL_EPSILON && (x) > -DBL_EPSILON)
#define ARE_ALMOST_EQUAL(x,y) IS_ALMOST_ZERO((x)-(y))

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

#endif

