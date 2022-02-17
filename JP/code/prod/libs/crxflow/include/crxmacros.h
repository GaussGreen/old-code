/*
***************************************************************************
** HEADER FILE: crxmacros.h
**
** Various useful macros. These are for use within the implementation files
** and not within the CX namespace.
**
** $Header$
***************************************************************************
*/

#ifndef CRX_MACROS_H
#define CRX_MACROS_H

/* Various macros */
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <alib/cgeneral.h>
#include <alib/cerror.h>
#include <alib/macros.h>
#include <alib/strutil.h>

#define ERROR_HANDLER GtoErrMsg
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

#undef STRDUP
#define STRDUP GtoStringDuplicate

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

