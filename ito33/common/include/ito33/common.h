/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/common.h
// Purpose:     common definitions for all ITO33 code
// Author:      Vadim Zeitlin
// Created:     19.12.02
// RCS-ID:      $Id: common.h,v 1.21 2005/05/26 19:29:47 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/common.h
    @brief  This file contains miscellaneous stuff which is used everywhere.

    This file should be always included, directly or indirectly, by all our
    source code.
 */

#ifndef _ITO33_COMMON_H_
#define _ITO33_COMMON_H_

// do this to disable some warnings (beforestd.h does it and afterstd.h is
// needed to restore the warning level -- but it doesn't reenable these
// particular warnings)
#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

// everyone should get ITO33_DLLDECL declaration
#include "ito33/dlldecl.h"

#ifdef _WIN32
  #define HAVE_WCHAR_T
#endif

// use file produced by configure if available
#ifdef HAVE_ITO33_CONFIG_H
    #include "ito33/config.h"
#endif // HAVE_ITO33_CONFIG_H

/**
    This namespace contains all our identifiers.

    Subnamespaces are used sometimes. Also, the macros we define are usually,
    but not always, prefixed with @c ITO33_ to avoid name clashes as the
    preprocessor knows nothing about namespaces.
 */
namespace ito33
{
} // namespace ito33

/**
    Macro returns the number of elements in the array.

    It can be used for any statically allocated array.
 */
#define SIZEOF(a)   (sizeof(a)/sizeof(a[0]))

/**
    Macro to be used to disable copy constructor and assignment operator.

    Often it is desirable to prevent the compiler from generating the default
    copy constructor and assignment operators because they wouldn't behave
    correctly for the objects of this class or simply because copying such
    objects doesn't make sense. The usual technique to do it is to declare (but
    not implement!) both of the methods as private -- this will catch most
    attempts to use them at compile-time and the remaining ones (from the class
    itself and its friends) during link-time.

    This works very well but is a bit verbose so we define a macro to make it
    simpler. Simply insert this macro somewhere in your class declaration
    (ideally at the end because it implicitly switches the scope to private and
    so the scope of fields and methods following it might be not what you
    expect).

    Final note: if you use this class, you should also have another ctor as the
    default ctor won't be generated automatically by the compiler neither.

    @param classname the name of the class this macro appears in
 */
#define NO_COPY_CLASS(classname)                                              \
private:                                                                      \
  classname(const classname&);                                                \
  classname& operator=(const classname&)

/**
    Same as NO_COPY_CLASS but can be used with template classes.

    Unfortunately this only works for template classes with a single template
    parameter so we also have to have NO_COPY_TEMPLATE_CLASS_2 and
    NO_COPY_TEMPLATE_CLASS_3 below.
 */
#define NO_COPY_TEMPLATE_CLASS(classname, arg)                                \
  private:                                                                    \
    classname(const classname<arg>&);                                         \
    classname<arg>& operator=(const classname<arg>&)

/// Same as NO_COPY_TEMPLATE_CLASS but for templates with 2 parameters.
#define NO_COPY_TEMPLATE_CLASS_2(classname, arg1, arg2)                       \
  private:                                                                    \
    classname(const classname<arg1, arg2>&);                                  \
    classname<arg1, arg2>& operator=(const classname<arg1, arg2>&)

/// Same as NO_COPY_TEMPLATE_CLASS but for templates with 3 parameters.
#define NO_COPY_TEMPLATE_CLASS_3(classname, arg1, arg2, arg3)                 \
  private:                                                                    \
    classname(const classname<arg1, arg2, arg3>&);                            \
    classname<arg1, arg2, arg3>& operator=(const classname<arg1, arg2, arg3>&)

/// Same as NO_COPY_TEMPLATE_CLASS but for templates with 4 parameters.
#define NO_COPY_TEMPLATE_CLASS_4(classname, arg1, arg2, arg3, arg4)           \
  private:                                                                    \
    classname(const classname<arg1, arg2, arg3, arg4>&);                      \
    classname<arg1, arg2, arg3, arg4>& operator=(                             \
                    const classname<arg1, arg2, arg3, arg4>&)

/**
    Macro which disables the assignment operator for a class.

    This is similar to but less restrictive then NO_COPY_CLASS() above: the
    copy ctor can still be explicitely defined (or left to be generated to the
    compiler) and only assignment is forbidden.

    This macro is useful with classes which contain reference or const fields
    as they can never be reassigned once initialized but copy ctor can usually
    be implemented without problems.
 */
#define NO_ASSIGN_CLASS(classname)                                            \
private:                                                                      \
  classname& operator=(const classname&)

/**
    Macro to suppress gcc warnings about class having private destructor and no
    friends.

    This warning is totally useless for the classes destroying themselves but
    there doesn't seem to be any way to disable just it so we have to work
    around it in the code.

    @param classname the name of the class this macro appears in
 */
#define NO_GCC_NO_FRIENDS_WARN(classname) friend class classname ## DummyFriend

/**
    @name printf() format checking macros

    ATTRIBUTE_PRINTF macros may be used in function declarations to enable
    checking of the printf() arguments.

    Note that only gcc supports this currently.
 */
//@{

#ifndef ATTRIBUTE_PRINTF

#ifdef __GNUC__
  #define ATTRIBUTE_PRINTF(m, n) __attribute__ ((__format__ (__printf__, m, n)))
#else
  #define ATTRIBUTE_PRINTF(m, n)
#endif /* ATTRIBUTE_PRINTF */

/// use this macro with a function whose first argument is the format string
#define ATTRIBUTE_PRINTF_1 ATTRIBUTE_PRINTF(1, 2)

/// use this macro with a function whose second argument is the format string
#define ATTRIBUTE_PRINTF_2 ATTRIBUTE_PRINTF(2, 3)

/// use this macro with a function whose third argument is the format string
#define ATTRIBUTE_PRINTF_3 ATTRIBUTE_PRINTF(3, 4)

/// use this macro with a function whose fourth argument is the format string
#define ATTRIBUTE_PRINTF_4 ATTRIBUTE_PRINTF(4, 5)

/// use this macro with a function whose fifth argument is the format string
#define ATTRIBUTE_PRINTF_5 ATTRIBUTE_PRINTF(5, 6)

#endif // !defined(ATTRIBUTE_PRINTF)

//@}

#if !defined(__GNUC__) && !(defined(_MSC_VER) && _MSC_VER >= 1300)
  /**
      Dummy __FUNCTION__ macro for the compilers which don't support the real
      thing.

      This allows us to show function names in the diagnostic messages for the
      compilers supporting it and fall back to __FILE__/__LINE__ for the older
      compilers.
   */
  #define __FUNCTION__ ""
#endif

#endif // _ITO33_COMMON_H_
