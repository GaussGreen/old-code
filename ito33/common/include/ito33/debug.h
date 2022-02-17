/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/debug.h
// Purpose:     various debugging helpers
// Author:      Vadim Zeitlin
// Created:     16.12.02
// RCS-ID:      $Id: debug.h,v 1.28 2006/05/30 12:38:54 wang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/debug.h
    @brief Debugging helpers: compile- and run-time asserts and more.

    This header defines the debug macros as well as MemoryUsageSnapshot class
    which may be useful for debugging memory leak/consumption problems.
 */

#ifndef _ITO33_DEBUG_H_
#define _ITO33_DEBUG_H_

#include "ito33/common.h"

// there is nothing exportable in this file and old versions of cpp2any die on 
// function pointer parameter to Assert()
#ifndef __CPP2ANY__

#ifndef NDEBUG
  #include "ito33/beforestd.h"
  #include <string>
  #include "ito33/afterstd.h"
#endif // NDEBUG

#include <cstdlib>             // for NULL
#include <climits>             // for CHAR_BIT used below

namespace ito33
{

// ----------------------------------------------------------------------------
// debugging assertions (these ones disappear in release build)
// ----------------------------------------------------------------------------

#ifndef NDEBUG

/// Type of the assert function
typedef void (*OnAssertFunc)(const char *, int, const char *, const char *);

/**
    This function may be redefined to do something non trivial and is called
    whenever one of debugging macros fails (ie condition is false in an
    assertion)

    @param filename file name of the ASSERT
    @param line line number of the ASSERT
    @param condition the condition which failed in text form
    @param message optional message explaining the reason
*/
extern void OnAssert(const char *filename,
                     int line,
                     const char *condition,
                     const char *message = NULL);

/**
    The global assert function pointer, can be changed for special asserts
    handling.

    Set to OnAssert() by default.
 */
extern OnAssertFunc AssertFunction;

/// Declaration of ITO333_ONASSERT_FUNC macros
#define ITO33_ONASSERT_FUNC ::ito33::AssertFunction

/**
    Call this function to break into the debugger unconditionally (assuming
    the program is running under debugger, of course)
 */
extern void Trap();

/**
    Helper function used to implement ASSERT and ASSERT_MSG

    @param cond the condition to check
    @param filename the name of the file which contains the assert
    @param line the line number in this file
    @param condition the condition being checked in text form
    @param message the message to show when assert fails
    @param onAssert function to be called whenever one of debugging macros
                        fails (ie condition is false in the assertion)
 */
inline void
Assert(int cond,
       const char *filename,
       int line,
       const char *condition,
       const char *message,
       OnAssertFunc onAssert)
{
  if ( !cond )
      (*onAssert)(filename, line, condition, message);
}

/**
    Overload of Assert() allowing to use std::string for the assert message.

    This is especially useful for messages constructed with String::Printf().
 */
inline void
Assert(int cond,
       const char *filename,
       int line,
       const char *condition,
       const std::string& message,
       OnAssertFunc onAssert)
{
  Assert(cond, filename, line, condition, message.c_str(), onAssert);
}

// note using "int" and not "bool" for cond to avoid VC++ warnings about
// implicit conversions when doing "Assert( pointer )" and also use of
// "!!cond" below to ensure that everything is converted to int

// macro name collision when ASSERT is defined for a mfc project
#ifndef ITO33MFC
/// generic assert macro
#define ASSERT(cond) \
  ::ito33::Assert(!!(cond), __FILE__, __LINE__, #cond, NULL, ITO33_ONASSERT_FUNC)
#endif

/// assert
#define ITO33ASSERT(cond) \
  ::ito33::Assert(!!(cond), __FILE__, __LINE__, #cond, NULL, ITO33_ONASSERT_FUNC)

/// assert with additional message explaining it's cause
#define ASSERT_MSG(cond, msg) \
  ::ito33::Assert(!!(cond), __FILE__, __LINE__, #cond, msg, ITO33_ONASSERT_FUNC)

/// helper macro to avoid warnings for things we useo nly in debug build
#define UNUSED_IN_RELEASE(arg)  arg

#else // NDEBUG

#define ITO33_ONASSERT_FUNC NULL

#ifndef ITO33MFC
#define ASSERT(cond)
#endif
#define ITO33ASSERT(cond)

#define ASSERT_MSG(cond, msg)

#define UNUSED_IN_RELEASE(arg)

#endif // !NDEBUG

/// special form of assert: always triggers it (in debug mode)
#define FAIL(msg) ASSERT_MSG(false, msg)

// ----------------------------------------------------------------------------
// the debugging macros which stay in release build as well
// ----------------------------------------------------------------------------

#ifndef NDEBUG

/**
    Helper for CHECK macros, similar to Assert() for ASSERT one.

    @param cond the condition to check
    @param filename the name of the file which contains the assert
    @param line the line number in this file
    @param condition the condition being checked in text form
    @param message the message to show when assert fails
    @param onAssert function to be called whenever one of debugging macros
                        fails (ie condition is false in the assertion)
 */
inline bool
Check(int cond,
      const char *filename,
      int line,
      const char *condition,
      const char *message,
      OnAssertFunc onAssert)
{
  return cond ? true : ((*onAssert)(filename, line, condition, message), false);
}

#else // NDEBUG

// Check() still exists in release build unlike Assert() but it simply checks
// the condition and doesn't do anything else then
inline bool
Check(int cond, const char *, int, const char *, const char *,
      void (*)(const char *, int, const char *, const char *))
{
    return !!cond;
}

#endif // !NDEBUG/NDEBUG

/**
    Check that the condition is true and return with the specified value if it
    isn't (also generating a FAIL(msg) then)
 */
#define CHECK(cond, rc, msg)    \
  if ( !::ito33::Check(!!(cond), __FILE__, __LINE__, \
                       #cond, msg, ITO33_ONASSERT_FUNC) ) \
    return rc

/**
    Same as CHECK but for the void functions.
 */
#define CHECK_VOID(cond, msg)    \
  if ( !::ito33::Check(!!(cond), __FILE__, __LINE__, \
                       #cond, msg, ITO33_ONASSERT_FUNC) ) \
    return

// ----------------------------------------------------------------------------
// compile time assertions
// ----------------------------------------------------------------------------

/*
  How this works (you don't have to understand it to be able to use the
  macros): we rely on the fact that it is invalid to define a named bit field
  in a struct of width 0. All the rest are just the hacks to minimize the
  possibility of the compiler warnings when compiling this macro: in
  particular, this is why we define a struct and not an object (which would
  result in a warning about unused variable) and a named struct (otherwise we'd
  get a warning about an unnamed struct not used to define an object!).
 */

/// MAKE_UNIQUE_ASSERT_NAME helper
#define MAKE_ASSERT_NAME(msg, line)         Assert_ ## msg ## line

/// COMPILE_TIME_ASSERT helper
#define MAKE_UNIQUE_ASSERT_NAME(msg)        MAKE_ASSERT_NAME(msg, __LINE__)

/**
  This macro generates a compile time error if the condition given to it is
  false.

  The second argument of this macro must be a valid C++ identifier and not a
  string. I.e. you should use it like this:

    COMPILE_TIME_ASSERT( sizeof(int) >= 2, YourIntsAreTooSmall );

  It may be used both within a function and in the global scope.
*/
#define COMPILE_TIME_ASSERT(expr, msg) \
    struct MAKE_UNIQUE_ASSERT_NAME(msg) { unsigned int msg: expr; }

/// helper for ASSERT_BITSIZE below
#define MAKE_BITSIZE_MSG(type, what, size) \
    type ## what ## size ## Bits

/**
    A special case of compile time assert: check that the size of the given
    type is at least the given number of bits
 */
#define ASSERT_MIN_BITSIZE(type, size) \
    COMPILE_TIME_ASSERT(sizeof(type) * CHAR_BIT >= size, \
                          MAKE_BITSIZE_MSG(type, SmallerThan, size))

/**
    Another special case of compile time assert: check that the size of the
    given type is exactly equal to the given number of bits
 */
#define ASSERT_BITSIZE(type, size) \
    COMPILE_TIME_ASSERT(sizeof(type) * CHAR_BIT == size, \
                          MAKE_BITSIZE_MSG(type, NotEqualTo, size))

// ----------------------------------------------------------------------------
// memory usage monitor class
// ----------------------------------------------------------------------------

/**
    This class allows to monitor memory consumption easily.

    Creating an object of this class takes the snapshot of the current memory
    usage and remembers the total number of objects allocated, total amount of
    memory consumed and also the memory consumption peak so far.

    Note that it is only available in debug build and only when using MS VC++
    for now. Otherwise all of the accessor functions simply always return 0.
 */
class MemoryUsageSnapshot
{
public:
    /// constructor takes memory snapshot
    MemoryUsageSnapshot();

    /// @name Accessors
    //@{

    /// get the number of blocks currently allocated
    unsigned long GetAllocCount() const { return m_blocksCur; }

    /// get memory, in bytes, currently in use
    unsigned long GetUsage() const { return m_bytesCur; }

    /// get max amount ever allocated so far
    unsigned long GetMaxUsage() const { return m_bytexMax; }

    //@}

private:
    // currently allocated blocks
    unsigned long m_blocksCur;

    // current memory consumption
    unsigned long m_bytesCur;

    // max memory consumption so far
    unsigned long m_bytexMax;
};


// ----------------------------------------------------------------------------
// debug messages
// ----------------------------------------------------------------------------

#ifndef NDEBUG
/**
    Prints the message to the Debug Output (win32).
    
    Uses OutputDebugString procedure to print the message, which will be
    visible in the Debug window of MSVC or in a standalone program like
    DbgView from sysinternals.com
    Has the same syntax as printf().
*/
extern void DebugPrintf(const char *format, ...);

/**
    Prints the message to the Debug Output (win32) if in Debug mode.

    Calls DebugPrintf() procedure if compiled in the Debug mode, disappears
    in the Release build.
    Because the number of arguments is not fixed, the additional pair of
    parenthesis should be supplied. Example:
    DEBUGPRINT(("Value is %d",value));
*/
#define DEBUGPRINT( exp ) DebugPrintf exp

#else  // NDEBUG

#define DEBUGPRINT( exp )

#endif  // ifndef NDEBUG


}  // namespace ito33

#endif // __CPP2ANY__

#endif // _ITO33_DEBUG_H_

