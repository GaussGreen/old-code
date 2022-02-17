///////////////////////////////////////////////////////////////////////////////
// Name:        atexit.h
// Purpose:     somewhat more flexible atexit() replacement
// Author:      Vadim Zeitlin
// Modified by:
// Created:     2004-10-22
// CVS-ID:      $Id: atexit.h,v 1.3 2006/05/16 01:16:54 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/atexit.h
   @brief RunAtExit class allows to run cleanup code automatically.

   This class allows to schedule things to run on program termination with
   minimum hassle.
 */

#ifndef _ITO33_ATEXIT_H_
#define _ITO33_ATEXIT_H_

#include "ito33/common.h"
#include "ito33/cpp.h"

namespace ito33
{

/**
   RunAtExit allows to run some cleanup code right before shutting down.

   This is the analog of standard atexit() but it will be called before
   static objects destruction and so you have more control on when exactly it
   is executed.

   To schedule something to be ran on program shutdown, create a static object
   of a class deriving from this one doing whatever you need to do in its Do().
 */
class RunAtExit
{
public:
   /**
      Constructor puts the object in the list of things to do on shutdown.
    */
   RunAtExit() : m_next(ms_first) { ms_first = this; }

   /**
      Virtual dtor just to suppress compiler warnings.

      The dtor might not be virtual because it must be trivial anyhow: it is
      ran too late to do anything!
    */
   virtual ~RunAtExit() { }

   /// Get the first object in the list
   static RunAtExit *GetFirst() { return ms_first; }

   /// Get the next one or NULL
   RunAtExit *GetNext() const { return m_next; }


   /**
      The operation to perform on shut down.

      This is called from MApplication dtor right before checking for memory
      leaks, all other subsystems can't be relied upon any longer!
    */
   virtual void Do() = 0;

private:
   // next object in the linked list
   RunAtExit * const m_next;

   // this is defined with ITO33_DEFINE_RUN_AT_EXIT
   static RunAtExit *ms_first;

   NO_COPY_CLASS(RunAtExit);
};

/**
   RunFunctionAtExit allows to run a simple function on exit.

   This just simplifies the life a little if the function to be executed on
   shutdown already exists.
 */
class RunFunctionAtExit : public RunAtExit
{
public:
   /**
      Constructor takes a pointer to the function to run.

      @param fn the function to execute on shutdown.
    */
   RunFunctionAtExit(void (*fn)()) : m_fn(fn) { }

   virtual void Do() { (*m_fn)(); }

private:
   void (*m_fn)();

   NO_COPY_CLASS(RunFunctionAtExit);
};

/**
    This macro can be used to schedule a function to be executed on program
    exit.

    For anything more complicated than a simple function call you will need to
    derive your own class from RunAtExit itself.

    @param func the function to call
 */
#define ITO33_AT_EXIT(func)                                                   \
    namespace                                                                 \
    {                                                                         \
      ::ito33::RunFunctionAtExit ITO33_MAKE_UNIQUE_NAME(ito33_AtExit_)(func); \
    }                                                                         \
    struct DummyStruct /* just to force semicolon after the macro */

/**
    Put this macro in some .cpp file to define static member of RunAtExit.

    This could be done in a separate file too, of course, but it is arguably
    simpler to use this macro in the same file where you use ITO33_RUN_AT_EXIT
    instead of remembering to add an extra almost empty file to makefile.

    Note that for all COM DLL projects there is no need to use it as it already
    appears in dllmain_impl.cpp.
 */
#define ITO33_DEFINE_RUN_AT_EXIT() \
  ::ito33::RunAtExit *::ito33::RunAtExit::ms_first = NULL

/**
    Use this macro to run the "at exit" code.

    Normally this should be placed as near to exit as possible but before
    main(), or equivalent, termination (otherwise there would be no advantage
    compared to using atexit() or static object dtors directly).
 */
#define ITO33_RUN_AT_EXIT()                                                   \
  do                                                                          \
  {                                                                           \
    using ::ito33::RunAtExit;                                                 \
                                                                              \
    /* this code is rather strange because we want it to be inside a do    */ \
    /* loop which executes exactly one time but writing standard while(0)  */ \
    /* at the end results in VC++ warnings so we work around it by putting */ \
    /* a "real" condition in the external loop and ensuring that it runs   */ \
    /* only once if there are no registered functions using the test below */ \
    RunAtExit *p = RunAtExit::GetFirst();                                     \
    if ( !p )                                                                 \
      break;                                                                  \
                                                                              \
    while ( p )                                                               \
    {                                                                         \
      p->Do();                                                                \
      p = p->GetNext();                                                       \
    }                                                                         \
  } while ( !::ito33::RunAtExit::GetFirst() )

} // namespace ito33

#endif // _ITO33_ATEXIT_H_

