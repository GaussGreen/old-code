/////////////////////////////////////////////////////////////////////////////
// Name:        src/unix/thread.cpp
// Purpose:     implementation of thread-related classes using POSIX API
// Author:      Vadim Zeitlin
// Created:     18.02.03
// RCS-ID:      $Id: thread.cpp,v 1.12 2005/04/02 14:03:08 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/string.h"

#include "ito33/thread.h"
#include "ito33/errnoexception.h"

#ifndef _NDEBUG
  #include <stdio.h>              // debug messages go to stderr
#endif

#include <errno.h>

// ============================================================================
// private classes
// ============================================================================

// ----------------------------------------------------------------------------
// PthreadAttr: wraps pthread_attr_t
// ----------------------------------------------------------------------------

class PthreadAttr
{
public:
  PthreadAttr()
  {
    const int rc = pthread_attr_init(&m_attr);
    if ( rc != 0 )
    {
      throw ERRNO_EXCEPTION_ERR("pthread_attr_init", rc);
    }
  }

  ~PthreadAttr()
  {
    const int rc = pthread_attr_destroy(&m_attr);
    if ( rc != 0 )
    {
#ifndef _NDEBUG
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      std::string msg =
        ERRNO_EXCEPTION_ERR("pthread_attr_destroy", rc).GetFullMessage();
      fputs(msg.c_str(), stderr);
#endif // _NDEBUG
    }
  }

  void SetScope(int scope)
  {
    ASSERT_MSG( scope == PTHREAD_SCOPE_SYSTEM ||
          scope == PTHREAD_SCOPE_PROCESS,
          "invalid scope argument to pthread_attr_setscope()" );

    const int rc = pthread_attr_setscope(&m_attr, scope);
    if ( rc != 0 )
    {
      throw ERRNO_EXCEPTION_ERR("pthread_attr_setscope", rc);
    }
  }

  void SetDetachState(int detach)
  {
    ASSERT_MSG( detach == PTHREAD_CREATE_DETACHED ||
          detach == PTHREAD_CREATE_JOINABLE,
          "invalid argument for pthread_attr_setdetachstate()" );

    const int rc = pthread_attr_setdetachstate(&m_attr, detach);
    if ( rc != 0 )
    {
      throw ERRNO_EXCEPTION_ERR("pthread_attr_setdetachstate", rc);
    }
  }

  // semi-implicit conversion to the real thing
  pthread_attr_t *operator()() { return &m_attr; }

private:
  pthread_attr_t m_attr;
};

// ----------------------------------------------------------------------------
// ThreadSpecificDataUntyped represents a Pthread key, see pthread_key_xxx()
//
// note that ThreadSpecificDataUntyped is not used directly, only
// ThreadSpecificData defined below is
// ----------------------------------------------------------------------------

class ThreadSpecificDataUntyped
{
public:
  ThreadSpecificDataUntyped()
  {
    const int rc = pthread_key_create(&m_key, NULL /* no dtor function */);
    if ( rc != 0 )
    {
      throw ERRNO_EXCEPTION_ERR("pthread_key_create", rc);
    }
  }

  ~ThreadSpecificDataUntyped()
  {
    const int rc = pthread_key_delete(m_key);
    if ( rc != 0 )
    {
#ifndef _NDEBUG
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      std::string msg =
        ERRNO_EXCEPTION_ERR("pthread_key_delete", rc).GetFullMessage();
      fputs(msg.c_str(), stderr);
#endif // _NDEBUG
    }
  }

  // these functions are only to be called by ThreadSpecificDataUntyped<>
protected:
  void Set(void *value)
  {
    const int rc = pthread_setspecific(m_key, value);
    if ( rc != 0)
    {
      throw ERRNO_EXCEPTION_ERR("pthread_setspecific", rc);
    }
  }

  void *Get() const
  {
    // there is no error return from pthread_getspecific(): errno is not set
    // byt it, we can only test whether the return value is NULL but it does
    // happen to be NULL for the main thread so we really can't do much here
    return pthread_getspecific(m_key);
  }

private:
  // the thread specific data key
  pthread_key_t m_key;
};

// a type safe wrapper around ThreadSpecificDataUntyped
template <typename T>
class ThreadSpecificData : public ThreadSpecificDataUntyped
{
public:
  // chained assignments don't make sense with TLS slots so operator=()
  // exceptionally returns nothing in this class
  void operator=(T value) { Set(reinterpret_cast<void *>(value)); }
  T operator()() const { return reinterpret_cast<T>(Get()); }
};

// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

// TSD key where we store the pointer to the current thread
//
// we use a function returning a reference to the static object to create it on
// first use only and not on the program startup
static inline ThreadSpecificData<ito33::Thread *>& KeyThisThread()
{
  static ThreadSpecificData<ito33::Thread *> s_keyThisThread;

  return s_keyThisThread;
}

// ============================================================================
// Thread implementation
// ============================================================================

namespace ito33
{

// ----------------------------------------------------------------------------
// Thread
// ----------------------------------------------------------------------------

/* static */
THREAD_RETTYPE THREAD_CALLCONV
Thread::ThreadStart(void *data)
{
  Thread * const thread = reinterpret_cast<Thread *>(data);

  // remember the thread pointer in the TLS to be able to retrieve it later
  KeyThisThread() = thread;

  // do run the thread now
  thread->m_run->Run();

  return 0;
}

/* static */
Thread *Thread::GetCurrent()
{
  return KeyThisThread()();
}

void Thread::Init(bool isJoinable)
{
  m_isJoinable = isJoinable;

  PthreadAttr attr;

  // we want to have system-wide threads everywhere and it's not always the
  // default state (it is under Linux but it isn't under Solaris)
  attr.SetScope(PTHREAD_SCOPE_SYSTEM);

  attr.SetDetachState(isJoinable ? PTHREAD_CREATE_JOINABLE
                 : PTHREAD_CREATE_DETACHED);

  int rc = pthread_create
      (
        &m_id,              // [out] thread id
        attr(),             // attributes
        ThreadStart,        // start routine
        this                // its parameter
      );

  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_create", rc);
  }
}

Thread::~Thread()
{
  // notify the thread that it has terminated
  m_run->OnThreadExit();
}

// ----------------------------------------------------------------------------
// JoinableThread
// ----------------------------------------------------------------------------

JoinableThread::~JoinableThread()
{
  if ( GetId() )
  {
    try
    {
      Join();
    }
#ifndef _NDEBUG
    catch ( ErrnoException& e )
    {
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      fputs(e.GetFullMessage().c_str(), stderr);
    }
#endif // _NDEBUG
    catch ( ... )
    {
      // do nothing -- dtor must not throw
    }
  }
}

void JoinableThread::Join()
{
  CHECK_VOID( GetId(), "Join() had been already called or invalid thread" );

  void *exitcode;

  const int rc = pthread_join(GetId(), &exitcode);
  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_join", rc);
  }

  // don't try to join it any more
  m_id = 0;
}

} // namespace ito33

