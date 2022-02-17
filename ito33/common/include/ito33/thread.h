/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread.h
// Purpose:     classes and functions for multithreaded programs
// Author:      Vadim Zeitlin
// Created:     18.12.02
// RCS-ID:      $Id: thread.h,v 1.28 2006/02/15 17:04:38 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_H_
#define _ITO33_THREAD_H_

/**
    @file   ito33/thread.h
    @brief  Various helpers for the multithreaded programs.

    This file provides various synchronization primitives to be used in the
    multithreaded programs. Currently only low-level classes are implemented
    but we may want to add more complicated things (read/write locks for
    example) later.
 */

#include "ito33/thread/defs.h"

namespace ito33
{

// forward declarations:
namespace thread
{
  class Condition;
};


/**
    This class represents a runnable task.

    Objects of this (or derived) class may be used with Thread to create tasks
    running in separate threads.

    This is an abstract base class, to really use it you have to either
    manually derive a class from it and implement Run() method in it or,
    maybe simpler, use RunnableAdaptor<> class and MakeRunnable() to create a
    Runnable object from an existing class or function.
 */
class Runnable
{
public:
  /**
      The entry point.

      The task represented by this object starts here.
   */
  virtual void Run() = 0;

  /**
      Method called when the associated thread is about to terminate.

      The default implementation deletes this object.
   */
  virtual void OnThreadExit() { delete this; }

  /**
      Virtual dtor as for any base class.

      Dtor doesn't do anything.
   */
  virtual ~Runnable() { }

private:
  // it doesn't make sense to copy the tasks
  Runnable& operator=(const Runnable&);
};

/**
    A specialization of a runnable task which can be requested to stop.

    A thread is often sitting in a loop waiting for something to happen or for
    the main thread to stop it. This class is useful for modelling this kind of
    behaviour.
 */
class RunnableLoop : public Runnable
{
public:
  RunnableLoop()
  {
    m_stop = false;
  }

  // default copy ctor is ok

  /**
      Overridden Run() uses CheckCondition() and DoIteration() to implement the
      loop.
   */
  virtual void Run()
  {
    while ( !StopRequested() && CheckCondition() )
    {
      DoIteration();
    }
  }

  /**
      Request the loop termination.

      The thread will terminate as soon as the call of DoIteration() in
      progress (if any) finishes.
   */
  virtual void Stop() { m_stop = true; }

protected:
  // check whether/when the thread should continue (independently of Stop())
  virtual bool CheckCondition() const { return true; }

  // do a slice of whatever the thread should be doing
  virtual void DoIteration() = 0;

  // check if we were asked to stop
  virtual bool StopRequested() const { return m_stop; }

private:
  // true if the thread should exit
  bool m_stop;

  RunnableLoop& operator=(const RunnableLoop&);
};


/**
    Adaptor allowing using free functions as Runnable objects.

    Using RunnableAdaptor<> it is possible to represent an existing function or an
    object of a class implementing operator() as a Runnable task.

    @sa MakeRunnable()
 */
template <class T>
class RunnableAdaptor : public Runnable
{
public:
  /// constructs a runnable object from any object with operator()
  RunnableAdaptor(const T& obj) : m_obj(obj) { }

  virtual void Run() { m_obj(); }

private:
  T m_obj;

  // can't be copied after creation
  RunnableAdaptor operator=(const RunnableAdaptor&);
};

// we must use this typedef or VC++ 6 can't grok the template specialization
// and we must use the specialization for the functions because it can't
// compile the general template for them for some strange reason
/// this is just an implementation detail, not to be used directly
typedef void RunnableFunction();

/// this specialization is just an implementation detail
template <>
class RunnableAdaptor<RunnableFunction> : public Runnable
{
public:
  /// constructs a runnable object from a function
  RunnableAdaptor(void (*func)()) : m_func(func) { }

  virtual void Run() { (*m_func)(); }

private:
  void (*m_func)();

  // can't be copied after creation
  RunnableAdaptor operator=(const RunnableAdaptor&);
};

/**
    Creates a Runnable object from anything which can run.

    This is just a convenience wrapper around RunnableAdaptor<> ctor which
    allows the compiler to deduce the template parameter automatically.

    @param func the function or object to run
    @return a runnable object which can be used with Thread class (it is
            allocated using new and must be deleted by the caller)
 */
template <typename T>
inline RunnableAdaptor<T> *MakeRunnable(T func)
{
  return new RunnableAdaptor<T>(func);
}



/**
    Thread class allows to create a new OS thread executing the given function.

    This class represents a detached thread, i.e. a thread of a "fire and
    forget" kind -- once it is created, it runs on its own and will terminate
    when Runnable::Run() exits.

    As a consequence of this, the Thread object shouldn't be destroyed until
    the thread terminates and the simplest way to achieve it is by letting it
    to auto destruct itself: this is the policy which this class uses. Hence
    all object of Thread type must be allocated on the heap using new!

    If you need to wait for the thread termination, consider using
    JoinableThread instead.
 */
class Thread
{
public:
  /**
      Ctor takes a Runnable object and starts it in the new thread.

      Note that both ctor and dtor of Runnable will be called in the calling
      thread context while Runnable::Run() will be executed in the new
      thread.

      @param run pointer to the object to run, we take ownership of it and
             will delete it
      @param name of the thread for debugging purposes (@sa SetName)
   */
  Thread(Runnable *run, const char *name = NULL);

  /**
      Dtor: it doesn't stop the thread.

      The OS thread may continue running after Thread object is destroyed.

      Dtor is virtual just to suppress gcc warnings, it doesn't mean that
      this class is supposed to be used polymorphically or even as a base
      class because it is not.
   */
  virtual ~Thread();

  /**
      Returns the thread the calling code is executing in.

      The function will only return non-NULL pointer if the current thread
      was created by this class. For the main thread and any thread created
      without using this class, the return value will be NULL.

      @return the thread object for the current thread or @c NULL
   */
  static Thread *GetCurrent();

  /**
      Return the id of the current thread.

      An id is unique for all threads in the system but can be reused as
      threads are created and destroyed so it should be used only for the
      threads currently being alive.

      Notice that this function works for any thread, including the main one,
      unlike GetCurrent().
   */
  static THREAD_IDTYPE GetCurrentId();

  /**
      Returns the thread id.

      Thread id is just some number uniquely representing the thread. It is
      mostly useful for diagnostic messages and so on.
   */
  THREAD_IDTYPE GetId() const { return m_id; }

  /**
      Sets the name of the thread shown in the debugger under Win32.

      This method does nothing under other platforms.

      @param name a static char string
   */
  void SetName(const char *name);

  /**
      Sleeps for the specified amount of microseconds.

      @param usec time to sleep in microseconds
   */
  static void Sleep(unsigned long usec);

protected:
  /// we need another version of the ctor for the joinable threads
  Thread(Runnable *run, const char *name, bool isJoinable);

  /// this is an implementation helper
  void Init(bool isJoinable);

  /// close the thread handle
  void Close();

#ifdef _WIN32
  /// thread handle
  HANDLE m_hThread;
#endif // _WIN32

  /// thread id
  THREAD_IDTYPE m_id;

private:
  // the real start point for the new threads
  static THREAD_RETTYPE THREAD_CALLCONV ThreadStart(void *data);

  // the object to run
  Runnable *m_run;

  // true if this is a joinable thread, false otherwise
  bool m_isJoinable;


  // threads can't be copied
  Thread(const Thread&);
  Thread& operator=(const Thread&);
};


/**
    This is a thread class for the non-daemon threads.

    Daemon threads are represented by Thread class and run entirely on their
    own. Joinable threads are represented by this class and can be joined, that
    is another can wait until their termination and retrieve their exit code.
    In fact, another thread @b must join each joinable thread and this must be
    done exactly once.
 */
class JoinableThread : public Thread
{
public:
  /**
      Ctor starts the thread.

      @param run pointer to the object to run, we take ownership of it and
             will delete it
      @param name for diagnostic purposes only
   */
  JoinableThread(Runnable *run, const char *name = NULL);

  /**
      Dtor implicitly joins the thread.

      If Join() hadn't been called, it is going to be done from the dtor.
      Note that this means that the dtor may block for an arbitrarily long
      amount of time!
   */
  ~JoinableThread();

  /**
      Join, or wait for, this thread.

      This function must be called in the context of another thread, i.e. not
      the thread this object represents. It also must be called exactly once
      for each JoinableThread and the dtor is going to call it if it hadn't
      been done yet.
   */
  void Join();

private:
  JoinableThread(const JoinableThread&);
  JoinableThread& operator=(const JoinableThread&);
};


/**
    This class keeps any synchronization object locked during its life time,
    ie it locks it in its ctor and unlocks in dtor.

    This is a template class which can be specialized for any object providing
    Lock() and Unlock() methods, for example CriticalSection or Mutex.
 */
template <class Lockable>
class Lock
{
public:
  /// ctor locks the resource
  Lock(Lockable& lock);

  /// dtor unlocks the resource
  ~Lock();

private:
  // Lock can't be copied
  Lock(const Lock&);
  Lock& operator=(const Lock&);

  Lockable& m_lock;
};



/**
    Critical section is the simplest synchronization primitive: it can be
    entered, or locked, by only one thread at a time.

    Note that CriticalSection satisfies Lockable requirments and so may -- and
    should! -- be used with Lock class. Calling Lock() and Unlock() directly is
    bad idea as it risks to leave the critical section in inconsistent state if
    an exception is thrown or even if the function simply returns
    "unexpectedly".

    CriticalSection never throws exceptions because Win32 critical section
    functions never fail at all and it would be unexpected to have methods
    throwing exceptions only under Unix. Also, even with pthreads the mutex
    functions usually don't fail and using throw would have severe performance
    implications.
 */
class CriticalSection
{
public:
  /**
      Ctor initializes the critical section.

      Note that it doesn't enter it.
   */
  CriticalSection();

  /**
      Dtor frees the resources associated with the critical section.

      Note that it doesn't leave the critical section and that destroying a
      locked section is an error because other threads waiting to enter it
      may be never waken up.
   */
  ~CriticalSection();

  /**
      Enters the critical section blocking until this becomes possible if
      needed.

      Reentering the critical section when it is already locked by the
      current thread is an error and leads to undefined behaviour.

      Do not use this method directly, use Lock class instead!
   */
  void Lock();

  /**
      Leaves the critical section.

      The critical section must have been previously entered by Lock().

      Do not use this method directly, use Lock class instead!
   */
  void Unlock();

  /**
      Check whether the calling thread is inside the critical section.

      This is only meant as a debugging helper and shouldn't be used for
      anything outside asserts.

      Currently only implemented for Win32, always returns -1 under Unix.

      @return 1 if the calling thread is inside the critical section,
              0 if it isn't and -1 if the function is not implemented on this
              platform.
   */
  int IsLockedByThisThread() const;

private:
#ifdef _WIN32
  DWORD GetOwningThread() const;

  // for efficiency, CriticalSection is implemented inline and doesn't
  // allocate memory dynamically under Win32 -- unfortunately this means that
  // we have to include <windows.h> from here to be able to use
  // CRITICAL_SECTION here
  CRITICAL_SECTION m_cs;
#else // Posix threads
  pthread_mutex_t m_mutex;
#endif // Win32/pthreads

  // critical sections can't be copied
  CriticalSection(const CriticalSection&);
  CriticalSection& operator=(const CriticalSection&);

  friend class thread::Condition;
};


/**
    This class provides a simple and fast way to execute some code only once.

    Often some code should only be executed once during the program life time.
    Normally this is easily achieved by putting it into a static object ctor,
    for example, but this is not MT-safe. To do it in a safe way, a global
    object of this class should be used and the "once only" code should be
    executed only if its IsFirstTime() method returns true.

    Please note that the object *must* be global, otherwise there would be a
    race condition during its initialization and more than one thread could
    sneak in the protected section!
 */
class RunOnce
{
public:
  RunOnce();

  /**
      Returns true exactly once, when it is called for the first time.

      All subsequent calls to this method (including "simultaneous" ones from
      the other threads) will return false.
   */
  bool IsFirstTime();

private:
#ifdef _WIN32
  // Win32 implementation uses InterlockedExchange() function which is much
  // faster than using critical sections

  // NB: this variable must be 4 bytes aligned but with default compiler
  //     options it should happen automatically
  LONG m_wasCalled;
#else // Posix threads
  // generic implementation using a critical section to protect the flag
  CriticalSection m_cs;
  bool m_wasCalled;
#endif // Win32/pthreads
};

// ============================================================================
// inline functions implementation
// ============================================================================

// ----------------------------------------------------------------------------
// CriticalSection
// ----------------------------------------------------------------------------

inline
CriticalSection::CriticalSection()
{
#ifdef _WIN32
  ::InitializeCriticalSection(&m_cs);
#else
  pthread_mutexattr_t attr;
  if ( pthread_mutexattr_init(&attr) != 0 )
  {
    FAIL( "pthread_mutexattr_init() failed" );
  }
  if ( pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE) != 0 )
  {
    FAIL( "pthread_mutexattr_settype() failed" );
  }
  if ( pthread_mutex_init(&m_mutex, &attr) != 0 )
  {
    FAIL( "pthread_mutex_init() failed" );
  }
  if ( pthread_mutexattr_destroy(&attr) != 0 )
  {
    FAIL( "pthread_mutexattr_destroy() failed" );
  }
#endif
}

inline
CriticalSection::~CriticalSection()
{
#ifdef _WIN32
  ::DeleteCriticalSection(&m_cs);
#else
  if ( pthread_mutex_destroy(&m_mutex) != 0 )
  {
    FAIL( "pthread_mutex_init() failed" );
  }
#endif
}

#ifdef _WIN32
inline
DWORD CriticalSection::GetOwningThread() const
{
  // NB: OwningThread is an id, not HANDLE, in spite of its type, and to add
  //     insult to injury VC++ complains about the cast
  #ifdef _MSC_VER
    #pragma warning(disable:4311)
  #endif

  return reinterpret_cast<DWORD>(m_cs.OwningThread);

  #ifdef _MSC_VER
    #pragma warning(default:4311)
  #endif
}
#endif // _WIN32

inline
int CriticalSection::IsLockedByThisThread() const
{
#ifdef _WIN32
  // thread ids could be reused but OwningThread handle wouldn't refer to a
  // thread which doesn't exist any longer unless it forgot to unlock m_cs in
  // which case we have worse problems than thread id reuse
  return GetOwningThread() == ::GetCurrentThreadId();
#else
  return -1;
#endif
}

inline
void CriticalSection::Lock()
{
#ifdef _WIN32
  // critical sections are recursive under Windows and so we can reenter the
  // same section again without even noticing it, but we shouldn't rely on
  // this behaviour, i.e. if it happens, it's a logical error in the code
  ASSERT_MSG( !IsLockedByThisThread(), "attempt to reenter a critical section" );

  ::EnterCriticalSection(&m_cs);
#else
  if ( pthread_mutex_lock(&m_mutex) != 0 )
  {
    FAIL( "pthread_mutex_lock() failed" );
  }
#endif
}

inline
void CriticalSection::Unlock()
{
  // only the thread which had locked us can unlock
  ASSERT_MSG( IsLockedByThisThread(), "unlocking not locked critical section" );

#ifdef _WIN32
  ::LeaveCriticalSection(&m_cs);
#else
  if ( pthread_mutex_unlock(&m_mutex) != 0 )
  {
    FAIL( "pthread_mutex_unlock() failed" );
  }
#endif
}

// ----------------------------------------------------------------------------
// Lock<>
// ----------------------------------------------------------------------------

template <class Lockable>
inline
Lock<Lockable>::Lock(Lockable& lock)
       : m_lock(lock)
{
  m_lock.Lock();
}

template <class Lockable>
inline
Lock<Lockable>::~Lock()
{
  m_lock.Unlock();
}

// ----------------------------------------------------------------------------
// RunOnce
// ----------------------------------------------------------------------------

inline
RunOnce::RunOnce()
{
  m_wasCalled = false;
}

inline
bool RunOnce::IsFirstTime()
{
#ifdef _WIN32
  // InterlockedExchange() returns the old value of m_wasCalled
  return !::InterlockedExchange(&m_wasCalled, true);
#else
  Lock<CriticalSection> lock(m_cs);

  if ( m_wasCalled )
    return false;

  m_wasCalled = true;

  return true;
#endif
}

// ----------------------------------------------------------------------------
// Thread
// ----------------------------------------------------------------------------

#if !defined(_WIN32) || !defined(_MSC_VER)
inline void Thread::SetName(const char * /* name */)
{
}
#endif // !VC++ under Win32

inline
Thread::Thread(Runnable *run, const char *name)
      : m_run(run)
{
  Init(false /* not joinable */);
  SetName(name);
}

inline
Thread::Thread(Runnable *run, const char *name, bool isJoinable)
      : m_run(run)
{
  Init(isJoinable);
  SetName(name);
}

/* static */ inline
void Thread::Sleep(unsigned long usec)
{
#ifdef _WIN32
  ::Sleep(usec / 1000);
#elif defined(HAVE_USLEEP)
  usleep(usec);
#else
  #error "Don't know how to implement Thread::Sleep() (-DHAVE_USLEEP?)"
#endif
}

/* static */ inline
THREAD_IDTYPE Thread::GetCurrentId()
{
#ifdef _WIN32
  return ::GetCurrentThreadId();
#else
  return pthread_self();
#endif
}

// ----------------------------------------------------------------------------
// JoinableThread
// ----------------------------------------------------------------------------

inline
JoinableThread::JoinableThread(Runnable *run, const char *name)
       : Thread(run, name, true /* joinable */)
{
}

} // namespace ito33

#endif // _ITO33_THREAD_H_

