/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/condition.h
// Purpose:     declares a Condition variable class
// Author:      Vadim Zeitlin
// Created:     2005-10-07
// RCS-ID:      $Id: condition.h,v 1.2 2005/12/15 14:05:45 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/thread/condition.h
    @brief Declaration of Condition class.

    Condition is used to block a thread until a certain condition becomes true,
    e.g. a resource becomes available (without polling for it).
 */

#ifndef _ITO33_THREAD_CONDITION_H_
#define _ITO33_THREAD_CONDITION_H_

#include "ito33/thread/mutex.h"

#ifdef _WIN32
  #include "ito33/thread/atomic.h"
  #include "ito33/win32/event.h"
#endif

namespace ito33
{

namespace thread
{

/**
    This is a portable equivalent of a POSIX condition variable.

    As pthread_cond_t, a Condition object is always used with a Mutex.
    The mutex must be locked before calling Wait() and will be atomically
    unlocked inside it and relocked before Wait() returns.

    Notice that spurious wake ups are possible so the thread which returns from
    Wait() @b must check the condition it is waiting on. As it should also
    check this condition before starting to wait (in case it's not necessary),
    the typical example of using a Condition is this:
      @code
        struct SharedData
        {
          Mutex mutex;
          Condition cond;
          int value;        // or whatever
        } g_data;
        ...
        // thread waiting for g_data.value to reach the threshold
        void WaitForValue()
        {
          Lock<Mutex> lock(g_data.mutex);

          // this loop is important!
          while ( g_data.value != 17 )
          {
            g_data.cond.Wait(g_data.mutex);
          }

          ... do whatever we have to do with data ...
          ... (while mutex is still locked) ...
        }
      @endcode
 */
class Condition
{
public:
  /**
      Constructor initializes a new condition variable.

      May throw if creation fails.
   */
  Condition();

  /**
      Destroys condition variable.

      No threads should be waiting on this condition!
   */
  ~Condition();

  /**
      Waits until the condition variable is signalled by another thread.

      This method has two preconditions:
        - the mutex passed to it @b must be locked
        - the condition for which we wait @b must be currently false

      When Wait() returns, the mutex is locked (again) and will have to be
      unlocked by the caller later.

      May throw if an error occurs.

      @todo Add timeout parameter here (and change return value to bool to
            indicate whether we timed out) or add separate WaitTimeout()?

      @param mutex the associated mutex, must be locked before the call
   */
  void Wait(Mutex& mutex);

  /**
      Signals the condition waking up one of the threads waiting for it.

      The mutex used to guard the shared data structure may be locked or not
      when Signal() is called. Typically it was locked recently to do the
      modification of this data structure which resulted in condition becoming
      true and it may be still locked; however it has to be unlocked before
      Wait() returns in the other thread (as it tries to reacquire the mutex as
      well).

      @todo Add Broadcast() which wakes up all the threads when we need it
   */
  void Signal();

private:
#ifdef _WIN32
  // number of waiting threads
  Atomic::IntType m_numWaiters;

  // a manual event which is set in Signal() if there are any waiters and reset
  // in Wait()
  Win32::Event m_evtSignal;
#else // pthreads
  pthread_cond_t m_cond;
#endif
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

#ifdef _WIN32

/*
    In general implementation of POSIX conditions in Win32 is *highly* non
    trivial, see

        http://www.cs.wustl.edu/~schmidt/win32-cv-1.html

    and

        http://sources.redhat.com/cgi-bin/cvsweb.cgi/pthreads/pthread_cond_wait.c?rev=1.10&cvsroot=pthreads-win32

    But we simplify our life here by not implementing Broadcast() at all.
 */

inline Condition::Condition()
                : m_numWaiters(0)
{
}

inline Condition::~Condition()
{
  ASSERT_MSG( !m_numWaiters, "waiters remain on condition being destroyed" );
}

inline void Condition::Wait(Mutex& mutex)
{
  Atomic::Inc(&m_numWaiters);

  // condition must currently be false as otherwise we wouldn't be here at all
  m_evtSignal.Reset();

  // it's ok to release the mutex: as m_numWaiters > 0 now, if Signal() is
  // called from now on, it will do SetEvent() and some thread will be woken up
  mutex.Unlock();

  m_evtSignal.Wait();

  Atomic::Dec(&m_numWaiters);

  // relock the mutex before returning
  mutex.Lock();
}

inline void Condition::Signal()
{
  // reading 32 bit value is atomic on all Windows platforms so no need to
  // protect it
  if ( m_numWaiters )
  {
    // at most one waiting thread will be woken up
    m_evtSignal.Set();
  }
}

#else // !_WIN32

inline Condition::Condition()
{
  if ( pthread_cond_init(&m_cond, NULL) != 0 )
  {
    FAIL( "pthread_cond_init() failed" );
  }
}

inline Condition::~Condition()
{
  if ( pthread_cond_destroy(&m_cond) != 0 )
  {
    FAIL( "pthread_cond_destroy() failed" );
  }
}

inline void Condition::Wait(Mutex& mutex)
{
  if ( pthread_cond_wait(&m_cond, &mutex.m_mutex) != 0 )
  {
    FAIL( "pthread_cond_wait() failed" );
  }
}

inline void Condition::Signal()
{
  if ( pthread_cond_signal(&m_cond) != 0 )
  {
    FAIL( "pthread_cond_signal() failed (threads waiting on condition?)" );
  }
}

#endif // _WIN32/!_WIN32

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_CONDITION_H_

