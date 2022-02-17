/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/queue.h
// Purpose:     Class scheduling works for worker threads
// Created:     2005/10/7
// RCS-ID:      $Id: queue.h,v 1.7 2006/02/21 22:07:02 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_QUEUE_H_
#define _ITO33_THREAD_QUEUE_H_

/**
    @file   ito33/thread/queue.h
    @brief  Class scheduling works for worker threads
 */

#include "ito33/beforestd.h"
#include <queue>
#include "ito33/afterstd.h"

#include "ito33/thread.h"
#include "ito33/thread/condition.h"

namespace ito33
{

namespace thread
{

/**
    Class scheduling works for worker threads.

    @todo Split this class in a non-template QueueBase and template Queue
          deriving from it to avoid template bloat.
 */
template<typename T>
class Queue
{
public:

  /// Ctor creates an empty queue.
  Queue() : m_isStopRequested(false) { }

  /// Dtor clears the queue, elements remaining in the queue are not handled.
  ~Queue() { }

  /**
      Get the mutex that needs to be locked when WaitForWork() and Get() are
      called.
   */
  Mutex& GetMutex() { return m_cs; }

  /**
      Add an element to this queue, also signal worker threads.

      @param element the element to be handled
   */
  void Add(const T& element)
  {
    Lock<Mutex> lock(m_cs);

    m_elements.push(element);

    m_condNotEmpty.Signal();
  }

  /**
      Register a worker thread on this queue.

      After calling RegisterWorker() you must call UnregisterWorker() some
      time later.
   */
  void RegisterWorker()
  {
    Lock<Mutex> lock(m_waitForStop.mutex);
    m_waitForStop.numRunning++;
  }

  /**
      Called by worker thread to wait until the queue becomes non empty.

      The mutex must be locked when this method is called. It is atomically
      unlocked here for the duration of the wait and relocked when the method
      returns.

      @return true if there is work to do
   */
  bool WaitForWork()
  {
    while ( m_elements.empty() && !m_isStopRequested )
    {
      m_condNotEmpty.Wait(m_cs);
    }

    return !m_isStopRequested;
  }

  /**
      Unregister a worker thread registered previously by calling
      RegisterWorker().
   */
  void UnregisterWorker()
  {
    Lock<Mutex> lock(m_waitForStop.mutex);
    m_waitForStop.numRunning--;
  }

  /**
      Stop all threads working on this queue.

      @todo use Condition::Broadcast() when available
   */
  void StopWorkerThreads()
  {
    // ensure that when the thread wakes up the next time from its wait on
    // m_condNotEmpty it will terminate
    {
      Lock<Mutex> lock(m_cs);
      m_isStopRequested = true;
    }

    // as we don't have Broadcast() yet, we busy wait for the threads to exit
    // and we have to wait for each one in turn as if we sent more than one
    // signal at once they could all wake up the same thread
    for ( ;; )
    {
      unsigned numRunningOld;
      {
        Lock<thread::Mutex> lock(m_waitForStop.mutex);
        numRunningOld = m_waitForStop.numRunning;
        if ( !numRunningOld )
          break;
      }

      m_condNotEmpty.Signal();

      for ( ;; )
      {
        Thread::Sleep(10);
        Lock<thread::Mutex> lock(m_waitForStop.mutex);
        if ( m_waitForStop.numRunning < numRunningOld )
          break;
      }
    }
  }

  /**
      Check if the queue is asked to stop scheduling.

      The caller must lock the mutex prior to calling this function.

      @return true if stop requested, false otherwise
   */
  bool StopRequested() const
  {
    ASSERT_MSG( m_cs.IsLockedByThisThread(), "must acquire lock before calling" );

    return m_isStopRequested;
  }

  /**
      Get the next scheduled element from this queue.

      Must be called with mutex locked by caller.

      @return the next scheduled element
   */
  T Get()
  {
    T element = m_elements.front();

    m_elements.pop();

    return element;
  }

private:

  struct WaitForStop
  {
    WaitForStop() : numRunning(0) { }

    Mutex mutex;
    unsigned numRunning;
  } m_waitForStop;

  std::queue<T> m_elements;

  mutable CriticalSection m_cs;
  Condition m_condNotEmpty;

  bool m_isStopRequested;

  NO_COPY_TEMPLATE_CLASS(Queue, T);
};

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_QUEUE_H_
