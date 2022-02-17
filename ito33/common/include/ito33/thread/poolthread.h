/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/poolthread.h
// Purpose:     Thread class monitoring a queue for task
// Created:     2005/10/7
// RCS-ID:      $Id: poolthread.h,v 1.4 2006/01/30 19:19:59 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_POOLTHREAD_H_
#define _ITO33_THREAD_POOLTHREAD_H_

/**
    @file   ito33/thread/poolthread.h
    @brief  Thread class monitoring a queue for task.
 */

#include "ito33/thread/atomic.h"
#include "ito33/thread/condition.h"
#include "ito33/thread/taskqueue.h"

#ifdef _MSC_VER
  // 'this' : used in base member initializer list
  #pragma warning(disable:4355)
#endif // _MSC_VER

namespace ito33
{

namespace thread
{

/**
    Thread class used to monitor a task queue, and works whenever there is
    request.

    This is a semi-private class used only by ThreadManager.

    Template parameter T is the type of objects stored in the queue.

    @see TaskPoolThread
 */
template<typename T>
class PoolThread : public Thread
{
public:
  /// Queue type for T
  typedef Queue<T> QueueType;

  /**
      Ctor takes a taskQueue object and starts monitoring it.

      @param taskQueue The task queue this thread is monitoring
      @param sem pointer to semaphore-like counter decremented by the thread
                 when it's ready to do some work; may be NULL if the user
                 doesn't use the counter
      @param name name of the thread for debugging purposes
   */
  PoolThread(QueueType& taskQueue, Atomic::IntType *sem,
             const char *name = NULL)
    : Thread(new Runnable(*this, taskQueue, sem), name)
  {
    taskQueue.RegisterWorker();
  }

  /**
      Trivial dtor.
   */
  ~PoolThread() { }

protected:
  /**
      Derived class must override this method to process a task fetched
      from the queue.

      @param task object fetched from the queue
   */
  virtual void HandleTask(const T& task) = 0;

private:
  // Main loop of this thread
  class Runnable : public ito33::Runnable
  {
  public:
    Runnable(PoolThread& thread, QueueType& taskQueue, Atomic::IntType *sem)
      : m_thread(thread), m_taskQueue(taskQueue), m_sem(sem)
    {
    }

  protected:
    virtual void Run()
    {
      for ( ; ; )
      {
        T task;

        {
          Lock<Mutex> lock( m_taskQueue.GetMutex() );

          // signal semaphor once to indicate that we've locked the mutex and
          // so are not going to lose any "work available" notifications now
          if ( m_sem )
          {
            Atomic::Dec(m_sem);
            m_sem = NULL;
          }

          if ( !m_taskQueue.WaitForWork() )
          {
            // WaitForWork returned because the queue was requested to stop
            m_taskQueue.UnregisterWorker();

            return;
          }

          task = m_taskQueue.Get();
        }
        // unlock the mutex before handling the queued item

        m_thread.HandleTask(task);
      }
    }

  private:
    PoolThread& m_thread;
    QueueType& m_taskQueue;
    Atomic::IntType *m_sem;

    NO_COPY_CLASS(Runnable);
  };

  friend class Runnable;

  NO_COPY_CLASS(PoolThread);
};

/**
    Specialization of PoolThread for running tasks from TaskQueue.
 */
class TaskPoolThread : public PoolThread<Task>
{
public:
  /**
      Ctor takes a taskQueue object and starts monitoring it.

      @param taskQueue The task queue this thread is monitoring
      @param sem pointer to semaphore-like counter decremented by the thread
                 when it's ready to do some work; may be NULL if the user
                 doesn't use the counter
      @param name name of the thread for debugging purposes
   */
  TaskPoolThread(TaskQueue& taskQueue, Atomic::IntType *sem,
                 const char *name = NULL)
    : PoolThread<Task>(taskQueue, sem, name)
  {
  }

protected:
  virtual void HandleTask(const Task& task)
  {
    task();
  }
};

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_POOLTHREAD_H_
