/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/threadmanager.h
// Purpose:     class for managing threads
// Created:     2005/10/6
// RCS-ID:      $Id: threadmanager.h,v 1.3 2006/01/30 19:19:59 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_MANAGER_H_
#define _ITO33_THREAD_MANAGER_H_

/**
    @file   ito33/thread/threadmanager.h
    @brief  class for managing threads.
 */

#include "ito33/vector.h"
#include "ito33/thread/poolthread.h"
#include "ito33/thread/taskqueue.h"

namespace ito33
{

namespace thread
{

class ThreadManager
{
public:

  /**
      Initialize the thread manager using up to the given number of threads.

      @param name the name of the threads created by this thread manager, for
                  debugging purposes
      @param maxThreads maximal number of threads we can use, default (0) means
                        determine the max number automatically to be the best 
                        one for the system (e.g. 2*number of CPUs)
   */
  ThreadManager(const char *name = NULL, unsigned maxThreads = 0);

  /**
      Shutdown the manager.

      This can block if there are any tasks still running.
   */
  ~ThreadManager();

  /**
      Queue another task for execution.

      @param task the function to execute, it must not throw any exceptions
   */
  void QueueTask(const Task& task);

private:

  std::vector<TaskPoolThread *> m_threads;

  unsigned m_maxThreads;

  TaskQueue m_taskQueue;

  NO_COPY_CLASS(ThreadManager);
};

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_THREADMANAGER_H_
