/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/threadmanager.cpp
// Purpose:     implementation of ThreadManager
// Created:     2005/10/7
// RCS-ID:      $Id: threadmanager.cpp,v 1.3 2006/01/30 19:19:59 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/thread.h"
#include "ito33/thread/threadmanager.h"

using namespace ito33;
using namespace ito33::thread;

ThreadManager::ThreadManager(const char *name, unsigned maxThreads)
             : m_maxThreads(maxThreads)
{
  // create 2 threads by default (TODO: use number of CPUs)
  if ( m_maxThreads == 0 )
    m_maxThreads = 2;

  Atomic::IntType sem(m_maxThreads);
  for ( unsigned i = 0; i < m_maxThreads; i++ )
  {
    m_threads.push_back(new TaskPoolThread(m_taskQueue, &sem, name));
  }

  // busy loop until the threads get really started and are ready to work
  //
  // we should be using a semaphore here but we don't have it yet and this busy
  // loop is normally very short and ThreadManager objects are created rarely
  // (i.e. once in program lifetime) so it shouldn't be a big deal
  for ( ;; )
  {
    Atomic::Inc(&sem);
    if ( !Atomic::Dec(&sem) )
      break;
  }
}

void ThreadManager::QueueTask(const Task& task)
{
  m_taskQueue.Add(task);
}

ThreadManager::~ThreadManager()
{
  m_taskQueue.StopWorkerThreads();
}
