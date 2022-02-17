/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/rwlock.h
// Purpose:     read-write lock implementation
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: rwlock.h,v 1.6 2006/02/18 23:44:33 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_RWLOCK_H_
#define _ITO33_RWLOCK_H_

/**
    @file   ito33/thread/rwlock.h
    @brief  Read-write lock and related helper classes.

    In multithreaded programs it is commonly ok for many threads to read some
    shared data structure but, of course, we can't be modifying and reading it
    at the same time -- this is where read-write lock comes in handy. It allows
    for an arbitrary number of readers as long as there are no writer threads
    but prevents the other threads from both reading and writing to the shared
    data while some thread is modifying it.
 */

#ifdef _WIN32
  #include "ito33/thread.h"
  #include "ito33/win32/event.h"

  using ito33::Win32::Event;
#else // !_WIN32
  #include "ito33/errnoexception.h"
#endif // _WIN32/!_WIN32

namespace ito33
{

/**
    RWLock is a read-write lock class.

    The class implemented here favours the writes over readers. The read locks
    are recursive, i.e. a thread which already has a read lock may acquire it
    again (and so may the other threads) but not upgradable, i.e. a thread
    holding a read lock must release it before trying to acquire a write lock.
    Write locks are, in general, not recursive even if on some platforms this
    could work -- but this behaviour is not portable.

    Under POSIX systems we use the native pthread_rwlock_xxx() functions to
    implement this class but under Win32 we have to roll our own.

    Consider using LockRead and LockWrite classes below instead of calling
    ReadLock/WriteLock() and ReadUnlock/WriteUnlock() directly.
 */
class RWLock
{
public:
  /**
      Constructor initializes a read write lock.

      It may throw an exception if any low level synchronization objects we
      need couldn't have been created.
   */
  RWLock();

  /**
      Acquire read lock.

      As long as there are no (active or pending) writers, this function
      succeeds immediately. Otherwise it blocks.

      It may throw an exception if something goes badly wrong at OS level.
   */
  void ReadLock();

  /**
      Acquire write lock.

      This function blocks until there are no more other readers nor writers.

      It may throw an exception if something goes badly wrong at OS level.
   */
  void WriteLock();

  /**
      Release read lock.

      This function never blocks. If the caller doesn't hold a read lock the
      behaviour is undefined but likely fatal.

      It may throw an exception if something goes badly wrong at OS level.
   */
  void ReadUnlock();

  /**
      Release write lock.

      This function never blocks. If the caller doesn't hold the write lock
      the behaviour is undefined but likely fatal.

      It may throw an exception if something goes badly wrong at OS level.
   */
  void WriteUnlock();

  /**
      Destructor is not virtual.

      This class is not meant to be used polymorphically.
   */
  ~RWLock();

private:
#ifdef _WIN32
  // this critical section is used for serializing the writers (hence the
  // name) but also for blocking any new readers from starting to read once a
  // writer is blocking waiting for the existing readers to finish
  CriticalSection m_csWrite;

  // this critical section is used to protect access to m_nReaders and to
  // ensure that m_eventCanWrite.Set() and Reset() calls are executed in the
  // right order (see comments in ReadLock/ReadUnlock)
  CriticalSection m_csRead;

  // this event is signaled when it's ok to start writing, i.e. when there
  // are no more readers (m_csWrite is used to ensure that there are never
  // more than one writer)
  Event m_eventCanWrite;

  // the number of readers, protected by m_csRead
  size_t m_nReaders;
#else // !_WIN32
  pthread_rwlock_t m_rwlock;
#endif

  // read-write locks can't be copied
  RWLock(const RWLock&);
  RWLock& operator=(const RWLock&);
};


/**
    LockRead automatically locks and unlocks a RWLock for reading.

    Calling RWLock::ReadLock() and ReadUnlock() manually is error prone as
    forgetting to call one of them or doing it one too many times is usually
    fatal. With this class you can simply create it on the stack and it will do
    everything for you. It will also ensure that RWLock is not left in a locked
    state if an exception is thrown.

    @sa Lock, LockWrite
 */
class LockRead
{
public:
  /// Constructor acquires the read lock
  LockRead(RWLock& rwlock) : m_rwlock(rwlock) { rwlock.ReadLock(); }

  /// Destructor releases the read lock
  ~LockRead() { m_rwlock.ReadUnlock(); }

private:
  RWLock& m_rwlock;

  // lockers can't be copied
  LockRead(const LockRead&);
  LockRead& operator=(const LockRead&);
};


/**
    LockWrite automatically locks and unlocks a RWLock for writing.

    Calling RWLock::WriteLock() and WriteUnlock() manually is error prone as
    forgetting to call one of them or doing it one too many times is usually
    fatal. With this class you can simply create it on the stack and it will do
    everything for you. It will also ensure that RWLock is not left in a locked
    state if an exception is thrown.

    @sa Lock, LockRead
 */
class LockWrite
{
public:
  /// Constructor acquires the write lock
  LockWrite(RWLock& rwlock) : m_rwlock(rwlock) { rwlock.WriteLock(); }

  /// Destructor releases the write lock
  ~LockWrite() { m_rwlock.WriteUnlock(); }

private:
  RWLock& m_rwlock;

  // lockers can't be copied
  LockWrite(const LockWrite&);
  LockWrite& operator=(const LockWrite&);
};

// ----------------------------------------------------------------------------
// RWLock inline methods implementation
// ----------------------------------------------------------------------------

#ifdef _WIN32

inline
RWLock::RWLock()
   : m_eventCanWrite(Event::Manual, Event::Signaled)
{
  m_nReaders = 0;
}

inline
void RWLock::ReadLock()
{
  // block if there are already any writers waiting, otherwise we'd have
  // a "writers starvation" situation
  Lock<CriticalSection> lockW(m_csWrite);

  // ensure that we don't call Reset() below before ReadUnlock() calls Set():
  // m_eventCanWrite state must correspond to whether the number of readers
  // is 0 or not and without this critical section we could break this
  // invariant when ReadUnlock() is called for the last reader but before it
  // has time to call m_eventCanWrite.Set() another reader arrives --
  // although it would do Reset() here, a later Set() in ReadUnlock() would
  // negate it
  Lock<CriticalSection> lockR(m_csRead);

  // incidentally, m_nReaders is now also protected by m_csRead so we can
  // modify it safely
  if ( !m_nReaders++ )
  {
    // we are the first reader, reset the event to indicate that we can't
    // write any more
    m_eventCanWrite.Reset();
  }
  //else: m_eventCanWrite should be already in non signaled state
}

inline
void RWLock::WriteLock()
{
  // entering this critical section serializes all writers and also prevents
  // new readers from continuing in ReadLock()
  m_csWrite.Lock();

  // now wait until the existing readers finish
  m_eventCanWrite.Wait();
}

inline
void RWLock::ReadUnlock()
{
  // see the long comment in ReadLock()
  Lock<CriticalSection> lockR(m_csRead);

  if ( !--m_nReaders )
  {
    // we were the last reader in existence, as there are no more readers
    // it becomes possible to acquire write lock again
    m_eventCanWrite.Set();
  }
}

inline
void RWLock::WriteUnlock()
{
  // let the other readers (and/or writers) blocking in Read/WriteLock()
  // continue
  m_csWrite.Unlock();
}

inline
RWLock::~RWLock()
{
  ASSERT_MSG( !m_nReaders, "RWLock: not all readers called ReadUnlock!" );
}

#else // !_WIN32

inline
RWLock::RWLock()
{
  // initialize with default attributes
  const int rc = pthread_rwlock_init(&m_rwlock, NULL);
  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_rwlock_init", rc);
  }
}

inline
void RWLock::ReadLock()
{
  const int rc = pthread_rwlock_rdlock(&m_rwlock);
  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_rwlock_rdlock", rc);
  }
}

inline
void RWLock::WriteLock()
{
  const int rc = pthread_rwlock_wrlock(&m_rwlock);
  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_rwlock_wrlock", rc);
  }
}

inline
void RWLock::ReadUnlock()
{
  const int rc = pthread_rwlock_unlock(&m_rwlock);
  if ( rc != 0 )
  {
    throw ERRNO_EXCEPTION_ERR("pthread_rwlock_unlock", rc);
  }
}

inline
void RWLock::WriteUnlock()
{
  ReadUnlock();
}

inline
RWLock::~RWLock()
{
  const int rc = pthread_rwlock_destroy(&m_rwlock);
  if ( rc != 0 )
  {
#ifndef _NDEBUG
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      std::string msg =
        ERRNO_EXCEPTION_ERR("pthread_rwlock_destroy", rc).GetFullMessage();
      fputs(msg.c_str(), stderr);
#endif // _NDEBUG
  }
}

#endif // _WIN32/!_WIN32

} // namespace ito33

#endif // _ITO33_RWLOCK_H_

