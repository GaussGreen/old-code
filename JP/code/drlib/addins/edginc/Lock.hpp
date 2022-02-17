// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 04/11/2006 Alex Fung
// (Adapted from driMagnet by Victor Paskhaver)
//
// $Header$
//

#ifndef EDG_LOCK_H
#define EDG_LOCK_H

#ifdef WIN32
#else
#include <pthread.h>
#endif

DRLIB_BEGIN_NAMESPACE

class Lock {
#ifdef WIN32
    void *  _mutex;
#else
    pthread_mutex_t _mutex;
#endif
public:
    Lock();
    void lock();
    void unlock();
    ~Lock();
};

class Guard {
    Lock* _lock;
  public:

    Guard(Lock* lock) : _lock(lock) {
        _lock->lock();
    }

    ~Guard() {
        _lock->unlock();
    }
};

DRLIB_END_NAMESPACE

#endif
