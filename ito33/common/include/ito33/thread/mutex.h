/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/mutex.h
// Purpose:     declares a Mutex class
// Author:      Vadim Zeitlin
// Created:     2005-10-07
// RCS-ID:      $Id: mutex.h,v 1.1 2005/10/10 12:22:16 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/thread/mutex.h
    @brief Declaration of Mutex class.

    Mutex is the same as CriticalSection for now, but it makes more sense to
    use mutexes than critical sections with Condition variables, so we add a
    synonym here.
 */

#ifndef _ITO33_THREAD_MUTEX_H_
#define _ITO33_THREAD_MUTEX_H_

#include "ito33/thread.h"

namespace ito33
{

namespace thread
{

/**
    Mutex is used for mutual exclusion between many threads and is useful for
    protecting a shared resource.
 */
typedef CriticalSection Mutex;

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_MUTEX_H_

