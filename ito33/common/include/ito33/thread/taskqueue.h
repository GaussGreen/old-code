/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/taskqueue.h
// Purpose:     class for task queuing
// Created:     2005/10/7
// RCS-ID:      $Id: taskqueue.h,v 1.1 2005/10/11 20:39:46 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_TASKQUEUE_H_
#define _ITO33_THREAD_TASKQUEUE_H_

/**
    @file   ito33/thread/taskqueue.h
    @brief  class for task queuing
 */

#include <boost/function.hpp>

#include "ito33/beforestd.h"
#include <queue>
#include "ito33/afterstd.h"

#include "ito33/thread/queue.h"

namespace ito33
{

namespace thread
{

typedef boost::function<void ()> Task;

typedef Queue<Task> TaskQueue;

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_TASKQUEUE_H_
