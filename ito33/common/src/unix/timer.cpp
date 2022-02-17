/////////////////////////////////////////////////////////////////////////////
// Name:        unix/timer.cpp
// Purpose:     implementation of timer-related functions
// Author:      Vadim Zeitlin
// Created:     25.05.03
// RCS-ID:      $Id: timer.cpp,v 1.5 2005/04/02 14:03:08 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <sys/time.h>

#include "ito33/common.h"
#include "ito33/timer.h"

#include "ito33/errnoexception.h"

using namespace ito33;

static const long MSEC_PER_SEC = 1000*1000;

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

/*
   We assume gettimeofday() is always available because it usually is in all
   modern day Unices, but if not we should add a test for it to configure.

   If we do this we could also use gethrtime() under Solaris which is more like
   Win32 ::QueryPerformanceCounter() than gettimeofday().
 */

Time::Ticks ito33::Time::GetCurrentTicks()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);

  Time::Ticks timeInUsec = tv.tv_sec;
  timeInUsec *= MSEC_PER_SEC;
  timeInUsec += tv.tv_usec;

  return timeInUsec;
}

Time::Ticks ito33::Time::GetTicksPerSecond()
{
  return MSEC_PER_SEC;
}

