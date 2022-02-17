/////////////////////////////////////////////////////////////////////////////
// Name:        win32/timer.cpp
// Purpose:     implementation of timer-related functions
// Author:      Vadim Zeitlin
// Created:     25.05.03
// RCS-ID:      $Id: timer.cpp,v 1.3 2004/10/05 09:13:47 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/timer.h"

#include "ito33/win32/winwrap.h"
#include "ito33/win32/exception.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

Time::Ticks ito33::Time::GetCurrentTicks()
{
  LARGE_INTEGER ticks;
  ::QueryPerformanceCounter(&ticks);

  return ticks.QuadPart;
}

Time::Ticks ito33::Time::GetTicksPerSecond()
{
  // we don't protect this static variable here from multiple threads because
  // doing it would skew the measurement significantly: cost of
  // entering/leaving a critical section may be high compared with thei
  // ntervals of time we want to measure, so sacrifice correctness (which may
  // be difficult to achieve when using QueryPerformanceCounter() on an SMP
  // system as it may return different values for different CPUs!) for
  // efficiency here, exceptionally
  static Time::Ticks s_freq = 0;

  if ( !s_freq )
  {
    LARGE_INTEGER freq;
    if ( !::QueryPerformanceFrequency(&freq) || !freq.QuadPart )
    {
      throw WIN32_EXCEPTION("QueryPerformanceFrequency");
    }

    s_freq = freq.QuadPart;
  }

  return s_freq;
}

