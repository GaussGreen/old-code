/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/timer.h
// Purpose:     functions and classes related to timing issues
// Author:      Vadim Zeitlin
// Created:     25.05.03
// RCS-ID:      $Id: timer.h,v 1.4 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/timer.h
    @brief  Functions to measure time intervals.

    Using only the standard C functions we can't measure time with precision
    better than 1 second which is woefully inadequate. We thus define functions
    to allow doing it with at least 1ms resolution and much better on some
    platforms.
 */

#ifndef _ITO33_TIMER_H_
#define _ITO33_TIMER_H_

#include "ito33/longlong.h"

namespace ito33
{

/**
    Time functions.
 */
namespace Time
{
  /// the type used for the timer ticks
  typedef ULongLong Ticks;

  /// get the current time in ticks
  Ticks GetCurrentTicks();

  /// get the number of ticks per second
  Ticks GetTicksPerSecond();
} // namespace Time


/**
    StopWatch allows measuring time intervals easily.

    To measure the elapsed time simply create a stop watch object as a local
    variable and then call its GetTime() function to get the elapsed time since
    its construction.

    Here is an example of using StopWatch for measuring its own overhead:
    @code
        for ( unsigned long N = 1; N <= 10000000; N *= 10 )
        {
            StopWatch sw;

            for ( unsigned long i = 0; i < N; ++i )
            {
                Time::GetCurrentTicks();
            }

            printf("%lu calls of GetCurrentTicks() took %lums\n", N, sw());
        }
    @endcode
 */
class StopWatch
{
public:
  /**
      Constructor immediately starts the stop watch.
   */
  StopWatch() { m_ticks0 = Time::GetCurrentTicks(); }

  /**
      Get the number of elapsed ticks.

      This is in general less useful than GetTime() but has two advantages:
          - you may get better resolution if needed (if the OS supports it)
          - the stop watch doesn't overflow after about 5 days (on 32b arch)

      You may use Time::GetTicksPerSecond() to convert ticks to abs time.
   */
  Time::Ticks GetTicks() const { return Time::GetCurrentTicks() - m_ticks0; }

  /**
      Get the time elapsed since the ctor was called.

      @return elapsed time in milliseconds
   */
  unsigned long GetTime() const
  {
    const Time::Ticks ticksPerSec = Time::GetTicksPerSecond();
    if ( !ticksPerSec )
      return (unsigned long)-1;

    // multiply by 1000 because we want the milliseconds
    return (unsigned long)((GetTicks() * 1000) / ticksPerSec);
  }

  /// implicit call of GetTime()
  unsigned long operator()() const { return GetTime(); }

private:
  /// the number of ticks when we started
  Time::Ticks m_ticks0;
};

} // namespace ito33

#endif // _ITO33_TIMER_H_

