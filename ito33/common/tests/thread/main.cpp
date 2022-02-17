/////////////////////////////////////////////////////////////////////////////
// Name:        test/thread/main.cpp
// Purpose:     main file of thread test program
// Author:      Vadim Zeitlin
// Created:     18.02.03
// RCS-ID:      $Id: main.cpp,v 1.10 2006/04/03 09:19:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/timer.h"
#include "ito33/cppunit.h"
#include "ito33/thread.h"
#include "ito33/thread/condition.h"

#include <stdio.h>

using namespace ito33;

// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

static long g_counter = 0;
static CriticalSection g_csIncDec;
static Time::Ticks g_timeIncDec = 0;

static const int NUM_ITERATIONS = 100000;

static struct SharedData
{
  thread::Mutex mutex;
  thread::Condition cond;
  int value;
} g_data;

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

class ThreadTestCase : public CppUnit::TestCase
{
public:
  ThreadTestCase() { }

  virtual void setUp() { g_counter = 0; }

protected:
  CPPUNIT_TEST_SUITE( ThreadTestCase );
    CPPUNIT_TEST( Condition );
    CPPUNIT_TEST( IncDecUnprotected );
    CPPUNIT_TEST( IncDecAtomic );
    CPPUNIT_TEST( IncDecProtected );
    CPPUNIT_TEST( OnceOnly );
  CPPUNIT_TEST_SUITE_END();

  void Condition();
  void IncDecUnprotected();
  void IncDecAtomic();
  void IncDecProtected();
  void OnceOnly();

  NO_COPY_CLASS(ThreadTestCase);
};

// ----------------------------------------------------------------------------
// various increment/decrement classes
// ----------------------------------------------------------------------------

// no protection
struct Increment
{
  static void Do(long& value) { ++value; }
};

struct Decrement
{
  static void Do(long& value) { --value; }
};

#ifdef _WIN32

// using InterlockedXXX() functions
struct IncrementInterlocked
{
  static void Do(long& value) { InterlockedIncrement(&value); }
};

struct DecrementInterlocked
{
  static void Do(long& value) { InterlockedDecrement(&value); }
};

#endif // _WIN32

// using critical section
struct IncrementCritSect
{
  static void Do(long& value)
  {
    Lock<CriticalSection> lock(g_csIncDec);

    ++value;
  }
};

struct DecrementCritSect
{
  static void Do(long& value)
  {
    Lock<CriticalSection> lock(g_csIncDec);

    --value;
  }
};

// ----------------------------------------------------------------------------
// example of a thread class
// ----------------------------------------------------------------------------

template <class T>
class IncDecThreadClass
{
public:
  IncDecThreadClass(int n) : m_n(n) { }

  void operator()()
  {
    for ( int i = 0; i < m_n; ++i )
    {
      T::Do(g_counter);
    }
  }

private:
  int m_n;
};

// ----------------------------------------------------------------------------
// test incrementing/decrementing a variable using the given policy
// ----------------------------------------------------------------------------

// we *must* use a struct and not a function because of a bug in MSVC6
// (explicit function template instantiation silently doesn't work)
template <class Inc, class Dec>
struct TestIncDec
{
  TestIncDec(const char *msg, Time::Ticks& g_timeIncDec)
  {
    StopWatch sw;

    {
      JoinableThread
        thr1(MakeRunnable(IncDecThreadClass<Inc>(NUM_ITERATIONS))),
        thr2(MakeRunnable(IncDecThreadClass<Dec>(NUM_ITERATIONS)));
    }

    const Time::Ticks diff = sw.GetTicks();
    DEBUGPRINT(("After %s:\t%ld (took %lums", msg, g_counter, sw()));
    if ( g_timeIncDec )
    {
      DEBUGPRINT((" or %g times more",
          ((double)(LongLong)diff) / (LongLong)g_timeIncDec));
    }
    else
    {
      g_timeIncDec = diff;
    }

    DEBUGPRINT((")\n"));
  }
};

// ----------------------------------------------------------------------------
// Atomic functions test
// ----------------------------------------------------------------------------

void ThreadTestCase::IncDecUnprotected()
{
  g_timeIncDec = 0;

  TestIncDec<Increment, Decrement>("using ++/--", g_timeIncDec);

  // it shouldn't be safe to manipulate the variable without any protection!
  CPPUNIT_ASSERT( g_counter != 0 );
}

void ThreadTestCase::IncDecAtomic()
{
  // IncrementInterlocked() is Win32-only
#ifdef _WIN32
  TestIncDec<IncrementInterlocked, DecrementInterlocked>("InterlockedXXX()", g_timeIncDec);
#endif // _WIN32

  CPPUNIT_ASSERT_EQUAL( 0L, g_counter );
}

void ThreadTestCase::IncDecProtected()
{
  TestIncDec<IncrementCritSect, DecrementCritSect>("crit sect", g_timeIncDec);

  CPPUNIT_ASSERT_EQUAL( 0L, g_counter );
}

// ----------------------------------------------------------------------------
// RunOnce test
// ----------------------------------------------------------------------------

static void DoItOnceOnly()
{
  static RunOnce once;
  if ( once.IsFirstTime() )
  {
    DEBUGPRINT(("Incrementing g_counter once from thread %lx\n",
       static_cast<unsigned long>(Thread::GetCurrent()->GetId())));

    g_counter++;
  }
}

void ThreadTestCase::OnceOnly()
{
  JoinableThread thr1(MakeRunnable(DoItOnceOnly)),
                 thr2(MakeRunnable(DoItOnceOnly)),
                 thr3(MakeRunnable(DoItOnceOnly)),
                 thr4(MakeRunnable(DoItOnceOnly)),
                 thr5(MakeRunnable(DoItOnceOnly)),
                 thr6(MakeRunnable(DoItOnceOnly)),
                 thr7(MakeRunnable(DoItOnceOnly)),
                 thr8(MakeRunnable(DoItOnceOnly)),
                 thr9(MakeRunnable(DoItOnceOnly));

  thr1.Join();
  thr2.Join();
  thr3.Join();
  thr4.Join();
  thr5.Join();
  thr6.Join();
  thr7.Join();
  thr8.Join();
  thr9.Join();

  CPPUNIT_ASSERT_EQUAL( 1L, g_counter );
}

// ----------------------------------------------------------------------------
// Condition variable test
// ----------------------------------------------------------------------------

static void ModifySharedDataAndSignal()
{
  for ( ;; )
  {
    Thread::Sleep(rand());

    Lock<thread::Mutex> lock(g_data.mutex);
    if ( !--g_data.value )
    {
      g_data.cond.Signal();
      return;
    }
  }
}

void ThreadTestCase::Condition()
{
  g_data.value = 10;

  JoinableThread thr(MakeRunnable(ModifySharedDataAndSignal));

  {
    Lock<thread::Mutex> lock(g_data.mutex);
    while ( g_data.value )
    {
      g_data.cond.Wait(g_data.mutex);
    }

    CPPUNIT_ASSERT_EQUAL( 0, g_data.value );
  }

  thr.Join();
}

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(ThreadTestCase::suite());

  return runner.run("") ? 0 : 1;
}
