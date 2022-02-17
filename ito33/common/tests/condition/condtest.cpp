/////////////////////////////////////////////////////////////////////////////
// Name:        test/condition/condtest.cpp
// Purpose:     main file of thread::Condition class test program
// Author:      Vadim Zeitlin
// Created:     2006-02-21
// RCS-ID:      $Id: condtest.cpp,v 1.2 2006/03/29 16:09:21 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/cppunit.h"
#include "ito33/thread.h"
#include "ito33/thread/condition.h"

#include <stdio.h>

using namespace ito33;
using namespace ito33::thread;

// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

static struct SharedData
{
  Mutex mutex;
  Condition cond;
  int value;
} g_data;

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

class ConditionTestCase : public CppUnit::TestCase
{
public:
  ConditionTestCase() { }

  virtual void setUp() { g_data.value = 0; }

protected:
  CPPUNIT_TEST_SUITE( ConditionTestCase );
    CPPUNIT_TEST( Signal );
    CPPUNIT_TEST( Broadcast );
  CPPUNIT_TEST_SUITE_END();

  void Signal();
  void Broadcast();

  NO_COPY_CLASS(ConditionTestCase);
};

void ConditionTestCase::Signal()
{
  struct Inc : Runnable
  {
    virtual void Run()
    {
      Lock<Mutex> lock(g_data.mutex);
      g_data.value++;
      g_data.cond.Signal();
    }
  };

  g_data.mutex.Lock();
  JoinableThread thr(new Inc());

  g_data.value = 16;

  Thread::Sleep(1000000);
  CPPUNIT_ASSERT_EQUAL( 16, g_data.value );

  g_data.cond.Wait(g_data.mutex);

  CPPUNIT_ASSERT_EQUAL( 17, g_data.value );

  g_data.mutex.Unlock();

  thr.Join();
}

void ConditionTestCase::Broadcast()
{
  static const unsigned N = 10;

  JoinableThread *threads[N];
  unsigned n;
  for ( n = 0; n < N; n++ )
  {
    struct Waiter : Runnable
    {
      virtual void Run()
      {
        Lock<Mutex> lock(g_data.mutex);
        g_data.value++;
        g_data.cond.Wait(g_data.mutex);
      }
    };

    threads[n] = new JoinableThread(new Waiter);
  }

  // wait until all threads are running
  for ( ;; )
  {
    Lock<Mutex> lock(g_data.mutex);
    if ( g_data.value == N )
    {
      //g_data.cond.Broadcast();
      for ( n = 0; n < N; n++ )
        g_data.cond.Signal();

      break;
    }
  }

  for ( n = 0; n < N; n++ )
  {
    threads[n]->Join();
    delete threads[n];
  }
}

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(ConditionTestCase::suite());

  return runner.run("") ? 0 : 1;
}
