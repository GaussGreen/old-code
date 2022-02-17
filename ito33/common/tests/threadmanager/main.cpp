/////////////////////////////////////////////////////////////////////////////
// Name:        test/threadmanager/main.cpp
// Purpose:     main file of ThreadManager test program
// Created:     2005/10/10
// RCS-ID:      $Id: main.cpp,v 1.3 2005/11/22 15:12:54 zhang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <iostream>

#include "ito33/debug.h"
#include "ito33/timer.h"
#include "ito33/cppunit.h"
#include "ito33/thread.h"
#include "ito33/thread/condition.h"
#include "ito33/thread/threadmanager.h"

using namespace ito33;

static CriticalSection g_csWork;
static int g_workDone = 0;
static CriticalSection g_csQuickWork;
static int g_quickWorkDone = 0;

// ----------------------------------------------------------------------------
// classes
// ----------------------------------------------------------------------------

class ThreadManagerTestCase : public CppUnit::TestCase
{
public:
  ThreadManagerTestCase() { }

protected:
  CPPUNIT_TEST_SUITE( ThreadManagerTestCase );

  CPPUNIT_TEST(SlowlyRun);
  CPPUNIT_TEST(Run);

  CPPUNIT_TEST_SUITE_END();

  void Run();
  void SlowlyRun();

  NO_COPY_CLASS(ThreadManagerTestCase);
};

// Just waste some time
void Work()
{
  double d;

  for ( int i = 0; i < rand(); i++)
    d = i * 100.;

  {
    Lock<CriticalSection> lock(g_csWork);

    g_workDone++;
  }
}

void QuickWork()
{
  Lock<CriticalSection> lock(g_csQuickWork);

  g_quickWorkDone++;
}

void ThreadManagerTestCase::SlowlyRun()
{
  thread::ThreadManager manager("test1", 2);
  int nTest = 100;
  for ( int i = 0; i < nTest; i++)
  {
    manager.QueueTask(Work);
    Thread::Sleep( 100 + 4000 * rand() / RAND_MAX );
  }

  CPPUNIT_ASSERT_EQUAL(nTest, g_workDone);
}

void ThreadManagerTestCase::Run()
{
  thread::ThreadManager manager("test2", 5);
  int nTest = 10000;
  for ( int i = 0; i < nTest; i++)
  {
    manager.QueueTask(QuickWork);
  }

  // wait a while
  Thread::Sleep( 4000 );
  CPPUNIT_ASSERT_EQUAL(nTest, g_quickWorkDone);
}

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(ThreadManagerTestCase::suite());

  return runner.run("") ? 0 : 1;
}

