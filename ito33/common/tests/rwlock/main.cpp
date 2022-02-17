/////////////////////////////////////////////////////////////////////////////
// Name:        test/rwlock/main.cpp
// Purpose:     main file of RWLock test program
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: main.cpp,v 1.6 2006/02/18 23:45:06 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/thread/rwlock.h"
#include "ito33/thread.h"
#include "ito33/timer.h"

// uncomment next line to verify that without locks the program doesn't work
// (the invariant can be broken)
//#define NO_LOCK

using namespace ito33;

// ----------------------------------------------------------------------------
// global variable
// ----------------------------------------------------------------------------

// this struct is protected by RWLock
struct Data
{
  // the lock protecting the global data
  RWLock lock;

  int x;
  int y;

  Data() { x = y = 0; }

private:
  // can't implement those because of RWLock
  Data(const Data&);
  Data& operator=(const Data&);
} g_data;

// ----------------------------------------------------------------------------
// local classes
// ----------------------------------------------------------------------------

class LoopThread
{
public:
  enum { NUM_ITERATIONS = 100 };

  LoopThread(const char *name, unsigned long delay) : m_name(name)
  {
    m_delay = delay;
  }

  void operator()()
  {
    for ( int i = 0; i < NUM_ITERATIONS; ++i )
    {
      Do();

      // lock for reading before accessing g_data
#ifndef NO_LOCK
      LockRead lockR(g_data.lock);
#endif // !NO_LOCK

      // if the data structure is in the consistent state the sum should
      // be always 0
      const int sum = g_data.x + g_data.y;
      if ( sum )
        printf("%s(%i):\tsum of values is %d.\n", m_name, i + 1, sum);

      Thread::Sleep(m_delay);
    }
  }

  virtual ~LoopThread() { }

protected:
  virtual void Do() const = 0;

  const char * const m_name;
  unsigned m_delay;

private:
  LoopThread& operator=(const LoopThread&);
};

class ReadThread : public LoopThread
{
public:
  ReadThread(const char *name, unsigned long delay) : LoopThread(name, delay)
    { }

  virtual ~ReadThread() { }

protected:
  virtual void Do() const
  {
  }

private:
  ReadThread& operator=(const ReadThread&);
};

class WriteThread : public LoopThread
{
public:
  WriteThread(const char *name, unsigned long delay) : LoopThread(name, delay)
    { }

  virtual ~WriteThread() { }

protected:
  virtual void Do() const
  {
    // try commenting out the line below: you'll see that the sum of the
    // values in the reader thread is going to be different from 0
    // sometimes then
#ifndef NO_LOCK
    LockWrite lockW(g_data.lock);
#endif // !NO_LOCK


    // change first one field
    ++g_data.x;

    // then wait a little... (1ms)
    Thread::Sleep(1000);

    // and change another one
    --g_data.y;
  }

private:
  WriteThread& operator=(const WriteThread&);
};

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  StopWatch sw;

  JoinableThread thrRead1(MakeRunnable(ReadThread("A", 15000))),
         thrRead2(MakeRunnable(ReadThread("B", 20000))),
         thrRead3(MakeRunnable(ReadThread("C", 10000))),
         thrWrite(MakeRunnable(WriteThread("W", 10000)));

  thrRead1.Join();
  thrRead2.Join();
  thrRead3.Join();
  thrWrite.Join();

  printf("Total time taken: %lums\n", sw());

  return 0;
}

