/////////////////////////////////////////////////////////////////////////////
// Name:        test/sharedptr/main.cpp
// Purpose:     test file of SharedPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: testsharedptr.cpp,v 1.10 2006/02/10 16:57:57 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/sharedptr.h"
#include "ito33/weakptr.h"
#include "ito33/weakalgo.h"
#include "ito33/selfptr.h"

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/thread.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testsharedptr.h"

using namespace ito33;

// check that we may declare a shared pointer for a forward declared class
class SharedPtrForwardTest;

SharedPtr<SharedPtrForwardTest> g_ptr;

class SharedPtrForwardTest
{
public:
  int n;
};

size_t Pointee2::ms_nObjects = 0;

ITO33_IMPLEMENT_SHAREDPTR(Pointee2);
ITO33_IMPLEMENT_SHAREDPTR(Derived);
ITO33_IMPLEMENT_SHAREDPTR(SharedPtrForwardTest);
ITO33_IMPLEMENT_SHAREDPTR(VoidSelfPtr);
ITO33_IMPLEMENT_SHAREDPTR(OneArgSelfPtr);
ITO33_IMPLEMENT_SHAREDPTR(TwoArgsSelfPtr);


void SharedPtrTestCase::Ctor()
{
  SharedPtr<Pointee2> ptr(new Pointee2);
}

void SharedPtrTestCase::CopyCtor()
{
  SharedPtr<Pointee2> ptr1(new Pointee2(17));
  SharedPtr<Pointee2> ptr2(ptr1);

  CPPUNIT_ASSERT( ptr1->GetValue() == ptr2->GetValue() );

  ptr2->SetValue(7);

  CPPUNIT_ASSERT( ptr1->GetValue() == ptr2->GetValue() );
}

void SharedPtrTestCase::AssignmentOperator()
{
  SharedPtr<Pointee2> ptr1(new Pointee2(17));
  SharedPtr<Pointee2> ptr2;
  ptr2 = ptr1;

  CPPUNIT_ASSERT( ptr1->GetValue() == ptr2->GetValue() );

  ptr2->SetValue(7);

  CPPUNIT_ASSERT( ptr1->GetValue() == ptr2->GetValue() );
}

void SharedPtrTestCase::AssignmentOperatorRaw()
{
  SharedPtr<Pointee2> ptr;
  ptr = new Pointee2(17);
}

void SharedPtrTestCase::Get()
{
  SharedPtr<Pointee2> ptr(new Pointee2(17));

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );

  ptr->SetValue(7);

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );

  ptr.get()->SetValue(9);

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );
}

void SharedPtrTestCase::OperatorBool()
{
  SharedPtr<Pointee2> ptr;

  CPPUNIT_ASSERT( !ptr );

  ptr = new Pointee2;

  CPPUNIT_ASSERT( ptr );
}

void SharedPtrTestCase::OperatorsStarAndArrow()
{
  SharedPtr<Pointee2> ptr(new Pointee2(17));

  CPPUNIT_ASSERT( (*ptr).GetValue() == ptr->GetValue() );

  (*ptr).SetValue(7);

  CPPUNIT_ASSERT( (*ptr).GetValue() == ptr->GetValue() );
}

void SharedPtrTestCase::ConversionFromDerived()
{
  SharedPtr<Derived> ptrDerived(new Derived(13));
  SharedPtr<Pointee2> ptr(ptrDerived);

  CPPUNIT_ASSERT( ptr->GetValue() == 13 );
}

template <class M>
void TestRelease()
{
  SharedPtr<Pointee2, M> ptr;
  CPPUNIT_ASSERT( !ptr.release() );

  Pointee2 *p = new Pointee2(17);
  SharedPtr<Pointee2, M> ptr2(p);
  ptr = ptr2;
  CPPUNIT_ASSERT( !ptr.release() );

  ptr2.reset();
  CPPUNIT_ASSERT_EQUAL( p, ptr.release() );
  CPPUNIT_ASSERT( !ptr.get() );

  CPPUNIT_ASSERT_EQUAL( 17, p->GetValue() );
  delete p;
}

void SharedPtrTestCase::Release()
{
  TestRelease<Policy::MTUnsafe>();
  TestRelease<Policy::MTSafe>();
}

void SharedPtrTestCase::Weak()
{
  typedef Policy::MTSafe M;
  typedef WeakPtr<Pointee2, M> WPtr;
  typedef SharedPtr<Pointee2, M> ShPtr;

  WPtr wp;
  CPPUNIT_ASSERT( !ShPtr(wp) );

  {
    ShPtr p1(new Pointee2(17));
    wp = p1;

    ShPtr p2(wp);
    CPPUNIT_ASSERT( p2 );
  }

  ShPtr p3(wp);
  CPPUNIT_ASSERT( !p3 );
}

// helper for WeakForEach test: this is used as ForEachWeak functor
class SumValues
{
public:
  SumValues() { m_sum = 0; }

  void operator()(const Pointee2& x) { m_sum += x.GetValue(); }

  int Result() const { return m_sum; }

private:
  int m_sum;
};

void SharedPtrTestCase::WeakForEach()
{
  typedef WeakPtr<Pointee2> WPtr;
  typedef SharedPtr<Pointee2> ShPtr;

  // first test iterating over a vector
  std::vector<WPtr> wv;
  ShPtr p1(new Pointee2(1));
  wv.push_back(p1);

  ShPtr p2(new Pointee2(2));
  wv.push_back(p2);

  ShPtr p3(new Pointee2(3));
  wv.push_back(p3);

  const std::vector<WPtr>::iterator bv = wv.begin(),
                                    ev = wv.end();

  // this should take into account all the pointers
  CPPUNIT_ASSERT_EQUAL(6, ForEachWeak(bv, ev, SumValues()).Result());

  // one pointer has expired, the corresponding weak pointer should be skipped
  p3.reset();
  CPPUNIT_ASSERT_EQUAL(3, ForEachWeak(bv, ev, SumValues()).Result());

  // another one...
  p1.reset();
  CPPUNIT_ASSERT_EQUAL(2, ForEachWeak(bv, ev, SumValues()).Result());


  // now test iterator over a map
  typedef std::map<int, WPtr> WMap;
  WMap wm;

  wm[1] = p1;
  wm[2] = p2;
  wm[3] = p3;

  const WMap::iterator bm = wm.begin(),
                       em = wm.end();
  CPPUNIT_ASSERT_EQUAL(2, ForEachWeak(bm, em, SumValues()).Result());
}

void SharedPtrTestCase::ThreadSafety()
{
  typedef SharedPtr<Pointee2, Policy::MTSafe> Ptr;
  typedef WeakPtr<Pointee2, Policy::MTSafe> WPtr;

  Ptr ptr(new Pointee2(1));

  {
    class TestRunnable : public Runnable
    {
    public:
      TestRunnable(const Ptr& ptr, unsigned sleep)
          : m_ptr(ptr), m_sleep(sleep) { }

      virtual void Run()
      {
        for ( size_t n = 0; n < 1000; ++n )
        {
          Ptr ptr = GetPtr();
          Ptr ptr2 = ptr;
          m_ptr = ptr;
          if ( m_sleep && (n % m_sleep) == 0 )
              Thread::Sleep(10);
        }
      }

      virtual void OnThreadExit() { }

    private:
      const Ptr& GetPtr() const { return m_ptr; }

      Ptr m_ptr;
      unsigned m_sleep;
    } run1(ptr, 0),
      run2(ptr, 0),
      run3(ptr, 2),
      run4(ptr, 7);

    JoinableThread thread1(&run1),
                   thread2(&run2);

    CPPUNIT_ASSERT( ptr->GetValue() == 1 );

    thread2.Join();
    thread1.Join();

    JoinableThread thread3(&run3),
                   thread4(&run4);

    CPPUNIT_ASSERT( ptr->GetValue() == 1 );

    thread4.Join();
    thread3.Join();
  }

  CPPUNIT_ASSERT( ptr->GetValue() == 1 );

  WPtr wptr(ptr);
  ptr.reset();
  ptr = wptr;

  CPPUNIT_ASSERT( !ptr );
}

void SharedPtrTestCase::SelfPtrClasses()
{
  VoidSelfPtr::Ptr ptr0(VoidSelfPtr::Create());
  CPPUNIT_ASSERT( ptr0->Check() );
  CPPUNIT_ASSERT( ptr0->DoCheck(ptr0->Self()) );

  OneArgSelfPtr::Ptr ptr1(OneArgSelfPtr::Create(17));
  CPPUNIT_ASSERT( ptr1->Check() );
  CPPUNIT_ASSERT( ptr1->DoCheck(ptr1->Self()) );

  TwoArgsSelfPtr::Ptr ptr2(TwoArgsSelfPtr::Create(17, "ok"));
  CPPUNIT_ASSERT( ptr2->Check() );
  CPPUNIT_ASSERT( ptr2->DoCheck(ptr2->Self()) );
}

// helper classes for PtrFromThis test
struct Monitor : public EnableSharedFromThis<Monitor, Policy::MTSafe>
{
  Monitor(int n) : m_n(n) { }
  ~Monitor()
  {
    // no shared pointers to us should exist any longer
    CPPUNIT_ASSERT( !SelfPtr() );
  }

  void OnWhatever() { m_n++; }

  int m_n;
};

ITO33_IMPLEMENT_SHAREDPTR(Monitor);

class Observer
{
public:
  Observer() { }

  void Register(const SharedPtr<Monitor, Policy::MTSafe>& monitor)
  {
    monitor->OnWhatever();
  }
};

void SharedPtrTestCase::PtrFromThis()
{
  SharedPtr<Monitor, Policy::MTSafe> pMon(new Monitor(17));
  CPPUNIT_ASSERT_EQUAL( 17, pMon->m_n );

  CPPUNIT_ASSERT_EQUAL( pMon, pMon->SelfPtr() );

  Observer obs;
  obs.Register(pMon);
  CPPUNIT_ASSERT_EQUAL( 18, pMon->m_n );

  obs.Register(pMon->SelfPtr());
  CPPUNIT_ASSERT_EQUAL( 19, pMon->m_n );
}

