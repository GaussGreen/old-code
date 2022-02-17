/////////////////////////////////////////////////////////////////////////////
// Name:        test/weakptr/main.cpp
// Purpose:     test file of WeakPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: testweakptr.cpp,v 1.2 2005/12/16 16:57:52 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/weakptr.h"

#include "ito33/list.h"
#include "ito33/vector.h"

#include "ito33/cppunit.h"

using namespace ito33;


// ----------------------------------------------------------------------------
// self-counting object
// 
// Called Pointee2, since the autoptr test uses the same class. While the
// individual tests are not mixed (and arguably should not be), both tests
// are included in other projects.  Hence, this class needs a different name
// ----------------------------------------------------------------------------

class Pointee2
{
public:
  Pointee2(int n = 0) : m_n(n) { ++ms_nObjects; }
  ~Pointee2() { CPPUNIT_ASSERT( ms_nObjects > 0 ); --ms_nObjects; }

  void SetValue(int n) { m_n = n; }
  int GetValue() const { return m_n; }

  static void CheckNoLeft() { CPPUNIT_ASSERT( ms_nObjects == 0 ); }

private:
  int m_n;

  static size_t ms_nObjects;
};

class Derived : public Pointee2
{
public:
  Derived(int n) : Pointee2(n) { }
};

// ----------------------------------------------------------------------------
// test class
// ----------------------------------------------------------------------------

class WeakPtrTestCase : public CppUnit::TestCase
{
public:
  WeakPtrTestCase() { }

  virtual void tearDown() { Pointee2::CheckNoLeft(); }

private:
  CPPUNIT_TEST_SUITE( WeakPtrTestCase );
    CPPUNIT_TEST( Ctor );
    CPPUNIT_TEST( CopyCtor );
    CPPUNIT_TEST( AssignmentOperator );
    CPPUNIT_TEST( AssignmentOperatorSharedPtr );
  CPPUNIT_TEST_SUITE_END();

  void Ctor();
  void CopyCtor();
  void AssignmentOperator();
  void AssignmentOperatorSharedPtr();

  NO_COPY_CLASS(WeakPtrTestCase);
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------

// check that we may declare a shared pointer for a forward declared class
class WeakPtrForwardTest;

WeakPtr<WeakPtrForwardTest> g_ptr;

class WeakPtrForwardTest
{
public:
  int n;
};

size_t Pointee2::ms_nObjects = 0;

ITO33_IMPLEMENT_SHAREDPTR(Pointee2);
ITO33_IMPLEMENT_SHAREDPTR(Derived);
ITO33_IMPLEMENT_SHAREDPTR(WeakPtrForwardTest);


void WeakPtrTestCase::Ctor()
{
  SharedPtr<Pointee2> shptr(new Pointee2);
  WeakPtr<Pointee2> ptr(shptr);
}

void WeakPtrTestCase::CopyCtor()
{
  SharedPtr<Pointee2> shptr(new Pointee2(17));
  WeakPtr<Pointee2> ptr1(shptr);
  WeakPtr<Pointee2> ptr2(ptr1);
}

void WeakPtrTestCase::AssignmentOperator()
{
  SharedPtr<Pointee2> shptr(new Pointee2(17));
  WeakPtr<Pointee2> ptr1(shptr);
  WeakPtr<Pointee2> ptr2;
  ptr2 = ptr1;

  CPPUNIT_ASSERT( SharedPtr<Pointee2>(ptr1)->GetValue() == shptr->GetValue() );

  SharedPtr<Pointee2>(ptr2)->SetValue(7);

  CPPUNIT_ASSERT( shptr->GetValue() == 7 );
}

void WeakPtrTestCase::AssignmentOperatorSharedPtr()
{
  SharedPtr<Pointee2> shptr(new Pointee2(17));
  WeakPtr<Pointee2> ptr;
  ptr = shptr;
}

CPPUNIT_TEST_SUITE_REGISTRATION(WeakPtrTestCase);
