/////////////////////////////////////////////////////////////////////////////
// Name:        test/sharedptr/main.cpp
// Purpose:     header file of SharedPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: testsharedptr.h,v 1.9 2006/02/09 15:43:10 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_TEST_SHAREDPTR_H_
#define _ITO33_TEST_SHAREDPTR_H_
// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/selfptr.h"

#include "ito33/cppunit.h"

// normally there should be no "using" in header but this is just a test...
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
// family of classes keeping a weak pointer to themselves
// ----------------------------------------------------------------------------

class VoidSelfPtr : public SelfPointer0<VoidSelfPtr>
{
public:
  bool Check() { return DoCheck(m_self); }

  static bool DoCheck(VoidSelfPtr::Ptr ptr)
  {
    return ptr != NULL;
  }

protected:
  friend SelfPointerBase;

  VoidSelfPtr() { }
};

class OneArgSelfPtr : public SelfPointer1<OneArgSelfPtr, int>
{
public:
  bool Check() { return DoCheck(m_self); }

  static bool DoCheck(OneArgSelfPtr::Ptr ptr)
  {
    return ptr != NULL && ptr->m_x == 17;
  }

protected:
  friend SelfPointerBase;

  OneArgSelfPtr(int x) : m_x(x) { }

  int m_x;
};

// NB: this shows that SelfPointer<> template parameters can differ from actual
//     class ctor parameters types as long as they're compatible
class TwoArgsSelfPtr : public SelfPointer2<TwoArgsSelfPtr, int, const char *>
{
public:
  bool Check() { return DoCheck(m_self); }

  static bool DoCheck(TwoArgsSelfPtr::Ptr ptr)
  {
    return ptr != NULL && ptr->m_x == 17 && ptr->m_s == "ok";
  }

protected:
  friend SelfPointerBase;

  TwoArgsSelfPtr(int x, const std::string& s) : m_x(x), m_s(s) { }

  int m_x;
  std::string m_s;
};

// ----------------------------------------------------------------------------
// test class
// ----------------------------------------------------------------------------

class SharedPtrTestCase : public CppUnit::TestCase
{
public:
  SharedPtrTestCase() { }

  virtual void tearDown() { Pointee2::CheckNoLeft(); }

private:
  CPPUNIT_TEST_SUITE( SharedPtrTestCase );
    CPPUNIT_TEST( Ctor );
    CPPUNIT_TEST( CopyCtor );
    CPPUNIT_TEST( AssignmentOperator );
    CPPUNIT_TEST( AssignmentOperatorRaw );
    CPPUNIT_TEST( Get );
    CPPUNIT_TEST( OperatorBool );
    CPPUNIT_TEST( OperatorsStarAndArrow );
    CPPUNIT_TEST( ConversionFromDerived );
    CPPUNIT_TEST( Release );
    CPPUNIT_TEST( Weak );
    CPPUNIT_TEST( WeakForEach );
    CPPUNIT_TEST( ThreadSafety );
    CPPUNIT_TEST( SelfPtrClasses );
    CPPUNIT_TEST( PtrFromThis );
  CPPUNIT_TEST_SUITE_END();

  void Ctor();
  void CopyCtor();
  void AssignmentOperator();
  void AssignmentOperatorRaw();
  void Get();
  void OperatorBool();
  void OperatorsStarAndArrow();
  void ConversionFromDerived();
  void Release();
  void Weak();
  void WeakForEach();
  void ThreadSafety();
  void SelfPtrClasses();
  void PtrFromThis();

  NO_COPY_CLASS(SharedPtrTestCase);
};

#endif // _ITO33_TEST_SHAREDPTR_H_

