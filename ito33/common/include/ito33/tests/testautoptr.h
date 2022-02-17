/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testautoptr.h
// Purpose:     header file of AutoPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: testautoptr.h,v 1.3 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_AUTOPTR_H_
#define _ITO33_TEST_AUTOPTR_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// self-counting object
// ----------------------------------------------------------------------------

class Pointee
{
public:
  Pointee(int n = 0) : m_n(n) { ++ms_nObjects; }
  ~Pointee() { CPPUNIT_ASSERT( ms_nObjects > 0 ); --ms_nObjects; }

  void SetValue(int n) { m_n = n; }
  int GetValue() const { return m_n; }

  static void CheckNoLeft() { CPPUNIT_ASSERT( ms_nObjects == 0 ); }

private:
  int m_n;

  static size_t ms_nObjects;
};


class AutoPtrTestCase : public CppUnit::TestCase
{
public:
  AutoPtrTestCase() { }

  virtual void tearDown() { Pointee::CheckNoLeft(); }

private:
  CPPUNIT_TEST_SUITE( AutoPtrTestCase );
    CPPUNIT_TEST( Ctor );
    CPPUNIT_TEST( CopyCtor );
    CPPUNIT_TEST( AssignmentOperator );
    CPPUNIT_TEST( Get );
    CPPUNIT_TEST( OperatorBool );
    CPPUNIT_TEST( OperatorsStarAndArrow );
  CPPUNIT_TEST_SUITE_END();

  void Ctor();
  void CopyCtor();
  void AssignmentOperator();
  void Get();
  void OperatorBool();
  void OperatorsStarAndArrow();

  NO_COPY_CLASS(AutoPtrTestCase);
};
#endif
