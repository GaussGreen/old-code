/////////////////////////////////////////////////////////////////////////////
// Name:        tests/error/main.cpp
// Purpose:     header file of error example/test program
// Author:      Vadim Zeitlin
// Created:     12.03.04
// RCS-ID:      $Id: testerror.h,v 1.4 2006/06/13 15:29:40 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_TEST_ERROR_H_
#define _ITO33_TEST_ERROR_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"

class ErrorCodeTestCase : public CppUnit::TestCase
{
public:
  ErrorCodeTestCase() { }

private:
  CPPUNIT_TEST_SUITE( ErrorCodeTestCase );
    CPPUNIT_TEST( GetCode );
    CPPUNIT_TEST( OperatorInt );
    CPPUNIT_TEST( GetMessage );
  CPPUNIT_TEST_SUITE_END();

  void GetCode();
  void OperatorInt();
  void GetMessage();

  NO_COPY_CLASS(ErrorCodeTestCase);
};
#endif
