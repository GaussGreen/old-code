/////////////////////////////////////////////////////////////////////////////
// Name:        testequity.h
// Purpose:     input validation for equity
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testequity.h,v 1.1 2004/11/19 16:16:02 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class EquityTest : public CppUnit::TestCase
{

public:
  EquityTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( EquityTest );
    CPPUNIT_TEST_EXCEPTION( NegativeSpotSharePrice, ito33::Exception );   
  CPPUNIT_TEST_SUITE_END();

  void NegativeSpotSharePrice();

  NO_COPY_CLASS(EquityTest);
};
