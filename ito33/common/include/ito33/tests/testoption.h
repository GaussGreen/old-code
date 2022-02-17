/////////////////////////////////////////////////////////////////////////////
// Name:        testoption.h
// Purpose:     acceptance test for option
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testoption.h,v 1.1 2004/11/18 20:14:17 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class OptionTest : public CppUnit::TestCase
{

public:
  OptionTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( OptionTest );
    CPPUNIT_TEST_EXCEPTION( NegativeStrike, ito33::Exception );   
    CPPUNIT_TEST( Dump );
  CPPUNIT_TEST_SUITE_END();

  void NegativeStrike();
  void Dump();

  NO_COPY_CLASS(OptionTest);
};
