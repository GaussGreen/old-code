/////////////////////////////////////////////////////////////////////////////
// Name:        test/normaldistrib/testND.cpp
// Purpose:     header file for normal distribution function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: testND.h,v 1.2 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"

#include "ito33/cppunit.h"

// Test class
class NormalDistTest : public CppUnit::TestCase { 
public: 
  NormalDistTest() {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( NormalDistTest );
    CPPUNIT_TEST ( normalof0 );
    CPPUNIT_TEST ( normalof1_96 );
    CPPUNIT_TEST ( cumulatednormalof0 );
    CPPUNIT_TEST ( cumulatednormalof1_96 );
    CPPUNIT_TEST ( cumulnormsymetry );
    CPPUNIT_TEST ( normalofd1 );
  CPPUNIT_TEST_SUITE_END();
  
  void normalof0();
  void normalof1_96();
  void cumulatednormalof0();
  void cumulatednormalof1_96();
  void cumulnormsymetry();
  void normalofd1();

  NO_COPY_CLASS( NormalDistTest );
};
