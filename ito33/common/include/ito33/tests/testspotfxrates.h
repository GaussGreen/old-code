/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/testspotfxrates.h
// Purpose:     test file for the class SpotFXRates
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: testspotfxrates.h,v 1.2 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004-  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

class SpotFXRatesTest : public CppUnit::TestCase 
{ 
public: 
  
  SpotFXRatesTest() {}

  void tearDown() {}


private:

  CPPUNIT_TEST_SUITE( SpotFXRatesTest );
  
    CPPUNIT_TEST_EXCEPTION ( GetWithoutSet,  ito33::Exception );
    
    CPPUNIT_TEST ( OneForSameNumeraires );
    CPPUNIT_TEST ( GetWhileSetInversely );
    CPPUNIT_TEST ( GetAfterUpdate );
    CPPUNIT_TEST ( MultiValues );

  CPPUNIT_TEST_SUITE_END();
  
  void GetWithoutSet();
  
  void OneForSameNumeraires();
  void GetWhileSetInversely();
  void GetAfterUpdate();
  void MultiValues();

  NO_COPY_CLASS( SpotFXRatesTest );

}; // Class SpotFXRatesTest
