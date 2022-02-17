/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testcrosscurrency.h
// Purpose:     testing cross-currency
// Author:      Nabil Ouachani
// Created:     2006/03/02
// RCS-ID:      $Id: testcrosscurrency.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/common.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class SessionData;
  class ConvertibleBond;
}

}

class CrossCurrencyTest : public CppUnit::TestCase
{
public:

  CrossCurrencyTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( CrossCurrencyTest );
    
    CPPUNIT_TEST_EXCEPTION( CC_TriggerInCurrencyOfUnderlyingButNoFixedFXRate, 
                            ito33::Exception );
    
    CPPUNIT_TEST_EXCEPTION( FQ_ButNotCrossCurrencyInstrument, 
                            ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( FQ_NeededFixedFXRateMissing, ito33::Exception );
    
    CPPUNIT_TEST_EXCEPTION( CBOption_With_FQ_NotSupported, 
                            ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( CBOption_And_CB_Not_Same_Currency, 
                            ito33::Exception );
  
  CPPUNIT_TEST_SUITE_END();

  // Test parameter validity for cross-currency (CC)
  void CC_TriggerInCurrencyOfUnderlyingButNoFixedFXRate();
  
  // Test parameter validity for fixed quanto (FQ)
  void FQ_ButNotCrossCurrencyInstrument();
  void FQ_NeededFixedFXRateMissing();

  // Test parameter validity for CBOption with cross-currency CB.
  void CBOption_With_FQ_NotSupported();
  void CBOption_And_CB_Not_Same_Currency();

  // Useful functions
  ito33::shared_ptr<ito33::finance::SessionData> Init_CC_SessionData();
  
  ito33::shared_ptr<ito33::finance::ConvertibleBond> 
  InitCrossCurrencyCB(const ito33::shared_ptr<ito33::finance::SessionData>& 
                      pSessionData);
  
  ito33::shared_ptr<ito33::finance::ConvertibleBond> 
  InitCB(const ito33::shared_ptr<ito33::finance::SessionData>& pSessionData);

  NO_COPY_CLASS(CrossCurrencyTest);
};
