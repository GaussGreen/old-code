/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volwithtimecomponent.h
// Purpose:     header file for testing parametrization_volwithtimecomponent.h
// Created:     2006/01/30
// RCS-ID:      $Id: testparametrization_volwithtimecomponent.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2006- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class TermStructureOption;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationVolWithTimeComponent;
}
}



class Parametrization_VolWithTimeComponentTest : public CppUnit::TestCase
{
public:
  Parametrization_VolWithTimeComponentTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolWithTimeComponentTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestHRs );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Test various hazard rate functions
  void TestHRs();

  // Generate some time only vols. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// The Option termstructure to fit
  ito33::shared_ptr<ito33::finance::TermStructureOption> m_ptsOption;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolWithTimeComponent> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolWithTimeComponentTest);
};

