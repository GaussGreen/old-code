/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_voltimeonlyhrwithtimecomponent.h
// Purpose:     header file for testing parametrization_voltimeonlthrwithtimecomponent.h
// Author:      ITO33
// Created:     2005/07/15
// RCS-ID:      $Id: testparametrization_voltimeonlyhrwithtimecomponent.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class Derivatives;
  class SessionData;
}
namespace ihg
{
  class ParametrizationVolTimeOnlyHRWithTimeComponent;
}
}



class Parametrization_VolTimeOnlyHRWithTimeComponentTest : public CppUnit::TestCase
{
public:
  Parametrization_VolTimeOnlyHRWithTimeComponentTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolTimeOnlyHRWithTimeComponentTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestSpotComponents );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Test various spot components
  void TestSpotComponents();

  // Generate some time only hazard rates. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// The general derivtive list to fit
  ito33::shared_ptr<ito33::finance::Derivatives> m_pDerivatives;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolTimeOnlyHRWithTimeComponent> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolTimeOnlyHRWithTimeComponentTest);
};

