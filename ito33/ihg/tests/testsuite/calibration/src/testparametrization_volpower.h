/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volpower.h
// Purpose:     header file for testing parametrization_volpower.h
// Author:      ITO33
// Created:     2004/12/29
// RCS-ID:      $Id: testparametrization_volpower.h,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class Option;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationVolPower;
}
}



class Parametrization_VolPowerTest : public CppUnit::TestCase
{
public:
  Parametrization_VolPowerTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolPowerTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Generate some vol powers. Should fit exactly
  void TestGeneratedValues();
  
  /// Setup the member variables by reading global xml file
  void Setup();

  /// Check if the given vol and hr fit the cds term structure
  void CheckFit( ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
                 ito33::shared_ptr<ito33::ihg::Volatility> pVol);

  /// The first option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption1;

  /// The second option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption2;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolPower> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolPowerTest);
};

