/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_hrwithspotcomponentpower.h
// Purpose:     header file for testing parametrization_hrwithspotcomponentpower.h
// Author:      ITO33
// Created:     2004/12/08
// RCS-ID:      $Id: testparametrization_hrwithspotcomponentpower.h,v 1.6 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class TermStructureCDS;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationHRWithSpotComponentPower;
}
}



class Parametrization_HRWithSpotComponentPowerTest : public CppUnit::TestCase
{
public:
  Parametrization_HRWithSpotComponentPowerTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_HRWithSpotComponentPowerTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestVols );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Test various vol functions
  void TestVols();

  // Generate some time only hazard rates. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// Check if the given vol and hr fit the cds term structure
  void CheckCDSFit( ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
                    ito33::shared_ptr<ito33::ihg::Volatility> pVol);

  /// The CDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationHRWithSpotComponentPower> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_HRWithSpotComponentPowerTest);
};

