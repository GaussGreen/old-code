/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volflathrwithspotcomponentpower.h
// Purpose:     header file for testing parametrization_volflathrwithspotcomponentpower.h
// Author:      ITO33
// Created:     2004/12/10
// RCS-ID:      $Id: testparametrization_volflathrwithspotcomponentpower.h,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class TermStructureCDS;
  class Option;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationVolFlatHRWithSpotComponentPower;
}
}



class Parametrization_VolFlatHRWithSpotComponentPowerTest : public CppUnit::TestCase
{
public:
  Parametrization_VolFlatHRWithSpotComponentPowerTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolFlatHRWithSpotComponentPowerTest );
    //CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Generate data by varying vol and power hazard rates
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// Check if the given vol and hr fit the cds term structure
  void CheckCDSFit( ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
                    ito33::shared_ptr<ito33::ihg::Volatility> pVol);

  /// Check if the given vol and hr fit the option
  void CheckOptionFit( ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
                       ito33::shared_ptr<ito33::ihg::Volatility> pVol);

  /// The CDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolFlatHRWithSpotComponentPower> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolFlatHRWithSpotComponentPowerTest);
};

