/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volpowerhrpower.h
// Purpose:     header file for testing parametrization_volpowerhr
// Author:      ITO33
// Created:     2005/01/03
// RCS-ID:      $Id: testparametrization_volpowerhrpower.h,v 1.4 2006/08/20 09:49:27 wang Exp $
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
  class ParametrizationVolPowerHRPower;
}
}



class Parametrization_VolPowerHRPowerTest : public CppUnit::TestCase
{
public:
  Parametrization_VolPowerHRPowerTest() 
  {
    m_dHRAlpha  = 0.01;
    m_dHRBeta   = 0.0;
    m_dVolAlpha = 0.001;
    m_dVolBeta  = -1.6;
  }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolPowerHRPowerTest );
    CPPUNIT_TEST( TestAsSpecified );
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

  /// Set generated test parameters. Return false if done testing.
  bool SetNextTest();

  /// The CDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The first option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption1;

  /// The second option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption2;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolPowerHRPower> 
    m_pParametrization;

  /// Test parameter values when generating data
  double m_dHRAlpha;
  double m_dHRBeta;
  double m_dVolAlpha;
  double m_dVolBeta;

  NO_COPY_CLASS(Parametrization_VolPowerHRPowerTest);
};

