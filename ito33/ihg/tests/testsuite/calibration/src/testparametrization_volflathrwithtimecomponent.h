/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volflathrwithtimecomponent.h
// Purpose:     header file for testing parametrization_volflathrwithtimecomponent.h
// Author:      ITO33
// Created:     2004/12/10
// RCS-ID:      $Id: testparametrization_volflathrwithtimecomponent.h,v 1.6 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class TermStructureParBond;
  class TermStructureCDS;
  class Option;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationVolFlatHRWithTimeComponent;
}
}



class Parametrization_VolFlatHRWithTimeComponentTest : public CppUnit::TestCase
{
public:
  Parametrization_VolFlatHRWithTimeComponentTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolFlatHRWithTimeComponentTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestGeneratedValues );
    CPPUNIT_TEST( TestGeneratedValuesCalibrateXXX );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Generate data by varying vol and power hazard rates
  void TestGeneratedValues();
  
  //test using specific calibration function
  void TestGeneratedValuesCalibrateXXX();


  /// Setup the member variables by reading global xml file
  void Setup();

  /// Check if the given vol and hr fit given instruments
  void CheckFitAll( ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
                    ito33::shared_ptr<ito33::ihg::Volatility> pVol);

  /// The ParBond termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureParBond> m_ptsParBond;

  /// The CDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolFlatHRWithTimeComponent> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolFlatHRWithTimeComponentTest);
};

