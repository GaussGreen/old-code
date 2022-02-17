/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_hrwithtimecomponent.h
// Purpose:     header file for testing parametrization_hrwithtimecomponent.h
// Author:      ITO33
// Created:     2004/11/29
// RCS-ID:      $Id: testparametrization_hrwithtimecomponent.h,v 1.7 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class TermStructureCDS;
  class TermStructureEDS;
  class TermStructureParBond;
  class SessionData;
}
namespace ihg
{
  class Volatility;
  class HazardRate;
  class ParametrizationHRWithTimeComponent;
}
}



class Parametrization_HRWithTimeComponentTest : public CppUnit::TestCase
{
public:
  Parametrization_HRWithTimeComponentTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_HRWithTimeComponentTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestVols );
    CPPUNIT_TEST( TestSpotComponents );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Test various vol functions
  void TestVols();

  // Test various spot components
  void TestSpotComponents();

  // Generate some time only hazard rates. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// The ParBond termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureParBond> m_ptsParBond;

  /// The CDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The EDS termstructure to fit (if defined)
  ito33::shared_ptr<ito33::finance::TermStructureEDS> m_ptsEDS;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationHRWithTimeComponent> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_HRWithTimeComponentTest);
};

