/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_volflathrpower.h
// Purpose:     header file for testing parametrizationVolFlatHRPower
// Author:      ITO33
// Created:     2005/07/29
// RCS-ID:      $Id: testparametrization_volflathrpower.h,v 1.3 2006/08/20 09:49:27 wang Exp $
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
  class ParametrizationVolFlatHRPower;
}
}



class Parametrization_VolFlatHRPowerTest : public CppUnit::TestCase
{
public:
  Parametrization_VolFlatHRPowerTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolFlatHRPowerTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Generate some time only hazard rates. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// The general derivtive list to fit
  ito33::shared_ptr<ito33::finance::Derivatives> m_pDerivatives;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolFlatHRPower> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolFlatHRPowerTest);
};

