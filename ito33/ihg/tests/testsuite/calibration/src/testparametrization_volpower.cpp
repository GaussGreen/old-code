#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"
#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/parametrization_volpower.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_volpower.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolPowerTest::Setup()
{

  // Reset the pricing variables before reading
 
  // Read from the currently active global input file
  ihg::XML::CalibrationReader reader(g_strInputFilename.c_str());

  // Get the session
  m_pSessionData = reader.ReadSessionData();

  CPPUNIT_ASSERT(m_pSessionData);

  // Get the parametrization
  ito33::ihg::ParametrizationVisitorGoodType param_visitor;
  ito33::finance::DerivativeVisitorGoodType deriv_visitor;
  ito33::finance::TermStructureEnumerator termstructures;
  shared_ptr<ito33::finance::Derivatives> 
    pDerivatives(new ito33::finance::Derivatives());
  reader.ReadCalibration(param_visitor, termstructures, deriv_visitor,
                         *pDerivatives);
  m_pParametrization = param_visitor.GetParametrizationVolPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect two optioms
  m_pOption1 = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption1);

  m_pOption2 = deriv_visitor.GetOption2();
  CPPUNIT_ASSERT(m_pOption2);
}


void Parametrization_VolPowerTest::CheckFit( 
  ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
  ito33::shared_ptr<ito33::ihg::Volatility> pVol)
{

  // Create and setup the pricing model
  ihg::TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  // compare the prices
  shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption1);
  double dScale = fabs( m_pOption1->GetMarketPrice() );
  if (dScale < 1.e-6)
    dScale = 1.0;

  double dError = fabs( output->GetPrice() - m_pOption1->GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);


  output = model.Compute(*m_pOption2);
  dScale = fabs( m_pOption2->GetMarketPrice() );
  if (dScale < 1.e-6)
    dScale = 1.0;

  dError = fabs( output->GetPrice() - m_pOption2->GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);

}


void Parametrization_VolPowerTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  shared_ptr<VolatilityPower>
    pVol = m_pParametrization->CalibrateWithOptions(*m_pOption1, *m_pOption2);

  // Check if it fits the data
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  
  CheckFit(pHR, pVol);
}


void Parametrization_VolPowerTest::TestGeneratedValues()
{
  
  // Read in the current data
  Setup();

  // Generate prices using vol powers.  Should be able to fit these prices 
  // exactly. However, local mins may be found, or the tolerance may be 
  // reached, so the parameters are not always the same
  double dS0 = m_pSessionData->GetSpotSharePrice();
  double dAlpha = 0.1;
  double dAlphaStep = 0.1;  
  while (dAlpha < 0.82)
  {
    double dBeta = -1.8;
    double dBetaStep = 0.2;

    while (dBeta < 1.81)
    {
      //dAlpha = 0.2;
      //dBeta = 1.6;

      //std::cout << "alpha = " << dAlpha
      //          << ", beta = " << dBeta
      //          << std::endl;

      shared_ptr<Volatility> pVol( new VolatilityPower(dAlpha, dBeta, dS0) );
      shared_ptr<HazardRate> pHR( new HazardRateFlat(0.1) );

      ihg::TheoreticalModel model;
      model.SetVolatility(pVol);
      model.SetHazardRate(pHR);

      // Take the computed prices as market price for testing
      shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption1);
      m_pOption1->SetMarketPrice( output->GetPrice() );

      output = model.Compute(*m_pOption2);
      m_pOption2->SetMarketPrice( output->GetPrice() );      
  
      // Calibrate
      m_pParametrization = make_ptr( new ParametrizationVolPower(pHR) );

      shared_ptr<VolatilityPower> pCalibratedVol = 
        m_pParametrization->CalibrateWithOptions(*m_pOption1, *m_pOption2);

      // Check if the data was reproduced exactly
      CheckFit(pHR, pCalibratedVol);

      // Check that the params were reproduced exactly
      // Local mins make this difficult/impossible
      /*
      double dCalibratedAlpha = pCalibratedVol->GetAlpha();
      double dCalibratedBeta = pCalibratedVol->GetBeta();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dAlpha, dCalibratedAlpha, 1.e-3 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dBeta, dCalibratedBeta, 1.e-3 );
      */

      dBeta += dBetaStep;
    } // loop over beta

    dAlpha += dAlphaStep;

  } // loop over alpha

}
