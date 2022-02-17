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

#include "ito33/ihg/parametrization_voltanh.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_voltanh.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolTanhTest::Setup()
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
  m_pParametrization = param_visitor.GetParametrizationVolTanh();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect two optioms
  m_pOption1 = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption1);

  m_pOption2 = deriv_visitor.GetOption2();
  CPPUNIT_ASSERT(m_pOption2);
}


void Parametrization_VolTanhTest::CheckFit( 
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


void Parametrization_VolTanhTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  shared_ptr<VolatilityTanh>
    pVol = m_pParametrization->CalibrateWithOptions(*m_pOption1, *m_pOption2);

  // Check if it fits the data
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  
  CheckFit(pHR, pVol);
}


void Parametrization_VolTanhTest::TestGeneratedValues()
{
  
  // Read in the current data
  Setup();

  
  // Generate prices using vol tanh.  Should be able to fit these prices 
  // exactly. However, local mins may be found, or the tolerance may be 
  // reached, so the parameters are not always the same
  double dStrike1 = m_pOption1->GetStrike();
  double dStrike2 = m_pOption2->GetStrike();
    
  double dScale, dS0 = 0.5 * (dStrike1 + dStrike2);;
  if ( dStrike1 < dStrike2 )
    dScale = 2. / (dStrike2 - dStrike1);
  else
    dScale = 2. / (dStrike1 - dStrike2);


  double dLeft = 0.1;
  double dLeftStep = 0.1;  
  while (dLeft < 0.81)
  {
    double dRight = 0.8;
    double dRightStep = -0.1;

    while (dRight > 0.09)
    {
      //dLeft = 0.3;
      //dRight = 0.1;

      //std::cout << "left = " << dLeft
      //          << ", right = " << dRight
      //          << std::endl;

      shared_ptr<Volatility> pVol( new VolatilityTanh(dLeft, dRight, dScale, dS0) );
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
      m_pParametrization = make_ptr( new ParametrizationVolTanh(pHR) );

      shared_ptr<VolatilityTanh> pCalibratedVol = 
        m_pParametrization->CalibrateWithOptions(*m_pOption1, *m_pOption2);

      // Check if the data was reproduced exactly
      CheckFit(pHR, pCalibratedVol);

      // Check that the params were reproduced exactly
      // Local mins make this difficult/impossible
      
      //double dCalibratedLeft = pCalibratedVol->GetLeft();
      //double dCalibratedRight = pCalibratedVol->GetRight();

      //CPPUNIT_ASSERT_DOUBLES_EQUAL( dLeft, dCalibratedLeft, 1.e-3 );
      //CPPUNIT_ASSERT_DOUBLES_EQUAL( dRight, dCalibratedRight, 1.e-3 );
      

      dRight += dRightStep;
    } // loop over beta

    dLeft += dLeftStep;

  } // loop over alpha

}
