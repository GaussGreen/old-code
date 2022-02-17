#include "ito33/beforestd.h"
#include <string>
#include <math.h>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/parametrization_volflat_hrpower.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_voltanhhrpower.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolTanhHRPowerTest::Setup()
{

  // Reset the pricing variables before reading
  m_pDerivatives = shared_ptr<Derivatives>( new Derivatives() );
  
  // Read from the currently active global input file
  ihg::XML::CalibrationReader reader(g_strInputFilename.c_str());

  // Get the session
  m_pSessionData = reader.ReadSessionData();

  // Get the parametrization
  ito33::ihg::ParametrizationVisitorGoodType param_visitor;
  ito33::finance::DerivativeVisitorGoodType deriv_visitor;
  ito33::finance::TermStructureEnumerator termstructures;
  reader.ReadCalibration(param_visitor, termstructures, 
                         deriv_visitor, *m_pDerivatives);
  m_pParametrization = 
    param_visitor.GetParametrizationVolTanhHRPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a general list of derivatives
  CPPUNIT_ASSERT(m_pDerivatives->GetAll().size() > 0);
}


void Parametrization_VolTanhHRPowerTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  m_pParametrization->CalibrateWithDerivatives(*m_pDerivatives);

  // Check the fit
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();
  CheckFit(pHR, pVol, *m_pDerivatives);
  
}


void Parametrization_VolTanhHRPowerTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // At the time this comment was written, the vol tanh, hr power calibrator
  // uses the general calibrator with ASA, which is slow.  For now, only 
  // enable one test.

  // Generate prices using flat volatility and power hazard rates. Should
  // be able to fit these prices exactly
  double dSpot = m_pSessionData->GetSpotSharePrice();

  //for (double dLeft = 0.1; dLeft < 0.91; dLeft += 0.4)
  for (double dLeft = 0.35; dLeft < 0.41; dLeft += 0.4)
  {
  //for (double dRight = 0.1; dRight < 0.91; dRight += 0.4)
  for (double dRight = 0.25; dRight < 0.31; dRight += 0.4)
  {
  //for (double dScale = 0.0; dScale < 21.0; dScale += 10.0)
  for (double dScale = 2.0; dScale < 3.0; dScale += 10.0)
  {

    shared_ptr<ihg::Volatility> 
      pVolatility( new ihg::VolatilityTanh(dLeft, dRight, dScale, dSpot) );

    //for (double dAlpha = 0.01; dAlpha < 0.82; dAlpha += 0.4)
    for (double dAlpha = 0.21; dAlpha < 0.22; dAlpha += 0.4)
    {

      //for (double dBeta = 0.0; dBeta < 1.61; dBeta += 0.8)
      for (double dBeta = 0.2; dBeta < 0.21; dBeta += 0.8)
      {
        //std::cout << "left = " << dLeft << ", right = " << dRight 
        //          << ", dScale = " << dScale << ", alpha = " << dAlpha
        //          << ", beta = " << dBeta 
        //          << std::endl;

        shared_ptr<ihg::HazardRate> 
          pHazardRate( new ihg::HazardRatePower(dAlpha, dBeta, dSpot) );

        ihg::TheoreticalModel model;
        model.SetVolatility( pVolatility );
        model.SetHazardRate( pHazardRate );

        // Reset the market prices in the derivative list
        Derivatives::Elements::const_iterator iter;
        const Derivatives::Elements& elements( m_pDerivatives->GetAll() );
        for (iter = elements.begin(); iter != elements.end(); ++iter)
        {
          Derivative& deriv = *iter->first;

          shared_ptr<finance::ModelOutput> output = model.Compute(deriv);

          // Take the computed prices as market price for testing
          deriv.SetMarketPrice( output->GetPrice() );
        }
 
        // Calibrate
        m_pParametrization->CalibrateWithDerivatives(*m_pDerivatives);

        CheckFit(m_pParametrization->GetHazardRate(), 
                 m_pParametrization->GetVolatility(), 
                 *m_pDerivatives);

      } // loop over beta

    } // loop over alpha

  } // loop over scale
  } // loop over right
  } // loop over left

}
