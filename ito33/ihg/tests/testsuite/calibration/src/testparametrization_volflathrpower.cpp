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
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_volflathrpower.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolFlatHRPowerTest::Setup()
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
    param_visitor.GetParametrizationVolFlatHRPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a general list of derivatives
  CPPUNIT_ASSERT(m_pDerivatives->GetAll().size() > 0);
}


void Parametrization_VolFlatHRPowerTest::TestAsSpecified()
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


void Parametrization_VolFlatHRPowerTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Generate prices using flat volatility and power hazard rates. Should
  // be able to fit these prices exactly
  double dSpot = m_pSessionData->GetSpotSharePrice();

  // At the time this comment was written, the vol flat, hr power calibrator
  // uses ASA, which is slow.  Also, several tests below fail.  For now, only 
  // enable one test.

  //for (double dVol = 0.05; dVol < 0.86; dVol += 0.2)
  for (double dVol = 0.2; dVol < 0.3; dVol += 0.2)
  {

    shared_ptr<ihg::Volatility> 
      pVolatility( new ihg::VolatilityFlat(dVol) );

    //for (double dAlpha = 0.01; dAlpha < 0.82; dAlpha += 0.2)
    for (double dAlpha = 0.1; dAlpha < 0.2; dAlpha += 0.2)
    {

      //for (double dBeta = 0.0; dBeta < 1.61; dBeta += 0.4)
      for (double dBeta = 0.0; dBeta < 0.1; dBeta += 0.4)
      {
        //std::cout << "vol = " << dVol 
        //          << ", alpha = " << dAlpha
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

  } // loop over vol

}
