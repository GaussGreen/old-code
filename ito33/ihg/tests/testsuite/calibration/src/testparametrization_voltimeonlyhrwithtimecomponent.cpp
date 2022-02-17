#include "ito33/beforestd.h"
#include <string>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_voltimeonlyhrwithtimecomponent.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolTimeOnlyHRWithTimeComponentTest::Setup()
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
    param_visitor.GetParametrizationVolTimeOnlyHRWithTimeComponent();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a general list of derivatives
  CPPUNIT_ASSERT(m_pDerivatives->GetAll().size() > 0);
}


void Parametrization_VolTimeOnlyHRWithTimeComponentTest::TestAsSpecified()
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


void Parametrization_VolTimeOnlyHRWithTimeComponentTest::TestSpotComponents()
{
  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolTimeOnlyHRWithTimeComponent() );

  double dS0 = m_pSessionData->GetSpotSharePrice();

  double dBeta = 0.0;
  double dBetaStep = 0.4;
  double dBetaMax = 1.61;

  while (dBeta < dBetaMax)
  {
    // Set a spot component
    shared_ptr<SpotComponent> pSpotComponent(new 
      HRSpotComponentPower(dBeta, dS0));

    m_pParametrization->SetSpotComponent(pSpotComponent);

    // Calibrate and check the fit
    m_pParametrization->CalibrateWithDerivatives(*m_pDerivatives);
    // m_pParametrization->Calibrate(*m_pDerivatives , *m_pDerivatives);

    CheckFit(m_pParametrization->GetHazardRate(), 
              m_pParametrization->GetVolatility(), 
              *m_pDerivatives);

    dBeta += dBetaStep;
  } // loop over beta
}


void Parametrization_VolTimeOnlyHRWithTimeComponentTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolTimeOnlyHRWithTimeComponent() );

  // Generate prices using power volatility and hazard rate. Should still
  // be able to fit these prices exactly
  double dSpot = m_pSessionData->GetSpotSharePrice();
  shared_ptr<ihg::Volatility> 
    pVolatility( new ihg::VolatilityPower(0.45, 0.25, dSpot) );

  shared_ptr<ihg::HazardRate> 
    pHazardRate( new ihg::HazardRatePower(0.2, 0.15, dSpot) );

  ihg::TheoreticalModel model;
  model.SetVolatility( pVolatility );
  model.SetHazardRate( pHazardRate );

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

}
