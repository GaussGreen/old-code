#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"
#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ito33/ihg/parametrization_volwithtimecomponent.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilitytimeonly.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"
#include "testparametrization_volwithtimecomponent.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolWithTimeComponentTest::Setup()
{

  // Reset the pricing variables before reading
  m_ptsOption = make_ptr( new TermStructureOption() );

  // Read from the currently active global input file
  ihg::XML::CalibrationReader reader(g_strInputFilename.c_str());

  // Get the session
  m_pSessionData = reader.ReadSessionData();

  // Get the parametrization
  ito33::ihg::ParametrizationVisitorGoodType param_visitor;
  ito33::finance::DerivativeVisitorGoodType deriv_visitor;
  ito33::finance::TermStructureEnumerator termstructures;
  shared_ptr<ito33::finance::Derivatives> 
    pDerivatives(new ito33::finance::Derivatives());
  reader.ReadCalibration(param_visitor, termstructures, deriv_visitor,
                         *pDerivatives);
  m_pParametrization = param_visitor.GetParametrizationVolWithTimeComponent();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect an option termstructure
  m_ptsOption = termstructures.GetTermStructureOption();
  CPPUNIT_ASSERT(m_ptsOption);
}


void Parametrization_VolWithTimeComponentTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // HR and calibrated vol
  shared_ptr<VolatilityWithTimeComponent> pVol;
  
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  if ( !pHR )
    pHR = make_ptr( new HazardRateFlat(0.05) ); 

  // Calibrate and check fit
  pVol = m_pParametrization->CalibrateWithOptions(*m_ptsOption);
  CheckFit(pHR, pVol, *m_ptsOption);

}


void Parametrization_VolWithTimeComponentTest::TestHRs()
{  

  // Read in the current data
  Setup();

  // Try some flat hazard rates
  size_t nIdx;
  double dHR = 0.01;
  double dStepHR = 0.05;
  double dMaxHR = 0.17;

  for (nIdx = 0; dHR <= dMaxHR; dHR += dStepHR, nIdx++)
  {
    shared_ptr<HazardRate> pHR( new HazardRateFlat(dHR) );

    // Make a new parametrization
    m_pParametrization = make_ptr( new ParametrizationVolWithTimeComponent(pHR) );

    // Calibrate and check the fit
    shared_ptr<VolatilityWithTimeComponent> pVol;

    pVol = m_pParametrization->CalibrateWithOptions(*m_ptsOption);
    CheckFit(pHR, pVol, *m_ptsOption);
  }

  // Test a non-flat hazardrate
  double dS0 = m_pSessionData->GetSpotSharePrice();
  shared_ptr<HazardRate> pHR(new HazardRatePower(0.1, 0.4, dS0));
    
  m_pParametrization = make_ptr( new ParametrizationVolWithTimeComponent(pHR) );

  // Calibrate and check the fit
  shared_ptr<VolatilityWithTimeComponent> pVol;
    
  pVol = m_pParametrization->CalibrateWithOptions(*m_ptsOption);

  CheckFit(pHR, pVol, *m_ptsOption);

}


void Parametrization_VolWithTimeComponentTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  TermStructureDerivative tsDeriv;
  tsDeriv = *m_ptsOption;

  // Generate prices using timeonly volatilities.  Should be able to
  // fit these prices exactly
  const TermStructureDerivative::Elements& pElements = tsDeriv.GetAll();
  size_t nNbElements = pElements.size();
  std::vector<Date> pDates(nNbElements);
  std::vector<double> pdValues(nNbElements);

  TermStructureDerivative::Elements::const_iterator ppElem;
  size_t nIdx;
  for ( ppElem = pElements.begin(), nIdx = 0; 
        ppElem != pElements.end(); 
        ++ppElem, nIdx++)
  {
    pDates[nIdx] = (*ppElem)->GetMaturityDate();
    pdValues[nIdx] = 0.1 + nIdx*0.04;
  }

  // Setup the model and generate prices
  shared_ptr<Volatility> pVol( new VolatilityTimeOnly(pDates, pdValues));
  shared_ptr<HazardRate> pHR( new HazardRateFlat(0.1));

  ihg::TheoreticalModel model;
  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);

  for ( ppElem = pElements.begin(); 
        ppElem != pElements.end(); 
        ++ppElem)
  {
    Derivative& deriv = *(*ppElem);

    shared_ptr<finance::ModelOutput> output = model.Compute(deriv);

    // Take the computed prices as market price for testing
    deriv.SetMarketPrice( output->GetPrice() );
  }

  // Calibrate
  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolWithTimeComponent(pHR) );  

  shared_ptr<VolatilityWithTimeComponent> pCalibratedVol;
    
  pCalibratedVol = m_pParametrization->CalibrateWithOptions(*m_ptsOption);

  CheckFit(pHR, pVol, *m_ptsOption);

  // Check that the params were reproduced exactly
  std::vector<double> pdCalValues = pCalibratedVol->GetTimeComponentValues();
  CPPUNIT_ASSERT( nNbElements == pdCalValues.size() );

  for (nIdx = 0; nIdx < nNbElements; nIdx++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdValues[nIdx], pdCalValues[nIdx], 1.e-3);   

}
