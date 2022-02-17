#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/parametrization_volpower_hrwithtimecomponent.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_volpowerhrwithtimecomponent.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolPowerHRWithTimeComponentTest::Setup()
{

  // Reset the pricing variables before reading
  m_ptsCDS = make_ptr( new TermStructureCDS() );

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
  m_pParametrization = 
    param_visitor.GetParametrizationVolPowerHRWithTimeComponent();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect two optioms
  m_pOption1 = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption1);

  m_pOption2 = deriv_visitor.GetOption2();
  CPPUNIT_ASSERT(m_pOption2);

  // Expect CDS term structure
  m_ptsCDS = termstructures.GetTermStructureCDS();
  CPPUNIT_ASSERT(m_ptsCDS);
}


void Parametrization_VolPowerHRWithTimeComponentTest::CheckFit( 
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


  const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    double dScale = fabs(cds.GetMarketPrice());
    if (dScale < 1.e-6)
      dScale = 1.0;

    double dError = fabs( output->GetPrice() - cds.GetMarketPrice() );
    CPPUNIT_ASSERT( dError/dScale < 1.e-3);
  }

}


void Parametrization_VolPowerHRWithTimeComponentTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  
  m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1, *m_pOption2, 
                                                  *m_ptsCDS);

  shared_ptr<VolatilityPower> pVol = m_pParametrization->GetVolatility();

  shared_ptr<HazardRateWithTimeComponent> pHR = m_pParametrization->GetHazardRate();

  
  // Check if it fits the data
  CheckFit(pHR, pVol);
}


void Parametrization_VolPowerHRWithTimeComponentTest::TestGeneratedValues()
{
  
  // Read in the current data
  Setup();


  // Generate prices using timeonly hazard rates and power vols. 
  // Should be able to fit these prices 
  // exactly. However, local mins may be found, or the tolerance may be 
  // reached, so the parameters are not always the same
  const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
  size_t nNbElements = pCDSElements.size();
  std::vector<Date> pDates(nNbElements);
  std::vector<double> pdValues(nNbElements);

  TermStructureCDS::Elements::const_iterator ppCDS;
  size_t nIdx;
  for ( ppCDS = pCDSElements.begin(), nIdx = 0; 
        ppCDS != pCDSElements.end(); 
        ++ppCDS, nIdx++)
  {
    const CDSLike& cds = *(*ppCDS);

    pDates[nIdx] = cds.GetMaturityDate();
    pdValues[nIdx] = 0.1 + nIdx*0.04;
  }

  double dS0 = m_pSessionData->GetSpotSharePrice();
  double dAlpha = 0.1;
  double dAlphaStep = 0.1;  
  while (dAlpha < 0.82)
  {
    double dBeta = -1.8;
    double dBetaStep = 0.2;

    while (dBeta < 1.81)
    {
      //std::cout << "alpha = " << dAlpha
      //          << ", beta = " << dBeta
      //          << std::endl;

      shared_ptr<Volatility> pVol( new VolatilityPower(dAlpha, dBeta, dS0) );
      shared_ptr<HazardRate> pHR( new HazardRateTimeOnly(pDates, pdValues));

      ihg::TheoreticalModel model;
      model.SetVolatility(pVol);
      model.SetHazardRate(pHR);

      // Take the computed prices as market price for testing
      shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption1);
      m_pOption1->SetMarketPrice( output->GetPrice() );

      output = model.Compute(*m_pOption2);
      m_pOption2->SetMarketPrice( output->GetPrice() );      
  
      for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
      {
        CDSLike& cds = *(*ppCDS);

        shared_ptr<finance::ModelOutput> output = model.Compute(cds);

        // Take the computed prices as market price for testing
        cds.SetMarketPrice( output->GetPrice() );
      }

      // Calibrate                               
      m_pParametrization = make_ptr( new ParametrizationVolPowerHRWithTimeComponent() );
      
      m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1, *m_pOption2, *m_ptsCDS);
     
      shared_ptr<VolatilityPower> pCalibratedVol = m_pParametrization->GetVolatility();
      shared_ptr<HazardRateWithTimeComponent> pCalibratedHR = m_pParametrization->GetHazardRate();

      // Check if the data was reproduced exactly
      CheckFit(pCalibratedHR, pCalibratedVol);

      dBeta += dBetaStep;
    } // loop over beta

    dAlpha += dAlphaStep;

  } // loop over alpha

} // TestGeneratedValues
