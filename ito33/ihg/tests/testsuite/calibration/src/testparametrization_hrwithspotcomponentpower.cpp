#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"
#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitycallback.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratepower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_hrwithspotcomponentpower.h"


extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;

extern void __stdcall 
parametricVol(double dTime, const double *pdS, double *pdValues, size_t nNbS,
              int i);


void Parametrization_HRWithSpotComponentPowerTest::Setup()
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
  m_pParametrization = param_visitor.GetParametrizationHRWithSpotComponentPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a cds termstructure 
  m_ptsCDS = termstructures.GetTermStructureCDS();
  CPPUNIT_ASSERT(m_ptsCDS);
}


void Parametrization_HRWithSpotComponentPowerTest::CheckCDSFit( 
  ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
  ito33::shared_ptr<ito33::ihg::Volatility> pVol)
{

  // Create and setup the pricing model
  ihg::TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  // compare the prices
  const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    //std::cout << "  Analytic (generated) price " << " = " 
    //          << cds.GetMarketPrice() << std::endl;
    
    //std::cout << "              Computed price " << " = " 
    //          << output->GetPrice() << std::endl;   

    double dScale = fabs( cds.GetMarketPrice() );
    if (dScale < 1.e-6)
      dScale = 1.0;

    double dError = fabs( output->GetPrice() - cds.GetMarketPrice() );
    CPPUNIT_ASSERT( dError/dScale < 1.e-3);

  }

}


void Parametrization_HRWithSpotComponentPowerTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  shared_ptr<HazardRateWithTimeComponent>
    pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);

  // Check if it fits the data
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();
  
  CheckCDSFit(pHR, pVol);
}


void Parametrization_HRWithSpotComponentPowerTest::TestVols()
{
  // Read in the current data
  Setup();

  // Try some flat vols
  size_t nIdx;
  double dVol = 0.0;
  double dStepVol = 0.2;
  for (nIdx = 0; dVol <= 0.8; dVol += dStepVol, nIdx++)
  {
    shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );

    // Make a new parametrization
    m_pParametrization = make_ptr( new ParametrizationHRWithSpotComponentPower(pVol) );

    // Calibrate
    shared_ptr<HazardRateWithTimeComponent>
      pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);

    // Check if it fits the data
    CheckCDSFit(pHR, pVol);
  }

  // Test a non-flat vol
  Date valuationDate = m_pSessionData->GetValuationDate();
  int iDate = valuationDate.GetExcel();
    
  shared_ptr<Volatility> pVol(new VolatilityCallBack(parametricVol, iDate, 0.0));
  m_pParametrization = make_ptr( new ParametrizationHRWithSpotComponentPower(pVol) );

  // Calibrate
  shared_ptr<HazardRateWithTimeComponent>
    pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);

  // Check if it fits the data
  CheckCDSFit(pHR, pVol);

}


void Parametrization_HRWithSpotComponentPowerTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Generate prices using spot component power hazard rates.  Should be 
  // able to fit these prices exactly. However, local mins may be found,
  // so the parameters are not always the same
  double dS0 = m_pSessionData->GetSpotSharePrice();
  double dAlpha = 0.1;
  double dAlphaStep = 0.1;  
  while (dAlpha < 0.82)
  {
    double dBeta = 0.0;
    double dBetaStep = 0.2;

    while (dBeta < 1.81)
    {
      shared_ptr<Volatility> pVol( new VolatilityFlat(0.2) );
      shared_ptr<HazardRate> pHR( new HazardRatePower(dAlpha, dBeta, dS0 ) );

      ihg::TheoreticalModel model;
      model.SetVolatility(pVol);
      model.SetHazardRate(pHR);

      const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
      TermStructureCDS::Elements::const_iterator ppCDS;
      for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
      {
        CDSLike& cds = *(*ppCDS);

        shared_ptr<finance::ModelOutput> output = model.Compute(cds);

        // Take the computed prices as market price for testing
        cds.SetMarketPrice( output->GetPrice() );
      }
  
      // Calibrate
      m_pParametrization = make_ptr( new ParametrizationHRWithSpotComponentPower(pVol) );

      shared_ptr<HazardRateWithTimeComponent>
        pCalibratedHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);

      // Check if the data was reproduced exactly
      CheckCDSFit(pCalibratedHR, pVol);

      // Check that the params were reproduced exactly
      // Local mins make this difficult/impossible
      /*
      shared_ptr<HRSpotComponentPower> pCalibratedPower = 
        m_pParametrization->GetHRSpotComponentPower();

      double dCalibratedAlpha = pCalibratedPower->GetAlpha();
      double dCalibratedBeta = pCalibratedPower->GetBeta();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dAlpha, dCalibratedAlpha, 1.e-3 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( dBeta, dCalibratedBeta, 1.e-3 );

      // The time components should all be 1.0
      std::vector<double> pdCalValues = pCalibratedHR->GetTimeComponentValues();
      size_t nNbElements = pCDSElements.size();
      CPPUNIT_ASSERT( nNbElements == pdCalValues.size() );

      size_t nIdx;
      for (nIdx = 0; nIdx < nNbElements; nIdx++)
      {
        //std::cout << nIdx 
        //          << ", " << pdValues[nIdx] 
        //          << ", " << pdCalValues[nIdx] 
        //          << std::endl;
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, pdCalValues[nIdx], 1.e-4);    
      }
      */

      dBeta += dBetaStep;
    } // loop over beta

    dAlpha += dAlphaStep;

  } // loop over alpha

}
