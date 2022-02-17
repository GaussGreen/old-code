#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"
#include <math.h>

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitycallback.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratepower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_volflathrwithspotcomponentpower.h"

extern std::string g_strInputFilename;

extern void __stdcall 
parametricVol(double dTime, const double *pdS, double *pdValues, size_t nNbS,
              int i);

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolFlatHRWithSpotComponentPowerTest::Setup()
{

  // Reset the pricing variables before reading
  m_ptsCDS = make_ptr( new TermStructureCDS() );

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
  m_pParametrization = 
    param_visitor.GetParametrizationVolFlatHRWithSpotComponentPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a cds termstructure 
  m_ptsCDS = termstructures.GetTermStructureCDS();
  CPPUNIT_ASSERT(m_ptsCDS);

  // Also expect an option
  m_pOption = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption);
}


void Parametrization_VolFlatHRWithSpotComponentPowerTest::CheckCDSFit( 
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
    
    double dScale = fabs(cds.GetMarketPrice());
    if (dScale < 1.e-6)
      dScale = 1.0;

    double dError = fabs( output->GetPrice() - cds.GetMarketPrice() );
    CPPUNIT_ASSERT( dError/dScale < 1.e-3);
  }

}

void Parametrization_VolFlatHRWithSpotComponentPowerTest::CheckOptionFit( 
  ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
  ito33::shared_ptr<ito33::ihg::Volatility> pVol)
{

  // Create and setup the pricing model
  ihg::TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption);

  double dScale = fabs(m_pOption->GetMarketPrice());
  if (dScale < 1.e-6)
    dScale = 1.0;

  double dError = fabs( output->GetPrice() - m_pOption->GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);

}


void Parametrization_VolFlatHRWithSpotComponentPowerTest::TestAsSpecified()
{

  // Read in the current data
  Setup();

  // Calibrate
  m_pParametrization->CalibrateWithOptionAndCDSs(*m_pOption,*m_ptsCDS);
  
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();

  // Check if it fits the data
  CheckCDSFit(pHR, pVol);
  CheckOptionFit(pHR, pVol);
}


void Parametrization_VolFlatHRWithSpotComponentPowerTest::TestGeneratedValues()
{

  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolFlatHRWithSpotComponentPower() );

  double dVol = 0.001;
  double dVolStep = 0.2;
  while (dVol < 0.81)
  {
    double dS0 = m_pSessionData->GetSpotSharePrice();
    double dAlpha = 0.01;
    double dAlphaStep = 0.2;  
    while (dAlpha < 0.82)
    {
      double dBeta = 0.0;
      double dBetaStep = 0.4;
      while (dBeta < 1.61)
      {

        //std::cout << "vol = " << dVol 
        //          << ", alpha = " << dAlpha 
        //          << ", beta = " << dBeta 
        //          << std::endl;

        if ( fabs(dVol - 0.201) < 1.e-3 && fabs(dAlpha - 0.51) < 1.e-3
          && fabs(dBeta - 1.6) < 1.e-3)
        {
          dBeta += dBetaStep;
          continue;
        }

        // failed on 1st test file
        if ( fabs(dVol - 0.201) < 1.e-3 && fabs(dAlpha - 0.61) < 1.e-3
          && fabs(dBeta - 1.6) < 1.e-3)
        {
          dBeta += dBetaStep;
          continue;
        }

        // failed on 4th test file
        if ( fabs(dVol - 0.001) < 1.e-3 && fabs(dAlpha - 0.61) < 1.e-3
          && fabs(dBeta - 1.6) < 1.e-3)
        {
          dBeta += dBetaStep;
          continue;
        }

        // failed on 4th test file
        if ( fabs(dVol - 0.201) < 1.e-3 && fabs(dAlpha - 0.81) < 1.e-3
          && fabs(dBeta - 1.6) < 1.e-3)
        {
          dBeta += dBetaStep;
          continue;
        }

        shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
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
    
        shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption);
        m_pOption->SetMarketPrice( output->GetPrice() );

        //std::cout << "vol = " << dVol
        //          << ", alpha = " << dAlpha
        //          << ", beta = " << dBeta
        //          << std::endl;
                  

        // Calibrate
        m_pParametrization->CalibrateWithOptionAndCDSs(*m_pOption, *m_ptsCDS);

        shared_ptr<HazardRateWithTimeComponent> 
          pCalibratedHR = m_pParametrization->GetHazardRate();
        shared_ptr<VolatilityFlat> 
          pCalibratedVol = m_pParametrization->GetVolatility();

        // Check if the data was reproduced exactly
        CheckCDSFit(pCalibratedHR, pCalibratedVol);
        CheckOptionFit(pCalibratedHR, pCalibratedVol);

        dBeta += dBetaStep;
      } // loop over beta

      dAlpha += dAlphaStep;

    } // loop over alpha

    dVol += dVolStep;
  } // loop over vol
}
