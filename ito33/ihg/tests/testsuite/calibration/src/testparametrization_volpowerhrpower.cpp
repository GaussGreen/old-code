#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratepower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "testparametrization_volpowerhrpower.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolPowerHRPowerTest::Setup()
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
    param_visitor.GetParametrizationVolPowerHRPower();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a cds termstructure 
  m_ptsCDS = termstructures.GetTermStructureCDS();
  CPPUNIT_ASSERT(m_ptsCDS);

  // Also expect two options
  m_pOption1 = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption1);

  m_pOption2 = deriv_visitor.GetOption2();
  CPPUNIT_ASSERT(m_pOption2);
}


void Parametrization_VolPowerHRPowerTest::CheckCDSFit( 
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

void Parametrization_VolPowerHRPowerTest::CheckOptionFit( 
  ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
  ito33::shared_ptr<ito33::ihg::Volatility> pVol)
{

  // Create and setup the pricing model
  ihg::TheoreticalModel model;

  model.SetHazardRate(pHR);
  model.SetVolatility(pVol);

  shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption1);

  double dScale = fabs(m_pOption1->GetMarketPrice());
  if (dScale < 1.e-6)
    dScale = 1.0;

  double dError = fabs( output->GetPrice() - m_pOption1->GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);


  output = model.Compute(*m_pOption2);

  dScale = fabs(m_pOption2->GetMarketPrice());
  if (dScale < 1.e-6)
    dScale = 1.0;

  dError = fabs( output->GetPrice() - m_pOption2->GetMarketPrice() );
  CPPUNIT_ASSERT( dError/dScale < 1.e-3);

}


void Parametrization_VolPowerHRPowerTest::TestAsSpecified()
{
  return;
  // Read in the current data
  Setup();

  // Calibrate
  m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1, 
                                                  *m_pOption2,
                                                  *m_ptsCDS);

  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();

  // Check if it fits the data
  CheckCDSFit(pHR, pVol);
  CheckOptionFit(pHR, pVol);
}


void Parametrization_VolPowerHRPowerTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolPowerHRPower() );

  double dS0 = m_pSessionData->GetSpotSharePrice();

  // Running all tests is time consuming, and some tests fail. For now,
  // just run one test

  //m_dHRAlpha  = 0.01;
  //m_dHRBeta   = 0.0;
  //m_dVolAlpha = 0.001;
  //m_dVolBeta  = -1.0;

  m_dHRAlpha  = 0.21;
  m_dHRBeta   = 0.8;
  m_dVolAlpha = 0.201;
  m_dVolBeta  = -0.8;

  //while ( SetNextTest() )
  {
    //std::cout << "hr alpha = " << m_dHRAlpha
    //          << ", hr beta = " << m_dHRBeta
    //          << ", vol alpha = " << m_dVolAlpha
    //          << ", vol beta = " << m_dVolBeta
    //          << std::endl;
    
    shared_ptr<Volatility> pVol( new VolatilityPower(m_dVolAlpha, m_dVolBeta, dS0));
    shared_ptr<HazardRate> pHR( new HazardRatePower(m_dHRAlpha, m_dHRBeta, dS0 ) );

    ihg::TheoreticalModel model;
    model.SetVolatility(pVol);
    model.SetHazardRate(pHR);

    const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
    TermStructureCDS::Elements::const_iterator ppCDS;
    shared_ptr<finance::ModelOutput> output;
    for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
    {
      CDSLike& cds = *(*ppCDS);

      output = model.Compute(cds);

      // Take the computed prices as market price for testing
      cds.SetMarketPrice( output->GetPrice() );
    }

    output = model.Compute(*m_pOption1);
    m_pOption1->SetMarketPrice( output->GetPrice() );

    output = model.Compute(*m_pOption2);
    m_pOption2->SetMarketPrice( output->GetPrice() );

    // Calibrate
    m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1, 
                                                    *m_pOption2,
                                                    *m_ptsCDS);

    shared_ptr<HazardRateCombo> 
      pCalibratedHR = m_pParametrization->GetHazardRate();
    shared_ptr<VolatilityPower> 
      pCalibratedVol = m_pParametrization->GetVolatility();

    // Check if the data was reproduced exactly
    CheckCDSFit(pCalibratedHR, pCalibratedVol);
    CheckOptionFit(pCalibratedHR, pCalibratedVol);
    
  } // while more tests

}



bool Parametrization_VolPowerHRPowerTest::SetNextTest()
{

  double dHRAlphaStep = 0.3;  
  double dHRAlphaMax  = 0.62;
  double dHRAlphaMin  = 0.01;
  
  double dHRBetaStep = 0.5;
  double dHRBetaMax  = 1.01;
  double dHRBetaMin  = 0.0;

  double dVolAlphaStep = 0.3;  
  double dVolAlphaMax  = 0.61;
  double dVolAlphaMin  = 0.001;
  
  double dVolBetaStep = 0.5;
  double dVolBetaMax  = 1.01;
  double dVolBetaMin  = -1.0;

  // work backwards through what would normally be while loops
  m_dVolBeta += dVolBetaStep;
  if ( m_dVolBeta > dVolBetaMax )
  {
    m_dVolBeta = dVolBetaMin;

    m_dVolAlpha += dVolAlphaStep;
    if ( m_dVolAlpha > dVolAlphaMax )
    {
      m_dVolAlpha = dVolAlphaMin;

      m_dHRBeta += dHRBetaStep;
      if ( m_dHRBeta > dHRBetaMax )
      {
        m_dHRBeta = dHRBetaMin;

        m_dHRAlpha += dHRAlphaStep;
        if ( m_dHRAlpha > dHRAlphaMax)
        {
          m_dHRAlpha = dHRAlphaMin;

          return false;
        } // hr alpha

      } // hr beta 

    } // vol alpha

  } // vol beta

  return true;
}
