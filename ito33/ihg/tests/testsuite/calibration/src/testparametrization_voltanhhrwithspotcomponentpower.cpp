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

#include "ito33/ihg/parametrization_volflat_hrpower.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_voltanhhrwithspotcomponentpower.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolTanhHRWithSpotComponentPowerTest::Setup()
{

  // Reset the pricing variables before reading
  m_ptsCDS = make_ptr( new TermStructureCDS() );
  
  // Read from the currently active global input file
  ihg::XML::CalibrationReader reader(g_strInputFilename.c_str());

  // Get the session
  m_pSessionData = reader.ReadSessionData();

  // Get the parametrization
  ParametrizationVisitorGoodType param_visitor;
  DerivativeVisitorGoodType deriv_visitor;
  TermStructureEnumerator termstructures;
  shared_ptr<Derivatives> pDerivatives(new Derivatives());
  reader.ReadCalibration(param_visitor, termstructures, 
                         deriv_visitor, *pDerivatives);
  m_pParametrization = 
    param_visitor.GetParametrizationVolTanhHRWithSpotComponentPower();
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


void Parametrization_VolTanhHRWithSpotComponentPowerTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1,
                                                  *m_pOption2,
                                                  *m_ptsCDS);
  
  // Check the fit
  shared_ptr<HazardRate> pHR = m_pParametrization->GetHazardRate();
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();

  CheckFit(pHR, pVol, *m_ptsCDS);

  CheckFit(pHR, pVol, *m_pOption1);

  CheckFit(pHR, pVol, *m_pOption2);

}


void Parametrization_VolTanhHRWithSpotComponentPowerTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Generate prices using flat volatility and power hazard rates. Should
  // be able to fit these prices exactly
  double dSpot = m_pSessionData->GetSpotSharePrice();

  for (double dLeft = 0.1; dLeft < 0.51; dLeft += 0.4)
  {
  for (double dRight = 0.2; dRight < 0.41; dRight += 0.2)
  {
  for (double dScale = 1.0; dScale < 12.0; dScale += 10.0)
  {

    shared_ptr<Volatility> 
      pVolatility( new VolatilityTanh(dLeft, dRight, dScale, dSpot) );

    for (double dAlpha = 0.01; dAlpha < 0.82; dAlpha += 0.4)
    {

      for (double dBeta = 0.0; dBeta < 0.81; dBeta += 0.8)
      {
        //std::cout << "left = " << dLeft << ", right = " << dRight 
        //          << ", dScale = " << dScale << ", alpha = " << dAlpha
        //          << ", beta = " << dBeta 
        //          << std::endl;

        shared_ptr<HazardRate> 
          pHazardRate( new HazardRatePower(dAlpha, dBeta, dSpot) );

        ihg::TheoreticalModel model;
        model.SetVolatility( pVolatility );
        model.SetHazardRate( pHazardRate );

        // Reset the market prices of the cds termstructure
        const TermStructureCDS::Elements& pCDSElements = m_ptsCDS->GetAll();
        TermStructureCDS::Elements::const_iterator ppCDS;
        for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
        {
          CDSLike& cds = *(*ppCDS);

          shared_ptr<finance::ModelOutput> output = model.Compute(cds);

          // Take the computed prices as market price for testing
          cds.SetMarketPrice( output->GetPrice() );
        }
    
        // Reset the market prices of the options
        shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption1);
        m_pOption1->SetMarketPrice( output->GetPrice() );

        output = model.Compute(*m_pOption2);
        m_pOption2->SetMarketPrice( output->GetPrice() );

        // Calibrate
        m_pParametrization->CalibrateWithOptionsAndCDSs(*m_pOption1,
                                                        *m_pOption2,
                                                        *m_ptsCDS);

        CheckFit(m_pParametrization->GetHazardRate(), 
                 m_pParametrization->GetVolatility(), 
                 *m_ptsCDS);

        CheckFit(m_pParametrization->GetHazardRate(), 
                 m_pParametrization->GetVolatility(), 
                 *m_pOption1);

        CheckFit(m_pParametrization->GetHazardRate(), 
                 m_pParametrization->GetVolatility(), 
                 *m_pOption2);

      } // loop over beta

    } // loop over alpha

  } // loop over scale
  } // loop over right
  } // loop over left

} // TestGenerated
