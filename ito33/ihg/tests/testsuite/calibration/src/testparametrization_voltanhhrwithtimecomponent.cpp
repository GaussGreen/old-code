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

#include "ito33/ihg/parametrization_voltanh_hrwithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_voltanhhrwithtimecomponent.h"

extern std::string g_strInputFilename;

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolTanhHRWithTimeComponentTest::Setup()
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
  reader.ReadCalibration(param_visitor, termstructures, 
                         deriv_visitor, *pDerivatives);
  m_pParametrization = 
    param_visitor.GetParametrizationVolTanhHRWithTimeComponent();
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


void Parametrization_VolTanhHRWithTimeComponentTest::TestAsSpecified()
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


void Parametrization_VolTanhHRWithTimeComponentTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Generate prices using timeonly hazard rates and tanh vols. 
  // Should be able to fit these prices exactly. However, local 
  // mins may be found, or the tolerance may be reached, so the 
  // parameters are not always the same
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

  double dSpot = m_pSessionData->GetSpotSharePrice();

  for (double dLeft = 0.1; dLeft < 0.51; dLeft += 0.4)
  {
    for (double dRight = 0.2; dRight < 0.41; dRight += 0.2)
    {
      for (double dScale = 1.0; dScale < 12.0; dScale += 10.0)
      {

        shared_ptr<ihg::Volatility> 
          pVolatility( new ihg::VolatilityTanh(dLeft, dRight, dScale, dSpot) );

        shared_ptr<HazardRate> pHR( new HazardRateTimeOnly(pDates, pdValues));

        //std::cout << "left = " << dLeft 
        //          << ", right = " << dRight 
        //          << ", dScale = " << dScale
        //          << std::endl;

        ihg::TheoreticalModel model;
        model.SetVolatility( pVolatility );
        model.SetHazardRate( pHR );

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

      } // loop over scale
    } // loop over right
  } // loop over left

} // TestGenerated
