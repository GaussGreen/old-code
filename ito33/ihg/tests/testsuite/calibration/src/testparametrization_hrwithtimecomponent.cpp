#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"
#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitycallback.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"
#include "testparametrization_hrwithtimecomponent.h"

extern std::string g_strInputFilename;

extern void __stdcall 
parametricVol(double dTime, const double *pdS, double *pdValues, size_t nNbS,
              int i);

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_HRWithTimeComponentTest::Setup()
{

  // Reset the pricing variables before reading
  m_ptsCDS = make_ptr( new TermStructureCDS() );
  m_ptsEDS = make_ptr( new TermStructureEDS() );

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
  m_pParametrization = param_visitor.GetParametrizationHRWithTimeComponent();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a cds, eds, or par bond termstructure
  m_ptsParBond = termstructures.GetTermStructureParBond();
  m_ptsCDS = termstructures.GetTermStructureCDS();
  m_ptsEDS = termstructures.GetTermStructureEDS();
  CPPUNIT_ASSERT(m_ptsParBond || m_ptsCDS || m_ptsEDS);
}


void Parametrization_HRWithTimeComponentTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Vol and calibrated hr
  shared_ptr<HazardRateWithTimeComponent> pHR;
  
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();
  if ( !pVol )
    pVol = make_ptr( new VolatilityFlat(0.1) ); 

  // Calibrate and check fit
  if (m_ptsParBond)
  {
    pHR = m_pParametrization->CalibrateWithParBonds(*m_ptsParBond);
    CheckFit(pHR, pVol, *m_ptsParBond);
  }

  if (m_ptsCDS)
  {
    pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);
    CheckFit(pHR, pVol, *m_ptsCDS);
  }

  if (m_ptsEDS)
  {
    pHR = m_pParametrization->CalibrateWithEDSs(*m_ptsEDS);
    CheckFit(pHR, pVol, *m_ptsEDS);
  }
}


void Parametrization_HRWithTimeComponentTest::TestVols()
{  

  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationHRWithTimeComponent() );

  // Try some flat vols
  size_t nIdx;
  double dVol = 0.0;
  double dStepVol = 0.2;
  double dMaxVol = 0.8;

  // EDS prices depend on vol. Can't go too high since we are comparing
  // to prices from the xml file, and these may have been generated
  // with small vol values
  if (m_ptsEDS)
  {
    dMaxVol = 0.4;
    dStepVol = 0.1;
  }

  for (nIdx = 0; dVol <= dMaxVol; dVol += dStepVol, nIdx++)
  {
    shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
    m_pParametrization->SetVolatility(pVol);

    // Calibrate and check the fit
    shared_ptr<HazardRateWithTimeComponent> pHR;

    if (m_ptsParBond)
    {
      pHR = m_pParametrization->CalibrateWithParBonds(*m_ptsParBond);
      CheckFit(pHR, pVol, *m_ptsParBond);
    }

    if (m_ptsCDS)
    {
      pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);
      CheckFit(pHR, pVol, *m_ptsCDS);
    }

    if (m_ptsEDS)
    {
      pHR = m_pParametrization->CalibrateWithEDSs(*m_ptsEDS);
      CheckFit(pHR, pVol, *m_ptsEDS);
    }
  }

  // Test a non-flat vol
  Date valuationDate = m_pSessionData->GetValuationDate();
  int iDate = valuationDate.GetExcel();
    
  shared_ptr<Volatility> pVol(new VolatilityCallBack(parametricVol,iDate,0.0));
  m_pParametrization->SetVolatility(pVol);

  // Calibrate and check the fit
  shared_ptr<HazardRateWithTimeComponent> pHR;

    
  if (m_ptsParBond)
  {
    pHR = m_pParametrization->CalibrateWithParBonds(*m_ptsParBond);
    CheckFit(pHR, pVol, *m_ptsParBond);
  }

  if (m_ptsCDS)
  {
    pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);
    CheckFit(pHR, pVol, *m_ptsCDS);
  }

  if (m_ptsEDS)
  {
    pHR = m_pParametrization->CalibrateWithEDSs(*m_ptsEDS);
    CheckFit(pHR, pVol, *m_ptsEDS);
  }
}


void Parametrization_HRWithTimeComponentTest::TestSpotComponents()
{
  // Read in the current data
  Setup();

  double dS0 = m_pSessionData->GetSpotSharePrice();

  double dBeta = 1.0;
  double dBetaStep = 0.2;
  double dBetaMax = 1.81;
  if( m_ptsParBond)
    dBetaMax = 1.21;

  while (dBeta < dBetaMax)
  {
    // Make a new parametrization
    m_pParametrization = make_ptr( new ParametrizationHRWithTimeComponent() );

    // Set a spot component
    shared_ptr<SpotComponent> 
      pSpotComponent(new HRSpotComponentPower(dBeta, dS0));

    m_pParametrization->SetSpotComponent(pSpotComponent);

    // Set a vol
    shared_ptr<Volatility> pVol( new VolatilityFlat(0.2) );
    m_pParametrization->SetVolatility(pVol);

    // Calibrate and check the fit
    shared_ptr<HazardRateWithTimeComponent> pHR;

    if (m_ptsParBond)
    {
      pHR = m_pParametrization->CalibrateWithParBonds(*m_ptsParBond);
      CheckFit(pHR, pVol, *m_ptsParBond);
    }

    if (m_ptsCDS)
    {
      pHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);
      CheckFit(pHR, pVol, *m_ptsCDS);
    }

    if (m_ptsEDS)
    {
      pHR = m_pParametrization->CalibrateWithEDSs(*m_ptsEDS);
      CheckFit(pHR, pVol, *m_ptsEDS);
    }

    dBeta += dBetaStep;
  } // loop over beta
}


void Parametrization_HRWithTimeComponentTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationHRWithTimeComponent() );

  TermStructureDerivative tsDeriv;
  if (m_ptsCDS)
    tsDeriv = *m_ptsCDS;
  else if (m_ptsEDS)
    tsDeriv = *m_ptsEDS;


  if (m_ptsCDS || m_ptsEDS)
  {
    // Generate prices using timeonly hazard rates.  Should be able to
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
      //const CDS& cds = *(*ppCDS);

      pDates[nIdx] = (*ppElem)->GetMaturityDate();
      pdValues[nIdx] = 0.1 + nIdx*0.04;
    }

    // Setup the model and generate prices
    shared_ptr<Volatility> pVol( new VolatilityFlat(0.2) );
      m_pParametrization->SetVolatility(pVol);

    shared_ptr<HazardRate> pHR( new HazardRateTimeOnly(pDates, pdValues));

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
    m_pParametrization->SetVolatility(pVol);

    shared_ptr<HazardRateWithTimeComponent> pCalibratedHR;
    
    if (m_ptsParBond)
    {
      pHR = m_pParametrization->CalibrateWithParBonds(*m_ptsParBond);
      CheckFit(pHR, pVol, *m_ptsParBond);
    }

    if (m_ptsCDS)
    {
      pCalibratedHR = m_pParametrization->CalibrateWithCDSs(*m_ptsCDS);
      CheckFit(pCalibratedHR, pVol, *m_ptsCDS);
    }
    else
    {
      pCalibratedHR = m_pParametrization->CalibrateWithEDSs(*m_ptsEDS);
      CheckFit(pCalibratedHR, pVol, *m_ptsEDS);
    }


    // Check that the params were reproduced exactly
    std::vector<double> pdCalValues = pCalibratedHR->GetTimeComponentValues();
    CPPUNIT_ASSERT( nNbElements == pdCalValues.size() );

    for (nIdx = 0; nIdx < nNbElements; nIdx++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pdValues[nIdx], pdCalValues[nIdx], 1.e-3);   
  }
}
