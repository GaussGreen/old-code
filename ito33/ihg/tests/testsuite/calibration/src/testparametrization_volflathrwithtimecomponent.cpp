#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include <math.h>
#include "ito33/afterstd.h"
#include "ito33/sharedptr.h"

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitycallback.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ihg/tests/checkfit.h"

#include "testparametrization_volflathrwithtimecomponent.h"

extern std::string g_strInputFilename;

extern void __stdcall 
parametricVol(double dTime, const double *pdS, double *pdValues, size_t nNbS,
              int i);

using namespace ito33;
using namespace ito33::ihg;
using namespace ito33::finance;


void Parametrization_VolFlatHRWithTimeComponentTest::Setup()
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
    param_visitor.GetParametrizationVolFlatHRWithTimeComponent();
  CPPUNIT_ASSERT(m_pParametrization);

  // Expect a cds termstructure 
  m_ptsCDS = termstructures.GetTermStructureCDS();
  m_ptsParBond = termstructures.GetTermStructureParBond();
  CPPUNIT_ASSERT(m_ptsCDS || m_ptsParBond );

  // Also expect an option
  m_pOption = deriv_visitor.GetOption();
  CPPUNIT_ASSERT(m_pOption);
}


void Parametrization_VolFlatHRWithTimeComponentTest::CheckFitAll( 
  ito33::shared_ptr<ito33::ihg::HazardRate> pHR,
  ito33::shared_ptr<ito33::ihg::Volatility> pVol)
{
  CheckFit(pHR, pVol, *m_pOption);
  if(m_ptsCDS)
    CheckFit(pHR, pVol, *m_ptsCDS);
  else
    CheckFit(pHR, pVol, *m_ptsParBond);
}


void Parametrization_VolFlatHRWithTimeComponentTest::TestAsSpecified()
{
  // Read in the current data
  Setup();

  // Calibrate
  if(m_ptsCDS)
   m_pParametrization->CalibrateWithOptionAndCDSs(*m_pOption,*m_ptsCDS);
  else
   m_pParametrization->CalibrateWithOptionAndParBonds(*m_pOption,*m_ptsParBond);
  
  shared_ptr<HazardRateWithTimeComponent> pHR = m_pParametrization->GetHazardRate();
  shared_ptr<Volatility> pVol = m_pParametrization->GetVolatility();

  // Check if it fits the data using the general calibration
  CheckFitAll(pHR, pVol);

    // Calibrate
  if(m_ptsCDS)
    m_pParametrization->Calibrate(*m_pOption,*m_ptsCDS);
  else
    m_pParametrization->Calibrate(*m_pOption,*m_ptsParBond);

  pHR = m_pParametrization->GetHazardRate();
  pVol = m_pParametrization->GetVolatility();


  CheckFitAll(pHR, pVol);

}


void Parametrization_VolFlatHRWithTimeComponentTest::TestGeneratedValues()
{
  // Read in the current data
  Setup();

  if(m_ptsParBond)
    return;

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolFlatHRWithTimeComponent() );

  // Generate prices using timeonly hazard rates and flat vols.  Should 
  // be able to fit these prices exactly
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

  double dVolStep = 0.1;
  double dVol = 0.2;
  while (dVol < 0.81)
  {
    // Setup the model and generate prices
    shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
    shared_ptr<HazardRate> pHR( new HazardRateTimeOnly(pDates, pdValues));

    ihg::TheoreticalModel model;
    model.SetVolatility(pVol);
    model.SetHazardRate(pHR);

    for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
    {
      CDSLike& cds = *(*ppCDS);

      shared_ptr<finance::ModelOutput> output = model.Compute(cds);

      // Take the computed prices as market price for testing
      cds.SetMarketPrice( output->GetPrice() );
    }

    shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption);
    m_pOption->SetMarketPrice( output->GetPrice() );
    
    //Calibrate    
    m_pParametrization->CalibrateWithOptionAndCDSs(*m_pOption, *m_ptsCDS);
    
    shared_ptr<HazardRateWithTimeComponent> 
      pCalibratedHR = m_pParametrization->GetHazardRate();
    shared_ptr<VolatilityFlat> 
      pCalibratedVol = m_pParametrization->GetVolatility();

    // Check if the data was reproduced exactly
    CheckFitAll(pHR, pVol);

    // Check that the params were reproduced exactly
    double dCalVol = pCalibratedVol->GetValue();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dCalVol, dVol, 1.e-4);   

    std::vector<double> pdCalValues = pCalibratedHR->GetTimeComponentValues();
    CPPUNIT_ASSERT( nNbElements == pdCalValues.size() );

    for (nIdx = 0; nIdx < nNbElements; nIdx++)
    {
      //std::cout << nIdx 
      //          << ", " << pdValues[nIdx] 
      //          << ", " << pdCalValues[nIdx] 
      //          << std::endl;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pdValues[nIdx], pdCalValues[nIdx], 1.e-4);   
    }

    dVol += dVolStep;
  } // loop over vols

}


void Parametrization_VolFlatHRWithTimeComponentTest::TestGeneratedValuesCalibrateXXX()
{
  // Read in the current data
  Setup();

  if(m_ptsParBond)
    return;

  // Make a new parametrization
  m_pParametrization = make_ptr( new ParametrizationVolFlatHRWithTimeComponent() );

  // Generate prices using timeonly hazard rates and flat vols.  Should 
  // be able to fit these prices exactly
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

  double dVolStep = 0.1;
  double dVol = 0.2;
  while (dVol < 0.81)
  {
    // Setup the model and generate prices
    shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
    shared_ptr<HazardRate> pHR( new HazardRateTimeOnly(pDates, pdValues));

    ihg::TheoreticalModel model;
    model.SetVolatility(pVol);
    model.SetHazardRate(pHR);

    for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
    {
      CDSLike& cds = *(*ppCDS);

      shared_ptr<finance::ModelOutput> output = model.Compute(cds);

      // Take the computed prices as market price for testing
      cds.SetMarketPrice( output->GetPrice() );
    }

    shared_ptr<finance::ModelOutput> output = model.Compute(*m_pOption);
    m_pOption->SetMarketPrice( output->GetPrice() );
    
    // Calibrate    
    m_pParametrization->CalibrateWithOptionAndCDSs(*m_pOption, *m_ptsCDS);

    shared_ptr<HazardRateWithTimeComponent> 
      pCalibratedHR = m_pParametrization->GetHazardRate();
    shared_ptr<VolatilityFlat> 
      pCalibratedVol = m_pParametrization->GetVolatility();

    // Check if the data was reproduced exactly
    CheckFitAll(pHR, pVol);

    // Check that the params were reproduced exactly
    double dCalVol = pCalibratedVol->GetValue();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dCalVol, dVol, 1.e-4);   

    std::vector<double> pdCalValues = pCalibratedHR->GetTimeComponentValues();
    CPPUNIT_ASSERT( nNbElements == pdCalValues.size() );

    for (nIdx = 0; nIdx < nNbElements; nIdx++)
    {
      //std::cout << nIdx 
      //          << ", " << pdValues[nIdx] 
      //          << ", " << pdCalValues[nIdx] 
      //          << std::endl;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pdValues[nIdx], pdCalValues[nIdx], 1.e-4);   
    }

    dVol += dVolStep;
  } // loop over vols

}
