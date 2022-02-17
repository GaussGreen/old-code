#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/option.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/exercisetype.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/payoffconstant.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/surfacedouble.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

#include "ihg/optionnumoutput.h"
#include "ihg/optionpricer.h"
#include "ihg/model.h"

#include "optiontester.h"

using namespace std;
using namespace ito33::finance;
using namespace ito33::pricing;

// Added for testing uniform mesh refinement
//extern size_t nGlobalRefine;

using ito33::ihg::OptionTester;

const double dTiny         = 1.e-13;
const double dRateCompare  = 0.5;
const double dPriceCompare = 1.e-4;
const double dGammaCompare = 1.e-3;

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{
  namespace finance
  {
    extern const char *OUTPUT_PRICE;
    extern const char *OUTPUT_DELTA;
    extern const char *OUTPUT_GAMMA;
    extern const char *OUTPUT_THETA;
    extern const char *OUTPUT_VEGA;
    extern const char *OUTPUT_RHO;
  }
}

void OptionTester::ReadContractParams(std::ifstream& sIn)
{

  std::string sMaturity;
  int iAmericanEuropean;
  double dStrike;
  sIn >> sMaturity >> dStrike >> m_iType >> iAmericanEuropean;

  Date maturityDate(sMaturity.c_str(), "%Y-%m-%d");

  finance::OptionType optionType = finance::Option_Call;
  if (m_iType == 0)
    optionType = finance::Option_Put;
  else if (m_iType == 1)
    optionType = finance::Option_Call;
  else if (m_iType == 2)
    optionType = finance::Option_Digital;
  else
    optionType = finance::Option_Other;

  finance::ExerciseType exerciseType = ExerciseType_American;
  if (!iAmericanEuropean)
    exerciseType = ExerciseType_European;

  // Save contract parameters in option object
  m_pOptionFinance = shared_ptr<finance::Option>
    ( new finance::Option(dStrike, maturityDate, optionType, exerciseType));

}

void OptionTester::ReadContractParams(std::string& sContractFile)
{
  ifstream sIn;
  sIn.open(sContractFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening option param file");
  ReadContractParams(sIn);
  sIn.close();
}


void OptionTester::RunPricer()
{

  // Specify which Greeks need to be computed
  m_bComputeDelta = true;
  m_bComputeGamma = true;
  m_bComputeVega = true;
  m_bComputeRho = true;
  m_bComputeTheta = true;
  m_bComputeArrays = true;
  m_bComputeSurfaces = true;

  shared_ptr<finance::Numeraire>
    pNumeraire( new finance::Numeraire("EUR") );

  shared_ptr<finance::RateData> 
    pRateData( new finance::RateData );

  pRateData->SetYieldCurve(pNumeraire, m_pYieldCurve);

  shared_ptr<finance::Equity>
    pEquity( new finance::Equity(m_dSpot, pNumeraire) );

  pEquity->SetBorrowCurve( m_pForeignCurve );
  pEquity->SetDividends( m_pDividends );

  shared_ptr<finance::SessionData> 
    pSessionData( new finance::SessionData(pRateData, pEquity, 
                                           m_ValuationDate) );

  m_pOptionFinance->SetSessionData( pSessionData );

  // Create the internal option and parameter objects
  pricing::Option option(*m_pOptionFinance);

  pricing::OptionParams optionParams(option);

  optionParams.SetValuationTime( GetDoubleFrom(m_ValuationDate) );

  optionParams.SetSpotSharePrice(m_dSpot);
  
  optionParams.SetYieldCurve( m_pYieldCurve );
  
  optionParams.SetYieldCurveForMesh( m_pYieldCurve );
  
  optionParams.SetForeignCurve( m_pForeignCurve );

  optionParams.SetDividends(m_pDividends);
 
  optionParams.SetNumParams(m_pNumParams);

  optionParams.SetMeshParams(m_pMeshParams);

  if (m_iType == 3)
  {
    shared_ptr<finance::Payoff> pPayoff(new finance::PayoffConstant(1.0));
    optionParams.SetPayoff(pPayoff);
  }
     
  shared_ptr<ComputationalFlags> pComputationalFlags(new ComputationalFlags);
 
  pComputationalFlags->SetComputeVega(m_bComputeVega);    
  pComputationalFlags->SetComputeRho(m_bComputeRho);
  pComputationalFlags->SetComputeSurface(m_bComputeSurfaces);
  

  ihg::Model model(m_pVolatility, m_pVolatility, m_pHazardRate);
  

  // The actual pricer
  OptionPricer pricer(optionParams, model, *pComputationalFlags);
   
  m_pOutput = pricer.Price()->GetModelOutput();

}

void OptionTester::OutputPlotData()
{
  // Run the pricer and create the output
  RunPricer();

  // TOFIX
  throw EXCEPTION_MSG( ITO33_UNEXPECTED,
    "This is not implemented after the recent modification of core design.");
 
}


void OptionTester::ReportConvergenceTest()
{
  ASSERT_MSG(m_nNbRuns > 0, 
        "Must call RunConvergenceTest before calling ReportConvergenceTest");

  cout << endl;
  cout << "Data at S=" << m_dSpot << endl << endl;

  string strText = "Price";
  ReportConvergence(m_nNbRuns, m_pdPrices.Get(), m_pdTimes.Get(), strText);

 
  if (m_bComputeDelta)
  {
    strText = "Delta";
    cout << endl;
    ReportConvergence(m_nNbRuns, m_pdDeltas.Get(), m_pdTimes.Get(), strText);
  }
  if (m_bComputeGamma)
  {
    strText = "Gamma";
    cout << endl;
    ReportConvergence(m_nNbRuns, m_pdGammas.Get(), m_pdTimes.Get(), strText);
  }
  if (m_bComputeVega)
  {
    strText = "Vega ";
    cout << endl;
    ReportConvergence(m_nNbRuns, m_pdVegas.Get(), m_pdTimes.Get(), strText);
  }
  if (m_bComputeRho)
  {
    strText = "Rho  ";
    cout << endl;
    ReportConvergence(m_nNbRuns, m_pdRhos.Get(), m_pdTimes.Get(), strText);
  }
  if (m_bComputeTheta)
  {
    strText = "Theta";
    cout << endl;
    ReportConvergence(m_nNbRuns, m_pdThetas.Get(), m_pdTimes.Get(), strText);
  }
  cout << endl;
  
}


double OptionTester::GetDelta()
{
  ASSERT_MSG(m_nNbRuns > 0, 
             "Must call RunConvergenceTest before GetDelta");

  return m_pdDeltas[m_nNbRuns-1];
}


double OptionTester::GetGamma()
{
  ASSERT_MSG(m_nNbRuns > 0, 
             "Must call RunConvergenceTest before GetGamma");

  return m_pdGammas[m_nNbRuns-1];
}


double OptionTester::GetVega()
{
  ASSERT_MSG(m_nNbRuns > 0, 
             "Must call RunConvergenceTest before GetVega");

  return m_pdVegas[m_nNbRuns-1];
}


bool OptionTester::RunFullTests(double dPrice, double dDelta, double dGamma, 
                                double dVega, size_t nNbTests)
{
  std::string sNumParamFile;
/*
  cout << "  Expected price = " << dPrice << std::endl; 
  cout << "  Expected delta = " << dDelta << std::endl; 
  cout << "  Expected gamma = " << dGamma << std::endl; 
  cout << "  Expected vega = "  << dVega  << std::endl; 
  cout << std::endl;
*/

  // Implicit testing
  sNumParamFile = "..\\txtfiles\\implicitnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  double dComputedPrice = GetPrice();
  bool bResult = ReportPass( (dPrice - dComputedPrice)/dPrice, 0.0, dPriceCompare*10.0);
  //cout << "  Implicit price = " << dComputedPrice << std::endl;

  if ( !bResult )
    return bResult;

 
  double dRate = GetConvergenceRate();
  bResult = ReportPass(2.0, dRate, dRateCompare);
  //cout << "  Implicit convergence rate = " << dRate << std::endl;
  
  if ( !bResult && !(dRate > 3.0) )
    return bResult;

  
  double dComputedDelta = GetDelta();
  bResult = ReportPass( (dDelta - dComputedDelta)/dDelta, 0.0, dPriceCompare*20.0);
  //cout << "  Implicit delta = " << dComputedDelta << std::endl;

  if ( !bResult )
    return bResult;

  
  double dComputedGamma = GetGamma();
  bResult = ReportPass( (dGamma - dComputedGamma)/dGamma, 0.0, dGammaCompare*10.0);
  //cout << "  Implicit gamma = " << dComputedGamma << std::endl;

  if ( !bResult )
    return bResult;

  double dComputedVega = GetVega();
  bResult = ReportPass( (dVega - dComputedVega)/dVega, 0.0, dPriceCompare*10.0);
  //cout << "  Implicit vega = " << dComputedVega << std::endl;

  if ( !bResult )
    return bResult;

  // Crank-Nicolson testing
  sNumParamFile = "..\\txtfiles\\cnnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  dComputedPrice = GetPrice();
  bResult = ReportPass( (dPrice - dComputedPrice)/dPrice, 0.0, dPriceCompare);
  //cout << "  Crank-Nicolson price = " << dComputedPrice << std::endl;

  if ( !bResult )
    return bResult;
    
  dRate = GetConvergenceRate();
  bResult = ReportPass(4.0, dRate, .9);
  //cout << "  Crank-Nicolson convergence rate = " << dRate << std::endl;

  if ( !bResult )
    return bResult;

  dComputedDelta = GetDelta();
  bResult = ReportPass( (dDelta - dComputedDelta)/dDelta, 0.0, dPriceCompare*20);
  //cout << "  Crank-Nicolson delta = " << dComputedDelta << std::endl;  

  if ( !bResult )
    return bResult;
  
  dComputedGamma = GetGamma();
  bResult = ReportPass( (dGamma - dComputedGamma)/dGamma, 0.0, dGammaCompare);
  //cout << "  Crank-Nicolson gamma = " << dComputedGamma << std::endl;

  if ( !bResult )
    return bResult;
  
  dComputedVega = GetVega();
  bResult = ReportPass( (dVega - dComputedVega)/dVega, 0.0, dPriceCompare);
  //cout << "  Crank-Nicolson vega = " << dComputedVega << std::endl;

  if ( !bResult )
    return bResult;

  // three-level (BDF) testing
  sNumParamFile = "..\\txtfiles\\bdfnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  dComputedPrice = GetPrice();
  bResult = ReportPass( (dPrice - dComputedPrice)/dPrice, 0.0, dPriceCompare);
  //cout << "  Three-level (BDF) price = " << dComputedPrice << std::endl;
    
  if ( !bResult )
    return bResult;

  dRate = GetConvergenceRate();
  bResult = ReportPass(4.0, dRate, dRateCompare);
  //cout << "  Three-level (BDF) convergence rate = " << dRate << std::endl;

  if ( !bResult )
    return bResult;
 
  dComputedDelta = GetDelta();
  bResult = ReportPass( (dDelta - dComputedDelta)/dDelta, 0.0, dPriceCompare*20);
  //cout << "  Three-level (BDF) delta = " << dComputedDelta << std::endl;  

  if ( !bResult )
    return bResult;

  dComputedGamma = GetGamma();
  bResult = ReportPass( (dGamma - dComputedGamma)/dGamma, 0.0, dGammaCompare);
  //cout << "  Three-level (BDF) gamma = " << dComputedGamma << std::endl;  
    
  if ( !bResult )
    return bResult;
 
  dComputedVega = GetVega();
  bResult = ReportPass( (dVega - dComputedVega)/dVega, 0.0, dPriceCompare);
  //cout << "  Three-level (BDF) vega = " << dComputedVega << std::endl;
  //cout << std::endl;

  return bResult;
}

