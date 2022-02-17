#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/constants.h"

#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/surfacedouble.h"
#include "ito33/finance/surfaceflag.h"
#include "ito33/finance/domain.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/schemetype.h"

#include "ihg/basetester.h"

// Added for testing uniform mesh refinement
//size_t nGlobalRefine;

using namespace std; 

using ito33::finance::YieldCurve;
using ito33::finance::Dividends;

const double dTiny = 1.e-13;
const double dRateCompare = 0.5;        // converged rate comparison
const double dPriceCompare = 2.*1.e-4;
const double dImplicitPriceCompare = 2.*1.e-3;



void __stdcall parametricHR(double /* dTime */, const double *pdS, double *pdValues, 
                            size_t nNbS, int /* i */)
{
  // Implement the standard HR formula hr(S,t) = alpha * (S0/S)^beta

  // Set default parameters
  const double dAlpha = 0.2;
  const double dBeta = 0.5;
  const double dS0 = 100.0;

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    if (pdS[nIdx] < 1.e-16)
      pdValues[nIdx] = dAlpha * pow((dS0 / 1.e-16), dBeta);
    else
      pdValues[nIdx] = dAlpha * pow((dS0 / pdS[nIdx]), dBeta);
  }
}


void __stdcall parametricVol(double dTime, const double *pdS, double *pdValues, 
                            size_t nNbS, int /* i */)
{
  // Implement the formula vol(S,t) = 0.2 + 0.00002*(100 - S)^2, 0.4 - t/10.0
  // Also, cap the values at 0.4, and make sure they are never negative 

  dTime -= ito33::Date(2003, ito33::Date::Jan, 1).GetExcel() * ito33::ONEDAY;

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    pdValues[nIdx] = 0.2 + 0.00002*(100.0 - pdS[nIdx])*(100.0 - pdS[nIdx]) - dTime/10.0;
    if (pdValues[nIdx] > 0.4) 
      pdValues[nIdx] = 0.4;
    if (pdValues[nIdx] < 0.0) 
      pdValues[nIdx] = 0.0;

  }
}


namespace ito33
{

namespace ihg
{


void BaseTester::ReadInputFile(string& sFilename)
{
  // Data file
  ifstream sIn(sFilename.c_str());
  
  ASSERT_MSG(sIn, 
             "Problem opening data file for testing");    

  // Read the common data
  ReadValuationDate(sIn);
  ReadYieldCurve(sIn);
  ReadForeignCurve(sIn);
  ReadDividends(sIn);
  ReadVolatility(sIn);
  ReadHazardRate(sIn);
  ReadMeshParams(sIn);
  ReadNumParams(sIn);

  sIn >> m_dSpot;

  // Read contract specific data  
  ReadContractParams(sIn);
  
  sIn.close();
}


void BaseTester::ReadFromMultipleFiles(std::string& sYieldCurveFile,
    std::string& sForeignCurveFile,
    std::string& sDividendFile,
    std::string& sVolatilityFile,
    std::string& sHazardRateFile,
    std::string& sMeshParamFile,
    std::string& sNumParamFile,
    std::string& sValuationDateFile,
    std::string& sContractFile)
{

  // No need for comments.
  ReadYieldCurve(sYieldCurveFile);
  ReadForeignCurve(sForeignCurveFile);
  ReadDividends(sDividendFile);
  ReadVolatility(sVolatilityFile);
  ReadHazardRate(sHazardRateFile);
  ReadMeshParams(sMeshParamFile);
  ReadNumParams(sNumParamFile);
  ReadValuationDate(sValuationDateFile);

  ReadContractParams(sContractFile);
}



void BaseTester::ReadValuationDate(ifstream& sIn)
{
  //int iDate;
  //sIn >> iDate;
  //m_ValuationDate.Set(iDate);

  std::string sDate;
  sIn >> sDate;
  m_ValuationDate.Set(sDate.c_str(), "%Y-%m-%d");
  
}

void BaseTester::ReadValuationDate(std::string& sValuationDateFile)
{
  ifstream sIn;
  sIn.open(sValuationDateFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening pricing date file");
  ReadValuationDate(sIn);
  sIn.close();
}

void BaseTester::ReadCurve(ifstream& sIn, int iWhich)
{
  // iWhich == 1 means the usual yield curve
  // iWhich == 2 mean the foreign curve
  ASSERT_MSG(iWhich == 1 || iWhich == 2, "Unknown curve in ReadCurve");
  int iType;
  sIn >> iType;

  std::string sDate;
  switch (iType)
  {
    case 0:
    {
      // Flat yield curve
      double dRate;

      sIn >> dRate;

      if (iWhich == 1)
        m_pYieldCurve = make_ptr( new finance::YieldCurveFlat(dRate) );
      else if (iWhich == 2)
        m_pForeignCurve = make_ptr( new finance::YieldCurveFlat(dRate) );

      break;
    }
    case 1:
    {
      // Zero-coupon yield curve
      size_t nIdxLeg, nNbLegs;
      
      sIn >> nNbLegs;
      sIn >> sDate;
      Date date(sDate.c_str(), "%Y-%m-%d");

      shared_ptr<finance::YieldCurveAnnuallyCompounded> 
        pYieldCurve(new finance::YieldCurveAnnuallyCompounded(date, nNbLegs));

      Array<double> pdMaturity = Array<double>(nNbLegs);
      Array<double> pdRate = Array<double>(nNbLegs);

      double dMaturity;
      double dRate;
      // loop, getting each maturity and rate
      for (nIdxLeg = 0; nIdxLeg < nNbLegs; nIdxLeg++)
      {
        sIn >> dMaturity >> dRate;

        pYieldCurve->AddLeg( size_t(dMaturity * 365. + 1.e-4), dRate );
      }
      if (iWhich == 1)
        m_pYieldCurve = pYieldCurve;
      else if (iWhich == 2)
        m_pForeignCurve = pYieldCurve;

      break;
    }
    default:
      ASSERT_MSG(false, 
             "Only flat and zero-coupon yield curves are supported in test files");
  }

}

void BaseTester::ReadYieldCurve(ifstream& sIn)
{
  ReadCurve(sIn, 1);
}

void BaseTester::ReadYieldCurve(std::string& sYieldCurveFile)
{
  ifstream sIn;
  sIn.open(sYieldCurveFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening yield curve file");
  ReadYieldCurve(sIn);
  sIn.close();
}


void BaseTester::ReadForeignCurve(ifstream& sIn)
{
  ReadCurve(sIn, 2);
}

void BaseTester::ReadForeignCurve(std::string& sForeignCurveFile)
{
  ifstream sIn;
  sIn.open(sForeignCurveFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening foreign curve file");
  ReadForeignCurve(sIn);
  sIn.close();
}

void BaseTester::ReadDividends(ifstream& sIn)
{ 
  m_pDividends = make_ptr( new Dividends() );

  size_t nNbDividends;

  sIn >> nNbDividends; 

  int iYieldCash;
  std::string sDate;
  double dValue;
  double dSita;

  // Read the dividends from the input file, and construct the dividends objet
  for (size_t nIdxDividend = 0; nIdxDividend < nNbDividends; nIdxDividend++)
  {
    sIn >> iYieldCash >> sDate >> dValue;
    Date date(sDate.c_str(), "%Y-%m-%d");

    if (iYieldCash)
      m_pDividends->AddYield(date, dValue);
    else
    {
      sIn >> dSita;

      m_pDividends->AddPseudoCash(date, dValue, dSita);
    }
  }
}


void BaseTester::ReadDividends(std::string& sDividendFile)
{
  ifstream sIn;
  sIn.open(sDividendFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening dividend file");
  ReadDividends(sIn);
  sIn.close();
}


void BaseTester::ReadMeshParams(ifstream& sIn)
{
  bool bUniformTimeGrid, bUniformSpaceGrid;
  sIn >> bUniformTimeGrid;
  sIn >> bUniformSpaceGrid;

  size_t nTimeAccumulation, nTimeCompression;
  size_t nSpaceAccumulation, nSpaceCompression;
  sIn >> nTimeAccumulation;
  sIn >> nTimeCompression;
  sIn >> nSpaceAccumulation;
  sIn >> nSpaceCompression;

  double dTimeStretch, dSpaceStretch, dGridSpan;

  sIn >> dTimeStretch;
  sIn >> dSpaceStretch;
  sIn >> dGridSpan;

  m_pMeshParams = make_ptr( new numeric::MeshParams() );
  m_pMeshParams->SetUniformTimeGrid(bUniformTimeGrid);
  m_pMeshParams->SetUniformSpaceGrid(bUniformSpaceGrid);
  m_pMeshParams->SetTimeAccumulation(nTimeAccumulation);
  m_pMeshParams->SetTimeCompression(nTimeCompression);
  m_pMeshParams->SetSpaceAccumulation(nTimeAccumulation);
  m_pMeshParams->SetSpaceCompression(nTimeCompression);
  m_pMeshParams->SetTimeStretch(dTimeStretch);
  m_pMeshParams->SetSpaceStretch(dSpaceStretch);
  m_pMeshParams->SetGridSpan(dGridSpan);

}


void BaseTester::ReadMeshParams(std::string& sMeshParamFile)
{
  ifstream sIn;
  sIn.open(sMeshParamFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening mesh param file");
  ReadMeshParams(sIn);
  sIn.close();
}


void BaseTester::ReadNumParams(ifstream& sIn)
{
  size_t nNbTimeSteps, nNbSpaceSteps, nSchemeType;
  sIn >> nNbTimeSteps;
  sIn >> nNbSpaceSteps;
  sIn >> nSchemeType;

  m_pNumParams = make_ptr( new numeric::NumParams() );

  m_pNumParams->SetNbTimeSteps(nNbTimeSteps);
  m_pNumParams->SetNbSpaceSteps(nNbSpaceSteps);

  switch (nSchemeType)
  {
  case 0:
    {
      m_pNumParams->SetSchemeType(numeric::SchemeType_Implicit);
     break;
    }
  case 1:
    {
      m_pNumParams->SetSchemeType(numeric::SchemeType_CrankNicolson);
      break;
    }
  case 2:
    {
      m_pNumParams->SetSchemeType(numeric::SchemeType_ThreeLevel);
      break;
    }
  default:
    ASSERT_MSG(false, 
             "Unknown scheme type in input file");
  }
  
}


void BaseTester::ReadNumParams(std::string& sNumParamFile)
{
  ifstream sIn;
  sIn.open(sNumParamFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening num param file");
  ReadNumParams(sIn);
  sIn.close();
}


void BaseTester::ReadVolatility(ifstream& sIn)
{
  // Read in the type of vol first (constant, surface, etc)
  int iVolType;
  sIn >> iVolType;

  switch (iVolType)
  {
  case 0:
    {
    // Constant/flat volatility
    double dVol;
    sIn >> dVol;
    m_pVolatility = make_ptr( new VolatilityFlat(dVol) ); 
    break;
    }
    /*
  case 1:
    {
    // polynomial volatility surface
    size_t nNbParams, nNbTimes;
    double pdC[5];
    double dTime;

    sIn >> nNbParams;
    sIn >> nNbTimes;
    PolynomialVol* pTmp = new PolynomialVol();
    m_pVolatility = pTmp;
    
    for (size_t nIdx = 0; nIdx < nNbTimes; nIdx++)
    { 
      sIn >> dTime;
      sIn >> pdC[0] >> pdC[1] >> pdC[2] >> pdC[3] >> pdC[4];
      //sIn >> pdC[0] >> pdC[1] >> pdC[2];
      pTmp->SetParams(dTime, pdC);
    }
    pTmp->Write(std::cout);
    break;
    }
    */
  default:
    ASSERT_MSG(false, 
             "Only flat volatility is currently supported");
  }
}

void BaseTester::ReadVolatility(std::string& sVolatilityFile)
{
  ifstream sIn;
  sIn.open(sVolatilityFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening volatility file");
  ReadVolatility(sIn);
  sIn.close();
}

void BaseTester::ReadHazardRate(ifstream& sIn)
{
  // Read in the type of hazard rate first (constant, surface, etc)
  int iHRType;
  sIn >> iHRType;

  switch (iHRType)
  {
  case 0:
    {
    // Constant/flat hazard rate
    double dHR;
    sIn >> dHR;
    // Fake a constant rate by having a constant rate for a very long period
    Array<Date> pDates(2);
    pDates[0] = Date("1970/01/01");
    pDates[1] = Date("2500/01/01");

    Array<double> pdRates(2);
    pdRates[0] = dHR;
    pdRates[1] = dHR;
    m_pHazardRate = make_ptr( new HazardRateTimeOnly
                                  ( pDates.Get(), pdRates.Get(), 2 ) ); 
    break;
    }
    
  case 1:
    {
    // expontially decaying hazard rate
    double dAlpha, dS0;
    sIn >> dAlpha;
    sIn >> dS0;
    m_pHazardRate = make_ptr( new HazardRateDecay(dAlpha, dS0) ); 
    break;
    }

  case 2:
    {
    // standard power hazard rate
    double dAlpha, dBeta, dS0;
    sIn >> dAlpha;
    sIn >> dBeta;
    sIn >> dS0;
    m_pHazardRate = make_ptr( new HazardRatePower(dAlpha, dBeta, dS0) ); 
    break;
    }
    
  default:
    ASSERT_MSG(false, 
             "Unknown hazard rate type");
  }
}

void BaseTester::ReadHazardRate(std::string& sHazardRateFile)
{
  ifstream sIn;
  sIn.open(sHazardRateFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening hazard rate file");
  ReadHazardRate(sIn);
  sIn.close();
}

double BaseTester::GetPrice()
{
  ASSERT_MSG(m_nNbRuns > 0, 
             "Must call RunConvergenceTest before GetPrice");

  return m_pdPrices[m_nNbRuns-1];
}


double BaseTester::GetConvergenceRate()
{
  ASSERT_MSG(m_nNbRuns > 2, 
             "Must have at least 3 runs to compute a convergence rate");

  return (m_pdPrices[m_nNbRuns-2] - m_pdPrices[m_nNbRuns-3]) /
         (m_pdPrices[m_nNbRuns-1] - m_pdPrices[m_nNbRuns-2]);
}


void BaseTester::ReportConvergence
(size_t nNbValues, double* pdValues, double* pdTimes, std::string& strHeading)
{

  ASSERT_MSG(nNbValues >= 3, 
             "Need at least 3 data points to compute convergence");

  cout.precision(8);

  cout << "Run   \t" << "Time  \t" 
       << strHeading << "    \tDifference\t" 
       << "Convergence" << endl;
  cout << "------\t" << "------\t" 
       << "--------\t" << "----------\t" 
       << "-----------" << endl;

  // The convergence rate is simply the ratio of successive differences
  // in the (price) data
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  {
    cout << nIdx << "\t" << pdTimes[nIdx] << "\t" << pdValues[nIdx];

    if (nIdx > 0)
      cout << "\t" << pdValues[nIdx] - pdValues[nIdx-1];
    if (nIdx > 1)
      cout << "\t" <<  (pdValues[nIdx - 1] - pdValues[nIdx - 2])
                     / (pdValues[nIdx] - pdValues[nIdx - 1]);

    cout << endl;
  }
}


bool BaseTester::ReportPass(double dExpected, double dComputed, double dTol)
{  
  if ( fabs(dComputed - dExpected ) < dTol )
    return true;

  return false;
}




void BaseTester::RunConvergenceTest(size_t nNbRuns)
{
  // Verify and sace the number of runs
  ASSERT_MSG(nNbRuns > 0, "Number of runs must be greater than zero");    

  m_nNbRuns = nNbRuns;

  // Allocate memory for convergence testing
  m_pdTimes = Array<double>(nNbRuns);
  m_pdPrices = Array<double>(nNbRuns);
  m_pdDeltas = Array<double>(nNbRuns);
  m_pdGammas = Array<double>(nNbRuns);
  m_pdVegas = Array<double>(nNbRuns);
  m_pdRhos = Array<double>(nNbRuns);
  m_pdThetas = Array<double>(nNbRuns);

  // Save the mesh sizes, since they need to be reset after being changed
  // during the convergence testing (in case this function is called
  // multiple times without re-reading mesh data)
  size_t nIdx = 0;
  size_t nNbSpacesteps = m_pNumParams->GetNbSpaceSteps();
  size_t nNbTimesteps = m_pNumParams->GetNbTimeSteps();
  size_t nNbSpacestepsSave = nNbSpacesteps;
  size_t nNbTimestepsSave = nNbTimesteps;

  for (nIdx = 0; nIdx < nNbRuns; nIdx++)
  {
    // Set the mesh sizes for this run 
    m_pNumParams->SetNbSpaceSteps(nNbSpacesteps);
    m_pNumParams->SetNbTimeSteps(nNbTimesteps);
    nNbSpacesteps = nNbSpacesteps * 2 - 1;
    nNbTimesteps = nNbTimesteps * 2 - 1;

    //nGlobalRefine = nIdx;
    //nGlobalRefine = 0;

    // Actually run the appropriate price
    clock_t t1 = clock();
    RunPricer();
    clock_t t2 = clock();

    // Save the output from this run
    m_pdTimes[nIdx] = (t2-t1)/1000.;

    m_pdPrices[nIdx] = m_pOutput->GetPrice();
    
    if (m_bComputeDelta)
      m_pdDeltas[nIdx] = m_pOutput->GetDelta();

    if (m_bComputeGamma)
      m_pdGammas[nIdx] = m_pOutput->GetGamma();

    if (m_bComputeVega)
      m_pdVegas[nIdx] = m_pOutput->GetVega();

    if (m_pOutput->HasRho())
      m_pdRhos[nIdx] = m_pOutput->GetRho();

    if (m_bComputeTheta)
      m_pdThetas[nIdx] = m_pOutput->GetTheta();
    
  }

  // reset the mesh params
  m_pNumParams->SetNbSpaceSteps(nNbSpacestepsSave);
  m_pNumParams->SetNbTimeSteps(nNbTimestepsSave);
}


bool BaseTester::RunPriceTests(double dExpected, double dTol)
{
  std::string sNumParamFile;

  //cout << "  Expected price = " << dExpected << std::endl; 
  //cout << std::endl;  

  // Implicit testing
  sNumParamFile = "..\\txtfiles\\implicitnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(1);
  double dPrice = m_pdPrices[0];
  bool bResult = ReportPass( (dExpected - dPrice)/dExpected, 0.0, dTol * 10.0);

  if ( !bResult )
    return bResult;

  //cout << "  Implicit price = " << dPrice << std::endl;  

  // Crank-Nicolson testing
  sNumParamFile = "..\\txtfiles\\cnnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(1);
  dPrice = m_pdPrices[0];
  bResult = ReportPass( (dExpected - dPrice)/dExpected, 0.0, dTol);

  if ( !bResult )
    return bResult;

  //cout << "  Crank-Nicolson price = " << dPrice << std::endl;
  
  // three-level (BDF) testing
  sNumParamFile = "..\\txtfiles\\bdfnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(1);
  dPrice = m_pdPrices[0];
  bResult = ReportPass( (dExpected - dPrice)/dExpected, 0.0, dTol);
  
  //cout << "  Three-level (BDF) price = " << dPrice << std::endl;
  
  //cout << std::endl;

  return bResult;
}


bool BaseTester::RunPriceConvergenceTests(double dExpected, size_t nNbTests)
{
  std::string sNumParamFile;

  //cout << "  Expected price = " << dExpected << std::endl; 
  //cout << std::endl;  

  // Implicit testing
  sNumParamFile = "..\\txtfiles\\implicitnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  double dPrice = GetPrice();  
  
  bool bResult = ReportPass( fabs(dExpected - dPrice)/fabs(dExpected), 0.0, dImplicitPriceCompare);
  
  if ( !bResult )
    return bResult;

  //cout << "  Implicit price = " << dPrice << std::endl;
  double dRate = GetConvergenceRate();
  bResult = ReportPass(2.0, dRate, dRateCompare);  
  
  if ( !bResult )
    return bResult;

  //cout << "  Implicit convergence rate = " << dRate << std::endl;  

  // Crank-Nicolson testing
  sNumParamFile = "..\\txtfiles\\cnnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  dPrice = GetPrice();
  bResult = ReportPass( (dExpected - dPrice)/dExpected, 0.0, dPriceCompare);
  
  if ( !bResult )
    return bResult;

  //cout << "  Crank-Nicolson price = " << dPrice << std::endl;
  dRate = GetConvergenceRate();
  bResult = ReportPass(4.0, dRate, dRateCompare);
  
  if ( !bResult )
    return bResult;

  //cout << "  Crank-Nicolson convergence rate = " << dRate << std::endl;
  
  
  // three-level (BDF) testing
  sNumParamFile = "..\\txtfiles\\bdfnumparams.txt";
  ReadNumParams(sNumParamFile);
  RunConvergenceTest(nNbTests);
  dPrice = GetPrice();
  bResult = ReportPass( (dExpected - dPrice)/dExpected, 0.0, dPriceCompare);
  
  if ( !bResult )
    return bResult;

  //cout << "  Three-level (BDF) price = " << dPrice << std::endl;
  dRate = GetConvergenceRate();
  bResult = ReportPass(4.0, dRate, dRateCompare);

 // cout << "  Three-level (BDF) convergence rate = " << dRate << std::endl;

  if ( !bResult )
    return bResult;

  //cout << std::endl;

  return bResult;
}


void BaseTester::OutputArray(const std::vector<double>& pdSpots, const std::vector<double>& pdValues,
    std::string sFileName)
{
  std::ofstream sOut;
  sOut.open(sFileName.c_str(), std::ios::out);
  ASSERT_MSG(sOut, "Problem opening output file");
  
  size_t nNb = pdSpots.size();
  for (size_t nIdx = 0; nIdx < nNb; nIdx++)
    sOut << pdSpots[nIdx] << " " << pdValues[nIdx] << std::endl;

  sOut.close();
}

void BaseTester::OutputSurface(const finance::SharedSurface& surface, 
    std::string sFileName, double dSLeft, double dSRight)
{
  
  std::ofstream sOut(sFileName.c_str(), std::ios::out);
  ASSERT_MSG(sOut, "Problem opening output file");

  shared_ptr<finance::Domain> domain = surface->GetDomain();
  finance::Domain::Spots spots = domain->GetUnderlyingSharePrices();
  finance::Domain::Dates dates = domain->GetDates();

  finance::Domain::Dates::const_iterator iterDates;

  size_t nIdxTime = 0;
  for (iterDates = dates.begin();
       iterDates != dates.end();
       ++iterDates)
  {
    // Get the values at this time
    finance::SurfaceDouble::Doubles values = surface->GetValuesAt(nIdxTime);
    nIdxTime++;
  
    for (size_t nIdx = 0; nIdx < spots.size(); nIdx++)
    {
      if (spots[nIdx] > dSLeft && spots[nIdx] < dSRight)
      {
        sOut << (*iterDates).GetExcel() * ONEDAY
              << " "
              << spots[nIdx]
              << " "
              << values[nIdx]
              << std::endl;
      }
    }
    sOut << std::endl;
  }

}


void BaseTester::OutputSurface(const shared_ptr<finance::SurfaceFlag>& surface, 
    std::string sFileName, double dSLeft, double dSRight)
{
  
  std::ofstream sOut(sFileName.c_str(), std::ios::out);
  ASSERT_MSG(sOut, "Problem opening output file");

  shared_ptr<finance::Domain> domain = surface->GetDomain();
  finance::Domain::Spots spots = domain->GetUnderlyingSharePrices();
  finance::Domain::Dates dates = domain->GetDates();

  finance::Domain::Dates::const_iterator iterDates;

  size_t nIdxTime = 0;
  for (iterDates = dates.begin();
       iterDates != dates.end();
       ++iterDates)
  {
    // Get the values at this time
    finance::SurfaceFlag::Flags values = surface->GetValuesAt(nIdxTime);
    nIdxTime++;
  
    for (size_t nIdx = 0; nIdx < spots.size(); nIdx++)
    {
      if (spots[nIdx] > dSLeft && spots[nIdx] < dSRight)
      {
        sOut << (*iterDates).GetExcel() * ONEDAY
              << " "
              << spots[nIdx]
              << " "
              << values[nIdx]
              << std::endl;
      }
    }
    sOut << std::endl;
  }

}



} // namespace ihg

} // namespace ito33
