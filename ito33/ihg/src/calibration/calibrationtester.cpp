#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/autoptr.h"

#include "ito33/finance/payoffoption.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/exercisetype.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/optionparams.h"

#include "ihg/optionpricer.h"

#include "ito33/ihg/hazardratedecay.h"
#include "ihg/polynomialvol.h"
#include "ihg/interpvol.h"
#include "ihg/polynomialhr.h"

#include "ihg/calibrator.h"
#include "ihg/calibrationtester.h"

using namespace std;

using ito33::ihg::CalibrationTester;

void CalibrationTester::ReadInputFile(string& sFilename)
{
  // Reset the number of runs (in case this object was called before)
  m_nNbRuns = 0;

  // Data file
  ifstream sIn(sFilename.c_str());
  
  ASSERT_MSG(sIn, 
             "Problem opening file");    

  // Read the common data (done by base class functions)
  ReadValuationDate(sIn);
  ReadYieldCurve(sIn);
  ReadForeignCurve(sIn);
  ReadDividends(sIn);

  // Read in the current spot
  sIn >> m_dSpot;

  // Read in the call data
  sIn >> m_iCallDataType;

  if (m_iCallDataType == 0)
    ReadPriceCallData(sIn);
  else if (m_iCallDataType == 1)
    ReadImpVolCallData(sIn);
  else
    ASSERT_MSG(false, "Unknown call data information type in input file");

  ReadMeshParams(sIn);

  if (m_iCallDataType == 1)
    ConvertImpVolsToPrices();


  sIn.close();
}

void CalibrationTester::ReadPriceCallData(ifstream& sIn)
{
  sIn >> m_nNbCalls;

  m_piCallMaturities = Array<int>(m_nNbCalls);
  m_pdCallStrikes = Array<double>(m_nNbCalls);
  m_pdCallPrices = Array<double>(m_nNbCalls);

  size_t nIdx;
  for (nIdx = 0; nIdx < m_nNbCalls; nIdx++)
    sIn >> m_piCallMaturities[nIdx] 
        >> m_pdCallStrikes[nIdx] 
        >> m_pdCallPrices[nIdx];
}

void CalibrationTester::ReadImpVolCallData(ifstream& sIn)
{
  sIn >> m_nNbCalls;

  m_piCallMaturities = Array<int>(m_nNbCalls);
  m_pdCallStrikes = Array<double>(m_nNbCalls);
  m_pdCallPrices = Array<double>(m_nNbCalls);
  m_pdImpliedVols = Array<double>(m_nNbCalls);

  size_t nIdx;
  for (nIdx = 0; nIdx < m_nNbCalls; nIdx++)
    sIn >> m_piCallMaturities[nIdx] 
        >> m_pdCallStrikes[nIdx] 
        >> m_pdImpliedVols[nIdx];


}

void CalibrationTester::ConvertImpVolsToPrices()
{
  /*
  // convert the implied vols to prices
  for (size_t nIdx = 0; nIdx < m_nNbCalls; nIdx++)
  {
    //ito33::pricing::Option option(;
    option.SetDividends(m_pDividends);
    option.SetForeignCurve(m_pForeignCurve);
    option.SetYieldCurve(m_pYieldCurve);
    option.SetOptionType(ito33::finance::Option_Call);
    option.SetExerciseType(ito33::finance::ExerciseType_European);
    option.SetValuationDate(m_iValuationDate);

    // Create the option params.  We really want an  OptionPricer, but the
    // OptionPricer needs an optionParams object
    OptionParams optionParams(option);
  
    shared_ptr<ComputationalFlags> pcf(new ComputationalFlags());
    pcf->SetComputeDelta(false);
    pcf->SetComputeGamma(false);
    pcf->SetComputeVega(false);
    pcf->SetComputeRho(false);
    pcf->SetComputeTheta(false);
    optionParams.SetComputationalFlags(pcf);
    // Need prices at all grid points
    optionParams.SetComputeArrays(false);  
    optionParams.SetComputeSurfaces(false);

    shared_ptr<ihg::HazardRate> pHazardRate(new HazardRateDecay(0.0, m_dSpot));
    optionParams.SetHazardRate(pHazardRate);

    optionParams.SetMaturity(m_piCallMaturities[nIdx] * ONEDAY);
    optionParams.SetNbRefine(0);
    optionParams.SetSpotSharePrice(m_dSpot);
    optionParams.SetStrike(m_pdCallStrikes[nIdx]);

    shared_ptr<finance::Volatility> pVolatility(new finance::FlatVolatility(m_pdImpliedVols[nIdx]/100.0));
    optionParams.SetVolatility(pVolatility);
  
    optionParams.SetNbS(m_nNbS);
    optionParams.SetNbTimesteps(m_nNbTimeStep);
  
    OptionPricer optionPricer(optionParams);

    double dTmp = optionPricer.Price()->GetPrice()->GetPrice();
    m_pdCallPrices[nIdx] = dTmp;
  }
  */
}


void CalibrationTester::Calibrate(size_t nNbRuns)
{

  // Select the number of runs for convergence testing
  m_pdTimes = Array<double>(nNbRuns);

  size_t nIdx = 0;
  for (nIdx = 0; nIdx < nNbRuns; nIdx++)
  {
//    clock_t t1 = clock();

    //shared_ptr<PolynomialVol> pVolatility( new PolynomialVol() );
    //shared_ptr<PolynomialVol2> pVolatility( new PolynomialVol2() );

    /*
    const nSize = 7;
    double pdPoints[nSize];
    pdPoints[0] = 0.3 * m_dSpot;
    pdPoints[1] = 0.7 * m_dSpot;
    pdPoints[2] = 0.9 * m_dSpot;
    pdPoints[3] = m_dSpot;
    pdPoints[4] = 1.1 * m_dSpot;
    pdPoints[5] = 1.5 * m_dSpot;
    pdPoints[6] = 2.0 * m_dSpot;
    */

    
    const nSize = 5;
    double pdPoints[nSize];
    pdPoints[0] = 0.5 * m_dSpot;
    pdPoints[1] = 0.9 * m_dSpot;
    pdPoints[2] = m_dSpot;
    pdPoints[3] = 1.1 * m_dSpot;
    pdPoints[4] = 1.5 * m_dSpot;

    shared_ptr<InterpVol> pVolatility( new InterpVol(pdPoints, nSize) );
    


    shared_ptr<PolynomialHR> pHazardRate( new PolynomialHR() );
    //Calibrator calibrator(pVolatility, pHazardRate);

    // Hack to get the code to compile. This hack causes the code to crash though
    shared_ptr<ParameterizedVol> pTmpVol( pVolatility.get() );
    shared_ptr<ParameterizedHR> pTmpHR( pHazardRate.get() );
    /*
    Calibrator calibrator(pTmpVol, pTmpHR);

    // Set calibration parameters
//    calibrator.SetValuationDate(m_iValuationDate);
//    calibrator.SetYieldCurve(m_pYieldCurve);
//    calibrator.SetForeignCurve(m_pForeignCurve);
//    calibrator.SetDidividends(m_pDividends);
//    calibrator.SetSpotSharePrice(m_dSpot);

    // Numerical parameters
//    calibrator.SetNbS(m_nNbS);
//    calibrator.SetNbTimesteps(m_nNbTimeStep);
//    calibrator.SetNbRefine(m_nRefine);

    // For convergence testing
//    m_nNbS *= 2;
//    m_nNbTimeStep *= 2;

    // Add the calibration data
    size_t nIdxCall;
    for (nIdxCall = 0; nIdxCall < m_nNbCalls; nIdxCall++)
      calibrator.AddCall(m_piCallMaturities[nIdxCall] * ONEDAY, 
                         m_pdCallStrikes[nIdxCall], 
                         m_pdCallPrices[nIdxCall]);

    // Actually calibrate
    calibrator.Calibrate();
  
    // Save the output from this run
    clock_t t2 = clock();

    m_pdTimes[nIdx] = (t2-t1)/1000.;
    */
  }
}

