/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/varianceswap.cpp
// Purpose:     Implementation of variance swap tests
// Created:     2006/03/08
// RCS-ID:      $Id: varianceswap.cpp,v 1.14 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/swaptype.h"
#include "ito33/finance/returntype.h"
#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/hg/theoreticalmodel.h"

#include "hg/tests/varianceswap.h"

using namespace ito33::finance;

namespace ito33
{

namespace hg
{

namespace VarianceSwapTest
{

// Global objects to work around CppUnit limitation
shared_ptr<TheoreticalModel> m_pModel;
shared_ptr<finance::VarianceSwap> m_pVarianceSwap;


void Setup(shared_ptr<TheoreticalModel> pModel,
           shared_ptr<finance::VarianceSwap> pVarianceSwap)
{
  m_pModel = pModel->Clone();
  m_pVarianceSwap = pVarianceSwap;
}


void VarianceSwapTest::BasicVarianceSwap()
{

  shared_ptr<finance::SessionData>
    pSessionData = m_pVarianceSwap->GetSessionData();

  // Create a new variance swap from the original (force variance)
  double dVolatilityStrike = m_pVarianceSwap->GetVolatilityStrike();
  ReturnType returnType = m_pVarianceSwap->GetReturnType();
  Date startOfSamplingPeriod = m_pVarianceSwap->GetStartOfSamplingPeriod();
  Date maturityDate = m_pVarianceSwap->GetMaturityDate();
  SwapType swapType = Swap_Variance;
  size_t nNbReturns = m_pVarianceSwap->GetNbSamplingReturns();

  shared_ptr<VarianceSwap> pVarianceSwap( 
    new VarianceSwap( dVolatilityStrike, maturityDate, swapType, 
                      startOfSamplingPeriod, nNbReturns) );

  pVarianceSwap->SetSessionData(pSessionData);
  pVarianceSwap->SetReturnType( returnType );

  // Price
  shared_ptr<finance::ModelOutput> pOutput(m_pModel->Compute(*pVarianceSwap));
  double dPrice = pOutput->GetPrice();

  // Compare to expected price
  // TODO: Get exected price from xml file, instead of hard-coding
  // Price found after 4 levels of refinement. Should be good for 2 digits.
  RelErrorCheck(dPrice,  0.00105192111860183, 3.e-2);

}


void VarianceSwapTest::BasicVolatilitySwap()
{

  shared_ptr<finance::SessionData>
    pSessionData = m_pVarianceSwap->GetSessionData();

  // Create a new volatility swap from the original (force volatility). 
  ReturnType returnType = m_pVarianceSwap->GetReturnType();
  Date startOfSamplingPeriod = m_pVarianceSwap->GetStartOfSamplingPeriod();
  Date maturityDate = m_pVarianceSwap->GetMaturityDate();
  size_t nNbReturns = m_pVarianceSwap->GetNbSamplingReturns();

  double dVolatilityStrike =  m_pVarianceSwap->GetVolatilityStrike();
  SwapType swapType = Swap_Volatility;

  shared_ptr<VarianceSwap> pVarianceSwap( 
    new VarianceSwap(dVolatilityStrike, maturityDate, swapType, 
                     startOfSamplingPeriod, nNbReturns) );

  pVarianceSwap->SetSessionData(pSessionData);
  pVarianceSwap->SetReturnType(returnType);

  // Price
  shared_ptr<finance::ModelOutput> pOutput(m_pModel->Compute(*pVarianceSwap));
  double dPrice = pOutput->GetPrice();

  // Compare to expected price
  // TODO: Get exected price from xml file instead of hard-coding
  // Price found after 4 levels of refinement. Should be good for 2 digits.  
  //RelErrorCheck( dPrice, 0.00137395566120967);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dPrice, 0.00137395566120967, 0.0005);

}


void VarianceSwapTest::ChangeParams()
{  

  // Save original params for later access
  double dVolatilityStrike   = m_pVarianceSwap->GetVolatilityStrike();
  SwapType swapType          = m_pVarianceSwap->GetSwapType();
  ReturnType returnType      = m_pVarianceSwap->GetReturnType();
  Date startOfSamplingPeriod = m_pVarianceSwap->GetStartOfSamplingPeriod();
  Date maturityDate          = m_pVarianceSwap->GetMaturityDate();
  size_t nNbReturns          = m_pVarianceSwap->GetNbSamplingReturns();

  // Price original contract
  shared_ptr<finance::ModelOutput> pOutput(m_pModel->Compute(*m_pVarianceSwap));
  double dPrice = pOutput->GetPrice();

  // Change the variance swap params

  // Price should scale with the nominal 
  double dNominal = 10.0;  
  shared_ptr<VarianceSwap> pVarianceSwap( 
    new VarianceSwap( dVolatilityStrike, maturityDate, swapType, 
                      startOfSamplingPeriod, nNbReturns) );
  pVarianceSwap->SetSessionData( m_pVarianceSwap->GetSessionData() );
  pVarianceSwap->SetReturnType( returnType );

  pOutput = m_pModel->Compute(*pVarianceSwap);
  double dPrice2 = dNominal * pOutput->GetPrice();

  RelErrorCheck( dPrice2 / dPrice, dNominal);

  // Price should decrease as strike increases
  double dStrikeNew = 2.0 * dVolatilityStrike;
  pVarianceSwap = new VarianceSwap(dStrikeNew, maturityDate, swapType,
                                   startOfSamplingPeriod, nNbReturns);
  pVarianceSwap->SetSessionData( m_pVarianceSwap->GetSessionData() );
  pVarianceSwap->SetReturnType( returnType );

  pOutput = m_pModel->Compute(*pVarianceSwap);
  dPrice2 = pOutput->GetPrice();

  CPPUNIT_ASSERT(dPrice2 < dPrice);

  // Price should increase as maturity increase, although the effect
  // is small.  Expect long term average to be hit quickly.
  Date maturityNew = maturityDate;
  maturityNew.AddMonths(1);
  pVarianceSwap = new VarianceSwap(dVolatilityStrike, maturityNew, swapType, 
                                   startOfSamplingPeriod, nNbReturns);
  pVarianceSwap->SetSessionData( m_pVarianceSwap->GetSessionData() );
  pVarianceSwap->SetReturnType( returnType );

  pOutput = m_pModel->Compute(*pVarianceSwap);
  dPrice2 = pOutput->GetPrice();

  CPPUNIT_ASSERT(dPrice2 > dPrice);

  // Price should increase as the start of the sampling period increases 
  // (the averaging period decreases, so less time to settle to long term
  // average.  Hence, the expected variance is higher)
  Date firstNew = startOfSamplingPeriod;
  firstNew.AddMonths(1);
  pVarianceSwap = new VarianceSwap( dVolatilityStrike, maturityDate, swapType,
                                    firstNew, nNbReturns);
  pVarianceSwap->SetSessionData( m_pVarianceSwap->GetSessionData() );
  pVarianceSwap->SetReturnType( returnType );

  pOutput = m_pModel->Compute(*pVarianceSwap);
  dPrice2 = pOutput->GetPrice();

  CPPUNIT_ASSERT(dPrice2 > dPrice);

  // Price should decrease with log returns.  The effect is small,
  // and could depend on parameters.  But this should be a safe
  // test for the default contract.
  CPPUNIT_ASSERT( returnType == Return_Actual );
  ReturnType returnTypeNew = Return_Log;

  pVarianceSwap = new VarianceSwap(dVolatilityStrike, maturityDate, swapType,
                                   startOfSamplingPeriod, nNbReturns);
  pVarianceSwap->SetSessionData( m_pVarianceSwap->GetSessionData() );
  pVarianceSwap->SetReturnType( returnTypeNew );

  pOutput = m_pModel->Compute(*pVarianceSwap);
  dPrice2 = pOutput->GetPrice();

  CPPUNIT_ASSERT(dPrice2 <= dPrice);
  CPPUNIT_ASSERT( fabs( (dPrice2 - dPrice) / dPrice) < 0.1);

}


void VarianceSwapTest::SetCurrentValues()
{  

  // Get the original session. Valuation date will be changed
  shared_ptr<SessionData> pSessionData = m_pVarianceSwap->GetSessionData();

  // Create a new variance swap from the original
  double dStrike = m_pVarianceSwap->GetVolatilityStrike();
  ReturnType returnType = m_pVarianceSwap->GetReturnType();
  Date startOfSamplingPeriod = m_pVarianceSwap->GetStartOfSamplingPeriod();
  Date maturityDate = m_pVarianceSwap->GetMaturityDate();
  SwapType swapType = m_pVarianceSwap->GetSwapType();
  size_t nNbReturns = m_pVarianceSwap->GetNbSamplingReturns();

  shared_ptr<VarianceSwap> pVarianceSwap( 
    new VarianceSwap(dStrike, maturityDate, swapType, 
                     startOfSamplingPeriod, nNbReturns) );  
  pVarianceSwap->SetReturnType( returnType );

  // Save original valuation date so it can be reset. Increase the new
  // valuation date
  Date valuationDate = pSessionData->GetValuationDate();
  Date valuationDateNew = valuationDate; 
  valuationDateNew.AddDays(30);   

  // Price original contract
  shared_ptr<finance::ModelOutput> pOutput;
  //pOutput = m_pModel->Compute(*m_pVarianceSwap);
  //double dPrice = pOutput->GetPrice();

  // Set the current volatility to the strike. Assume 21 trading days in
  // 30 calender days
  pVarianceSwap->SetCurrentValues(dStrike, 21);
  pSessionData->SetValuationDate(valuationDateNew);
  pVarianceSwap->SetSessionData(pSessionData);

  pOutput = m_pModel->Compute(*pVarianceSwap);
  double dPriceStrike = pOutput->GetPrice();

  // Set the current volatility less than the strike. Should decrease price.
  pVarianceSwap->SetCurrentValues(dStrike / 1.1, 21);
  pSessionData->SetValuationDate(valuationDateNew);
  pVarianceSwap->SetSessionData(pSessionData);

  pOutput = m_pModel->Compute(*pVarianceSwap);
  double dPriceLow = pOutput->GetPrice();

  // Set current volatility higher than the strike. Should increase price.
  pVarianceSwap->SetCurrentValues(dStrike * 1.1, 21);
  pSessionData->SetValuationDate(valuationDateNew);
  pVarianceSwap->SetSessionData(pSessionData);

  pOutput = m_pModel->Compute(*pVarianceSwap);
  double dPriceHigh = pOutput->GetPrice();

  // Undo the valuation date change
  pSessionData->SetValuationDate(valuationDate);

  // Check results
  CPPUNIT_ASSERT(dPriceHigh > dPriceStrike);
  CPPUNIT_ASSERT(dPriceStrike > dPriceLow);

}

void VarianceSwapTest::Compare2Dand3D()
{  
  // Make sure session has no dividends, so similarity reduction is used
  shared_ptr<SessionData> pSessionData = m_pVarianceSwap->GetSessionData();
  shared_ptr<Dividends> pDividendsOrig = pSessionData->GetEquity()->GetDividends();  
  
  shared_ptr<Dividends> pDividendsNew(new Dividends());  
  pSessionData->GetEquity()->SetDividends(pDividendsNew);

  shared_ptr<finance::ModelOutput> pOutput(m_pModel->Compute(*m_pVarianceSwap));
  double dPrice2D = pOutput->GetPrice();

  // Now add dividends so full 3D pricing is used
  Date dividendDate = pSessionData->GetValuationDate();
  dividendDate.AddDays(1);
  pDividendsNew->AddCash(dividendDate, 1.e-10);
  pSessionData->GetEquity()->SetDividends(pDividendsNew);

  pOutput = m_pModel->Compute(*m_pVarianceSwap);
  double dPrice3D = pOutput->GetPrice();

  // Compare
  CPPUNIT_ASSERT( fabs(dPrice3D - dPrice2D) < 0.1 );
  //RelErrorCheck( dPrice2D, dPrice3D);

  // Undo changes to session data
  pSessionData->GetEquity()->SetDividends(pDividendsOrig);
}


void VarianceSwapTest::RelErrorCheck(double dVal1, double dVal2, 
                                     double dPrecision)
{
  // precision defaults to 1.e-3
  double dRelError = (dVal1 - dVal2) / dVal2;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(dRelError, 0.0, dPrecision);
}

void VarianceSwapTest::ClosedForm()
{
  // Create an underlying process with one regime that doesn't default
  size_t nNbRegimes = 1;
  double dVol = 0.02;
  double dIntensity = 0.;
  std::vector<double> pdVols(nNbRegimes, dVol);

  std::vector<double> pdIntensities(nNbRegimes, dIntensity);
  
  shared_ptr<UnderlyingProcess> 
    pUP( new UnderlyingProcess(nNbRegimes, pdVols, pdIntensities) );

  Jumps jumps;
  jumps.push_back(Jump(0.2, -0.2));

  pUP->SetJumps(0, 0, jumps);

  TheoreticalModel tm(pUP);

  // Create a variance swap
  Date valuationDate(2002, Date::Jan, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  shared_ptr<RateData> pRateData(new RateData);

  double dContinuousRate = 0.;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  pRateData->SetYieldCurve(pCurrency, pyc);

  double dS0 = 100.0;
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  pEquity->SetPreviousSharePrice(dS0);

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  // Create the variance swap
  Date firstObservationDate = valuationDate;
  Date maturityDate = valuationDate;

  // This value here doesn't really matter, according to the closed form
  size_t nNbObservations = 100;
  maturityDate.AddDays(nNbObservations);

  // Just use the same bs vol so that the resulted price should be small
  double dVolatilityStrike = dVol;

  // Close form exists only for variance swap
  SwapType swapType = Swap_Variance;
  
  shared_ptr<VarianceSwap> pVarianceSwap( 
    new VarianceSwap( dVolatilityStrike, maturityDate, swapType, 
                      firstObservationDate, nNbObservations) );
  
  pVarianceSwap->SetSessionData(pSessionData);

  shared_ptr<finance::ModelOutput> pOutput = tm.Compute(*pVarianceSwap);

  double dComputed = pOutput->GetPrice();
  
  double dT = GetDoubleFrom(maturityDate) - GetDoubleFrom(valuationDate);
  double dDeltaT = dT / (nNbObservations);

  const Jumps& internalJumps( pUP->GetJumps(0, 0) );
  Jumps::const_iterator internalJump = internalJumps.begin();

  double dSum1 = 0, dSum2 = 0;
  for ( ; internalJump != internalJumps.end(); internalJump++)
  { 
    const double dIntensity = internalJump->GetIntensity();
    const double dAmplitude = internalJump->GetAmplitude();
    const double dTmp = log(1 + dAmplitude);

    dSum1 += dIntensity * dTmp * dTmp;
    dSum2 += dIntensity * (dTmp - dAmplitude);
  }

  double dTmp = dContinuousRate - 0.5 * dVol * dVol + dSum2;

  // Compute the expected theoretical price
  double 
    dExpected = exp(- dContinuousRate * dT)
              * (     INVERSEONETRADINGDAY * dDeltaT 
                    * (dVol * dVol + dSum1 + dDeltaT * dTmp * dTmp)
                  - dVolatilityStrike * dVolatilityStrike);

  // TODO: add a relative check for doubles
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dExpected, dComputed, 1.e-4);
}

} // namespace VarianceSwapTest

} // namespace hg

} // namespace ito33
