/////////////////////////////////////////////////////////////////////////////
// Name:        test/ycurve/main.cpp
// Purpose:     implementation file of RWLock test program
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: testycurve.cpp,v 1.8 2006/08/19 23:22:41 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"

#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_swap.h"

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testycurve.h"

using namespace ito33;
using namespace ito33::finance;


// ----------------------------------------------------------------------------
// YCurveFlatTestCase test classes
// ----------------------------------------------------------------------------

void YCurveFlatTestCase::GetFlatRate()
{
  YieldCurveFlat yc(0.07);

  ITO33_ASSERT_DOUBLES_EQUAL( yc.GetRate(), 0.07 );
}

void YCurveFlatTestCase::NegativeRate()
{
  YieldCurveFlat yc(- 0.07);
}

// ----------------------------------------------------------------------------
// YCurveAnnuallyCompoundedTestCase test classes
// ----------------------------------------------------------------------------

void YCurveAnnuallyCompoundedTestCase::AddLegDuplicate()
{
  // check that adding duplicates doesn't work
  m_yc.AddLeg(1, 0.04);
  m_yc.AddLeg(1, 0.08);
  m_yc.Validate();
}

void YCurveAnnuallyCompoundedTestCase::AddInvalidLeg()
{
  // check that creating an invalid leg object doesn't work
  m_yc.AddLeg(10, -0.2);
}

void YCurveAnnuallyCompoundedTestCase::Perturb()
{
  m_yc.AddLeg(1, 0.02);
  m_yc.AddLeg(10, 0.02);
  m_yc.AddLeg(30, 0.025);
  m_yc.AddLeg(70, 0.03);
  m_yc.AddLeg(500, 0.05);

  double pdDays[] = {0, 0.1, 0.2, 0.5, 0.7, 0.8, 1, 2};
  size_t nNb = SIZEOF(pdDays);

  Array<double>
    pdRatesOld(nNb),
    pdRatesNew(nNb);
  double
    dShift = 0.001;

  m_yc.GetZeroRates(pdDays, pdRatesOld.Get(), nNb);

  m_yc.Perturb(dShift)->GetZeroRates(pdDays, pdRatesNew.Get(), nNb);

  for ( size_t n = 0; n < nNb; ++n )
    ITO33_ASSERT_DOUBLES_EQUAL( pdRatesOld[n] + dShift, pdRatesNew[n] );

}


void YCurveAnnuallyCompoundedTestCase::Clear()
{
  m_yc.AddLeg(0, 0.04);
  m_yc.AddLeg(30, 0.04);

  CPPUNIT_ASSERT( m_yc.GetAll().size() == 2 );

  m_yc.Clear();

  // shouldn't throw
  m_yc.AddLeg(0, 0.04);
}

void YCurveAnnuallyCompoundedTestCase::NoLeg()
{
  m_yc.AddLeg(0, 0.04);
  m_yc.AddLeg(30, 0.04);

  m_yc.Clear();

  m_yc.GetAll();
}

void YCurveAnnuallyCompoundedTestCase::SingleLeg()
{
  m_yc.AddLeg(30, 0.04);

  m_yc.GetAll();
}

void YCurveAnnuallyCompoundedTestCase::GetLegs()
{
  static const int days[] = { 7, 30, 365 };
  static const double rates[] = { 0.1, 0.2, 0.3 };

  static const size_t nLegs = SIZEOF(days);

  for(size_t n = 0; n < nLegs; n++)
    m_yc.AddLeg(days[n], rates[n]);

  const YieldCurveAnnuallyCompounded::Legs& legs = m_yc.GetAll();

  CPPUNIT_ASSERT( legs.size() == nLegs);

  YieldCurveAnnuallyCompounded::Legs::const_iterator i;

  for ( i = legs.begin(), n = 0; i != legs.end(); ++i, ++n )
  {
    ITO33_ASSERT_DOUBLES_EQUAL( days[n], i->GetDuration() );
    ITO33_ASSERT_DOUBLES_EQUAL( rates[n], i->GetRate() );
  }
}

void YCurveAnnuallyCompoundedTestCase::AddLegsReverse()
{
  // check that the legs are sorted in the right order
  m_yc.Clear();

  m_yc.AddLeg(30, 0.04);
  m_yc.AddLeg(0, 0.04);

  const YieldCurveAnnuallyCompounded::Legs& legs = m_yc.GetAll();

  ITO33_ASSERT_DOUBLES_EQUAL( 0, legs.begin()->GetDuration() );
  ITO33_ASSERT_DOUBLES_EQUAL( 30.0, legs.rbegin()->GetDuration() );
}


// ----------------------------------------------------------------------------
// YCurveSwapTestCase test classes
// ----------------------------------------------------------------------------

void YCurveSwapTestCase::AddCashDuplicate()
{
  // check that adding duplicates doesn't work
  CashRates cashRates;
  cashRates.AddLeg(0.05, 1, TimeUnit_Month);
  cashRates.AddLeg(0.05, 1, TimeUnit_Day);
  cashRates.AddLeg(0.05, 1, TimeUnit_Month);
  m_yc.SetCashRates(cashRates);
  m_yc.Validate();
}

void YCurveSwapTestCase::AddSwapDuplicate()
{
  // check that adding duplicates doesn't work
  Frequency f = Frequency_Annual;

  SwapRates swapRates;
  swapRates.AddLeg(0.05, 1, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 3, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 2, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 1, TimeUnit_Year, f);
  m_yc.SetSwapRates(swapRates);
  m_yc.Validate();
}

void YCurveSwapTestCase::AddSwapTermBeforeCashTerm()
{
  CashRates cashRates;
  cashRates.AddLeg(0.05, 20, TimeUnit_Month);
  m_yc.SetCashRates(cashRates);

  Frequency f = Frequency_Annual;
  SwapRates swapRates;
  swapRates.AddLeg(0.05, 1, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 3, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 2, TimeUnit_Year, f);
  swapRates.AddLeg(0.05, 1, TimeUnit_Year, f);
  m_yc.SetSwapRates(swapRates);
  m_yc.Validate();
}

void YCurveSwapTestCase::Clear()
{
  CashRates cashRates;
  cashRates.AddLeg(0.05, 2, TimeUnit_Month);
  m_yc.SetCashRates(cashRates);

  Frequency f = Frequency_Annual;
  SwapRates swapRates;
  swapRates.AddLeg(0.06, 3, TimeUnit_Year, f);
  swapRates.AddLeg(0.07, 1, TimeUnit_Year, f);
  m_yc.SetSwapRates(swapRates);

  m_yc.Validate();

  m_yc.Clear();

  m_yc.Validate();
}

void YCurveSwapTestCase::Perturb()
{
  CashRates cashRates;
  cashRates.AddLeg(0.05, 2, TimeUnit_Month);
  m_yc.SetCashRates(cashRates);

  Frequency f = Frequency_Annual;
  SwapRates swapRates;
  swapRates.AddLeg(0.06, 3, TimeUnit_Year, f);
  swapRates.AddLeg(0.07, 1, TimeUnit_Year, f);
  m_yc.SetSwapRates(swapRates);

  m_yc.Validate();

  double dShift = 0.0001;
  shared_ptr<YieldCurve> pNewYcTmp( m_yc.Perturb(dShift) );

  shared_ptr<YieldCurveSwap> 
    pNewYc( dynamic_pointer_cast<YieldCurveSwap>(pNewYcTmp) );

  const std::vector<SwapRate>& pSwapRates(pNewYc->GetSwapRates().GetAll());

  ITO33_ASSERT_DOUBLES_EQUAL( 0.06 + dShift, pSwapRates[0].GetRate() );
  ITO33_ASSERT_DOUBLES_EQUAL( 0.07 + dShift, pSwapRates[1].GetRate() );


  const std::vector<CashRate>& pCashRates(pNewYc->GetCashRates().GetAll());

  ITO33_ASSERT_DOUBLES_EQUAL( 0.05 + dShift, pCashRates[0].GetRate() );
  CPPUNIT_ASSERT_EQUAL( 2, int(pCashRates[0].GetMaturityDuration()) );

}

