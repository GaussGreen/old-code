/////////////////////////////////////////////////////////////////////////////
// Name:        testoption.cpp
// Purpose:     Acceptance test for option 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testparbond.cpp,v 1.3 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>
#include "ito33/sharedptr.h"

#include "ito33/finance/parbond.h"
#include "ito33/finance/frequency.h"
#include "ito33/finance/cashflowstream.h"

#include "ito33/tests/testparbond.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Conversion schedule reset
// ----------------------------------------------------------------------------

static Date g_valuationDate = Date("2003/02/01");
static double g_dYTM = 0.03;
static double g_dSpread = 0.005;
static Frequency g_frequency = Frequency_Quarterly;
static Date::DayCountConvention g_dcc = Date::DayCountConvention_Act365;
static size_t g_nMaturity = 36;
static double g_dRecoveryRate = 0.2;

void ParBondTest::BadYTM()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               -0.01,
               g_dSpread,
               g_frequency,
               g_dcc,
               g_dRecoveryRate);                                 
}

void ParBondTest::BadSpread()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               g_dYTM,
               -0.004,
               g_frequency,
               g_dcc,
               g_dRecoveryRate);
                                 
} 

void ParBondTest::BadRecoveryRate()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               g_dYTM,
               g_dSpread,
               g_frequency,
               g_dcc,
               1.1);
                                 
} 

void ParBondTest::BadFrequency()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               g_dYTM,
               g_dSpread,
               Frequency(7),
               g_dcc,
               g_dRecoveryRate);
                                 
} 

void ParBondTest::BadDayCountConvention()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               g_dYTM,
               g_dSpread,
               g_frequency,
               Date::DayCountConvention(-1),
               g_dRecoveryRate);
                                 
} 

void ParBondTest::BadFrequencyMaturity()
{
  ParBond bond(g_valuationDate, 26,
               g_dYTM,
               g_dSpread,
               g_frequency,
               g_dcc,
               g_dRecoveryRate);
                                 
} 

void ParBondTest::BadMaturity()
{
  ParBond bond(g_valuationDate, 0,
               g_dYTM,
               g_dSpread,
               g_frequency,
               g_dcc,
               g_dRecoveryRate);
                                 
} 

void ParBondTest::GetCouponTrivial()
{
  ParBond bond(g_valuationDate, 12,
               g_dYTM,
               g_dSpread,
               Frequency_Annual,
               g_dcc,
               g_dRecoveryRate);
  
//  std::cout << bond.GetCouponRate() << std::endl;

  // the coupon can't be a little different due to day count convention
  CPPUNIT_ASSERT_DOUBLES_EQUAL(g_dYTM + g_dSpread,
                               bond.GetCouponRate(),
                               1.e-4);

}

void ParBondTest::GetCouponAndCashFlowStream()
{
  ParBond bond(g_valuationDate, g_nMaturity,
               g_dYTM,
               g_dSpread,
               g_frequency,
               g_dcc,
               g_dRecoveryRate);
  
  shared_ptr<CashFlowStream> pcfs = bond.GetCashFlowStream();

  double
    dCoupon = bond.GetCouponRate() / g_frequency,
    dSum = 0,
    dDiscount = 0;

  CashFlowStream::const_iterator iter = pcfs->begin();
  for(;iter != pcfs->end(); iter++)
  {
    dDiscount = pow(1 + (g_dYTM + g_dSpread) / g_frequency,
                    -g_frequency * 
                    Date::YearsDiff(bond.GetContractingDate(),
                                    iter->first,
                                    g_dcc));
    dSum += dDiscount * dCoupon;
  }
  dSum += dDiscount;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(dSum, 1, 1.e-6);
}
