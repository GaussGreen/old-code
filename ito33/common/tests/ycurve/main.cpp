/////////////////////////////////////////////////////////////////////////////
// Name:        test/ycurve/main.cpp
// Purpose:     main file of yield curve test program
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: main.cpp,v 1.16 2005/04/18 14:06:01 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "ito33/cppunit.h"

#include "ito33/tests/testycurve.h"

#include "ito33/finance/yieldcurve_swap.h"

using namespace ito33;
using namespace ito33::finance;

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(YCurveFlatTestCase::suite());
  runner.addTest(YCurveAnnuallyCompoundedTestCase::suite());
  runner.addTest(YCurveSwapTestCase::suite());

  

  return runner.run("") ? 0 : 1;
}

/*
void PrintRates(const YieldCurveSwap& yc)
{ 
  std::vector<Date> pdateTerms;
  size_t nNb = 2000;
  for(size_t i = 0; i < nNb; i++)
  {
    Date t("2005/01/01");
    t.AddDays(i+1);
    pdateTerms.push_back(t);
  }
  std::ofstream sOut("c:\\ito33\\output\\rate\\rate.txt");
  for(size_t n = 0; n < nNb; n++)
  {
    double dFrac = Date::YearsDiff( Date("2005/01/01"),
                                    pdateTerms[n],
                                    Date::DayCountConvention_Act365);

    double dFactD = Date::YearsDiff( Date("2005/01/01"),
                                    pdateTerms[n],
                                    Date::DayCountConvention_30360);

    sOut << dFrac << " " << yc.GetZeroRate(dFrac) << "\n";
  }
}

void testswap()
{
  YieldCurveSwap yc(Date("2000/09/05"));

  yc.AddCashRate(0.05, 1, TimeUnit_Day);
  yc.AddCashRate(0.0502, 7, TimeUnit_Day);
  yc.AddCashRate(0.0513, 1, TimeUnit_Month);
  yc.AddCashRate(0.0523, 3, TimeUnit_Month);
  yc.AddCashRate(0.054, 6, TimeUnit_Month);
  yc.AddCashRate(0.0562, 1, TimeUnit_Year);
  

  yc.AddSwapRate(0.0571, 2, TimeUnit_Year, Frequency_SemiAnnual);
  yc.AddSwapRate(0.058, 3, TimeUnit_Year, Frequency_SemiAnnual);
  yc.AddSwapRate(0.061, 4, TimeUnit_Year, Frequency_SemiAnnual);

  yc.Validate();

  PrintRates(yc);

}
*/