/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/tests/showconvergence.h
// Purpose:     header file for showing convergence
// Author:      Wang
// Created:     2004/09/13
// RCS-ID:      $Id: showconvergence.h,v 1.12 2005/06/27 15:09:43 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TESTS_SHOWCONVERGENCE_H_
#define _ITO33_TESTS_SHOWCONVERGENCE_H_

// has to check this macro so the setters on NumParamsReference are defined
#ifdef ITO33_TEST_MODENV 

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/array.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"

#endif // #ifdef ITO33_TEST_MODENV 

#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

namespace ito33
{

#ifdef ITO33_TEST_MODENV 

class InitParams
{
public:
  
  InitParams()
  {
    nNbTimesInit = numeric::NumParamsReference::GetNbTimeStepsFor5Years();
    nNumberSpotsInit = numeric::NumParamsReference::GetMinNbSpaceSteps();
    dDeltaXInit = numeric::NumParamsReference::GetMaxDeltaLogS();
    schemeTypeInit = numeric::NumParamsReference::GetSchemeType();
  }

  ~InitParams()
  {
  
    numeric::NumParamsReference::SetNbTimeStepsFor5Years(nNbTimesInit);
      numeric::NumParamsReference::SetMinNbSpaceSteps(nNumberSpotsInit);
      numeric::NumParamsReference::SetMaxDeltaLogS(dDeltaXInit);
      numeric::NumParamsReference::SetSchemeType(schemeTypeInit);
  }

  size_t GetNbTimes() { return nNbTimesInit; }
  size_t GetNbSpots() { return nNumberSpotsInit; }
  double GetDeltaX() { return dDeltaXInit; }

private:

  size_t nNbTimesInit;
  size_t nNumberSpotsInit;
  double dDeltaXInit;
  numeric::SchemeType schemeTypeInit;
};

#endif // #ifdef ITO33_TEST_MODENV 

/**
   A helper fucntion to show the convergence of numeric result.
   
   it will do nothing if the macro ITO33_TEST_MODENV is not defined so that 
   the code to call it can always compile.

   Although it doesn't make much sense to inline this function, it's just
   a convenient way to avoid having to add one more source file manually to a
   projet. It would be included inside a single test file anyway.
 */
inline void 
ShowConvergence(const finance::TheoreticalModel& tm, 
                const finance::Derivative& deriv,
                std::vector<double>& pdPrices,
                std::vector<double>& pdDiffs,
                std::vector<double>& pdRatings,
                bool bQuiet = true
                )
{
#ifdef ITO33_TEST_MODENV 

  InitParams numParams;

  size_t nNbTimes = numParams.GetNbTimes();
  size_t nNumberSpots = numParams.GetNbSpots();
  double dX = numParams.GetDeltaX() * nNumberSpots;

  size_t n;
  size_t nNbTests = pdPrices.size();

  pdDiffs.resize(nNbTests);
  pdRatings.resize(nNbTests);

  for (n = 0; n < nNbTests; n++, nNbTimes *=2 , nNumberSpots*=2 )
  {
    numeric::NumParamsReference::SetNbTimeStepsFor5Years(nNbTimes);
    numeric::NumParamsReference::SetMinNbSpaceSteps(nNumberSpots);
    numeric::NumParamsReference::SetMaxDeltaLogS(dX / nNumberSpots);

    pdPrices[n] = tm.Compute(deriv)->GetPrice();

    if(n > 0)
      pdDiffs[n] = pdPrices[n] - pdPrices[n-1];
    if(n > 1)
      pdRatings[n] = pdDiffs[n - 1] / pdDiffs[n];
    
    if(!bQuiet)
    {
      if (n == 0)
      std::cout << pdPrices[n] << std::endl;
      else if (n == 1)
        std::cout << pdPrices[n] << "  " << pdDiffs[n] << std::endl;
      else
        std::cout << pdPrices[n] << "  " << pdDiffs[n] << "  "
                  << pdRatings[n] << std::endl;
    }
  }

#endif // #ifdef ITO33_TEST_MODENV 

}

/**
   A helper fucntion to show the convergence of numeric result.
   
   it will do nothing if the macro ITO33_TEST_MODENV is not defined so that 
   the code to call it can always compile.

   Although it doesn't make much sense to inline this fucntion, it's just
   a convenient way to avoid having add manually one more source file to a
   projet. It would be included inside a single test file anyway.
 */
inline void 
ShowConvergence(const finance::TheoreticalModel& tm, 
                const finance::Derivative& deriv,
                size_t nNbTests = 5)
{
  std::cout.precision(15);

  std::cout << "\n\n\nconvergence test.\n\n";

  std::vector<double> pdPrices(nNbTests);
  std::vector<double> pdDiffs;
  std::vector<double> pdRatings;

  const bool bNotQuiet = false;
  ShowConvergence(tm, deriv, pdPrices, pdDiffs, pdRatings, bNotQuiet);
}


} // namespace ito33

#endif // _ITO33_TESTS_SHOWCONVERGENCE_H_
