/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/showconvergence.h
// Purpose:     header file for showing convergence
// Author:      Zhang
// Created:     2004/09/26
// RCS-ID:      $Id: showconvergence.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_TESTS_SHOWCONVERGENCE_H_
#define _IHG_TESTS_SHOWCONVERGENCE_H_

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

#include "ito33/tests/showconvergence.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/finance/modeloutput.h"

namespace ito33
{

namespace ihg
{

#ifdef ITO33_TEST_MODENV 

/** calculate difference and rating of "name" if it has been the theorectical
    value of "name" has been calculated. disply the results if necessary.
  */
#define IHG_SHOW_CONVERGENCE_HANDLE_POSSIBLE(name)  \
  if( pResults[0]->Has ##name() )                    \
     IHG_SHOW_CONVERGENCE_HANDLE(name)


/** calculate difference and rating of "name" and disply the results
    if necessary.
  */
#define IHG_SHOW_CONVERGENCE_HANDLE(name)                                     \
{                                                                             \
  {                                                                           \
    for( n = 0; n < nNbTests; n++)                                            \
    {                                                                         \
      if(n > 0)                                                               \
      {                                                                       \
        pDiffs[n]->Set ##name                                                 \
                    (pResults[n]->Get ##name() - pResults[n-1]->Get ##name());\
      }                                                                       \
      if(n > 1)                                                               \
      {                                                                       \
        if( pDiffs[n]->Get ##name() != 0)                                     \
          pRatings[n]->Set ##name                                             \
                      (pDiffs[n - 1]->Get ##name() / pDiffs[n]->Get ##name());\
        else                                                                  \
          pRatings[n]->Set ##name(10);                                        \
      }                                                                       \
    }                                                                         \
                                                                              \
    if( !bQuiet )                                                             \
    {                                                                         \
    std::cout.precision(15);                                                  \
    std::cout << "___________Convergence_on_"                                 \
              <<  #name  << "__________________________\n";                   \
      for( n = 0; n < nNbTests; n++)                                          \
      {                                                                       \
        if(n == 0)                                                            \
          std::cout << pResults[n]->Get ##name () << "\n";                    \
        else if (n == 1)                                                      \
          std::cout << pResults[n]->Get ##name () << "  "                     \
                    << pDiffs[n]->Get ##name() << std::endl;                  \
        else                                                                  \
          std::cout << pResults[n]->Get ##name () << "  "                     \
                    << pDiffs[n]->Get ##name()    << "  "                     \
                    << pRatings[n]->Get ##name() << std::endl;                \
      }                                                                       \
    }                                                                         \
  }                                                                           \
} 


#endif // #ifdef ITO33_TEST_MODENV 

/**
   A helper fucntion to show the convergence of numeric result.
   
   it will do nothing if the macro ITO33_TEST_MODENV is not defined so that 
   the code to call it can always compile.

   Although it doesn't make much sense to inline this function, it's just
   a convenient way to avoid having add one more source file manually to a
   projet. It would be included inside a single test file anyway.
 */
inline void 
ShowConvergence(const ihg::TheoreticalModel& tm, 
                const finance::Derivative& deriv,
                std::vector< shared_ptr<finance::ModelOutput> >& pResults,
                std::vector< shared_ptr<finance::ModelOutput> >& pDiffs,
                std::vector< shared_ptr<finance::ModelOutput> >& pRatings,
                bool bQuiet = true
                )
{
#ifdef ITO33_TEST_MODENV 

  InitParams numParams;
 
  size_t nNbTimes = numParams.GetNbTimes();
  size_t nNumberSpots = numParams.GetNbSpots();
  double dX = numParams.GetDeltaX() * nNumberSpots;

  size_t n;
  size_t nNbTests = pResults.size();

  pDiffs.resize(nNbTests);
  pRatings.resize(nNbTests);

  for (n = 0; n < nNbTests; n++, nNbTimes *=2 , nNumberSpots*=2 )
  {
    numeric::NumParamsReference::SetNbTimeStepsFor5Years(nNbTimes);
    numeric::NumParamsReference::SetMinNbSpaceSteps(nNumberSpots);
    numeric::NumParamsReference::SetMaxDeltaLogS(dX / nNumberSpots);

    pResults[n] = tm.Compute(deriv);
    pDiffs[n] = make_ptr( new finance::ModelOutput );
    pRatings[n] = make_ptr( new finance::ModelOutput );
  }

  IHG_SHOW_CONVERGENCE_HANDLE(Price);
  IHG_SHOW_CONVERGENCE_HANDLE(Delta);
  IHG_SHOW_CONVERGENCE_HANDLE(Gamma);
  IHG_SHOW_CONVERGENCE_HANDLE(Theta);
  IHG_SHOW_CONVERGENCE_HANDLE_POSSIBLE(Vega);
  IHG_SHOW_CONVERGENCE_HANDLE_POSSIBLE(Rho);
  IHG_SHOW_CONVERGENCE_HANDLE_POSSIBLE(Fugit);

#endif // #ifdef ITO33_TEST_MODENV 

}

} // namespace ihg

} // namespace ito33

#endif // _IHG_TESTS_SHOWCONVERGENCE_H_
