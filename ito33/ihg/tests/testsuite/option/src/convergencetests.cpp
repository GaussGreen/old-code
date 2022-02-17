/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/optionconvergencetests.cpp
// Purpose:     test that the numerical convergence property are valid
// Author:      Ito 33 Canada
// Created:     2005/06/13
// RCS-ID:      $Id: convergencetests.cpp,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////



#include "ito33/beforestd.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/derivative.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "convergencetests.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/showconvergence.h"

namespace ito33 
{

namespace ihg 
{

namespace test
{

bool DerivativeCheckConvergence(shared_ptr<TheoreticalModel> pModel, 
                      shared_ptr<finance::Derivative> pDeriv, TestParam testParam,
                      ito33::XML::Tag &tag, size_t nNbTests)
{
  std::vector< shared_ptr<finance::ModelOutput> > pResults(nNbTests);
  std::vector< shared_ptr<finance::ModelOutput> > pDiffs;
  std::vector< shared_ptr<finance::ModelOutput> > pRatings;

  ihg::ShowConvergence(*pModel, *pDeriv, pResults, pDiffs, pRatings);

  shared_ptr<finance::ModelOutput> pOutput = pResults[0];

  std::vector<double> pdPrices(nNbTests, 0.);
  std::vector<double> pdDiffs(nNbTests, 0.);
  std::vector<double> pdRatings(nNbTests, 0.);

  size_t n;

  for(n = 0; n < nNbTests; n++)
  {
    pdPrices[n] = pResults[n]->GetPrice();
    if(n > 0)
      pdDiffs[n] = pDiffs[n]->GetPrice();
    if(n > 1)
      pdRatings[n] = pRatings[n]->GetPrice();
  }        
    
  XML::Tag localtest("test", tag);
  localtest.Element("name")("convergence_test");
  localtest.Element("file")(testParam.GetFileName());
  localtest.Element("result")("pass");

  ito33::XML::ConvergenceParameterValue paramvalue;
   paramvalue.init("price","difference","ratio");

  for( n = 0; n < nNbTests; n++)
  {
    paramvalue.m_prices.push_back(  pdPrices[n] );
    paramvalue.m_diffs.push_back( pdDiffs[n] ) ;
    paramvalue.m_ratings.push_back( pdRatings[n] ) ;
  }

  //output details
  localtest.Element("summary",paramvalue);

  return true;
}

} //end test namespace
} //end ihg namespace
}//end ito33 namespace
