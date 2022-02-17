/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedhr.cpp
// Purpose:     Paramaterized hazard rate surface
// Author:      David 
// Created:     04/01/26
// RCS-ID:      $Id: parameterizedhr.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ihg/parameterizedhr.h"

namespace ito33
{

namespace ihg
{

void ParameterizedHR::WriteGnuPlot(std::ostream &out) 
{
  const size_t nSize = 17;
  double pdS[nSize] = {40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 
    115, 120, 125, 130, 160};
  double pdHRs[nSize];

  std::list<double>::iterator iterTimes;
  for (iterTimes = m_listOfTimes.begin();
       iterTimes != m_listOfTimes.end();
       ++iterTimes)
  {

    GetHazardRates(*iterTimes, pdS, pdHRs, nSize);

    for (size_t nIdx = 0; nIdx < nSize; nIdx++)
    {
      out << pdS[nIdx] << " " << *iterTimes << " " << pdHRs[nIdx] << std::endl;
    }
    out << std::endl;

    std::list<double>::iterator iterTmp = iterTimes;
    ++iterTmp;
    if (iterTmp != m_listOfTimes.end())
    {
      double dTime = 0.5 * (*iterTmp + *iterTimes);
      GetHazardRates(dTime, pdS, pdHRs, nSize); 

      for (size_t nIdx = 0; nIdx < nSize; nIdx++)
      {
        out << pdS[nIdx] << " " << dTime << " " << pdHRs[nIdx] << std::endl;
      }
      out << std::endl;
    }

  }

}

 
} // namespace ihg

} // namespace ito33
