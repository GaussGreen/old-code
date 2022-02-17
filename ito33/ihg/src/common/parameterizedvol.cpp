/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedvol.cpp
// Purpose:     Parameterized volatility surface class
// Author:      David
// Created:     04.01.24
// RCS-ID:      $Id: parameterizedvol.cpp,v 1.12 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parameterizedvol.cpp
*/

#include "ito33/beforestd.h"
#include <math.h>
#include <iostream>
#include "ito33/afterstd.h"

#include "ihg/parameterizedvol.h"

using namespace ito33::ihg;

void ParameterizedVol::GetVolsSquared(double dTime, const double *Spots,
                                      double *Vols, size_t nNbSpots) const
{
  GetVols(dTime, Spots, Vols, nNbSpots);
  for (size_t nIdx = 0; nIdx < nNbSpots; nIdx++)
    Vols[nIdx] *= Vols[nIdx];
}

void ParameterizedVol::WriteGnuPlot(std::ostream &out) 
{
  const size_t nSize = 17;
  double pdS[nSize] = {40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 
    115, 120, 125, 130, 160};
  double pdVols[nSize];

  std::list<double>::iterator iterTimes;
  for (iterTimes = m_listOfTimes.begin();
       iterTimes != m_listOfTimes.end();
       ++iterTimes)
  {

    GetVols(*iterTimes, pdS, pdVols, nSize); 

    for (size_t nIdx = 0; nIdx < nSize; nIdx++)
    {
      out << pdS[nIdx] << " " << *iterTimes << " " << pdVols[nIdx] << std::endl;
    }
    out << std::endl;

    std::list<double>::iterator iterTmp = iterTimes;
    ++iterTmp;
    if (iterTmp != m_listOfTimes.end())
    {
      double dTime = 0.5 * (*iterTmp + *iterTimes);
      GetVols(dTime, pdS, pdVols, nSize); 

      for (size_t nIdx = 0; nIdx < nSize; nIdx++)
      {
        out << pdS[nIdx] << " " << dTime << " " << pdVols[nIdx] << std::endl;
      }
      out << std::endl;
    }

  }

}
