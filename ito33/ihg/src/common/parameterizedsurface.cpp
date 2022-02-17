/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedsurface.cpp
// Purpose:     Paramterized surface class
// Author:      David
// Created:     04.01.24
// RCS-ID:      $Id: parameterizedsurface.cpp,v 1.4 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parameterizedsurface.cpp
*/

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include <math.h>

#include "ihg/parameterizedsurface.h"


namespace ito33
{

namespace ihg
{

void ParameterizedSurface::Write(std::ostream& out) const
{

  // Just write the number of parameters, the number of times, and 
  // then the data
  out << m_nNbParams << std::endl;
  out << m_listOfTimes.size() << std::endl;

  std::list< Array<double> >::iterator iterParams;
  std::list< double >::iterator iterTimes;

  for (iterTimes = m_listOfTimes.begin(), iterParams = m_listOfParams.begin();
    iterTimes != m_listOfTimes.end();
    ++iterTimes, ++iterParams)
  {
    out << *iterTimes << std::endl;
    for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
      out << (*iterParams)[nIdx] << " ";
    out << std::endl;
  }
}


void ParameterizedSurface::SetParams(double dTime, double* pdParamValues)
{

  // Construct an array for the new data
  Array<double> pdData(m_nNbParams);
  for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
    pdData[nIdx] = pdParamValues[nIdx];

  // Special case of empty list
  if (m_listOfTimes.size() == 0)
  {
    m_listOfTimes.push_back(dTime);
    m_listOfParams.push_back(pdData);

    // Set the interval iterators to valid data (it doesn't matter
    // what they are, so long as they can be dereferenced)
    m_iterParamsLeft = m_listOfParams.begin();
    m_iterParamsRight = m_listOfParams.begin();
    m_iterTimesLeft = m_listOfTimes.begin();
    m_iterTimesRight = m_listOfTimes.begin();
    return;
  }

  // Add to list, but first check for duplication
  std::list<double>::iterator iterTimes;
  std::list< Array<double> >::iterator iterParams;
  for (iterTimes = m_listOfTimes.begin(), iterParams = m_listOfParams.begin();
       iterTimes != m_listOfTimes.end();
       ++iterTimes, ++iterParams)
  {
    // Duplicate time. Replace the data
    if ( fabs(dTime - *iterTimes) < 1.e-14 )
    {
      *iterParams = pdData;
      break;
    } 
    else if ( dTime < *iterTimes)
    {
      m_listOfTimes.insert(iterTimes, dTime);
      m_listOfParams.insert(iterParams, pdData);
      break;
    }

  }

  // Check if we need to add to the end
  if ( iterTimes == m_listOfTimes.end() )
  {
    m_listOfTimes.push_back(dTime);
    m_listOfParams.push_back(pdData);
  }

}

void ParameterizedSurface::FindInterval(double dTime) const
{

  ASSERT_MSG(m_listOfTimes.size() >= 2, "Need more data in paramvolsurface");

  // Check if we are alreay in the right interval.  Values are initialized
  // in SetParams
  if (*m_iterTimesLeft <= dTime && *m_iterTimesRight >= dTime)
    return;

  // Find the place in the list so that TimeLeft <= dTime <= TimeRight
  m_iterTimesRight = m_listOfTimes.begin();
  m_iterParamsRight = m_listOfParams.begin();
  
  while (m_iterTimesRight != m_listOfTimes.end() && dTime > *m_iterTimesRight)
  {
    ++m_iterTimesRight;
    ++m_iterParamsRight;
  }

  // Check boundary conditions, and reset pointers
  if (m_iterTimesRight == m_listOfTimes.begin())
  {
    // Left boundary
    m_iterTimesLeft = m_iterTimesRight;
    ++m_iterTimesRight;

    m_iterParamsLeft = m_iterParamsRight;
    ++m_iterParamsRight;
  }
  else if (m_iterTimesRight == m_listOfTimes.end())
  {
    // Right boundary
    --m_iterTimesRight;
    m_iterTimesLeft = m_iterTimesRight;
    --m_iterTimesLeft;

    --m_iterParamsRight;
    m_iterParamsLeft = m_iterParamsRight;
    --m_iterParamsLeft;
  }
  else
  {
    // Interior
    m_iterTimesLeft = m_iterTimesRight;
    --m_iterTimesLeft;

    m_iterParamsLeft = m_iterParamsRight;
    --m_iterParamsLeft;
  }
}

}  // namespace ihg

} // namespace ito33
