/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/timestruct.cpp
// Purpose:     implementation of time array classes
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: timestruct.cpp,v 1.8 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/timestruct.h"
#include "ito33/autoptr.h"
#include "ito33/arraycheckers.h"

namespace ito33
{

// implement autoptr
ITO33_IMPLEMENT_AUTOPTR(numeric::TimeStructConstant);


namespace numeric
{

TimeStruct::TimeStruct(const double *pdTimes, size_t nNumber)
{
  CheckIncreasingOrder(pdTimes, nNumber);
  m_nNumber = nNumber;
  Array<double> pTmp(nNumber);
  m_pdTimes = pTmp;
  for(size_t n = 0; n < nNumber; n++)
    m_pdTimes[n] = pdTimes[n];
}

} // namespace numeric

} // namespace ito33
