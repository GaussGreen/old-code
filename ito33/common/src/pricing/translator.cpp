/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/translator.cpp
// Purpose:     translate between an array of parameter into model params
// Created:     2005/07/04
// RCS-ID:      $Id: translator.cpp,v 1.3 2006/06/01 21:13:51 dave Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file common/src/pricing/translator.h
    @brief translate between an array of parameter into model params
 */

#include "ito33/vector.h"
#include "ito33/debug.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/translator.h"

namespace ito33
{

namespace pricing
{

void Translator::SetFlags(const std::vector<bool>& flags)
{
  ASSERT( flags.size() >= m_nNbFlags );

  for (size_t n = 0; n < m_nNbFlags; n++)
    m_flags[n] = flags[n];
}


size_t Translator::GetNbParameters() const
{
  size_t nNbParams = 0;

  for (size_t n = 0; n < m_nNbFlags; n++)
    if ( m_flags[n] )
      nNbParams++;

  return nNbParams;
}


std::vector<double> Translator::GetParameters() const
{
  std::vector<double> parameters( GetNbParameters() );
  
  GetParameters(&parameters[0]);

  return parameters;
}


std::vector<double> Translator::ApplyBounds(const double* pdX)
{
  // The return vector
  size_t nNbX = GetNbParameters();
  std::vector<double> pdXNew(nNbX);

  // Avoid code duplication by just getting the bounds in vector form. Not
  // the most efficient approach though.
  const std::vector<double> pdLowerBounds = GetLowerBounds();
  const std::vector<double> pdUpperBounds = GetUpperBounds();

  ASSERT_MSG(nNbX == pdLowerBounds.size(), "Size mismatch in ApplyBounds");
  ASSERT_MSG(nNbX == pdUpperBounds.size(), "Size mismatch in ApplyBounds");

  // Enforce the bounds
  for (size_t nIdx = 0; nIdx < nNbX; nIdx++)
  {
    if ( pdX[nIdx] < pdLowerBounds[nIdx] )
      pdXNew[nIdx] = pdLowerBounds[nIdx];
    else if ( pdX[nIdx] > pdUpperBounds[nIdx] )
      pdXNew[nIdx] = pdUpperBounds[nIdx];
    else
      pdXNew[nIdx] = pdX[nIdx];
  }

  return pdXNew;
}


} //namespace pricing

} //namespace ito33
