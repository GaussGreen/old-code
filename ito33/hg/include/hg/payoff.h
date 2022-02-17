/////////////////////////////////////////////////////////////////////////////
// Name:        hg/payoff.h
// Purpose:     Class holds the prices using HG model at a given time
// Created:     2005/06/13
// RCS-ID:      $Id: payoff.h,v 1.1 2005/06/13 17:56:23 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/payoff.h
    @brief Class holds the prices using HG model at a given time
 */

#ifndef _HG_PAYOFF_H_
#define _HG_PAYOFF_H_

#include "ito33/array.h"
#include "ito33/numeric/interpolation.h"


namespace ito33 
{

namespace hg 
{

/**
   Class holds the prices using HG model at a given time.
 */
class Payoff
{

public:

  /**
     Constructor sets the discrete points and the corresponding values.
   */
  Payoff(const double* pdS, const double* pdValues, size_t nNbS, 
         size_t nNbRegimes)
  {
    m_nNbS = nNbS;
    m_nNbRegimes = nNbRegimes;

    m_pdS = Array<double>(m_nNbS);
    m_pdValues = Array<double>(m_nNbS * m_nNbRegimes);

    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
      m_pdS[nIdxS] = pdS[nIdxS];

    for (size_t nIdx = 0; nIdx < m_nNbS * m_nNbRegimes; nIdx++)
      m_pdValues[nIdx] = pdValues[nIdx];
  }
 
  // Default dtor is ok

  /**
     Get the option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double *pdS, double *pdPrices, size_t nNbS) const
  {
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      numeric::QuadraticInterpolate(m_pdS.Get(), 
                                    m_pdValues.Get() + nIdxR * m_nNbS, 
                                    m_nNbS,
                                    pdS, 
                                    pdPrices + nIdxR * nNbS,
                                    nNbS);
    }
  }


protected:

  /// The number of points
  size_t m_nNbS;

  /// The number of regime 
  size_t m_nNbRegimes;

  /// The discrete points
  Array<double> m_pdS;

  /// The values at the discrete points
  Array<double> m_pdValues;
};


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PAYOFF_H_
