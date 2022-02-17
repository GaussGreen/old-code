/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/payoffdiscrete.h
// Purpose:     Payoff class defined by values at discrete points
// Author:      WANG Xuewen
// Created:     2003/09/25
// RCS-ID:      $Id: payoffdiscrete.h,v 1.6 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/payoffdiscrete.h
    @brief Payoff class defined by values at discrete points
 */

#ifndef _ITO33_FINANCE_PAYOFFDISCRETE_H_
#define _ITO33_FINANCE_PAYOFFDISCRETE_H_

#include "ito33/array.h"
#include "ito33/finance/payoff.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"


namespace ito33 
{

namespace finance 
{

/**
  Payoff defined by values at discrete points.

  Interpolation is used to return values at arbitrary points (that may or may 
  not correspond to the discrete points defining the payoff).
*/
class PayoffDiscrete : public Payoff
{
public:

  /**
    constructor, set the discrete points and the corresponding values
  */
  PayoffDiscrete(const double* pdPoints, const double* pdValues, 
                 size_t nNbPoints)
  {
    m_nNbPoints = nNbPoints;
    m_pdPoints = Array<double>(m_nNbPoints);
    m_pdValues = Array<double>(m_nNbPoints);

    for (size_t nIdx = 0; nIdx < m_nNbPoints; nIdx++)
    {
      m_pdPoints[nIdx] = pdPoints[nIdx];
      m_pdValues[nIdx] = pdValues[nIdx];
    }
  }
 
  /// dummy virtual destructor
  virtual ~PayoffDiscrete() { }

  /**
     Get the option payoff values at an array of spots

     @param pdS the spots
     @param pdPrices the payoff values at the spots
     @param nNbS the number of spots
   */
  void Get(const double *pdS, double *pdPrices, size_t nNbS) const
  {
    numeric::QuadraticInterpolate(m_pdPoints.Get(), 
                                  m_pdValues.Get(), 
                                  m_nNbPoints,
                                  pdS, pdPrices, nNbS,
                                  numeric::ExtrapolationMode_Linear,
                                  numeric::ExtrapolationMode_Linear);
  }

protected:

  /// The number of points
  size_t m_nNbPoints;

  /// The discrete points
  Array<double> m_pdPoints;

  /// The values at the discrete points
  Array<double> m_pdValues;
  

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PAYOFFDISCRETE_H_
