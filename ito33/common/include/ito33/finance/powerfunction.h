/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ops/powerfunction.h
// Purpose:     PowerFunction class
// Created:     2006/07/05
// RCS-ID:      $Id: powerfunction.h,v 1.3 2006/08/16 13:52:34 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ops/powerfunction.h
    @brief Declaration of PowerFunction class
 */

#ifndef _ITO33_FINANCE_POWERFUNCTION_H_
#define _ITO33_FINANCE_POWERFUNCTION_H_

#include "ito33/common.h"

namespace ito33
{

namespace finance
{

/**
    A power function. A power function is defined by the following
    formula
      \f[ power (S) = \alpha * (S_0 / S)^\beta \f]
 */
class PowerFunction
{
public:
  /**
      Creates PowerFunction object.

      @param alpha value of \f$\alpha\f$ variable in the power function.
      @param beta value of \f$\beta\f$ variable in the power function.
      @param S0 value of S0 variable in the power function.
   */
  PowerFunction(double alpha, double beta, double S0)
    : m_dAlpha(alpha),
      m_dBeta(beta),
      m_dS0(S0)
  {
  }

  /**
      Gets the value of Alpha variable of the power function.

      @return the value of \f$\alpha\f$ variable of the power function.
   */
  double GetAlpha() const
  {
    return m_dAlpha;
  }

  /**
      Gets the value of Beta variable of the power function.

      @return the value of \f$\beta\f$ variable of the power function.
   */
  double GetBeta() const
  {
    return m_dBeta;
  }

  /**
      Gets the value of S0 variable of the power function.

      @return the value of S0 variable of the power function.
   */
  double GetS0() const
  {
    return m_dS0;
  }

private:
  double m_dAlpha;

  double m_dBeta;

  double m_dS0;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_OPS_HAZARDRATEPOWER_H_
