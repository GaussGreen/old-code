///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/normaldistribution.h
// Purpose:          useful functions related to normal distributions
// Author:           laurence
// Created:          18/09/2003
// RCS-ID:           $Id: normaldistribution.h,v 1.5 2006/03/23 09:27:14 yann Exp $
// Copyright:        (c) 2003 Trilemma LLP all rights reserved
///////////////////////////////////////////////////////////////////////////////

/**

  @file ito33/numeric/normaldistribution.h

  @brief Defines useful functions related to normal distributions. See 
  "handbook of mathematical functions" by Milton Abramowitz and 
  Irene A. Stegun, p931-932

 */

#ifndef _ITO33_NUMERIC_NORMALDISTRIBUTION_H_
#define _ITO33_NUMERIC_NORMALDISTRIBUTION_H_

namespace ito33
{

namespace numeric
{

/// constant
const double INVSQRT2PI = 0.3989422804014326779;

/**
   Normal law density probability
   @brief calculate \f$ e^{-x^2/2}/ \sqrt{2 \pi} \f$
  
 */
double Normal(double dX);

/**
   Cumulated normal distribution

   @brief Calculate the cumulated normal 
   \f$ \int_{-\infty}^x e^{-t^2/2}/\sqrt{2 \pi}\f$
   using the formulae found in "handbook of mathematical functions" by Milton
   Abramowitz and Irene A. Stegun, p931-932, according to expressions 26.2.2
   and 26.2.17

 */
double CumulatedNormal(double dX);

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_NORMALDISTRIBUTION_H_
