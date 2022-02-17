/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/piecewiseconstantfunctor.h
// Purpose:     Class for piecewise constant function
// Author:      Wang
// Created:     2004/06/03
// RCS-ID:      $Id: piecewiseconstantfunctor.h,v 1.5 2005/04/08 13:30:11 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/piecewiseconstantfunctor.h
    @brief Class for piecewise constant function.

    The name PiecewiseConstantFunctor is used instead of 
    PiecewiseconstantFunction to indicate that it's not a function in c++,
    but a class represents a piecewise constant function.
*/

#ifndef _ITO33_NUMERIC_PIECEWISECONSTANTFUNCTOR_H_
#define _ITO33_NUMERIC_PIECEWISECONSTANTFUNCTOR_H_

#include "ito33/vector.h"

namespace ito33
{

namespace numeric
{


/// Class for piecewise constant function
class PiecewiseConstantFunctor
{

public:

  /**
     Ctor takes two pointers to construct the step function

     Note that the element pdX[nNbX - 1] doesn't make much sense with the 
     actual convention.
     
     @param pdX the abscissas
     @param pdY the values of the step function
     @param nNbX the number of intervals
   */
  PiecewiseConstantFunctor(const double *pdX, const double *pdY, size_t nNbX);

  /**
     @internal
     @brief Returns the value of the function at a given x coordinate

     Please note that the function will search manually the interval, so won't 
     be efficient if number of intervals is big.

     @param dX the point at which the function is evalued
     @return the value of the step function at dX

     @noexport
   */
  double operator()(double dX) const;

  /**
     Gets the X coodinates

     @return a reference to the internal coordinate vectors
   */
  const std::vector<double>& GetX() const { return m_pdX; }

  /**
     Gets the Y coodinates

     @return a reference to the function values vector 
   */
  const std::vector<double>& GetY() const { return m_pdY; }


protected:

  std::vector<double> m_pdX;

  std::vector<double> m_pdY;

}; // class StepFunctor


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_PIECEWISECONSTANTFUNCTOR_H_
