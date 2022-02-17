/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/deltagamma.h
// Purpose:     some functions to compute the delta and gamma of a curve
// Author:      ICARE 
// Created:     2003/28/11
// RCS-ID:      $Id: deltagamma.h,v 1.8 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/deltagamma.h

    Functions to help the calculation of the delta and gamma of a curve.
*/

#ifndef _ITO33_NUMERIC_DELTAGAMMA_H_
#define _ITO33_NUMERIC_DELTAGAMMA_H_

#include <cstddef> 

namespace ito33
{

namespace numeric
{

/**
   Compute the delta according using spline.

   @param pdS the grid points
   @param pdValues the function values at the grid points
   @param nNbS, number of points   
   @param pdDeltas the array to be filled
 */
void ComputeDelta(const double *pdS, const double *pdValues, size_t nNbS, 
                  double *pdDeltas);

/**
   Compute the delta according to the finite difference scheme used during
   discretization.

   @param pdS the grid points
   @param piFDconst the finite difference scheme
   @param pdValues the function values at the grid points
   @param nNbS, number of points   
   @param pdDeltas the array to be filled 
 */
void ComputeDelta(const double *pdS, 
                  const int *piFDconst, 
                  const double *pdValues, 
                  size_t nNbS, 
                  double *pdDeltas);
/**
   Compute the gamma according using spline.

   @param pdS the grid points
   @param pdValues the function values at the grid points
   @param nNbS, number of points   
   @param pdGammas the array to be filled 
 */
void ComputeGamma(const double *pdS, const double *pdValues, size_t nNbS, 
                  double *pdGammas);

/**
   Compute the gamma according to the finite difference scheme used during
   discretization.

   @param pdS the grid points
   @param pdValues the function values at the grid points
   @param nNbS, number of points   
   @param pdGammas the array to be filled by the computed Gamma.
 */
void ComputeGammaFD(const double *pdS, const double *pdValues, size_t nNbS, 
                  double *pdGammas);
/**
   Compute the second derivative according using a central difference scheme, 
   ignoring the boundary.

   @param pdX the grid points
   @param pdPrices the function values at the grid points
   @param nNbS, number of points   
   @param pdValues the array to be filled by the computed 1st.
 */
void Compute2nd(const double* pdX, const double* pdPrices, size_t nNbS, 
                double* pdValues);

/**
   Compute the first derivative according to a central difference scheme.
   Assuming equally spaced grid ignoring the boundary.

   @param pdX the grid points
   @param pdPrices the function values at the grid points
   @param nNbS, number of points   
   @param pdValues the array to be filled by the computed 1st.
 */
void Compute1st(const double* pdX, const double* pdPrices, size_t nNbS, 
                double* pdValues);


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_DELTAGAMMA_H_
