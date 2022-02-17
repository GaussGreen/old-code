/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/interpolationmatrix.h
// Purpose:     Matrix helps for interpolation
// Created:     2005/06/04
// RCS-ID:      $Id: interpolationmatrix.h,v 1.10 2005/09/20 12:25:36 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_INTERPOLATIONMATRIX_H_
#define _ITO33_NUMERIC_INTERPOLATIONMATRIX_H_

#include "ito33/vector.h"

#include "ito33/numeric/interpolation.h"

namespace ito33
{

namespace numeric
{


/// Struct as an element of the matrix
struct WeightedPoint
{
  /// Ctor by using the index and the corresponding weight
  WeightedPoint(size_t nIdx, double dWeight)
  {
    m_nIdx = nIdx;
    m_dWeight = dWeight;
  }

  /// The index of the point used to do the interpolation
  size_t m_nIdx;

  /// The weight of the point used to do the interpolation
  double m_dWeight;
};

/// A collection of weighted points as a row of the matrix
class WeightedPoints : public std::vector<WeightedPoint>
{
public:

  /**
     Does the interpolation using the given values at the interpolation
     points.
    
     @param pdValues The value of the points to be used for interpolation
   */
  double InterpolateWith(const double* pdValues) const;

  /**
     Reverses the interpolation, as transpose matrix vector product.

     @param dValue The value at the interpolated point
     @param pdValues The result of the transpose matrix vector product,
                     on the points used to do the interpolation
   */
  void TransposeInterpolate(double dValue, double* pdValues) const;
};


/// A class helps for interpolation. 
class InterpolationMatrix
{

public:

  /**
     Ctor takes the points and the space mesh.

     @param pdX The points at which we want to get the interpolated value
     @param nNbPoints The number of points at which we want to get the value
     @param pdS The points to be used for interpolation.
     @param nNbS The number of points to be used for interpolation.
     @param nNbSubSystem The number of subsystem. We are supposing that the 
                         sub systems have the same mesh.
     @param bIsBackard true if the interpolation will be from each sub system
                       to itself, false if values interpolated at all sub 
                       systems will be summed up.
     @param interpMethod The interpolation method(linear, quadratic etc)
   */
  InterpolationMatrix
  (const double* pdX, size_t nNbPoints, 
   const double* pdS, size_t nNbS,
   size_t nNbSubSystem,
   bool bIsBackard = true,
   ExtrapolationMode emLeft = ExtrapolationMode_Linear, 
   ExtrapolationMode emRight = ExtrapolationMode_Linear,
   InterpolationMethod interpMethod = InterpolationMethod_Linear);

  /**
     Multiplies the matrix by a scalar.

     @param dScalar The scalar that the matrix will be multiplied by
   */
  void MultiplyBy(double dScalar);

  /**
     Add interpolation weight to matrix to particular entry

     @param dWeight weight to add 
     @param nIdxP1 row number
     @param nIdxP2 column number
   */
  void Add(size_t nIdxP1, size_t nIdxP2, double dWeight);

  /**
     Product matrix vector, actually is just the interpolation.

     @param pdX The values at the points to be used for interpolation
     @param pdF The interpolated values
     @param bAddon true if the interpolated values are already initialized,
                   false otherwise
   */
  void 
  ProductMatrixVector
  (const double* pdX, double* pdF, bool bAddon = false) const;

  /**
     Product transpose matrix vector.

     @param pdX The interpolated value
     @param pdF The contribution of the interpolated values on the points
                used for interpolation.
     @param true if the values are already initialized, false otherwise
   */
  void
  ProductTransposeMatrixVector
  (const double* pdX, double* pdF, bool bAddon = false) const;


private:

  /// The number of points at which we want to get the values
  size_t m_nNbPoints;

  /// The number of points the interpolation will be based on
  size_t m_nNbS;

  /// The number of sub systems. 
  size_t m_nNbSubSystem;

  /// if the interpolation will be from each sub system to itself
  bool m_bIsBackward;

  /// The storage of the matrix
  std::vector<WeightedPoints> m_points;

}; // class InterpolationMatrix


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_INTERPOLATIONMATRIX_H_
