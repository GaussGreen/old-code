/**************************************************************************
 * file      : ito33/numeric/interpolation.h
 * Purpose   : interpolation function header
 * Author    : ZHANG Yunzhi
 * Created   : 3/10/2001
 * RCS-ID    : $Id: interpolation.h,v 1.21 2005/12/02 21:08:52 yann Exp $
 * Copyright : (c) 2003 - 2005 Trilemma LLP all rights reserved
 **************************************************************************/

#ifndef _ITO33_NUMERIC_INTERPOLATION_H_
#define _ITO33_NUMERIC_INTERPOLATION_H_

#include "ito33/debug.h"

#include "ito33/numeric/extrapolationmode.h"

namespace ito33
{

namespace numeric
{


/**
    Helper for doing linear interpolation: given (x1, y1) and (x2, y2) it gives
    the value of y for any x.

    This class is a functor, i.e. it redefines operator() to do the
    interpolation -- this makes using it similar to calling an interpolation
    function even though this is in fact a method call.

    Typical usage is:
    @code
      // construct the interpolator outside the loop
      numeric::LinearInterpolator valueAt(dateStart, valueStart,
                                          dateEnd, valueEnd);
      while ( ... )
      {
        // use it here
        double value = valueAt(date);
      }
    @endcode
 */
class LinearInterpolator
{
public:
  /**
      Constructor takes the points between which we interpolate.

      Degenerate case when x1 == x2 is supported but in this case y1 should be
      equal to y2 as well, otherwise the data is inconsistent. Also,
      interpolation can only be done for the same point x == x1 == x2 then.
      This is still convenient as it allows to avoid checks for special cases
      in the user code.

      @param x1 abscissa of the first point
      @param y1 ordinate of the first point
      @param x2 abscissa of the second point
      @param y2 ordinate of the second point
   */
  LinearInterpolator(double x1, double y1, double x2, double y2)
  {
    m_x1 = x1;
    m_y1 = y1;

    if ( !(m_interpolate = x1 != x2) )
    {
      ASSERT_MSG( y1 == y2, "can't interpolate between (x, y1) and (x, y2)" );
      m_dx =
      m_dy = 0.;
    }
    else // real interpolation
    {
      m_dx = x2 - x1;
      m_dy = y2 - y1;
    }
  }

  /**
      Interpolation function.

      @param x point for which to do the interpolation.
      @return the interpolated value
   */
  double operator()(double x) const
  {
    if ( !m_interpolate )
    {
      ASSERT_MSG( x == m_x1, "can't extrapolate from a single point" );

      return m_y1;
    }

    return (m_dy*(x - m_x1))/m_dx + m_y1;
  }

private:
  double m_x1,
         m_dx,
         m_y1,
         m_dy;

  bool m_interpolate;
};


/// Enum for different interpolation methods
enum InterpolationMethod
{
  /// Linear interpolation
  InterpolationMethod_Linear,

  /// Quadratic interpolation
  InterpolationMethod_Quadratic,

  InterpolationMethod_Max
};

/** 
   Search the intervals on a given mesh of a set of new, ordered points.

   Note that if a point is equal exactly to one of the point of the mesh,
   the interval of this point can be the left or the right, so there is 
   no default choice(it shoudn't matter).

   @param pdS The given mesh where the points will be searched on
   @param nNbS The size of the given mesh
   @param pdNewS The points that need to find the intervals on the given mesh
   @param piIntervals Output, the intervals of the points on the given mesh
   @param nNbNewS The number of points to be searched
 */
void SearchIntervals(const double* pdS, size_t nNbS,
                     const double* pdNewS, int* piIntervals, size_t nNbNewS);

/**
   Interpolation function uses an additional parameter to switch between
   linear and quadratic method.

   @param pdX the mesh on which will be done the interpolation
   @param pdY the known values at pdX  
   @param nNbX the number of the mesh
   @param pdXNew the points that we need to get the interpolation values
   @param pdYNew the interpolated values
   @param nNbNew the number of the points to be interpolated
   @param emLeft the extrapolation type on the left
   @param emRight the extrapolation type on the right
   @param interpMethod The interpolation method(linear or quadratic)
 */
void Interpolate
    (const double *pdX, const double *pdY, size_t nNbX, 
     const double *pdXnew, double *pdYnew, size_t nNbNew,
     ExtrapolationMode emLeft = ExtrapolationMode_Linear,
     ExtrapolationMode emRight = ExtrapolationMode_Linear,
     InterpolationMethod interpMethod = InterpolationMethod_Linear);

/**
   Linear interpolation function. 

   @param pdX the mesh on which will be done the interpolation
   @param pdY the known values at pdX  
   @param nNbX the number of the mesh
   @param pdXNew the points that we need to get the interpolation values
   @param pdYNew the interpolated values
   @param nNbNew the number of the points to be interpolated
   @param emLeft the extrapolation type on the left
   @param emRight the extrapolation type on the right
 */
void LinearInterpolate
    (const double *pdX, const double *pdY, size_t nNbX, 
     const double *pdXnew, double *pdYnew, size_t nNbNew,
     ExtrapolationMode emLeft = ExtrapolationMode_Linear,
     ExtrapolationMode emRight = ExtrapolationMode_Linear);

/**
   Quadratically interpolates a function known on a set of values on a 
   new set of value.

   @param pdX array of the set of initial values
   @param pdY the known values at pdX 
   @param nNbX the number of the mesh
   @param pdXNew the points that we need to get the interpolation values
   @param pdYNew the interpolated values
   @param nNbNew the number of the points to be interpolated
   @param emLeft the extrapolation type on the left
   @param emRight the extrapolation type on the right
 */
void QuadraticInterpolate
    (const double *pdX, const double *pdY, size_t nNbX, 
     const double *pdXnew, double *pdYnew, size_t nNbNew,
     ExtrapolationMode emLeft = ExtrapolationMode_Linear, 
     ExtrapolationMode emRight = ExtrapolationMode_Linear);

/**
   Computes the second derivative of the function given on 
   the set of points ppdX[i] by ppdY[i]. 
  
   See Numerical Recipes in C p113-166 for details. 
   
   The result is to be used in a third order interpolation. ppdY2[] 
   verifies A ppdY2[] = B where A is a tridiagonal matrix using values of
   pdX[] and B a vector using values of pdX[] and pdY[]. The system is solved
   thanks to the tri diagonal algorithm that can be found in Numerical 
   Recipes p51
  
   @param pdX contains the set of points on which the function is known  
   @param pdY contains the vectors of values of the function at points pdX[] 
   @param nN is the size of vectors pointed by pdX and pdY
   @param dYP0 is the first derivative of the function at point ppdX[0]
   @param dYPNm1 is the first derivative of the function at point ppdX[N-1]
   @param *pdY2 is the final vector containing second derivative

   @todo when tested with a third order polynom, it returns 0 when
         dY0 -> infinity  
 */
void Spline(const double *pdX, const double *pdY, size_t nNbX, 
            double dYP0, double dYPNm1, 
            double *pdY2);

/**
   Computes the cubic spline interpolation on a point,
   given a set of points, of values of the interpolated function and of 
   its second derivative at these points, the latter having been calculated
   with spline function. 
   
   See Numerical Recipes in C p 113-116.
    
   @param pdX is the set of points on which the function is known
   @param pdY is the vectors of values of the function at points pdX[]
   @param dY2 is the vector containing second derivative
   @param nNbX is the number of given points
   @param dX is the point where the function given on *pdX by *pdY is
             to be calculated  
 */
double Splint(const double *pdX, const double *pdY, const double *pdY2, 
              size_t nNbX, 
              double dX);

/**
   Computs the second derivatives in dX2 direction dor a multi dimensional
   function. 
   
   See Numerical Recipes in C, p123-128

   @param pdX2 given points in dX2 direction
   @param ppdY values of a function of iM variables given on points *pdX2 of the dX2
               direction
   @param nNbY number of variables on which depends dYa
   @param nNbX number of points in dX2 direction
   @param ppdY2 arrays containing second derivatives of dYa in dX2 direction 
 */
void Splie2(const double *pdX2, const double **ppdY, 
            size_t nNbY, size_t nNbX,
            double **ppdY2);


/**
   Computes the cubic spline interpolation for a multivariable function.
  
   See Numerical Recipes in C p 113-116.

   @param pdX1 set of points in first direction 
   @param pdX2 set of points in second direction
   @param ppdY set of values of the function known at all points
               (pdX1[i], pdX2[j])
   @param ppdY2 set of values for the second derivative of the function in the
                dX2 direction, previously computed with Splie2 function
   @param nNbY number of variables on which depends the function (here 2)
   @param nNbX number of values in the dX1 direction
   @param dX1 the value in the first direction on which the interpolation
              has to be made (point of the interpolation is (dX1, dX2))
   @param dX2 the value in the second direction on which the interpolation has 
              to be made
 */
double Splin2(const double *pdX1, const double *pdX2, 
              const double **ppdY, const double **ppdY2, 
              size_t nNbY, size_t nNbX,
              double dX1, double dX2);


/**
   Computes a linear 2 D interpolation.

      
           | f1 (x1,y1)     f2,(x2,y1) |
        ---+---------------------------+-- 
           |         f (x,y)           |                            
        ---+---------------------------+--     the x values can all be different
           | f3 (x3,y2)     f4 (x4,y2) | 


   @param dx1,dx2,dx3,dx4 x-coordinates of the four corners of the box containing
                          the point (x,y)
          

   @param dy1,dy2 y coordinates

   @param x,y coordinates of the point to interpolate

   @param df_i, i=1,..4 values at each corners
   @return interpolated value f at (x,y)
 */
double LinearInterpolate2D(double dx,double dy,double dx1,double dx2,
                           double dx3,double dx4,double dy1,double dy2,
                           double df1,double df2,double df3,double df4);


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_INTERPOLATION_H_
