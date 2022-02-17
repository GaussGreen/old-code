/**************************************************************************
 * file      : common/src/numeric/interpolation.cpp
 * Purpose   : interpolation functions
 * Author    : ZHANG Yunzhi
 * Created   : 3/10/2001
 * RCS-ID    : $Id: interpolation.cpp,v 1.16 2005/06/14 16:47:23 wang Exp $
 * Copyright : (c) 2001 - 3005 Trilemma LLP all rights reserved
 **************************************************************************/

/**
   @todo The Spline function (that computes the second derivative
         of the function to be interpolated) returns 0 when 
         dY0 -> infinity with a third order polynom  
 */

#include "ito33/debug.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/morsematrix.h"

namespace ito33
{

namespace numeric
{


void SearchIntervals(const double* pdS, size_t nNbS,
                     const double* pdNewS, int* piIntervals, size_t nNbNewS)
{
  ASSERT_MSG( nNbS != 0,
              "The number of points on the mesh can't be null in "
              "SearchIntervals Function"
            );

  ASSERT_MSG( pdS != 0 && pdNewS != 0,
              "Null pointer in SearchIntervals Function"
            );

  int iIdxNewBegin = 0;
  while ( iIdxNewBegin < int(nNbNewS) && pdNewS[iIdxNewBegin] < pdS[0] )
    piIntervals[iIdxNewBegin++] = -1;

  int iIdxNewEnd = (int)nNbNewS - 1;
  while ( iIdxNewEnd >= 0 && pdNewS[iIdxNewEnd] >= pdS[nNbS - 1] )
    piIntervals[iIdxNewEnd--] = int(nNbS) - 1;

  int iIdx = 0;
  for (int iIdxNew = iIdxNewBegin; iIdxNew <= iIdxNewEnd; iIdxNew++)
  {
    while (pdS[iIdx] < pdNewS[iIdxNew])
      iIdx++;
    
    piIntervals[iIdxNew] = iIdx - 1;
  }
}

void Interpolate
    (const double *pdX, const double *pdY, size_t nNbX, 
     const double *pdXnew, double *pdYnew, size_t nNbNew,
     ExtrapolationMode emLeft,
     ExtrapolationMode emRight,
     InterpolationMethod interpMethod)
{
  ASSERT(    interpMethod == InterpolationMethod_Linear
          || interpMethod == InterpolationMethod_Quadratic );

  if ( interpMethod == InterpolationMethod_Linear )
    LinearInterpolate(pdX, pdY, nNbX, pdXnew, pdYnew, nNbNew, 
                      emLeft, emRight);
  else if ( interpMethod == InterpolationMethod_Quadratic )
    QuadraticInterpolate(pdX, pdY, nNbX, pdXnew, pdYnew, nNbNew, 
                         emLeft, emRight);
}

/*
  Caution should be taken on extrapolation. Either 0 or linear, it is usualy
  not very precise, especially for not so small step.
*/
void LinearInterpolate
     (const double *pdX, const double *pdY, size_t nNbX, 
      const double *pdXnew, double *pdYnew, size_t nNbNew,
      ExtrapolationMode emLeft, ExtrapolationMode emRight)
{
  size_t
    nX;

  ASSERT_MSG( nNbX != 0,
              "The number of points on the mesh can't be null in "
              "Interpolate Function"
            );

  ASSERT_MSG( pdX != 0 && pdY != 0,
              "Null pointer in Interpolate Function"
            );

  if (nNbX == 1)
  {
    for (nX = 0; nX < nNbNew; nX++)
      pdYnew[nX] = pdY[0];

    return;
  }

  double
    dTang;

  nX = 0;
  dTang = (pdY[1] - pdY[0]) / (pdX[1] - pdX[0]);
 
  size_t
    nNew = 0;

  // find XNew such that Xnew < X[0] and then extrapolate
  switch (emLeft)
  {
    case ExtrapolationMode_Constant :
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; nNew++)
        pdYnew[nNew] = pdY[0];
      break;

    case ExtrapolationMode_Linear : 
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; ++nNew)
        pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];
      break;

    case ExtrapolationMode_Zero :
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; nNew++)
        pdYnew[nNew] = 0.; 
  }

  if (nNbX > 2)
    for ( ;; )
    {
      // calculate pdYnew for all points between [x[ix], x[ix+1]]
      while (nNew < nNbNew && pdXnew[nNew] <= pdX[nX + 1])
      {
        pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];
        nNew++;
      }
    
      if(nNew == nNbNew)
        break;

      // here Xn[iN] > X[ix]
      nX++;
      // find nx such that x[ix] < Xn[iN] <= x[ix + 1]
      while(nX < nNbX - 2 && pdXnew[nNew] > pdX[nX + 1])
        nX++;
    
      dTang = (pdY[nX + 1] - pdY[nX]) / (pdX[nX + 1] - pdX[nX]);
    
      if (nX == nNbX - 2)
        break;
    }

  for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX + 1]; nNew++)
    pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];

  if (nNew == nNbNew)
    return;

  // nX should be nNbX - 2 now
  switch (emRight)
  {
    case ExtrapolationMode_Constant :
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = pdY[nX + 1]; 
      break;

    case ExtrapolationMode_Linear :
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];
      break;

    case ExtrapolationMode_Zero : 
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = 0.; 
  }
}

/*
  quadratically interpolates a function known on a set of values on a 
  new set of value
*/
void QuadraticInterpolate(const double *pdX, const double *pdY, size_t nNbX, 
                          const double *pdXnew, double *pdYnew, 
                          size_t nNbNew, ExtrapolationMode emLeft, 
                          ExtrapolationMode emRight)
{
  size_t
    nX;

  ASSERT_MSG( nNbX != 0,
              "The number of points on the mesh can't be null in "
              "QuadraticInterpolate function"
            );

  ASSERT_MSG( pdX != 0 && pdY != 0,
              "Null pointer in Quadratic Interpolate Function"
            );

  if (nNbX == 1)
  {
    for (nX = 0; nX < nNbNew; nX++)
      pdYnew[nX] = pdY[0];

    return;
  }

  double
    dTang;

  nX = 0;
   
  size_t
    nNew = 0;

  // find XNew such that Xnew < X[0] and then interpolate
  switch (emLeft)
  {
    case ExtrapolationMode_Constant :
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; nNew++)
        pdYnew[nNew] = pdY[0];
      break;

    case ExtrapolationMode_Linear : 
      dTang = (pdY[1] - pdY[0]) / (pdX[1] - pdX[0]);
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; nNew++)
        pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];
      break;

    case ExtrapolationMode_Zero :
      for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX]; nNew++)
        pdYnew[nNew] = 0.; 
  }

  if (nNbX > 2)
    for ( ;; )
    {
      // Look to the right for the extra point.  
      double x1 = pdX[nX];
      double x2 = pdX[nX+1];
      double x3 = pdX[nX+2];

      double y1 = pdY[nX];
      double y2 = pdY[nX+1];
      double y3 = pdY[nX+2];

      // calculate pdYnew for all points between [x[nX], x[nX+1]]
      while (nNew < nNbNew && pdXnew[nNew] <= pdX[nX + 1])
      {
        // Setup Lagrange interpolant
        double x = pdXnew[nNew];
        double l1 = (x-x2)*(x-x3) / ( (x1-x2)*(x1-x3) );
        double l2 = (x-x1)*(x-x3) / ( (x2-x1)*(x2-x3) );
        double l3 = (x-x1)*(x-x2) / ( (x3-x1)*(x3-x2) );

        pdYnew[nNew] = y1*l1 + y2*l2 + y3*l3;
        nNew++; 
      }
    
      if(nNew == nNbNew)
        break;

      // for the moment Xn[iN] > X[ix]
      nX++;
      // find nx such that x[ix] < Xn[iN] <= x[ix + 1]
      while(nX < nNbX - 2 && pdXnew[nNew] > pdX[nX + 1])
        nX++;
    
      if (nX == nNbX - 2)
        break;
    }

  // Look to the left for the extra point
  double x1 = pdX[nX-1];
  double x2 = pdX[nX];
  double x3 = pdX[nX+1];

  double y1 = pdY[nX-1];
  double y2 = pdY[nX];
  double y3 = pdY[nX+1];

  for (; nNew < nNbNew && pdXnew[nNew] <= pdX[nX + 1]; nNew++)
  {
    double x = pdXnew[nNew];
    double l1 = (x-x2)*(x-x3) / ( (x1-x2)*(x1-x3) );
    double l2 = (x-x1)*(x-x3) / ( (x2-x1)*(x2-x3) );
    double l3 = (x-x1)*(x-x2) / ( (x3-x1)*(x3-x2) );

    pdYnew[nNew] = y1*l1 + y2*l2 + y3*l3;
  }

  if (nNew == nNbNew)
    return;

  // nX should be nNbX - 2 now
  switch (emRight)
  {
    case ExtrapolationMode_Constant :
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = pdY[nX + 1]; 
      break;

    case ExtrapolationMode_Linear :
      dTang = (pdY[nX + 1] - pdY[nX]) / (pdX[nX + 1] - pdX[nX]);
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = dTang * (pdXnew[nNew] - pdX[nX]) + pdY[nX];
      break;

    case ExtrapolationMode_Zero : 
      for (; nNew < nNbNew; nNew++)
        pdYnew[nNew] = 0.; 
  }
}

/*
  This function computes the second derivative of the function given on 
  the set of points pdX[i] by pdY[i]. See Numerical Recipes in C p113-166 for
  details. The result is to be used in a third order interpolation. pdY2[] 
  verifies A pdY2[] = B where A is a tridiagonal matrix using values of
  pdX[] and B a vector using values of pdX[] and pdY[]. The system is solved
  thanks to the tri diagonal algorithm that can be found in Numerical 
  Recipes p51

  @todo when tested with a third order polynom, it returns 0 when
        dY0 -> infinity 
*/
void Spline(const double *pdX, const double *pdY, size_t nNbX, 
            double dYP0, double dYPNm1, 
            double *pdY2)
{
  size_t
    nX;
  
  double
    dP,
    dTmp,
    dSig,
    dUn;
    
  ASSERT_MSG( nNbX != 0,
              "The number of points on the mesh can't be null in "
              "function Spline"
            );

  ASSERT_MSG( pdX != 0 && pdY != 0,
              "Null pointer in function Spline"
            );

  double *pdU = new double [nNbX];

  for ( nX = 1; nX < nNbX; nX++ ) //B vector
  {
    pdU[nX] = (pdY[nX] - pdY[nX - 1]) / (pdX[nX] - pdX[nX - 1]);
  }

  if (dYP0 > 0.99e10 || dYP0 < -0.99e10) // f'(x0) -> +/- infinity
    pdY2[0] = pdU[0] = 0.0;               
  else // lower boundary condition have a specified first derivative
  {
    pdY2[0] = - 0.5;
    pdU[0] = (3.0 / (pdX[1] - pdX[0])) * ((pdY[1] - pdY[0]) 
            / (pdX[1] - pdX[0]) - dYP0);
  }

  //decomposition loop for tridiagonal algorithm
  for (nX = 1; nX < nNbX - 1; nX++) 
  {
    dTmp = (pdX[nX + 1] - pdX[nX - 1]);
    dSig = (pdX[nX] - pdX[nX - 1]) / dTmp;
    dP = dSig * pdY2[nX - 1] + 2.0;
    pdY2[nX] = (dSig - 1.0) / dP;
    pdU[nX] = pdU[nX+1] - pdU[nX];  
    pdU[nX] = (6.0 * pdU[nX] / dTmp - dSig * pdU[nX - 1]) / dP;
  }

  if (dYPNm1 > 0.99e10 || dYPNm1 < -0.99e10) // f'(X[n-1]) -> +/- infinity
    dUn = (6.0 / (pdX[nNbX - 1] - pdX[nNbX - 2]))
        * ( - (pdY[nNbX - 1] - pdY[nNbX - 2]) 
            / (pdX[nNbX - 1] - pdX[nNbX - 2]) );
  else           
  {
    dUn = (3.0 / (pdX[nNbX - 1] - pdX[nNbX - 2])) 
        * ( dYPNm1 - (pdY[nNbX - 1] - pdY[nNbX - 2])
                   / (pdX[nNbX - 1] - pdX[nNbX - 2]) );
  }

  pdY2[nNbX - 1] = (dUn - 0.5 * pdU[nNbX - 2]) / ( 1. + 0.5 * pdY2[nNbX - 2] );
  
  for(nX = nNbX - 2; nX > 0; nX--)
    pdY2[nX] = pdY2[nX] * pdY2[nX + 1] + pdU[nX];

  pdY2[0] = pdY2[0] * pdY2[1] + pdU[0];
  
  delete [] pdU;

}

/*
  This function computes the cubic spline interpolation on a point,
  given a set of points, of values of the interpolated function and of 
  its second derivative at these points, the latter having been calculated
  with spline function. See Numerical Recipes in C p 113-116.
*/
double Splint(const double *pdX, const double *pdY, const double *pdY2,
              size_t nNbX, 
              double dX)
{ 
  ASSERT_MSG
  (
     nNbX != 0,
     "The number of points on the mesh can't be null in function Splint"
  );

  ASSERT_MSG( pdX != 0 && pdY != 0 && pdY2 != 0,
              "Null pointer in function Splint" 
            ) ;

  int
    iIMin = 0,
    iIMax = (int)nNbX-1,
    iIMilieu;

  double
    dH,
    dB,
    dA;

  while(iIMax > iIMin)
  {
    iIMilieu = (iIMin + iIMax) >> 1;
    if( pdX[iIMilieu] < dX)
      iIMin = iIMilieu + 1;
    else
      iIMax = iIMilieu;
  }

  if(iIMin > 0)
    iIMax = iIMin--;
  else    // ajout 27/04/2000 (z)
    {
    iIMin = 0;
    iIMax = 1;
    }
    
  dH = pdX[iIMax] - pdX[iIMin];
  dA = (pdX[iIMax] - dX) / dH;
  dB = (dX - pdX[iIMin]) / dH;

  return   dA * pdY[iIMin] + dB * pdY[iIMax] 
         + (  (dA * dA * dA - dA) * pdY2[iIMin] 
            + (dB * dB * dB - dB) * pdY2[iIMax]  ) * dH * dH / 6.0;
}


/*
  computs the second derivatives in dX2 direction for a multi dimensional
  function. See Numerical Recipes in C, p123-128
*/
void Splie2(const double *pdX2, const double **ppdY, size_t nNbY, size_t nNbX,
            double **ppdY2)
{
  ASSERT_MSG( nNbX != 0 && nNbY !=0,
               "One of the array sizes is null in function Splie2"
            );

  ASSERT_MSG( pdX2 != 0 && ppdY != 0,
               "Null pointer in function Splie2"
            );

  for (size_t nY = 0; nY < nNbY; nY++)
    Spline(pdX2, ppdY[nY], nNbX, 1.0e30, 1.0e30, ppdY2[nY]);
}

/*
  computes the cubic spline interpolation for a multivariable function.
  See Numerical Recipes in C p 113-116.
*/
double Splin2(const double *pdX1, const double *pdX2, 
              const double **ppdY, const double **ppdY2, 
              size_t nNbY, size_t nNbX, 
              double dX1, double dX2)
{
  double
    dY;
    
  ASSERT_MSG( nNbX != 0 && nNbY !=0,
              "One of the array sizes is null in function Splin2"
            );

  ASSERT_MSG( pdX1 != 0 && pdX2 != 0 && ppdY != 0 && ppdY2 != 0,
              "Null pointer in function Splin2"
            );

  double 
    *pdYTmp = new double [nNbY],
    *pdYYTmp = new double [nNbY];

  for (size_t nY = 0; nY < nNbY; nY++)
    pdYYTmp[nY] = Splint(pdX2, ppdY[nY], ppdY2[nY], nNbX, dX2);

  Spline(pdX1, pdYYTmp, nNbY, 1.0e30, 1.0e30, pdYTmp);

  dY = Splint(pdX1, pdYYTmp, pdYTmp, nNbY, dX1);

  delete [] pdYTmp;
  delete [] pdYYTmp;

  return dY;
}


double LinearInterpolate2D(double dx,double dy,double dx1,double dx2,
                           double dx3,double dx4,double dy1,double dy2,
                           double df1,double df2,double df3,double df4)
{ 
  // horizontal interpolation
  double dWeight   = (dx - dx1) / (dx2 - dx1); 
  double dValAbove = (1.0 - dWeight) * df1 + dWeight * df2;
 
  dWeight          = (dx - dx3) / (dx4 - dx3);
  double dValBelow = (1.0 - dWeight) * df3 + dWeight * df4;

  // vertical interpolation between dValAbove and dValBelow 
  dWeight = (dy - dy1) / (dy2 - dy1);

  return  (1.0 - dWeight) * dValAbove + dWeight * dValBelow;

} // LinearInterpolate2D


} // namespace numeric

} // namespace ito33
