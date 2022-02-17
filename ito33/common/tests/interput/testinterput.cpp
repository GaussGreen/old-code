/**
  @file tests/interput/testIU.cpp

  @brief test functions of interputil

 */

#include <cmath>

#include "ito33/array.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"

#include "ito33/tests/testinterput.h"

using namespace ito33;
using namespace ito33::numeric;

// Global test function
double polynom(double x)
{
  return ((x - 1.) * x + 1.) * x - 1.;
}

double fstderpol(double x)
{
  return (3.*x - 2.)* x + 1.;
}

double scdderpol(double x)
{
  return 6.*x - 2.;
}

void InterpUtTest::LinearSearch()
{
  // Create a to be searched array
  size_t nNbS = 10;
  Array<double> pdS(nNbS);
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
    pdS[nIdx] = double(nIdx);

  // Create a new array that we want to search the intervals
  size_t nNbS2 = nNbS + 1;
  Array<double> pdS2(nNbS2);
  for (nIdx = 0; nIdx < nNbS2; nIdx++)
    pdS2[nIdx] = double(nIdx) - 0.1;

  Array<int> piIntervals(nNbS2);

  SearchIntervals(pdS.Get(), nNbS, pdS2.Get(), piIntervals.Get(), nNbS2);

  for (nIdx = 0; nIdx < nNbS2; nIdx++)
    CPPUNIT_ASSERT(piIntervals[nIdx] == int(nIdx) - 1);
}

void InterpUtTest::LinearInterp()
{
  const size_t
    nX = 2001;

  double 
    dX[nX],
    dY[nX],
    dRes[6],
    dXmin = -2.,
    dXmax = 3.,
    dStep,
    //-1.8937645
    dP[6] = {-2.001, -1.0173398, -0.4377286, 
             0.8937645, 1.0173398, 2.3572946},
    dYnew[6];

  ExtrapolationMode
    EMLeft = ExtrapolationMode_Linear,
    EMRight = ExtrapolationMode_Linear;

  dStep = fabs(dXmax - dXmin)/(nX - 1);

  for (int i=0; i <= nX-1; i++)
  {
    dX[i] = dXmin + dStep * i;
    dY[i] = polynom(dX[i]);
  }

  Interpolate(dX, dY, nX, dP, dYnew, 6, EMLeft, EMRight);
    
  for (int j=0; j< 6; j++)
  {
    dRes[j] = polynom(dP[j]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dYnew[j], dRes[j], 1e-4  );
  }
}

void InterpUtTest::QuadraInt()
{
  const size_t
    nX = 2001;

  double 
    dX[nX],
    dY[nX],
    dRes[6],
    dXmin = -2.,
    dXmax = 3.,
    dStep,
    //-1.8937645
    dP[6] = {-2.001, -1.0173398, -0.4377286, 
             0.8937645, 1.0173398, 2.3572946},
    dYnew[6];

  ExtrapolationMode
    EMLeft = ExtrapolationMode_Linear,
    EMRight = ExtrapolationMode_Linear;

  dStep = fabs(dXmax - dXmin)/(nX - 1);

  for (int i=0; i <= nX-1; i++)
  {
    dX[i] = dXmin + dStep * i;
    dY[i] = polynom(dX[i]);
  }

  QuadraticInterpolate(dX, dY, nX, dP, dYnew, 6, EMLeft, EMRight);
    
  for (int j=0; j< 6; j++)
  {
    dRes[j] = polynom(dP[j]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dYnew[j], dRes[j], 1e-4  );
  }
}


void InterpUtTest::Spline()
{
  const size_t
    nX = 21;

  double 
    dX[nX],
    dY[nX],
    dY2[nX],
    dRes[nX],
    dXmin = -2.,
    dXmax = 3.,
    dStep,
    dYP0,
    dYPN;

  dStep = fabs(dXmax - dXmin)/(nX - 1);

  for (int i=0; i<=nX-1; i++)
  {
    dX[i] = dXmin + dStep * i;
    dY[i] = polynom(dX[i]);
    dRes[i] = scdderpol(dX[i]);
  }

  dYP0 = fstderpol(dXmin);
  dYPN = fstderpol(dXmax);

  numeric::Spline(dX, dY, nX, dYP0, dYPN, dY2);
  
  // 
  for (int j=0; j<= nX-1; j++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dY2[j], dRes[j], 1e-5  );
}

void InterpUtTest::Splint()
{
  const size_t
    nX = 21;

  double 
    dX[nX],
    dY[nX],
    dY2[nX],
    dRes[6],
    dXmin = -2.,
    dXmax = 3.,
    dStep,
    dYP0,
    dYPN,
    dP[6] = {-2.351, -1.0173398, -0.4377286, 
             0.8937645, 1.0173398, 2.3572946};

  dStep = fabs(dXmax - dXmin)/(nX - 1);

  for (int i=0; i< nX; i++)
  {
    dX[i] = dXmin + dStep * i;
    dY[i] = polynom(dX[i]);
  }

  dYP0 = fstderpol(dXmin);
  dYPN = fstderpol(dXmax);

  numeric::Spline(dX, dY, nX, dYP0, dYPN, dY2);
    
  for (int j=0; j< 6; j++)
  {
    dRes[j] = polynom(dP[j]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(numeric::Splint(dX, dY, dY2, nX, dP[j]),
                                 dRes[j], 1e-5  );
  }
}


void InterpUtTest::LinearInterpolateUsingMatrix()
{
  const size_t nNbX = 2001;
  const size_t nNbNewX = 6;

  double 
    pdX[nNbX],
    pdY[nNbX],
    dRes,
    dXmin = -2.,
    dXmax = 3.,
    dStep,
    //-1.8937645
    pdNewX[nNbNewX] = {-2.001, -1.0173398, -0.4377286, 
                     0.8937645, 1.0173398, 2.3572946},
    pdnewY[nNbNewX];

  dStep = fabs(dXmax - dXmin) / (nNbX - 1);

  for (int i = 0; i <= nNbX - 1; i++)
  {
    pdX[i] = dXmin + dStep * i;
    pdY[i] = polynom(pdX[i]);
  }

  InterpolationMatrix matrix(pdNewX, nNbNewX, pdX, nNbX, 1);
    
  matrix.ProductMatrixVector(pdY, pdnewY);

  for (int j = 0; j < 6; j++)
  {
    dRes = polynom(pdNewX[j]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdnewY[j], dRes, 1e-4  );
  }
}

void InterpUtTest::QuadraticInterpolateUsingMatrix()
{
  const size_t nNbX = 200;
  double pdX[nNbX], pdY[nNbX]; 
  const size_t nNbNewX = nNbX + 1;
  double pdNewX[nNbNewX], pdNewY[nNbNewX], pdNewYI[nNbNewX];

  // Setup the given function and the points that we want to interpolate at
  for (size_t n = 0; n < nNbX; n++)
  {
    pdX[n] = double(n);
    pdY[n] = sin(pdX[n]);

    pdNewX[n] = pdX[n] - 0.004 * (n + 1);
  }
  pdNewX[nNbNewX - 1] = double(nNbX) + 0.01;

  // Test with ExtrapolationMode_Linear
  InterpolationMatrix 
    matrix(pdNewX, nNbNewX, pdX, nNbX, 1, true, 
           ExtrapolationMode_Linear, ExtrapolationMode_Linear,
           InterpolationMethod_Quadratic);

  matrix.ProductMatrixVector(pdY, pdNewY);

  QuadraticInterpolate(pdX, pdY, nNbX, pdNewX, pdNewYI, nNbNewX);

  for (size_t n = 0; n < nNbNewX; n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdNewYI[n], pdNewY[n], 1.e-12);
  
  // Test with ExtrapolationMode_Constant
  InterpolationMatrix 
    matrix1(pdNewX, nNbNewX, pdX, nNbX, 1, true, 
            ExtrapolationMode_Constant, ExtrapolationMode_Constant,
            InterpolationMethod_Quadratic);

  matrix1.ProductMatrixVector(pdY, pdNewY);

  QuadraticInterpolate(pdX, pdY, nNbX, pdNewX, pdNewYI, nNbNewX,
                       ExtrapolationMode_Constant, ExtrapolationMode_Constant);

  for (size_t n = 0; n < nNbNewX; n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdNewYI[n], pdNewY[n], 1.e-12);
}
