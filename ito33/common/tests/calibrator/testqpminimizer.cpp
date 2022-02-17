/////////////////////////////////////////////////////////////////////////////
// Name:        common/tests/calibrator/testqpminimizer.cpp
// Purpose:     unit tests for wrapper of NAG QP minimizer
// Created:     2005/06/21
// RCS-ID:      $Id: testqpminimizer.cpp,v 1.1 2005/06/23 08:29:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/numeric/qp_nag.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testqpminimizer.h"

#include <cmath>

using namespace std;
using namespace ito33;

void QPMinimizerTest::Identity()
{
  const size_t nDim = 2;
  double elements[nDim * nDim];
  Array<double *> matrix(nDim);
  double pdC[nDim];
  double pdX[nDim];

  for (size_t n = 0; n < nDim; n++)
    matrix[n] = elements + n * nDim;

  for (size_t n1 = 0; n1 < nDim; n1++)
    for (size_t n2 = 0; n2 < nDim; n2++)
      if ( n1 == n2 )
        matrix[n1][n2] = 1.;
      else
        matrix[n1][n2] = 0.;

  for (size_t n = 0; n < nDim; n++)
    pdC[n] = 0.5 * ( n + 1);

  for (size_t n = 0; n < nDim; n++)
    pdX[n] = 0.;

  numeric::QPMinimizerNAG minimizer(nDim);

  minimizer(matrix.Get(), pdC, pdX);

  for (size_t n = 0; n < nDim; n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( - pdC[n], pdX[n], 1.e-6);
}

void QPMinimizerTest::Appended()
{
  const size_t nDim = 2;
  double elements[nDim * (nDim + 1)];
  Array<double *> matrix(nDim + 1);
  double pdX[nDim];

  for (size_t n = 0; n < nDim + 1; n++)
    matrix[n] = elements + n * nDim;

  matrix[0][0] = 2.;
  matrix[0][1] = 1.;
  matrix[1][0] = 1.;
  matrix[1][1] = 3.;

  matrix[2][0] = 4.;
  matrix[2][1] = 7.;

  for (size_t n = 0; n < nDim; n++)
    pdX[n] = 0.;

  numeric::QPMinimizerNAG minimizer(nDim);

  minimizer(matrix.Get(), matrix[nDim], pdX);

  for (size_t n = 0; n < nDim; n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( - (double(n) + 1), pdX[n], 1.e-6);
}