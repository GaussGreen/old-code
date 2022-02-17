/////////////////////////////////////////////////////////////////////////////
// Name:        common/tests/densematrix/testdensematrix.cpp
// Purpose:     unit tests for dense matrix
// Created:     2006/07/06
// RCS-ID:      $Id: testdensematrix.cpp,v 1.1 2006/07/06 16:02:26 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/sharedptr.h"
#include "ito33/cppunit.h"

#include "ito33/numeric/densematrix.h"
#include "ito33/numeric/computeexpmatrix.h"

using namespace ito33;
using namespace ito33::numeric;

class DenseMatrixTest : public CppUnit::TestCase
{

public:
 
  DenseMatrixTest() { }
  
private:

  CPPUNIT_TEST_SUITE( DenseMatrixTest );
    CPPUNIT_TEST( ComputeExpScalar );
    CPPUNIT_TEST( ComputeExpIdentity );
  CPPUNIT_TEST_SUITE_END();

  void ComputeExpScalar();
  void ComputeExpIdentity();

  NO_COPY_CLASS(DenseMatrixTest);
};

CPPUNIT_TEST_SUITE_REGISTRATION(DenseMatrixTest);

void DenseMatrixTest::ComputeExpScalar()
{
  double d = 0.8;
  DenseMatrix matrix(1);
  matrix[0][0] = d;
  DenseMatrix expMatrix(1);

  ComputeExpMatrix( 1, matrix.Get(), expMatrix.Get() );

  CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(d), expMatrix[0][0], 1.e-5);
}

void DenseMatrixTest::ComputeExpIdentity()
{
  size_t nSize = 3;
  DenseMatrix matrix(nSize);
  for (size_t nI = 0; nI < nSize; nI++)
    for (size_t nJ = 0; nJ < nSize; nJ++)
      matrix[nI][nJ] = 0;

  for (size_t nI = 0; nI < nSize; nI++)
    matrix[nI][nI] = 1.;

  DenseMatrix expMatrix(nSize);
  ComputeExpMatrix( nSize, matrix.Get(), expMatrix.Get() );
  for (size_t nI = 0; nI < nSize; nI++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL
    (exp(matrix[nI][nI]), expMatrix[nI][nI], 1.e-5);
}
