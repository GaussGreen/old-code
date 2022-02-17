/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/tridiagonalmatrix.h
// Purpose:     tridiagonal matrix sclass
// Author:      ICARE
// Created:     2004/1/28
// RCS-ID:      $Id: tridiagonalmatrix.h,v 1.13 2005/06/01 12:52:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/tridiagonalmatrix.h
    @brief Tridiagonal matrix class
 */

#ifndef _ITO33_NUMERIC_TRIDIAGONALMATRIX_H_
#define _ITO33_NUMERIC_TRIDIAGONALMATRIX_H_

#include "ito33/common.h"
#include "ito33/vector.h"
#include "ito33/debug.h"


namespace ito33
{

namespace numeric
{


/// A simple tridiagonal matrix
class TridiagonalMatrix
{

public:

  /**
     ctor, allocate the memory for the arrays

     @param nNb the estimated initial dimension of the tridiagonal matrix

     REQUIRE: nNb can't be 0
   */
  TridiagonalMatrix(size_t nNb = 10) :
      m_nNb(nNb),
      m_nMaxSize(nNb),
      m_pdA(nNb),
      m_pdB(nNb),
      m_pdC(nNb)
  {
    ASSERT_MSG(nNb != 0, "The matrix dimension can't be 0.");
  }
 
  // Default dtor is ok.

  /// Get the dimension of the matrix
  size_t Dimension() const { return m_nNb; }

  /// Get a pointer to the lower diagonal part
  double* GetA() { return &m_pdA[0]; }

  /// Get a const pointer to the lower diagonal part
  const double* GetA() const { return &m_pdA[0]; }

  /// Get a pointer to the diagonal part
  double* GetB() { return &m_pdB[0]; }

  /// Get a const pointer to the diagonal part
  const double* GetB() const { return &m_pdB[0]; }

  /// Get a pointer to the upper diagonal part
  double* GetC() { return &m_pdC[0]; }

  /// Get a const pointer to the upper diagonal part
  const double* GetC() const { return &m_pdC[0]; }

  /** 
     Do the matrix-vector product.

     @param pdX the vector to be multiplied
     @param pdF the product of the matrix-vector product
     @param bAddon Add the product to pdF or not
   */
  void 
  ProductMatrixVector
  (const double *pdX, double *pdF, bool bAddon = false) const;

  /** 
     Do the transpose matrix vector product.

     @param pdX the vector to be multiplied
     @param pdF the product of the transpose matrix-vector product
     @param bAddon Add the product to pdF or not
   */
  void 
  ProductTransposeMatrixVector
  (const double* pdX, double* pdF, bool bAddon = false) const;

  /**
     Define the dimension of the matrix.

     @param nDimension new dimension of the matrix

     REQUIRE: nDimension can't be 0
   */
  void SetDimension(size_t nDimension)
  {
    ASSERT_MSG(nDimension != 0, "The matrix dimension can't be 0.");

    m_nNb = nDimension;

    // resize the arrays
    if (nDimension > m_nMaxSize)
    {
      m_nMaxSize = nDimension;
      m_pdA.resize(m_nMaxSize);
      m_pdB.resize(m_nMaxSize);
      m_pdC.resize(m_nMaxSize);
    }
  }

  TridiagonalMatrix* GetTranspose() const;
  

private:

  /// The dimension of the matrix
  size_t m_nNb;
  
  /// The maximum size of the matrix
  size_t m_nMaxSize;

  /// The lower diagonal part 
  std::vector<double> m_pdA;

  /// The diagonal
  std::vector<double> m_pdB;
  
  /// The upper diagonal part
  std::vector<double> m_pdC;

}; // class TridiagonalMatrix


} // namespace numeric

} // namespace ito33

#endif // _ITO33_NUMERIC_TRIDIAGONALMATRIX_H_
