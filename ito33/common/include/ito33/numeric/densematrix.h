/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/densematrix.h
// Purpose:     Functions for dealing with dense matrices
// Created:     2005/10/18
// RCS-ID:      $Id: densematrix.h,v 1.4 2006/07/06 16:01:31 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/densematrix.h
    @brief Functions for dealing with dense matrices
 */

#ifndef _ITO33_NUMERIC_DENSEMATRIX_H_
#define _ITO33_NUMERIC_DENSEMATRIX_H_

#include "ito33/array.h"

namespace ito33
{

namespace numeric
{

/// Class manages memory for a dense matrix represented as a double pointer
class DenseMatrix
{
public:

  /**
      Default ctor intializes an empty matrix, Init() must be called
      before being used.
   */
  DenseMatrix() { }

  /**
      Ctor takes the two dimensions and allocates the memory.
      
      @param nDim1 The number of line
      @param nDim2 The number of column
   */
  DenseMatrix(size_t nDim1, size_t nDim2)
  {
    Init(nDim1, nDim2);
  }

  /**
      Ctor takes the dimension for a square matrix and allocates the memory. 
      
      @param nDim The number of line
   */
  DenseMatrix(size_t nDim)
  {
    Init(nDim, nDim);
  }

  /**
      Allocates the memory for the matrix.

      @param nDim1 The number of line
      @param nDim2 The number of column
   */
  void Init(size_t nDim1, size_t nDim2);
  
  /**
      Returns the ith line.

      @param i The line number
   */
  double* operator[](size_t i) { return m_ppd[i]; }
  
  /// Returns as a double pointer, const version.
  const double* const * Get() const { return m_ppd.Get(); }

  /// Returns as a double pointer, non const version.
  double** Get() { return m_ppd.Get(); }

private:

  /// The matrix, can be used as a double pointer
  Array<double *> m_ppd;

  /// Storage
  Array<double> m_pd;
};
 
/**
    Euclidean norm of a square matrix.

    @param nSize The size of the matrix
    @param ppdA The matrix
    @return The Euclidean norm of the matrix
 */
double EuclideanNorm(size_t nSize, double const* const* ppdA);

/**
    LU factorization of dense matrix.

    The factors are stored in the original matrix.

    @param ppdA  The matrix to factorize. Stores LU on output.
    @param nSize The size of the matrix
    @param nIndices Stores the permutations
 */
void LU_Factor(double **ppdA, size_t nSize, size_t *nIndices);

/**
    Solves system A x = B using previously computed LU factorization.

    This function should be called LU_Factor.

    @param ppdA  The LU factorization of original matrix A
    @param nSize The size of the matrix
    @param nIndices The permutations done during the factorization
    @param pdB The right hand side
    @param pdX (output) The solution
 */
void LU_Solve(double** ppdA, size_t nSize, size_t* nIndices, 
              double* pdB, double* pdX);

/**
    Singular value decomposition (SVD) of a dense matrix.

    The 'U' part of the factorization is stored in the original matrix.
    A = U W V^T

    @param a  (input/output) The matrix to factor. Stores U on output
    @param m The number of rows
    @param n The number of columns
    @param w (output) The singular values
    @param v (output) The V part of the factorization

    @return 0 if something went wrong, 1 otherwise.
 */
int SVD(double **a, int m, int n, double *w, double **v);

/**
    Solves system A x = b using previously computed SVD.

    @param ppdU The U matrix
    @param pdW The singular values
    @param ppdV The v matrix
    @param m The number of rows
    @param n The number of columns
    @param b The RHS
    @param x (output) The computed solution
 */
void SVD_Solve(const double* const* ppdU, 
               const double* pdW, 
               const double* const* ppdV, 
               size_t nRows, size_t nColumns, 
               const double *pdB, 
               double *pdX);

} // namespace numeric

} // namespace ito33

#endif // #ifdef _ITO33_NUMERIC_DENSEMATRIX_H_
