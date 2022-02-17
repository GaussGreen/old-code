/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/densesolver_direct.h
// Purpose:     Some functions to solve linear systems for dense matrix
// Author:      Nabil
// Created:     2006/07/04
// RCS-ID:      $Id: densesolver_direct.h,v 1.1 2006/07/11 14:28:11 nabil Exp $
// Copyright:   (c) 2003-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/densesolver_direct.h
    @brief Some functions to solve linear systems for dense matrix.
    
 */

#ifndef _ITO33_NUMERIC_DENSESOLVER_DIRECT_H_
#define _ITO33_NUMERIC_DENSESOLVER_DIRECT_H_

namespace ito33
{

namespace numeric
{

/**
   Solves a two dimensional linear system.

   Solves the two dimensional linear system A.X = B.
   A is a 2x2 matrix and Aij is the element of the ith row and jth column of A.

   @param dA11 A(1,1)
   @param dA12 A(1,2)
   @param dB1 B(1)
   @param dA21 A(2,1)
   @param dA22 A(2,2)
   @param dB2 B(2)
   @param dX1 (output) X(1)
   @param dX2 (output) X(2)

   @return true if successful, false otherwise.
 */
bool Solve2DLinearSystem(double dA11, double dA12, double dB1, double dA21, 
                        double dA22, double dB2, double& dX1, double& dX2);

/**
   Solves the linear system A.X = B using the Gauss-Jordan elimination.

   A is a NxN matrix and B is a NxM matrix. 
   See Numerical recipes in C p36 for more details.

   @param ppdA (Input and output) The matrix A. It becomes the inverse of A.
   @param ppdB (Input and output) The matrix B. It becomes the unknown X.
   @param dX1 (output) X(1)
   @param dX2 (output) X(2)

   @return true if successful, false otherwise.
 */
bool SolveLinearSystem(double** ppdA, double** ppdB, int iN, int iM);

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_DENSESOLVER_DIRECT_H_
