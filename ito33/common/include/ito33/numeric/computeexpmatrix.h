/////////////////////////////////////////////////////////////////////////////
// Name:      ito33/numeric/computeexpmatrix.h
// Purpose:   Compute the exp of a matrix
// Created:   2006/07/05
// RCS-ID:    $Id: computeexpmatrix.h,v 1.1 2006/07/06 16:01:31 wang Exp $
// Copyright: (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/computeexpmatrix.h
    @brief Compute the exp of a matrix 
 */

#ifndef _ITO33_NUMERIC_COMPUTEEXPMATRIX_H_
#define _ITO33_NUMERIC_COMPUTEEXPMATRIX_H_

namespace ito33
{

namespace numeric
{

/**
    Computes the exp of a matrix (small and dense) using taylor series.

    @param nNb The dimension of the matrix
    @param ppdA The matrix
    @param ppdExpA The output matrix
 */
void ComputeExpMatrix(size_t nNb, double const* const* ppdA, double** ppdExpA);

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_COMPUTEEXPMATRIX_H_
