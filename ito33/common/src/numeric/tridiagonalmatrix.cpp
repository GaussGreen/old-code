/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/tridiagonalmatrix.cpp
// Purpose:     implementation for the tridiagonal matrix class
// Created:     2004/1/28
// RCS-ID:      $Id: tridiagonalmatrix.cpp,v 1.12 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalmatrix.h"

using ito33::numeric::TridiagonalMatrix;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::TridiagonalMatrix);
} 

void TridiagonalMatrix::ProductMatrixVector
                        ( const double *pdX, double *pdF, bool bAddon ) const
{
  if ( bAddon )
  {
    if ( m_nNb == 1 )
    {
      pdF[0] += pdX[0] * m_pdB[0];
      return;
    }

    size_t nIdx;

    nIdx = 0;
    pdF[nIdx] += m_pdB[nIdx] * pdX[nIdx] + m_pdC[nIdx] * pdX[nIdx + 1];

    for (nIdx = 1; nIdx < m_nNb - 1; nIdx++)
      pdF[nIdx] += m_pdA[nIdx] * pdX[nIdx - 1]
                 + m_pdB[nIdx] * pdX[nIdx]
                 + m_pdC[nIdx] * pdX[nIdx + 1];

    pdF[nIdx] += m_pdA[nIdx] * pdX[nIdx - 1] + m_pdB[nIdx] * pdX[nIdx];
  }
  else
  {
    if ( m_nNb == 1 )
    {
      pdF[0] = pdX[0] * m_pdB[0];
      return;
    }

    size_t nIdx;

    nIdx = 0;
    pdF[nIdx] = m_pdB[nIdx] * pdX[nIdx] + m_pdC[nIdx] * pdX[nIdx + 1];

    for (nIdx = 1; nIdx < m_nNb - 1; nIdx++)
      pdF[nIdx] = m_pdA[nIdx] * pdX[nIdx - 1]
                + m_pdB[nIdx] * pdX[nIdx]
                + m_pdC[nIdx] * pdX[nIdx + 1];

    pdF[nIdx] = m_pdA[nIdx] * pdX[nIdx - 1] + m_pdB[nIdx] * pdX[nIdx];
  }
}

void TridiagonalMatrix::ProductTransposeMatrixVector
                        (const double* pdX, double* pdF, bool bAddon) const
{
  if ( bAddon )
  {
    if ( m_nNb == 1 )
    {
      pdF[0] += pdX[0] * m_pdB[0];
      return;
    }

    size_t nIdx;

    nIdx = 0;
    pdF[nIdx] += m_pdB[nIdx] * pdX[nIdx] + m_pdA[nIdx + 1] * pdX[nIdx + 1];

    for (nIdx = 1; nIdx < m_nNb - 1; nIdx++)
      pdF[nIdx] += m_pdC[nIdx - 1] * pdX[nIdx - 1]
                 + m_pdB[nIdx] * pdX[nIdx]
                 + m_pdA[nIdx + 1] * pdX[nIdx + 1];

    pdF[nIdx] += m_pdC[nIdx - 1] * pdX[nIdx - 1] + m_pdB[nIdx] * pdX[nIdx];
  }
  else
  {
    if ( m_nNb == 1 )
    {
      pdF[0] = pdX[0] * m_pdB[0];
      return;
    }

    size_t nIdx;

    nIdx = 0;
    pdF[nIdx] = m_pdB[nIdx] * pdX[nIdx] + m_pdA[nIdx + 1] * pdX[nIdx + 1];

    for (nIdx = 1; nIdx < m_nNb - 1; nIdx++)
      pdF[nIdx] = m_pdC[nIdx - 1] * pdX[nIdx - 1]
                + m_pdB[nIdx] * pdX[nIdx]
                + m_pdA[nIdx + 1] * pdX[nIdx + 1];

    pdF[nIdx] = m_pdC[nIdx - 1] * pdX[nIdx - 1] + m_pdB[nIdx] * pdX[nIdx];
  }
}

TridiagonalMatrix* TridiagonalMatrix::GetTranspose() const
{
  TridiagonalMatrix* pMatrix = new TridiagonalMatrix(m_nNb);
  
  for (size_t n = 0; n < m_nNb; n++)
    pMatrix->m_pdB[n] = m_pdB[n];

  pMatrix->m_pdA[0] = pMatrix->m_pdC[m_nNb - 1] = 0;

  for (size_t n = 0; n < m_nNb - 1; n++)
  {
    pMatrix->m_pdC[n] = m_pdA[n + 1];
    pMatrix->m_pdA[n + 1] = m_pdC[n];
  }

  return pMatrix;
}
