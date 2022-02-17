/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/densesolver_direct.cpp
// Purpose:     Some functions to solve linear systems for dense matrix
// Author:      Nabil
// Created:     2006/07/04
// RCS-ID:      $Id: densesolver_direct.cpp,v 1.2 2006/07/12 09:45:19 nabil Exp $
// Copyright:   (c) 2003-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/vector.h"
#include "ito33/numeric/densesolver_direct.h"

using namespace std;

namespace ito33
{

namespace numeric
{

bool Solve2DLinearSystem(double dA11, double dA12, double dB1, double dA21, 
                         double dA22, double dB2, double& dX1, double& dX2)
{
  double
    dDeterminant = dA22 * dA11 - dA21 * dA12;

  if (dDeterminant != 0.)
  {
    dX1 = (dA22 * dB1 - dA12 * dB2) / dDeterminant;
    dX2 = (dA11 * dB2 - dA21 * dB1) / dDeterminant;

    return true;
  }
  else
    return false;
}

bool SolveLinearSystem(double** ppdA, double** ppdB, int iN, int iM)
{
  int
    I,
    J,
    K,
    ICol = 0,
    IRow = 0;

  vector<int> 
    IPiv (iN),
    IndxC (iN),
    IndxR (iN);

  double
    Big,
    PivInv,
    Temp;

  for (I = 0; I < iN; I++)
  {
    Big = 0.0;
    for (J = 0; J < iN; J++)
      if (IPiv[J] != 1)
        for (K = 0; K < iN; K++)
          if (IPiv[K] == 0)
          {
            if (fabs(ppdA[J][K]) >= Big)
              Big = fabs(ppdA[IRow = J][ICol = K]);
          }
          else if (IPiv[K] > 1)
          {
            return false;
          }
    IPiv[ICol]++;
    if (IRow != ICol)
    {
      for (J = 0; J < iN; J++)
      {
        Temp = ppdA[IRow][J];
        ppdA[IRow][J] =	ppdA[ICol][J];
        ppdA[ICol][J] = Temp;
      }
      for (J = 0; J < iM; J++)
      {
        Temp = ppdB[IRow][J];
        ppdB[IRow][J] =	ppdB[ICol][J];
        ppdB[ICol][J] = Temp;
      }
    }
    IndxR[I] = IRow;
    IndxC[I] = ICol;
    if (ppdA[ICol][ICol] == 0.0)
    {
      return false;
    }
    PivInv = 1.0 / ppdA[ICol][ICol];
    ppdA[ICol][ICol] = 1.0;
    for (J = 0; J < iN; J++)
      ppdA[ICol][J] *= PivInv;
    for (J = 0; J < iM; J++)
      ppdB[ICol][J] *= PivInv;
    for (J = 0; J < iN; J++)
      if (J != ICol)
      {
        Temp = ppdA[J][ICol];
        ppdA[J][ICol] = 0.0;
        for (K = 0; K < iN; K++)
          ppdA[J][K] -=	ppdA[ICol][K] * Temp;
        for (K = 0; K < iM; K++)
          ppdB[J][K] -=	ppdB[ICol][K] * Temp;
      }
  }
  for (I = iN - 1; I >= 0; I--)
    if (IndxR[I] != IndxC[I])
      for (J = 0; J < iN; J++)
      {
        Temp = ppdA[J][IndxR[I]];
        ppdA[J][IndxR[I]] =	ppdA[J][IndxC[I]];
        ppdA[J][IndxC[I]] = Temp;
      }

  return true;
}

} // namespace numeric

} // namespce ito33
