/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/densematrix.cpp
// Purpose:     Functions for dealing with dense matrices
// Created:     2005/10/18
// RCS-ID:      $Id: densematrix.cpp,v 1.6 2006/07/06 16:01:12 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/vector.h"
#include "ito33/exception.h"
#include "ito33/numeric/exception.h"
#include "ito33/error.h"
#include "ito33/gettext.h"

#include "ito33/numeric/densematrix.h"

extern const ito33::Error ITO33_BAD_DATA, ITO33_MAX_ITER;

namespace ito33
{

namespace numeric
{

void DenseMatrix::Init(size_t nDim1, size_t nDim2)
{
  ASSERT(nDim1 > 0 && nDim2 > 0);

  m_pd = Array<double>(nDim1 * nDim2);
  m_ppd = Array<double *>(nDim1);
  for (size_t nIdx = 0; nIdx < nDim1; nIdx++)
    m_ppd[nIdx] = m_pd.Get() + nIdx * nDim2;
}

double EuclideanNorm(size_t nSize, double const* const* ppdA)
{
  double dNorm = 0;
  for (size_t nI = 0; nI < nSize; nI++)
    for (size_t nJ = 0; nJ < nSize; nJ++)
      dNorm += ppdA[nI][nJ] * ppdA[nI][nJ];

  return sqrt(dNorm);
}

void LU_Factor(double **ppdA, size_t nSize, size_t *nIndices)
{
  double dTINY = 1.e-14;


  // Find the max element in each row. Used for scaling below.
  std::vector<double> pdVV(nSize);
  size_t nI, nJ, nK;
  for(nI = 0; nI < nSize; nI++)
  {
    double dMax = 0.0;
    for (nJ = 0; nJ < nSize; nJ++)
    {
      double dTemp = fabs( ppdA[nI][nJ] );
      if ( dTemp > dMax ) 
        dMax = dTemp;
    }

    if ( dMax == 0.0) 
      throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Singular matrix.")
          );
    
    pdVV[nI] = 1.0 / dMax;
  }
 
  for (nJ = 0; nJ < nSize; nJ++)
  {
    for (nI = 0; nI < nJ; nI++)
    {
	    double dSum = ppdA[nI][nJ];
	    for (nK = 0; nK < nI; nK++)
	      dSum -= ppdA[nI][nK] * ppdA[nK][nJ];
	    ppdA[nI][nJ] = dSum;
    }

    double dMax = 0.0;
    size_t nImax = 0;
    for (nI = nJ; nI < nSize; nI++)
    {
	    double dSum = ppdA[nI][nJ];
	    for (nK = 0; nK < nJ; nK++)
	      dSum -= ppdA[nI][nK] * ppdA[nK][nJ];

	    ppdA[nI][nJ] = dSum;

      double dTmp = pdVV[nI] * fabs(dSum);
	    if ( dTmp >= dMax )
	    {
	      dMax = dTmp;
	      nImax = nI;
	    }
    }
     
    if (nJ != nImax)
    {
	    for (nK = 0; nK < nSize; nK++)
	    {
	      double dTmp = ppdA[nImax][nK];
	      ppdA[nImax][nK] = ppdA[nJ][nK];
	      ppdA[nJ][nK] = dTmp;
	    }

	    pdVV[nImax] = pdVV[nJ];
    }

    nIndices[nJ] = nImax;
    if (ppdA[nJ][nJ] == 0.0)
    {
	    //printf("zero a[j][j]\n");
	    ppdA[nJ][nJ] = dTINY;
    }

    if (nJ != nSize)
    {
	    double dTmp = 1.0/ppdA[nJ][nJ];
	    for (nI = nJ+1; nI < nSize; nI++)
	      ppdA[nI][nJ] *= dTmp;
    }
  }

}

void LU_Solve(double** ppdA, size_t nSize, size_t* nIndices, double* pdB, double* pdX)
{

  int nI, nJ;
  for (nI = 0; nI < (int) nSize; nI++)
    pdX[nI] = pdB[nI];
  
  for (nI = 0; nI < (int) nSize; nI++)
  {
    size_t nip;
    nip = nIndices[nI];

    double dSum;
    dSum = pdX[nip];
    pdX[nip] = pdX[nI];

    for (nJ = 0; nJ <= nI-1; nJ++) 
      dSum -= ppdA[nI][nJ] * pdX[nJ];

    pdX[nI] = dSum;
  }

  for(nI = nSize-1; nI >= 0; nI--)
  {
    double dSum = pdX[nI];
    for ( nJ = nI+1; nJ < (int) nSize; nJ++)
	    dSum -= ppdA[nI][nJ] * pdX[nJ];
    
    if (ppdA[nI][nI] <= 1.e-14)
      pdX[nI] = 0.0;
    else
      pdX[nI] = dSum/ppdA[nI][nI];
  }
}


static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int SVD(double **a, int m, int n, double *w, double **v)
{
  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  
  if (m < n) 
    return(0);

  std::vector<double> rv1(n);

  /*  avoid warnings */
  l = 0;
  nm = 0;

  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++) 
  {
    /* left-hand reduction */
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m) 
    {
      for (k = i; k < m; k++) 
        scale += fabs(a[k][i]);
      if (scale) 
      {
        for (k = i; k < m; k++) 
        {
          a[k][i] = a[k][i]/scale;
          s += (a[k][i] * a[k][i]);
        }
        f = a[i][i];
        g = (f < 0) ? sqrt(s) : -sqrt(s);
        h = f * g - s;
        a[i][i] = (f - g);
        if (i != n - 1) 
        {
          for (j = l; j < n; j++) 
          {
            for (s = 0.0, k = i; k < m; k++) 
              s += (a[k][i] * a[k][j]);
            f = s / h;
            for (k = i; k < m; k++) 
              a[k][j] += (f * a[k][i]);
          }
        }
        for (k = i; k < m; k++) 
          a[k][i] = (a[k][i]*scale);
      }
    }
    w[i] = (scale * g);

    /* right-hand reduction */
    g = s = scale = 0.0;
    if (i < m && i != n - 1) 
    {
      for (k = l; k < n; k++) 
        scale += fabs(a[i][k]);
      if (scale) 
      {
        for (k = l; k < n; k++) 
        {
          a[i][k] = a[i][k] / scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = (f < 0) ? sqrt(s) : -sqrt(s);
        h = f * g - s;
        a[i][l] = f - g;
        for (k = l; k < n; k++) 
          rv1[k] = a[i][k] / h;
        if (i != m - 1) 
        {
          for (j = l; j < m; j++) 
          {
            for (s = 0.0, k = l; k < n; k++) 
              s += a[j][k] * a[i][k];
            for (k = l; k < n; k++) 
              a[j][k] += s * rv1[k];
          }
        }
        for (k = l; k < n; k++) 
          a[i][k] = a[i][k] * scale;
      }
    }
    double tmp = (fabs(w[i]) + fabs(rv1[i]));
    if ( tmp > anorm)
      anorm = tmp;
  }
  
  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--) 
  {
    if (i < n - 1) 
    {
      if (g) 
      {
        for (j = l; j < n; j++)
          v[j][i] = ((a[i][j] / a[i][l]) / g);
        /* double division to avoid underflow */
        for (j = l; j < n; j++) 
        {
          for (s = 0.0, k = l; k < n; k++) 
            s += (a[i][k] * v[k][j]);
          for (k = l; k < n; k++) 
            v[k][j] += (s * v[k][i]);
        }
      }
      for (j = l; j < n; j++) 
        v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  
  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--) 
  {
    l = i + 1;
    g = w[i];
    if (i < n - 1) 
      for (j = l; j < n; j++) 
        a[i][j] = 0.0;
    if (g) 
    {
      g = 1.0 / g;
      if (i != n - 1) 
      {
        for (j = l; j < n; j++) 
        {
          for (s = 0.0, k = l; k < m; k++) 
            s += (a[k][i] * a[k][j]);
          f = (s / a[i][i]) * g;
          for (k = i; k < m; k++) 
            a[k][j] += f * a[k][i];
        }
      }
      for (j = i; j < m; j++) 
        a[j][i] = a[j][i] * g;
    }
    else 
    {
      for (j = i; j < m; j++) 
        a[j][i] = 0.0;
    }
    ++a[i][i];
  }

  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--) 
  {                             /* loop over singular values */
    for (its = 0; its < 30; its++) 
    {                         /* loop over allowed iterations */
      flag = 1;
      for (l = k; l >= 0; l--) 
      {                     /* test for splitting */
        nm = l - 1;
        if (fabs(rv1[l]) + anorm == anorm) 
        {
          flag = 0;
          break;
        }
        if (fabs(w[nm]) + anorm == anorm) 
          break;
      }
      if (flag) 
      {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) 
        {
          f = s * rv1[i];
          if (fabs(f) + anorm != anorm) 
          {
            g = w[i];
            h = PYTHAG(f, g);
            w[i] = h; 
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < m; j++) 
            {
              y = a[j][nm];
              z = a[j][i];
              a[j][nm] = (y * c + z * s);
              a[j][i] = (z * c - y * s);
            }
          }
        }
      }
      z = w[k];
      if (l == k) 
      {                  /* convergence */
        if (z < 0.0) 
        {              /* make singular value nonnegative */
          w[k] = (-z);
          for (j = 0; j < n; j++) 
            v[j][k] = (-v[j][k]);
        }
        break;
      }
      if (its >= 30) 
      {
        throw EXCEPTION_MSG
        (
          ITO33_MAX_ITER,
          TRANS("Maximum number of iterations exceeded in SVD.")
        );
        
      }
    
      /* shift from bottom 2 x 2 minor */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / ((f<0)?(f-g):(f+g)) ) - h)) / x;
            
      /* next QR transformation */
      c = s = 1.0;
      for (j = l; j <= nm; j++) 
      {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++) 
        {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = (x * c + z * s);
          v[jj][i] = (z * c - x * s);
        }
        z = PYTHAG(f, h);
        w[j] = z;
        if (z) 
        {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = (c * g) + (s * y);
        x = (c * y) - (s * g);
        for (jj = 0; jj < m; jj++) 
        {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = (y * c + z * s);
          a[jj][i] = (z * c - y * s);
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  return(1);
}


void SVD_Solve(const double* const* ppdU, 
               const double* pdW, 
               const double* const* ppdV, 
               size_t nRows, size_t nColumns, 
               const double *pdB, 
               double *pdX)
{
  // A x = b  => U W V^T x = b
  //                     x = (V^T)^-1 W^-1 U^-1 b
  //                     x = V W^-1 U^T b
  // If the matrix is singular (or nearly so) some of the singular
  // values will be zero (small).
  size_t nJJ, nJ, nI;

  std::vector<double> pdTmp(nColumns);
  std::vector<double> pdWTmp(nColumns);

  // Find the largest singular value
  double dWMax = 0.0;
  size_t nIdx;
  for (nIdx = 0; nIdx < nColumns; nIdx++)
  {
    if ( pdW[nIdx] > dWMax )
      dWMax = pdW[nIdx];
  }

  // Zero the small (relative to max element) singular values. Also
  // copy W into WTmp since W is const, and could perhaps be reused.
  double dTol = dWMax * 1.e-10;
  if (dTol < 1.e-24)
    dTol = 1.e-24;
  for (nIdx = 0; nIdx < nColumns; nIdx++)
  {
    if ( pdW[nIdx] < dTol )
      pdWTmp[nIdx] = 0.0;
    else
      pdWTmp[nIdx] = pdW[nIdx];
  }

  // Calculate W^-1 U^T B
  double dSum;
  for (nJ = 0; nJ < nColumns; nJ++)
  {
    dSum = 0.0;
    if ( pdWTmp[nJ] > 0.0 )
    {
      for (nI = 0; nI < nRows; nI++)
        dSum += ppdU[nI][nJ] * pdB[nI];

      dSum /= pdWTmp[nJ];
    }
    pdTmp[nJ] = dSum;
  }

  // matrix multiply by V to get the solution
  for ( nJ = 0; nJ < nColumns; nJ++ )
  {
    dSum = 0.0;
    for ( nJJ = 0; nJJ < nColumns; nJJ++)
      dSum += ppdV[nJ][nJJ] * pdTmp[nJJ];

    pdX[nJ] = dSum;
  }

}


} // namespace numeric

} // namespace ito33
