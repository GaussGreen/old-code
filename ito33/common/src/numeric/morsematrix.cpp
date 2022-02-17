/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/morsematrix.cpp
// Purpose:     Implementation for the sparse matrix using morse storage
// Created:     2005/01/15
// RCS-ID:      $Id: morsematrix.cpp,v 1.16 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/debugparameters.h"

#include "ito33/numeric/morsematrix.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(numeric::MorseStruct);
  ITO33_IMPLEMENT_AUTOPTR(numeric::MorseMatrix);
}

namespace ito33
{

namespace numeric
{


MorseStruct::MorseStruct(std::list<size_t> *ppnColumnLists, size_t nNbX)
{
  // The dimension of the underlying matrix
  m_nDimension = nNbX;

  // Get the number of nonzero entries by getting the unique,
  // sorted column numbers for the lists for the rows
  m_pnColumnCounts = Array<size_t>(nNbX);

  m_nNbElems = 0;
  for (size_t nIdx = 0; nIdx < nNbX; nIdx++)
  {
    ppnColumnLists[nIdx].sort();
    ppnColumnLists[nIdx].unique();

    m_pnColumnCounts[nIdx] = ppnColumnLists[nIdx].size();
    m_nNbElems += m_pnColumnCounts[nIdx];
  }

  // Allocate the matrix that holds the column entries.  It is allocated
  // as an array, but accessed as a matrix
  m_pnColumns = Array<size_t>(m_nNbElems);
  m_ppnColumns = Array<size_t *>(nNbX);

  // Setup the morse struct
  m_nNbElems = 0;
  for (size_t nIdx = 0; nIdx < nNbX; nIdx++)   
  {
    m_ppnColumns[nIdx] = m_pnColumns.Get() + m_nNbElems;

    std::list<size_t>::iterator listIter;
    size_t nCount = 0;
    for (listIter = ppnColumnLists[nIdx].begin();
         listIter != ppnColumnLists[nIdx].end();
         ++listIter)
    {
      m_ppnColumns[nIdx][nCount++] = *listIter;
    }
    
    m_nNbElems += m_pnColumnCounts[nIdx];
  }
}

MorseMatrix::MorseMatrix(const shared_ptr<MorseStruct>& pMS) : m_pMS(pMS)
{
  Init(pMS);
}

void MorseMatrix::Init(const shared_ptr<MorseStruct>& pMS) 
{
  m_pMS = pMS; 

  m_pdElems = Array<double>(pMS->m_nNbElems);
  m_ppdElems = Array<double *>(m_pMS->m_nDimension);

  size_t nNbElements = 0;
  for (size_t nIdx = 0; nIdx < m_pMS->m_nDimension; nIdx++)
  {
    m_ppdElems[nIdx] = m_pdElems.Get() + nNbElements;
    nNbElements += pMS->m_pnColumnCounts[nIdx];
  }
}

void MorseMatrix::Init(const shared_ptr<MorseStruct>& pMS, double dInit) 
{
  Init(pMS);

  for (size_t nIdx = 0; nIdx < m_pMS->m_nNbElems; nIdx++)
    m_pdElems[nIdx] = dInit;
}

double& MorseMatrix::operator()(size_t nIdxI, size_t nIdxJ) 
{
  size_t nIdx;

  for (nIdx = 0; 
          nIdx < m_pMS->m_pnColumnCounts[nIdxI]
       && m_pMS->m_ppnColumns[nIdxI][nIdx] < nIdxJ; 
       nIdx++)
    ;

  ASSERT(nIdxJ == m_pMS->m_ppnColumns[nIdxI][nIdx]);
    
  return m_ppdElems[nIdxI][nIdx];
}

double MorseMatrix::operator()(size_t nIdxI, size_t nIdxJ) const
{
  size_t nIdx;

  for (nIdx = 0; 
          nIdx < m_pMS->m_pnColumnCounts[nIdxI]
       && m_pMS->m_ppnColumns[nIdxI][nIdx] < nIdxJ; 
       nIdx++)
    ;

  if ( nIdxJ == m_pMS->m_ppnColumns[nIdxI][nIdx])
    return m_ppdElems[nIdxI][nIdx];
  else
    return 0.0;
}

void MorseMatrix::ProductMatrixVector
                  (const double* pdX, double* pdF, bool bAddon) const
{
  if ( bAddon )
  {
    for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    {      
      for (size_t nIdxJ = 0; nIdxJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJ++)
        pdF[nIdxI] += m_ppdElems[nIdxI][nIdxJ]
                    * pdX[m_pMS->m_ppnColumns[nIdxI][nIdxJ]];
    }
  }
  else
  {
    for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    {
      pdF[nIdxI] = 0;
      
      for (size_t nIdxJ = 0; nIdxJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJ++)
        pdF[nIdxI] += m_ppdElems[nIdxI][nIdxJ]
                    * pdX[m_pMS->m_ppnColumns[nIdxI][nIdxJ]];
    }
  }
}

void MorseMatrix::ProductMatrixVector
                  (const CVecteurDouble& X, CVecteurDouble& F) const
{
  ProductTransposeMatrixVector( X.get(), F.get() );
}

CVecteurDouble MorseMatrix::operator*(const CVecteurDouble& X) const
{
  CVecteurDouble F(X.size()); 

  ProductMatrixVector(X, F);
  
  return F;
}

void MorseMatrix::ProductTransposeMatrixVector
                  (const double* pdX, double* pdF) const
{
  size_t nIdxI;

  for (nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    pdF[nIdxI] = 0;

  for (nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++)
    for (size_t nIdxJJ = 0; nIdxJJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJJ++)
      pdF[m_pMS->m_ppnColumns[nIdxI][nIdxJJ]] += m_ppdElems[nIdxI][nIdxJJ]
                                               * pdX[nIdxI];
}

void MorseMatrix::ProductMatrixVector
                    (const double* pdX, double* pdF, size_t nIdxColumn) const
{
  for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++)
    for (size_t nIdxJJ = 0; nIdxJJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJJ++) 
    {
      size_t nIdxJ = m_pMS->m_ppnColumns[nIdxI][nIdxJJ];
      if ( nIdxJ == nIdxColumn )
      {
        pdF[nIdxI] += m_ppdElems[nIdxI][nIdxJJ] * pdX[nIdxJ];
        break;
      }
    }
}

void MorseMatrix::SetColumn(size_t nIdxColumn)
{
  for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    for (size_t nIdxJJ = 0; nIdxJJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJJ++)
      if (m_pMS->m_ppnColumns[nIdxI][nIdxJJ] == nIdxColumn)
      {
        m_ppdElems[nIdxI][nIdxJJ] = 0.;
        break;
      }
}

void MorseMatrix::SetColumn(const MorseMatrix& matrix, size_t nIdxColumn)
{
  ASSERT(m_pMS == matrix.m_pMS);

  for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    for (size_t nIdxJJ = 0; nIdxJJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJJ++)
      if (m_pMS->m_ppnColumns[nIdxI][nIdxJJ] == nIdxColumn)
      {
        m_ppdElems[nIdxI][nIdxJJ] = matrix.m_ppdElems[nIdxI][nIdxJJ];
        break;
      }
}

void MorseMatrix::SetRow(const MorseMatrix& matrix, size_t nIdxRow)
{
  ASSERT(m_pMS == matrix.m_pMS);
  
  for (size_t nIdxJJ = 0; nIdxJJ < m_pMS->m_pnColumnCounts[nIdxRow]; nIdxJJ++)
    m_ppdElems[nIdxRow][nIdxJJ] = matrix.m_ppdElems[nIdxRow][nIdxJJ];
}

void MorseMatrix::MultiplyBy(double dScalar)
{
  for (size_t nIdx = 0; nIdx < m_pMS->m_nNbElems; nIdx++)
    m_pdElems[nIdx] *= dScalar;
}

void MorseMatrix::Dump()
{
  std::string filename(DebugParameters::GetDebugDir());
  
  filename +=  "morsematrix.txt";

  std::ofstream of(filename.c_str());

  for (size_t nIdxI = 0; nIdxI < m_pMS->m_nDimension; nIdxI++) 
    for (size_t nIdxJ = 0; nIdxJ < m_pMS->m_pnColumnCounts[nIdxI]; nIdxJ++)
    {
      if (m_ppdElems[nIdxI][nIdxJ] != 0)
        of << nIdxI << '\t' << m_pMS->m_ppnColumns[nIdxI][nIdxJ] << '\t'
           << m_ppdElems[nIdxI][nIdxJ] << '\n';
    }
}


} // namespace numeric

} // namespace ito33
