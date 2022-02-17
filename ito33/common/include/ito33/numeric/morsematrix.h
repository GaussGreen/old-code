/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/morsematrix.h
// Purpose:     sparse matrix class using morse storage
// Created:     2005/01/15
// RCS-ID:      $Id: morsematrix.h,v 1.17 2006/08/19 22:26:42 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/numeric/morsematrix.h
   @brief sparse matrix class using morse storage

   @todo It might speed up the frozen method if the storage is column-major.
   @todo use boost::ublas and rip out completely this manaully written code!
 */
#ifndef _ITO33_NUMERIC_MORSEMATRIX_H_
#define _ITO33_NUMERIC_MORSEMATRIX_H_

#include "ito33/array.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"

#include "ito33/mv/mvvd.h"

namespace ito33
{

namespace numeric
{


/// Class for morse storage for sparse matrix
class MorseStruct 
{
public:
  
  /**
     Ctor initializes the morse storage by using an external list that contains
     the index of non zero elements.

     @param ppnColumnLists the list that contains the morse struct, it may 
                           contains duplicated elements
     @param nNbX the dimension of the underlying matrix
   */
  MorseStruct(std::list<size_t> *ppnColumnLists, size_t nNbX);

  /**
     Ctor intializes the morse storage, it takes also the ownership of the 
     memory. 

     @param nDimension The dimension of the square matrix
     @param pnColumnCounts Array of number of non zero elements for each row
     @param pnColumns Column indexs of non zero elements
   */
  MorseStruct(size_t nDimension, size_t* pnColumnCounts, size_t* pnColumns) 
  {
    m_nDimension = nDimension;
    m_pnColumnCounts = Array<size_t>(pnColumnCounts);
    m_pnColumns = Array<size_t>(pnColumns);
    
    m_ppnColumns = Array<size_t *>(m_nDimension);
    m_nNbElems = 0;
    for (size_t nIdx = 0; nIdx < nDimension; nIdx++)
    {
      m_ppnColumns[nIdx] = pnColumns + m_nNbElems;
      m_nNbElems += m_pnColumnCounts[nIdx];
    }
  }
  
  /// The dimension of the morse matrix
  size_t m_nDimension;

  /// Array of non zero element counts for each line
  Array<size_t> m_pnColumnCounts;
 
  /// Number of non zero elements 
  size_t m_nNbElems;

  /// Array for all column index of non zero elements 
  Array<size_t> m_pnColumns;

  /// Array of position at m_pnColumns for first element of each line  
  Array<size_t *> m_ppnColumns;
  
  friend class MorseMatrix;

}; // class MorseStruct 


/// Sparse matrix class using morse storage
class MorseMatrix 
{
public:

  /// Default ctor initalizes an empty matrix, unusable.
  MorseMatrix()
  {
  }
  
  /**
     Ctor initializes the matrix using a well defined morse struct

     @param pMS A well defined morse struct
   */
  MorseMatrix(const shared_ptr<MorseStruct>& pMS);
  
  /// Accessor to the morse struct
  const shared_ptr<MorseStruct>& GetMorseStruct() const { return m_pMS; }

  /// Accessor directly to the non zero elements
  double& DirectElement(size_t nIdxI, size_t nIdxJ)
  {
    return m_ppdElems[nIdxI][nIdxJ];
  }

  /// initialize the storage of the morse matrix
  void Init(const shared_ptr<MorseStruct>& pMS);

  /// initialize the storage of the morse matrix and fill the elements
  void Init(const shared_ptr<MorseStruct>& pMS, double dInit);
  
  /// operator which can be used to fill the matrix
  double& operator()(size_t, size_t);

  /// operator which can be used to query the matrix
  double operator()(size_t, size_t) const;

  /// operator to be used with the MV library
  CVecteurDouble operator*(const CVecteurDouble& X) const;

  ///{ @name Helper functions to operate on a given row or column

  void ProductMatrixVector
       (const double* pdX, double* pdF, size_t nIdxColumn) const;

  /**
     Sets the elements of the given column to zero.

     @param nIdxColumn The given column.
   */
  void SetColumn(size_t nIdxColumn);

  /**
     Sets the elements of the given column by copying the given matrix.
     
     @param matrix The given matrix from which the elements will be copied
     @param nIdxColumn The given column.
   */
  void SetColumn(const MorseMatrix& matrix, size_t nIdxColumn);

  /**
     Sets the elements of the given row by copying the given matrix.
     
     @param matrix The given matrix from which the elements will be copied
     @param nIdxRow The given row.
   */
  void SetRow(const MorseMatrix& matrix, size_t nIdxRow);

  ///@} 

  /** 
     Do the matrix-vector product.

     @param pdX the vector to be multiplied
     @param pdF the product of the matrix-vector product
     @param bAddon Add the product to pdF or not
   */
  void 
  ProductMatrixVector
  (const double* pdX, double* pdF, bool bAddon = false) const;

  /// Helper function to be used with the MV library for matrix-vector product
  void ProductMatrixVector(const CVecteurDouble& b, CVecteurDouble& res) const; 
 
  /// Function doing transport matrix-vector product
  void ProductTransposeMatrixVector(const double*, double*) const;

  /// Dump the matrix elements into a file to help debug
  void Dump();

  /// Multiply the matrix by a scalar
  void MultiplyBy(double dScalar);


private:

  /// The morse struct
  shared_ptr<MorseStruct> m_pMS;

  /// The storage of the non zero elements
  Array<double> m_pdElems;

  /// Non zero elements corresponding to the morse struct
  Array<double *> m_ppdElems;

  NO_COPY_CLASS(MorseMatrix);

}; // class MorseMatrix


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MORSEMATRIX_H_
