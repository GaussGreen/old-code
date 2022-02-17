/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mempoolmatrixs.h,v $
 * Revision 1.1  2003/12/22 16:47:59  ebenhamou
 * Initial revision
 */

/*! \file mempoolmatrix.h
 *
 *  \brief memory pool type matrix
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _INGPINFRA_MEMPOOLMATRIX_H
#define _INGPINFRA_MEMPOOLMATRIX_H

#include "port.h"
#include "env.h"

#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)

/// kernel
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_MemPool_Matrix
/// \brief 
///		compared to normal matrix classes
///	ARM_MemPool_Matrix is a non contiguous memory 
///	matrix class that allows very fast removal and 
///	insertion of rows (not column!) 
/// since it is implemented  with a memory 
/// pool using vector (using a vector of vector )
/// and a vector of used and unused elements
//////////////////////////////////////////////

class ARM_MemPool_Matrix
{
public:
	typedef vector< double > singleVector;
	typedef singleVector::iterator iterator;

private:
	vector< singleVector > itsMemoryPool;
	vector< size_t > itsUsedElements;
	vector< size_t > itsUnusedElements;

public:
    explicit ARM_MemPool_Matrix(size_t rowsNb=0,size_t colsNb=0);

	/// operator for matrix like accessor!
	/// inline for fast access!
    inline size_t GetRowsNb() const {return itsUsedElements.size();}
    ///inline size_t GetColsNb() const 
    inline size_t GetColsNb() const { return itsUsedElements.empty()? 0 : itsMemoryPool[0].size();}

    inline double operator() (size_t i,size_t j) const;
    inline double& operator() (size_t i,size_t j);

	/// this function clears out all the memory locked by the memory pool!
	void FreeMemoryPool();

	/// function to remove a vector of elements
	/// assumed that the vector indexVec is sorted!
	void Remove_NRowsSorted( const vector<size_t>& indexVec );

	/// function to remove a vector of elements
	/// assumed that the vector indexVec is sorted and with no duplicate
	void Remove_NRows( vector<size_t>& indexVec );

	/// function to remove a vector of elements
	/// does not assumed that the vector indexVec is sorted, 
	///	does not either assume that the vector is with unique index!
	void Remove_NRowsWithDuplicates( vector<size_t>& indexVec );

	/// function to remove one column with index j
// FIXMEFRED: mig.vc8 (22/05/2007 15:58:33):return type missing
	inline void Remove_OneRow( size_t j ) { Remove_NRows( vector<size_t>(1,j) ); }
	/// function to add N columns
	void Add_NRows( size_t newRowNb,size_t colsNb );
	/// function to add one column
	inline void Add_OneRow(size_t colsNb) { Add_NRows(1,colsNb); }

	/// resizeRow using previous function Remove and Add Rows!
    void resizeRows(size_t rowsNb,size_t colsNb );
	/// resize function (very cheap only if colsNb is unchanged!)
    void resize(size_t rowsNb,size_t colsNb);

	/// GetRow non const to allow to give by reference to avoid copy!
	inline singleVector& GetRow( size_t i);

	/// Gets Iterator at the beginning of the i^th row
	inline iterator GetIthRowIterator( size_t i);
};

inline double ARM_MemPool_Matrix::operator() (size_t i,size_t j) const
{
#ifdef _DEBUG
	if( i>= GetRowsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside row range!" );
	if( j>= GetColsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside col range!" );
#endif
	return itsMemoryPool[itsUsedElements[i]][j];
}

inline double& ARM_MemPool_Matrix::operator() (size_t i,size_t j)
{
#ifdef _DEBUG
	size_t rowsNB = GetRowsNb();
	size_t colsNb = GetColsNb();

	if( i>= GetRowsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside row range!" );
	if( j>= GetColsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside col range!" );
#endif
	return itsMemoryPool[itsUsedElements[i]][j];
}


inline ARM_MemPool_Matrix::singleVector& ARM_MemPool_Matrix::GetRow( size_t i)
{
#ifdef _DEBUG
	if( i>= GetRowsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside row range!" );
#endif
	return itsMemoryPool[itsUsedElements[i]];
}

inline ARM_MemPool_Matrix::iterator ARM_MemPool_Matrix::GetIthRowIterator( size_t i)
{
#ifdef _DEBUG
	if( i>= GetRowsNb() )
       	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"trying to access element outside row range!" );
#endif
	return (itsMemoryPool[itsUsedElements[i]]).begin();
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

