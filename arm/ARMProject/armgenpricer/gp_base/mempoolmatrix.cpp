/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mempoolmatrix.cpp,v $
 * Revision 1.1  2003/12/08 16:44:43  ebenhamou
 * Initial revision
 *
 */

/*! \file mempoolmatrix.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/mempoolmatrix.h"
#include "gpbase/ostringstream.h"
#include <algorithm>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MemPool_Matrix
///	Routine: Constructor
///	Returns: 
///	Action : Constructor allocating memory, the memory pool is
///			a vector of size rowsNb of vector of size colsNb
///			the unusedElements is empty since everything is used!
///			there is one little exception which is if the rowsNb == 0
///			but colsNB = 0... in that case, we put 1 vector in memoryPool
///			and put it as unused!
////////////////////////////////////////////////////
ARM_MemPool_Matrix::ARM_MemPool_Matrix(size_t rowsNb,size_t colsNb)
:	itsMemoryPool(rowsNb, singleVector(colsNb) ), itsUnusedElements()
{
	/// at start we use all elements from the memory pool!
	/// in the increasing order
	size_t i;
	itsUsedElements.reserve( rowsNb );
	for(i=0;i<rowsNb; ++i )
		itsUsedElements.push_back(i);

	/// in case of empty rowsNb
	if(rowsNb==0 && colsNb!=0)
	{
		itsMemoryPool.push_back(singleVector(colsNb));
		itsUnusedElements.push_back(0);
	}
};



////////////////////////////////////////////////////
///	Class  : ARM_MemPool_Matrix
///	Routine: FreeMemoryPool
///	Returns: 
///	Action : the freeMemoryPool function removes
///			all the unused elements and copy to a 
///			compact and contiguous memory block the used
///			elements
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::FreeMemoryPool()
{
	if( !itsUnusedElements.empty() )
	{
		/// tmpMemoryPool to receive all the used elements!
		size_t rowsNb = itsUsedElements.size();
		vector< singleVector > tmpMemoryPool( rowsNb, singleVector( GetRowsNb() ) );
		
		/// first copies the value to a temp and swap with the memory pool
		size_t i;
		for(i=0;i<rowsNb;++i)
			tmpMemoryPool[i] = itsMemoryPool[itsUsedElements[i]];
		itsMemoryPool.swap(tmpMemoryPool);
		
		/// second remove all unused elements
		itsUnusedElements = vector< size_t >();

		/// third update the used elements
		/// since we copy elements exactly in the order
		/// specified by used Elements
		/// the used elements are in order 1..2...n
		itsUsedElements = vector< size_t >();
		itsUsedElements.reserve(rowsNb);
		for(i=0; i<rowsNb; ++i)
			itsUsedElements[i] = i;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MemPool_Matrix
///	Routines: Remove_NRowsSorted
///	Returns :
///	Action  : remove all the columns with index contained in indexVec!
///				assumed that the indexes to remove are sorted
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::Remove_NRowsSorted( const vector<size_t>& indexVec )
{
	size_t i;

#ifdef __GP_STRICT_VALIDATION
	if( !indexVec.empty() )
	{
		if( indexVec.size() > itsUsedElements.size() )
       		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"trying to remove more than what is used!" );

		if( indexVec[0] >= GetRowsNb() )
       		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"trying to access element outside row range in indexVec[0]!" );
		
		for(i=1; i<indexVec.size(); ++i )
		{
			if( indexVec[i] <= indexVec[i-1] )
			{
				CC_Ostringstream os;
				os << "the vector is not with strictly increasing elements for " << i;
	       		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
		}

		if( indexVec[indexVec.size()-1] >= GetRowsNb() )
		{
			CC_Ostringstream os;
			os << "trying to access element outside row range in indexVec[" << indexVec.size()-1 <<"]";
       		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
#endif

	if( !indexVec.empty() )
	{
		size_t 
			j,
			end,
			index,
            indexVecSize = indexVec.size();

        /// forecase the new size
        itsUnusedElements.reserve(itsUnusedElements.size()+indexVecSize);

		for(i=0; i<indexVecSize-1; ++i)
		{
			/// get the index corresponding to the column to remove
			/// shift all the used element by 1
			///	and push the index to the unused elements
			end		= indexVec[i+1];
			index	= itsUsedElements[indexVec[i]];

			for(j=indexVec[i]-i; j<end-i; ++j)
				itsUsedElements[j]=itsUsedElements[j+i+1];

			/// this is now an unused index
			itsUnusedElements.push_back(index);
		}

		/// last one is specific on the end
		index	= itsUsedElements[indexVec[i]];
		end		= itsUsedElements.size()-1;

		for(j=indexVec[i]-i; j<end-i; ++j)
			itsUsedElements[j]=itsUsedElements[j+i+1];

		/// move to unused index
		itsUnusedElements.push_back(index);

		/// resize the used elements
		itsUsedElements.resize(itsUsedElements.size()-indexVecSize );
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MemPool_Matrix
///	Routines: Remove_NRows
///	Returns :
///	Action  : remove all the columns with index contained in indexVec!
///				indexVex should not be sorted... but it will be after..
///				hence the non const argument! Warning it assumes that the
///				indexes are with no duplicates! Check in strict validation mode
///				but is at you own risk in release mode
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::Remove_NRows( vector<size_t>& indexVec )
{
	if( indexVec.size() > 1 )
	{
		bool needToBeSorted = false;

		for(size_t i=1; i<indexVec.size(); ++i )
		{

#ifdef __GP_STRICT_VALIDATION
			if( indexVec[i] == indexVec[i-1] )
			{
				CC_Ostringstream os;
				os << "indexes [" << i << "] and [" << i-1 << "] equal with value = " << indexVec[i-1];
       			ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
			}
#endif
			if( indexVec[i] < indexVec[i-1] )
			{
				needToBeSorted= true;
				break;
			}
		}

		if( needToBeSorted )
			CC_NS(std,sort)(indexVec.begin(), indexVec.end() );
	}
	
	Remove_NRowsSorted( indexVec );
}

////////////////////////////////////////////////////
///	Class   : Remove_NRowsWithDuplicates
///	Routines: Remove_NRows
///	Returns :
///	Action  : remove all the columns with index contained in indexVec!
///				indexVex should not be sorted... and may contain duplicates
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::Remove_NRowsWithDuplicates( vector<size_t>& indexVec )
{
	if( indexVec.size() > 1 )
	{
		CC_NS(std,sort)( indexVec.begin(), indexVec.end() );
		CC_NS(std,vector)<size_t>::iterator pos  = CC_NS(std,unique)( indexVec.begin(), indexVec.end() );
		indexVec.resize(pos-indexVec.begin());
	}

	Remove_NRowsSorted( indexVec );
}

////////////////////////////////////////////////////
///	Class   : ARM_MemPool_Matrix
///	Routines: Add_NRows
///	Returns :
///	Action  : remove column with index j
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::Add_NRows( size_t newRowsNb,size_t colsNb )
{
	if( newRowsNb > 0 )
	{
		/// can we use the unused elements to add the n cols?
		if( newRowsNb <= itsUnusedElements.size() )
		{
			size_t i;
			size_t usedElmSize = itsUsedElements.size();
			
			/// reserve the correct memory
			itsUsedElements.reserve(usedElmSize+newRowsNb);
			for(i=0;i<newRowsNb; ++i)
				itsUsedElements.push_back( itsUnusedElements[itsUnusedElements.size()-i-1] );
			
			/// remove the last elements from the unused elements!
			itsUnusedElements.resize(itsUnusedElements.size()-newRowsNb);
		}
		
		/// we first need to used all the unused elements 
		/// we do it this way since we can benefit from the 
		/// dynamic growth of vector and hence may not need
		/// to create a new vector!
		else
		{
			/// first get all the unused elements
			if( itsUnusedElements.size() )
			{
				newRowsNb -= itsUnusedElements.size();
				Add_NRows(itsUnusedElements.size(),colsNb);
			}

			size_t i;

			/// now do the real job of allocating more memory!
			/// reserve the necessary memory
			/// this may be a no effect operation
			/// if the memory is already reserved!
			itsMemoryPool.reserve(GetRowsNb()+newRowsNb);
			itsUsedElements.reserve(GetRowsNb()+newRowsNb);
			size_t currentIndex = itsMemoryPool.size();
			for( i=0; i<newRowsNb; ++i)
			{
				itsMemoryPool.push_back( singleVector( colsNb ) );
				itsUsedElements.push_back( currentIndex );
				++currentIndex;
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MemPool_Matrix
///	Routines: resizeRows
///	Returns :
///	Action  : smart resize using the property of the memory pool
///				use the cheap functions Add, Remove Rows... 
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::resizeRows(size_t rowsNb,size_t colsNb)
{
	if(rowsNb<GetRowsNb())
	{
		/// to be STL compliant, we need to remvoe the latest
		///	indexes, hence the need to remove
		///	itsUsedElements.size()-1 ... to itsUsedElements.size()-indexToRemoveVec

		size_t elmToRemoveNb = GetRowsNb()-rowsNb;
		vector<size_t> indexToRemoveVec(elmToRemoveNb);
		size_t lastUsedIndex = itsUsedElements.size()-1;
		size_t i;
		for(i=0;i<elmToRemoveNb; ++i)
			indexToRemoveVec[i] = lastUsedIndex-i;

		/// just now call remove!
		Remove_NRows( indexToRemoveVec );
	}
	else
	{
		size_t elmToAddNb = rowsNb - GetRowsNb();
		Add_NRows(elmToAddNb,colsNb);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MemPool_Matrix
///	Routines: resize
///	Returns : void
///	Action  :	1) if the colsNb is unchanged, use the cheap resizeRows function
///		which uses underneath functions Add, Remove Rows... 
///				2) else resize first columns (using STL resize) to be in the first case
///		and use then the resizeRows
///
///		resize is STL compliant in the sense that it 
///		keeps previous values!
////////////////////////////////////////////////////
void ARM_MemPool_Matrix::resize(size_t rowsNb,size_t colsNb)
{
	if( colsNb == GetColsNb() )
		resizeRows(rowsNb,colsNb);
	else
	{
		for(size_t i=0;i<itsMemoryPool.size(); ++i)
			itsMemoryPool[i].resize(colsNb);
		resizeRows(rowsNb,colsNb);
	}
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
