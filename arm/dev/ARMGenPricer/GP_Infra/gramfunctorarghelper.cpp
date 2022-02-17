/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctorarghelper.cpp,v $
 * Revision 1.1  2004/05/01 07:52:11  ebenhamou
 * Initial revision
 *
 */

/*! \file gramfunctorarghelper.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */


#include "gpinfra/gramfunctorarghelper.h"
#include "gpinfra/gramfunctorarg.h"
#include "gpbase/gpmatrix.h"

CC_BEGIN_NAMESPACE( ARM )

const int ARM_GramFctorArgVectorHelper::IndexNotFound = -1;

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgVectorHelper
///	Routine: FindFirstVectorIndex
///	Returns: int
///	Action : returns the index of the first vector looking at all args
///				with position given in position 
///				if no vector is found returns
///				ARM_GramFctorArgVectorHelper::IndexNotFound(-1)!
////////////////////////////////////////////////////
int ARM_GramFctorArgVectorHelper::FindFirstVectorIndex( const ARM_GramFctorArgVector& arg, const vector<size_t>& position )
{
	int result = ARM_GramFctorArgVectorHelper::IndexNotFound;
	for( size_t i=0; i<position.size(); ++i )
	{
		if( arg[position[i] ].GetType() == GFAT_VECTOR_TYPE )
		{
			result = position[i];
			break;
		}
	}
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgVectorHelper
///	Routine: FindFirstVectorIndex (version with iterator flavor)
///	Returns: int
///	Action : returns the index of the first vector looking at all args
///				with position given in position 
///				if no vector is found returns
///				ARM_GramFctorArgVectorHelper::IndexNotFound(-1)!
////////////////////////////////////////////////////

int ARM_GramFctorArgVectorHelper::FindFirstVectorIndex( const ARM_GramFctorArgVector& arg, const size_t* positionIteratorStart, const size_t* positionIteratorEnd )
{
	int result = ARM_GramFctorArgVectorHelper::IndexNotFound;
	
	for( const size_t* iter=positionIteratorStart; iter<positionIteratorEnd; ++iter )
	{
		if( arg[*iter ].GetType() == GFAT_VECTOR_TYPE )
		{
			result = *iter;
			break;
		}
	}
	return result;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



