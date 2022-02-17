/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file env.cpp
 *
 *  \brief various routines to manipulate vectors
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/vectormanip.h"
#include "gpbase/checkinputs.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"
#include <algorithm>
#include <vector>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Routine: MergeSortedVectorNoDuplicates
///	Action : function to merge two sorted vector in increasing order
///			remove duplicates as well!
////////////////////////////////////////////////////
ARM_GP_Vector* MergeSortedVectorNoDuplicates( const ARM_GP_Vector& vec1, const ARM_GP_Vector& vec2 )
{
#if defined(__GP_CHECK_SORTED_VECTOR)
	CheckVectorIncreasing( vec1, "vec1", "MergeSortedVector",__LINE__,__FILE__);
	CheckVectorIncreasing( vec2, "vec2", "MergeSortedVector",__LINE__,__FILE__);
#endif

	CC_NS(std,vector)<double> mergeVec(vec1.size() + vec2.size() );

	/// merge model and eventTime
	CC_NS(std,merge)(vec1.begin(), vec1.end(),
		vec2.begin(), vec2.end(),
		mergeVec.begin() );

	CC_NS(std,vector)<double>::iterator iter =
		( CC_NS(std,unique) ( mergeVec.begin(), mergeVec.end() ) );
	
	return new ARM_GP_Vector( iter-mergeVec.begin(), &mergeVec[0] );
}



////////////////////////////////////////////////////
///	Routine: MergeSortedVector
///	Action : function to merge two sorted vector in increasing order
////////////////////////////////////////////////////
ARM_GP_Vector* MergeSortedVectorWithDuplicates( const ARM_GP_Vector& vec1, const ARM_GP_Vector& vec2 )
{
#if defined(__GP_CHECK_SORTED_VECTOR)
	CheckVectorIncreasing( vec1, "vec1", "MergeSortedVector",__LINE__,__FILE__);
	CheckVectorIncreasing( vec2, "vec2", "MergeSortedVector",__LINE__,__FILE__);
#endif

	CC_NS(std,vector)<double> mergeVec(vec1.size() + vec2.size());

	/// merge model and eventTime
	CC_NS(std,merge)(vec1.begin(), vec1.end(),
		vec2.begin(), vec2.end(),
		mergeVec.begin() );

	return new ARM_GP_Vector( mergeVec.size(), &mergeVec[0] );
}

////////////////////////////////////////////////////
///	Routine: SortSTLBased
///	Action : sort using the efficient STL sort algorithm
////////////////////////////////////////////////////
ARM_GP_Vector* SortSTLBased( const ARM_GP_Vector& vec )
{
	ARM_GP_Vector* result = new ARM_GP_Vector(vec);
	CC_NS(std,sort)( result->begin(),result->end() );
	return result;
}


////////////////////////////////////////////////////
///	Routine: SortSTLBased
///	Action : remove duplicates
////////////////////////////////////////////////////
ARM_GP_Vector* VectorUnique(const ARM_GP_Vector& vec )
{
#if defined(__GP_CHECK_SORTED_VECTOR)
	CheckVectorIncreasing( vec, "vec", "VectorUnique",__LINE__,__FILE__);
#endif

	CC_NS(std,vector)<double> result( vec.size() );
	CC_NS(std,copy)( vec.begin(), vec.end(), result.begin() );
	CC_NS(std,vector)<double>::iterator 
		vbegin = result.begin(),
		vend = result.end(),
		pos = CC_NS(std,unique)( vbegin, vend );
	return new ARM_GP_Vector(pos-result.begin(), &result[0] );
}


////////////////////////////////////////////////////
///	Routine: VectorToString
///	Action : stringify a vector
////////////////////////////////////////////////////
CC_NS(std,string) VectorToString( const ARM_GP_Vector& vec )
{
	CC_Ostringstream os;
	os << "\nVector of size =" << vec.size() << "\n";
	for( size_t i=0; i<vec.size(); ++i )
		os << "\t" << vec[i];
	os << "\n";
	return os.str();
}


bool ExistsInVector( const ARM_GP_Vector& dateVector, double date )
{
	size_t k=0;
	while( k<dateVector.size() && dateVector.Elt(k) < date )
		++k;
	return (k < dateVector.size() && abs( date - dateVector.Elt(k) ) < K_DOUBLE_TOL ) ? true : false;
}


size_t IdxFromValue( const ARM_GP_Vector& dateVector, double date )
{
	size_t k=0;
	while( k<dateVector.size() && dateVector.Elt(k) < date )
		++k;

	if( k >= dateVector.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_ProcessBuilderPDE::IdxFromDate: date not found in the map!");

	return k;
}


size_t IdxFromValueWithTol( const ARM_GP_Vector& dateVector, double date, double daysNb , bool NoThrow)
{
	int i=-1, N = dateVector.size();

	while( (abs( date - dateVector.Elt(++i) ) > daysNb ) && i<N ) {}

	if( i == N && NoThrow == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_ProcessBuilderPDE::IdxFromDate: date not found in the map!");

	return i == N ? -1 : i;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

