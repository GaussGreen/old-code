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
std::vector<double>* MergeSortedVectorNoDuplicates( const std::vector<double>& vec1, const std::vector<double>& vec2 )
{
#if defined(__GP_CHECK_SORTED_VECTOR)
	CheckVectorIncreasing( vec1, "vec1", "MergeSortedVector",__LINE__,__FILE__);
	CheckVectorIncreasing( vec2, "vec2", "MergeSortedVector",__LINE__,__FILE__);
#endif

	std::vector<double> mergeVec(vec1.size() + vec2.size() );

	/// merge model and eventTime
	std::merge(vec1.begin(), vec1.end(),vec2.begin(), vec2.end(), mergeVec.begin() );

	std::vector<double>::iterator it = std::unique(mergeVec.begin(), mergeVec.end());

	mergeVec.resize( std::distance(mergeVec.begin(),it) );

	return new std::vector<double>( mergeVec.begin(), mergeVec.end() );
}



////////////////////////////////////////////////////
///	Routine: MergeSortedVector
///	Action : function to merge two sorted vector in increasing order
////////////////////////////////////////////////////
std::vector<double>* MergeSortedVectorWithDuplicates( const std::vector<double>& vec1, const std::vector<double>& vec2 )
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

	return new std::vector<double>( mergeVec );
}

////////////////////////////////////////////////////
///	Routine: SortSTLBased
///	Action : sort using the efficient STL sort algorithm
////////////////////////////////////////////////////
std::vector<double>* SortSTLBased( const std::vector<double>& vec )
{
	std::vector<double>* result = new std::vector<double>(vec);
	CC_NS(std,sort)( result->begin(),result->end() );
	return result;
}


////////////////////////////////////////////////////
///	Routine: SortSTLBased
///	Action : remove duplicates
////////////////////////////////////////////////////
std::vector<double>* VectorUnique(const std::vector<double>& vec )
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
	return new std::vector<double>(result.begin(), result.end() );
}


////////////////////////////////////////////////////
///	Routine: VectorToString
///	Action : stringify a vector
////////////////////////////////////////////////////
CC_NS(std,string) VectorToString( const std::vector<double>& vec )
{
	CC_Ostringstream os;
	os << "\nVector of size =" << vec.size() << "\n";
	for( size_t i=0; i<vec.size(); ++i )
		os << "\t" << vec[i];
	os << "\n";
	return os.str();
}


bool ExistsInVector( const std::vector<double>& dateVector, double date )
{
	size_t k=0;
	while( k<dateVector.size() && dateVector[k] < date )
		++k;
	return (k < dateVector.size() && abs( date - dateVector[k] ) < K_DOUBLE_TOL ) ? true : false;
}


size_t IdxFromValue( const std::vector<double>& dateVector, double date )
{
	size_t k=0;
	while( k<dateVector.size() && dateVector[k] < date )
		++k;

	if( k >= dateVector.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_ProcessBuilderPDE::IdxFromDate: date not found in the map!");

	return k;
}


size_t IdxFromValueWithTol( const std::vector<double>& dateVector, double date, double daysNb , bool NoThrow)
{
	int i=-1, N = dateVector.size();

	while( (abs( date - dateVector[++i] ) > daysNb ) && i<N ) {}

	if( i == N && NoThrow == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_ProcessBuilderPDE::IdxFromDate: date not found in the map!");

	return i == N ? -1 : i;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

