/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorargdict.cpp
 *
 *  \brief gramfunctorarg are arguments for the functor 
 *	corresponding to grammar functions
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/gramfunctorargdict.h"

/// ARM Kernel
#include <glob/expt.h>

// STL
#include <iomanip>
#include <algorithm>

CC_USING_NS(std,pair)
CC_USING_NS(std,find)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: Constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_GramFctorArgDict::ARM_GramFctorArgDict( const string& key, const ARM_GramFctorArg& arg )
:	itsContents()
{
	SetName(ARM_GRAMFCTORARGDICT);
	if( key != "" )
		itsContents.insert( pair< const string, gramFctorArgStringPair >( stringGetUpper( key), gramFctorArgStringPair( arg, key ) ) );
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: Copy Constructor
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_GramFctorArgDict::ARM_GramFctorArgDict( const ARM_GramFctorArgDict & rhs )
:	ARM_RootObject( rhs ), itsContents( rhs.itsContents )
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GramFctorArgDict& ARM_GramFctorArgDict::operator=( const ARM_GramFctorArgDict & rhs )
{
	if( this != &rhs )
	{
		ARM_RootObject::operator =(rhs);
		itsContents = rhs.itsContents;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: Destructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////
ARM_GramFctorArgDict::~ARM_GramFctorArgDict ()
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: Clone, View, toString
///	Returns: built object
///	Action : 
////////////////////////////////////////////////////
ARM_Object* ARM_GramFctorArgDict::Clone() const
{
	return new ARM_GramFctorArgDict( *this );
}


string ARM_GramFctorArgDict::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "ARM_GramFunctor Arguments Dictionary:\n\n";

	if( !itsContents.empty() )
	{
		os << CC_NS(std,setw)(35) << CC_NS( std, left ) << " keys :" << "\t" << "Value type and type\n";
		strGramFctorArgMap::const_iterator  iter = itsContents.begin(),
			end = itsContents.end();

		while( iter != end )
		{
			os << CC_NS(std,setw)(35) << CC_NS( std, left ) << (*iter).second.second << "\t" << (*iter).second.first.toString() << "\n";
			++iter;
		}

		os << CC_NS(std,endl);
	}

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: InsertData
///	Returns: a reference to the element with key key
///	Action : if the element does not exist, creates one!
////////////////////////////////////////////////////
ARM_GramFctorArg& ARM_GramFctorArgDict::operator[]( const string& key )
{
	string upperCaseKey	= stringGetUpper( key);
	itsContents[upperCaseKey]= gramFctorArgStringPair( ARM_GramFctorArg(), key );
	return itsContents[upperCaseKey].first;
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: InsertData
///	Returns: Says whether an object has been inserted or not
///	Action : 
////////////////////////////////////////////////////
bool ARM_GramFctorArgDict::InsertData( const string& key, const ARM_GramFctorArg& arg )
{
	return itsContents.insert( pair< const string, gramFctorArgStringPair  >( key, gramFctorArgStringPair( arg, key ) ) ).second;
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: GetKeys
///	Returns: Returns all the keys
///	Action : the original ones and not the upper case ones!
////////////////////////////////////////////////////
vector<string> ARM_GramFctorArgDict::GetKeys() const
{
	vector<string> result;
	result.reserve(itsContents.size());

	strGramFctorArgMap::const_iterator iter, end = itsContents.end();
	for( iter = itsContents.begin(); iter != end; ++iter )
		result.push_back( (*iter).second.second );
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: IsDataExist
///	Returns: bool
///	Action : Check if the key exist in the map
////////////////////////////////////////////////////
bool ARM_GramFctorArgDict::IsDataExist(const string& key) const
{
    strGramFctorArgMap::const_iterator iter = itsContents.find(key);

    return (iter!=itsContents.end());
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArgDict
///	Routine: GetData
///	Returns: Throws an exception if the key does not exist!
///	Action : case insensitive....
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GramFctorArgDict::GetData( const string& key ) const
{
	strGramFctorArgMap::const_iterator  iter,
		end = itsContents.end();
	 
	iter = itsContents.find( stringGetUpper( key ) );

	if( iter == end )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + ": could not find key " + key );
	else
		return (*iter).second.first;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

