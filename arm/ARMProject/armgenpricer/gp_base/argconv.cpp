/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconv.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/argconv.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "expt.h"

CC_USING_NS( std, pair )
CC_USING_NS( ARM, stringGetUpper)
CC_USING_NS( ARM, StrUpper )

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ArgConv
///	Routine: Constructor
///	Returns: built the object
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////
void ARM_ArgConv::InsertOneLine(const char* name, int value, size_t i )
{
	bool InsertedSucces = itsConvention.insert( pair< const string, int >( StrUpper( string(name) ), value ) ).second;
	if( !InsertedSucces  )
	{
#if defined( _DEBUG )
		/// to catch this bug at run time, we use a special assembler instruction to break on error
#if defined( WIN32 )
		__asm int 3
#endif
#endif
			CC_Ostringstream os;
		os << "Problem with inserting arg " << i << " in argconv table "<< name;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ArgConv
///	Routine: Constructor
///	Returns: built the object
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////

ARM_ArgConv::ARM_ArgConv( const ARGConvTable* tableInit, const string& name )
: itsName( name )
{
	/// MAXSIZE to avoid exploding table
	/// table can have maximum MAXISIZE rows!
	unsigned long MAXSIZE = 100;
	size_t i=0;
	while( i<MAXSIZE && strcmp( tableInit[i].itsName, ENDOFLINE_CHAR ) != 0 )
	{
		InsertOneLine( tableInit[i].itsName, tableInit[i].itsValue, i );
		i++;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ArgConv
///	Routine: Constructor with ARGConvRevTable
///	Returns: built the object
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////

ARM_ArgConv::ARM_ArgConv( const ARGConvRevTable* tableInit, const string& name )
: itsName( name )
{
	/// MAXSIZE to avoid exploding table
	/// table can have maximum MAXISIZE rows!
	unsigned long MAXSIZE = 100;
	size_t i=0;
	while( i<MAXSIZE && strcmp( tableInit[i].itsName, ENDOFLINE_CHAR ) != 0 )
	{
		InsertOneLine( tableInit[i].itsName, tableInit[i].itsValue, i );
		i++;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ArgConv
///	Routine: GetNumber
///	Returns: corresponding number for a given string
///	Action : similar to the convention conversion function
///		in interglob.cpp
////////////////////////////////////////////////////

int ARM_ArgConv::GetNumber( const string& s ) const
{
	CC_NS( std, map)< string, int, less< string > >::const_iterator 
		iter, end;
	
	iter= itsConvention.find( stringGetUpper( s ) );
	end	= itsConvention.end();
	
	if( iter == end )
	{
		CC_Ostringstream os;
		os << "Invalid argument for " << itsName << " Valid are (case insensitive) ";
		
		CC_NS( std, map)< string,int, less<string> >::const_iterator names = itsConvention.begin();
		for( ; names != end; ++names )
			os <<  (*names).first<< ", ";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	
	return (*iter).second;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_ArgConvReverse
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_ArgConvReverse
///	Routine: Constructor
///	Returns: built the object
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////
void ARM_ArgConvReverse::InsertOneLine( int value, const char* name, size_t i )
{
	/// check that we can safely insert!
	bool InsertedSucces = itsReverseConvention.insert( pair< const int, string>( value, StrUpper( string(name) ) ) ).second;
	if( !InsertedSucces  )
	{
		/// to catch this bug at run time, we use a special assembler instruction to break on error
#if defined( _DEBUG )
		/// to catch this bug at run time, we use a special assembler instruction to break on error
#if defined( WIN32 )
		__asm int 3
#endif
#endif
			CC_Ostringstream os;
		os << "Problem with inserting arg " << i << " in argconv reverse table "<< name;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ArgConvReverse
///	Routine: Constructor
///	Returns: void
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////

ARM_ArgConvReverse::ARM_ArgConvReverse( const ARGConvTable* tableInit, const string& name )
: itsName( name )
{
	/// MAXSIZE to avoid exploding table
	/// table can have maximum MAXISIZE rows!
	unsigned long MAXSIZE = 100;
	size_t i=0;
	while( i<MAXSIZE && strcmp( tableInit[i].itsName, ENDOFLINE_CHAR ) != 0 )
	{
		InsertOneLine( tableInit[i].itsValue, tableInit[i].itsName, i );
		++i;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_ArgConvReverse
///	Routine: Constructor with ARGConvRevTable
///	Returns: void
///	Action : Set the name and usee the tableInit to 
///				fill the map
////////////////////////////////////////////////////

ARM_ArgConvReverse::ARM_ArgConvReverse( const ARGConvRevTable* tableInit, const string& name )
:	itsName( name )
{
	/// MAXSIZE to avoid exploding table
	/// table can have maximum MAXISIZE rows!
	unsigned long MAXSIZE = 100;
	size_t i=0;
	while( i<MAXSIZE && strcmp( tableInit[i].itsName, ENDOFLINE_CHAR ) != 0 )
	{
		InsertOneLine( tableInit[i].itsValue, tableInit[i].itsName, i );
		++i;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_ArgConvReverse
///	Routine: GetString
///	Returns: corresponding string for a given number
///	Action : similar to the convention conversion function
///		in interglob.cpp
////////////////////////////////////////////////////

string ARM_ArgConvReverse::GetString( int d ) const
{
	CC_NS( std, map)< int, string, less< int > >::const_iterator 
		iter, end;
	
	iter= itsReverseConvention.find( d );
	end	= itsReverseConvention.end();
	
	if( iter == end )
	{
		CC_Ostringstream os;
		os << "Invalid argument for" << d << " Valid are (case insensitive) ";
		
		CC_NS( std, map)< int, string, less<int> >::const_iterator names = itsReverseConvention.begin();

		/// flag for blank
		os.unsetf(std::ios::skipws);
		for( ; names != end; ++names )
			os <<  (*names).first << ' ';
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	
	return (*iter).second;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

