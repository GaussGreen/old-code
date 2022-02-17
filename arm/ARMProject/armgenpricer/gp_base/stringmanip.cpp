/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file stringmanip.cpp
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/stringmanip.h"

/// ARM Kernel
#include "expt.h"			/// neceassary for exception throwing

#include <string>
#include <algorithm>
#include <vector>

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// function to convert a string to uppercase
/////////////////////////////////////////////////////////////////
void stringToUpper( string& s )
{
	/// be aware that we are forced to use transform
	/// because of function ambiguity!
	std::transform( s.begin(), s.end(), s.begin(), toupper );
}


/////////////////////////////////////////////////////////////////
///	Routine: routine that converts a string to lower case
/////////////////////////////////////////////////////////////////
string& StrUpper( string& s )
{
	std::transform( s.begin(), s.end(), s.begin(), toupper );
	return s;
}



/////////////////////////////////////////////////////////////////
/// function to returns the upper case version of a string
/////////////////////////////////////////////////////////////////
string stringGetUpper( const string& s )
{
	string a = s;
	stringToUpper( a );
	return a;
}


/////////////////////////////////////////////////////////////////
/// function to returns the version of a string without blank
/////////////////////////////////////////////////////////////////
void stringTrim(string& str )
{
	string str2;
	for (int i = 0; i < str.size(); ++i)
		if ((str[i] != ' ') && (str[i] != '\t'))
			str2 += str[i];

	str = str2;
}

/////////////////////////////////////////////////////////////////
//// function to check that a string is really in uppercase
/////////////////////////////////////////////////////////////////
void StringUpperCaseValidation( const string& s, const string& stringName )
{
	if( stringGetUpper( s ) != s ) 
	{
		char msg[255];
		sprintf( msg, "the %s tag %s, is not in upper case! Please advise", stringName.c_str(), s.c_str() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
}


/////////////////////////////////////////////////////////////////
/// split a given string into a vector of string according to a given delimiter
/////////////////////////////////////////////////////////////////
vector< string> splitString( const string& s, const string& delimiter )
{
	vector<string> subTags;
	string::size_type pos = 0, prev_pos = 0;
	while( (pos = s.find_first_of(delimiter, pos ) ) != string::npos )
	{
		subTags.push_back( s.substr( prev_pos,pos-prev_pos ) );
		prev_pos = ++pos;
	}

	/// last one not to be forgotten
	subTags.push_back( s.substr( prev_pos,pos-prev_pos ) );
	return subTags;
}


/////////////////////////////////////////////////////////////////
/// string equal case insensitive
/////////////////////////////////////////////////////////////////
bool StrIEqual( string& s1, string& s2 )
{
	return StrUpper( s1 ) == StrUpper( s2 ); 
}


/////////////////////////////////////////////////////////////////
/// string less case insensitive
/////////////////////////////////////////////////////////////////
bool StrILessOp( string& s1, string& s2 )
{
	return StrUpper( s1 ) < StrUpper( s2 ); 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

