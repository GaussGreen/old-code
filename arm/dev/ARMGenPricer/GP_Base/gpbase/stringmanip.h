/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: stringmanip.h,v $
 * Revision 1.5  2003/11/14 15:05:43  ebenhamou
 *  added functor for less case insensitive
 *
 * Revision 1.4  2003/10/16 16:12:00  ebenhamou
 * change name because old name was really confusing
 *
 * Revision 1.3  2003/10/09 19:04:34  ebenhamou
 * added StringFctor
 *
 * Revision 1.2  2003/10/03 09:34:52  ebenhamou
 * added function StrUpper
 *
 * Revision 1.1  2003/09/17 18:06:22  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file stringmanip.h
 *
 *  \brief Various string manipulations
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _STRINGMANIP_H
#define _STRINGMANIP_H

#include "port.h"
#include <string>
#include <functional>
#include <vector>

CC_USING_NS( std, string )
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

/// conversion to upper case
/// directly on the string
void stringToUpper( string& s );

/// similar but return the string
string& StrUpper( string& s );

/// returns the uppercase version of the string
/// on a const string
string stringGetUpper( const string& s );

/// function to returns the version of a string 
/// without blank
void stringTrim( string& str );

/// checks that something is a real upper case string
void StringUpperCaseValidation( const string& s, const string& stringName = "" );

/// split a given string into a vector of string according to a given delimiter
vector<string> splitString( const string& s, const string& delimiter );

/// compare case insensitive two strings!
bool StrIEqual( string& s1, string& s2 );

/// less<string> case insensitive oeprator
bool StrILessOp( string& s1, string& s2 );


//// check case insensitive
//// StrEqualFunctor predicate
struct StrEqualFunctor : std::unary_function< string, bool >
{
	string s1;
	StrEqualFunctor( const string& s ) : s1( s ){};
	
	bool operator()( string& s2 ){
		return StrIEqual( s1, s2 ); 
	}
};

struct StrILess : std::binary_function< string, string, bool >
{
	bool operator()( string& s1, string& s2 ) const {
		return StrILessOp( s1, s2 );
	}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

