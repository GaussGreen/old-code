/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconv.h
 *
 *  \brief files to handle string argument ... similar
 *			to the interface interglob.cpp file
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPBASE_ARGCONV_H
#define _INGPBASE_ARGCONV_H

/// this header comes firts as it includes some preprocessor constants!
#include "removeidentifiedwarning.h"

/// use our macro for namespace
#include "port.h"

#include <string>
CC_USING_NS( std, string )

#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )

CC_BEGIN_NAMESPACE( ARM )

#define ENDOFLINE_CHAR	"ENDOFLINE"
#define ENDOFLINE_INT	-1111111111


//// struct for simple initialisation
struct ARGConvTable
{
	const char*	itsName;
	int			itsValue;
};

struct ARGConvRevTable
{
	int			itsValue;
	const char*	itsName;
};

/// class for conversion of string to the corresponding enum!
class ARM_ArgConv
{
private:
	map< string, int, less<string> > itsConvention;
	string itsName;
	void InsertOneLine(const char* name, int value, size_t i );
public:
	ARM_ArgConv( const ARGConvTable* tableInit, const string& name );
	ARM_ArgConv( const ARGConvRevTable* tableInit, const string& name );
	int GetNumber( const string& s ) const;
};


/// class to reverse arg convention
/// warning a table can be reversed if and only if
/// there is a unique one to one mapping between string and
///	int!
class ARM_ArgConvReverse
{
private:
	map< int, string, less<int> > itsReverseConvention;
	string itsName;
	void InsertOneLine(int value, const char* name, size_t i);
public:
	ARM_ArgConvReverse( const ARGConvTable* tableInit, const string& name );
	ARM_ArgConvReverse( const ARGConvRevTable* tableInit, const string& name );
	string GetString( int d ) const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
