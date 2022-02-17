/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctorargdict.h,v $
 * Revision 1.1  2004/15/02 18:53:24  ebenhamou
 * initial version
 *
 *
 */

/*! \file gramfunctorargdict.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */


#ifndef _INGPINFRA_GRAMFUNCTORARGDICT_H
#define _INGPINFRA_GRAMFUNCTORARGDICT_H

#include "gpbase/port.h"
#include "typedef.h"
#include "gramfunctorarg.h"
#include "gpbase/rootobject.h"
#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )
#include <vector>
CC_USING_NS( std, vector )

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// \class ARM_GramFctorArgDict
/// \brief class to implement the concept of a generic dictionary
///	 that support multi-type... the idea underlying is to map
///  a string to something of the type double, date, string or vector of double
///	 leveraging on the STL map concept with string and ARM_GramFctorArgDict 
///  the design decision was to provide a minimalistic interface and to provide
///	 iterator access for user that want more advanced functionalities.
///
///	 beware that this is case insensitive, hence the key is in upper case
///	 while the original string is kept as the second element of the pair  
/////////////////////////////////////////////////////////////////

class ARM_GramFctorArgDict : public ARM_RootObject
{
public:
	ARM_GramFctorArgDict( const string& key = "", const ARM_GramFctorArg& arg = ARM_GramFctorArg() );
	bool InsertData( const string& key, const ARM_GramFctorArg& arg );

	/// accessors!
	vector<string> GetKeys() const;
	ARM_GramFctorArg GetData( const string& key ) const;
    bool IsDataExist(const string& key) const;
	ARM_GramFctorArg& operator[]( const string& key );				/// if the element does not exist, creates one!
	const ARM_GramFctorArg& operator[] ( const string& key ) const; /// if the element does not exist, creates one!

	/// copy constructor, assignment and destructor
	ARM_GramFctorArgDict( const ARM_GramFctorArgDict & rhs );
	ARM_GramFctorArgDict & operator=( const ARM_GramFctorArgDict & rhs );
	virtual ~ARM_GramFctorArgDict ();

	/// standard ARM Object support
	/// because of the clone implementation in terms of the copy constructor
	/// there is no need to define BitwiseCopy and Copy!
	virtual ARM_Object* Clone() const;
	/// for debugging and exception and so on
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// iterators are the ones of the STL
	/// this class only wraps up the STL iterator
	/// we use a pair of string and ARM_GramFctorArg to implement
	/// case insensitive key!
	typedef map< string, CC_NS( std, pair)< ARM_GramFctorArg, string >, less< string > > strGramFctorArgMap;
	typedef CC_NS( std, pair)< ARM_GramFctorArg, string > gramFctorArgStringPair;
	typedef strGramFctorArgMap::const_iterator const_iterator;
	typedef strGramFctorArgMap::iterator iterator;

	/// iterator support const version
	const_iterator begin() const { return itsContents.begin(); }
	const_iterator end() const { return itsContents.end(); }

	/// non const version
	iterator begin() { return itsContents.begin(); }
	iterator end() { return itsContents.end(); }

private:
	strGramFctorArgMap itsContents;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

