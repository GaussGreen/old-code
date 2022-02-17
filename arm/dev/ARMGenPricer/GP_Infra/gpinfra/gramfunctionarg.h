/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctionarg.h,v $
 * Revision 1.1  2003/10/08 16:46:45  ebenhamou
 * Initial revision
 *
 *
 */


/*! \file gramfunctionarg.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_GRAMFUNCTIONARG_H
#define _INGPINFRA_GRAMFUNCTIONARG_H



#include "gpbase/port.h"
#include "gramargs.h"
#include "gpbase/rootobject.h"
#include "gpbase/assignop.h"

#include <string>
CC_USING_NS( std, string )


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
// Class	   : ARM_GramFunctionArg 
// Description : Declaration of an argument to a tableau function.
/////////////////////////////////////////////////////////////////

class ARM_GramFunctionArg : public ARM_RootObject
{
public:
	/// default constructor
	ARM_GramFunctionArg();	

	/// real constructor
	ARM_GramFunctionArg( GramFuncArgType Type, bool Required, bool IsVaArg, const string& Name,
		const string& Description, bool IsVisible = true, const string& DefaultValue = "", bool AddRefNode = false );
	
	/// copy constructor
	ARM_GramFunctionArg( const ARM_GramFunctionArg& rhs );

	/// destructor
	virtual ~ARM_GramFunctionArg();

	/// assignment operator
	ASSIGN_OPERATOR(ARM_GramFunctionArg)

	/// accessors all inline
	inline bool Required() const { return itsRequired; }
	inline bool IsVaArg() const { return itsIsVaArg; }
	inline GramFuncArgType  Type() const { return itsType; }
	inline string DefaultValue() const { return itsDefaultValue; }
	inline bool IsVisible() const { return itsIsVisible; }
	inline bool AddRefNode() const { return itsAddRefNode; }

	/// comparison function
	bool equal( ARM_GramFunctionArg const& rhs ) const	
	{ 
		return ( ( itsName  ==  rhs.itsName ) && ( itsType == rhs.itsType ) );
	}

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;


private:
	GramFuncArgType  itsType;

	bool itsRequired;

	bool itsIsVaArg;

	string itsName;

	string itsDescription;

    bool itsIsVisible;

	string itsDefaultValue;

	bool itsAddRefNode;
};


inline bool operator==( ARM_GramFunctionArg const& lhs, ARM_GramFunctionArg const& rhs )
{
	return lhs.equal( rhs );
}

inline bool operator!=( ARM_GramFunctionArg const& lhs, ARM_GramFunctionArg const& rhs )
{
	return !lhs.equal( rhs );
}

/// for output purposes!
CC_NS(std,ostream)& operator<<( CC_NS(std,ostream)& s, const ARM_GramFunctionArg & );


CC_END_NAMESPACE()

#endif 



