/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctionarg.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpinfra/gramfunctionarg.h"
#include "gpbase/ostringstream.h"

/// for easy printing
#include <ostream>
CC_USING( std::ostream )

#include <iomanip>

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunctionArg 
///	Routine: ARM_GramFunctionArg
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunctionArg::ARM_GramFunctionArg()
:	itsType( GFAT_UNKNOWN_TYPE)
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunctionArg 
///	Routine: ARM_GramFunctionArg
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunctionArg::ARM_GramFunctionArg(
	GramFuncArgType Type, 
	bool Required, 
	bool IsVaArg,
	const string& Name,
	const string& Description,
    bool IsVisible,
	const string& DefaultValue,
	bool AddRefNode
) :
	itsType( Type ), 
	itsRequired( Required ), 
	itsIsVaArg(IsVaArg),
	itsName( Name ),
	itsDescription( Description ),
    itsIsVisible( IsVisible ),
	itsDefaultValue( DefaultValue ),
	itsAddRefNode( AddRefNode )
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunctionArg 
///	Routine: ARM_GramFunctionArg
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunctionArg::ARM_GramFunctionArg( const ARM_GramFunctionArg &Arg )
:
	itsType( Arg.itsType ), 
	itsRequired( Arg.itsRequired ),
	itsIsVaArg( Arg.itsIsVaArg ),
	itsName( Arg.itsName ),
	itsDescription( Arg.itsDescription ),
    itsIsVisible( Arg.itsIsVisible ),
	itsDefaultValue( Arg.itsDefaultValue ),
	itsAddRefNode( Arg.itsAddRefNode ) 
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunctionArg 
///	Routine: ~ARM_GramFunctionArg
///	Returns: nothing
///	Action : Destructor
/////////////////////////////////////////////////////////////////

ARM_GramFunctionArg::~ARM_GramFunctionArg()
{
	/// nothing
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunctionArg 
///	Routine: toString() 
///	Returns: The function argument as a string.
///	Action : 
/////////////////////////////////////////////////////////////////

string ARM_GramFunctionArg::toString(const string& indent, const string& nextIndent) const
{
    /// If internal argument, dump nothing !
    if(!itsIsVisible)
        return "";

	CC_Ostringstream os;
	const int ARGSIZE = 20;
	os << indent;

	/// get old flags
	CC_NS( std, ios )::fmtflags oldFlags = os.flags();
	os << CC_NS( std, setiosflags )( CC_NS( std, ios)::left );

	size_t argTextLength = itsName.size() + 2;
	if( argTextLength < ARGSIZE+2)
		argTextLength = ARGSIZE+2;

	if( !itsRequired )
	{
		os  << CC_NS( std, setw )(argTextLength)  
			<< string("[") + itsName + string("]");
	}
	else
	{
		os  << CC_NS( std, setw )(argTextLength)  
			<< string(" ") + itsName + string(" ");
	}
	
	os << " /// ";
	os << itsDescription;

	if( !itsRequired )
		os << " (optional DEFAULT = " << itsDefaultValue << " )";

	if (itsIsVaArg)
		os << " ...";

	/// restore default
	os.flags( oldFlags  );

	return os.str();
}


/////////////////////////////////////////////////////////////////
///	Routine: operator<<
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ostream& operator<<( ostream& s,  const ARM_GramFunctionArg &Arg )
{
	s << Arg.toString("");
	return s;
}


/////////////////////////////////////////////////////////////////
///	Routine: Clone
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_GramFunctionArg::Clone() const
{
	return new ARM_GramFunctionArg( * this );
}



CC_END_NAMESPACE()
