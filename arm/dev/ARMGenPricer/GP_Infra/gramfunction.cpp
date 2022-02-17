/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunction.cpp,v $
 * Revision 1.1  2003/10/08 16:41:04  ebenhamou
 * Initial revision
 *
 *
 *
 */


/*! \file gramfunction.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramfunction.h"
#include "gpinfra/gramfunctiontable.h"
#include "gpbase/ostringstream.h"
#include "gpbase/env.h"
#include "gpbase/gpmatrix.h"

#ifndef unix
	#include <iomanip>
#endif

#include <algorithm>
CC_USING_NS( std, pair )

#include "gpbase/stringmanip.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction 
///	Routine: IsAFunctionName
///	Returns: bool
///	Action : tells wheter a string is a function name
/////////////////////////////////////////////////////////////////
bool ARM_GramFunction::IsAFunctionName( const string& funcName )
{
	/// case insensitive search!
	MapStringToVoid* FTable			= ARM_GramFunction::TheFuncTable;
	MapStringToVoid::iterator iter	= FTable->find( stringGetUpper( funcName ) );
	return iter != FTable->end();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: Declaration
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::Declaration::Declaration( 
	const string& Name, 
	GramFuncArgType ReturnType,
	const string& Description,
	const ArgVector& Args,
	const string& Category,
	ARM_GramFctor* Fctor,
	GramNodeBuilderFunc gramNodeBuilderFunc )
:
	itsFuncName( Name ), 
	itsReturnType( ReturnType ), 
	itsDescription( Description ), 
	itsArgs( Args ),
	itsFctor( Fctor ),
	itsCategory( Category ),
	itsGramNodeBuilderFunc(gramNodeBuilderFunc )
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: Declaration
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::Declaration::Declaration() 
:	itsReturnType( GFAT_UNKNOWN_TYPE ), itsFctor( NULL ) 
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: Declaration
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::Declaration::Declaration(
	const ARM_GramFunction::Declaration &FDecl)
:
	itsFuncName( FDecl.itsFuncName ), 
	itsReturnType( FDecl.itsReturnType ), 
	itsDescription( FDecl.itsDescription ), 
	itsArgs( FDecl.itsArgs ),
	itsFctor( FDecl.itsFctor ),
	itsCategory( FDecl.itsCategory ),
	itsGramNodeBuilderFunc( FDecl.itsGramNodeBuilderFunc )
{}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: ~Declaration
///	Returns: nothing
///	Action : Destructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::Declaration::~Declaration()
{

}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: operator= 
///	Returns: *this with new values
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_GramFunction::Declaration& ARM_GramFunction::Declaration::operator=(
	const ARM_GramFunction::Declaration& rhs )
{
	if( this != &rhs )
	{
		itsFuncName 			= rhs.itsFuncName; 
		itsReturnType 			= rhs.itsReturnType; 
		itsDescription 			= rhs.itsDescription; 
		itsArgs 				= rhs.itsArgs;
		itsFctor				= rhs.itsFctor;
		itsCategory				= rhs.itsCategory;
		itsGramNodeBuilderFunc	= rhs.itsGramNodeBuilderFunc;
	}
	return *this;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction::Declaration 
///	Routine: toString 
///	Returns: The function declaration as a string.
///	Action : 
/////////////////////////////////////////////////////////////////

string ARM_GramFunction::Declaration::toString(const string& indent, const string& nextIndent) const
{
	const int FUNCNAMESIZE = 14;
	const string SEPARATOR("  ");
	CC_Ostringstream os;
	os << "ARM_GramFunction \n";
	os << "*********** Function    : ";

#ifndef unix
	/// get old flags
	CC_NS( std, ios )::fmtflags oldFlags = os.flags();
	os << CC_NS( std, setiosflags )( CC_NS( std, ios)::left );
	os << CC_NS( std, setw)(FUNCNAMESIZE);
	os << itsFuncName;
	/// restore default
	os.flags( oldFlags  );
#else
	os << itsFuncName;
#endif
	
	os << " ****************\n";
	os << "Category	    \t: " + itsCategory << "\n";
	os << "Description	\t: " + itsDescription << "\n";
	os << "Syntax		\t:\n";
	os << SEPARATOR << itsFuncName << "(\n";

	ArgVector::const_iterator
			Arg    = itsArgs.begin(),
			ArgEnd = itsArgs.end();

	for( ;Arg != ArgEnd; ++Arg )
    {
        if(Arg->IsVisible())
		    os <<  Arg->toString(SEPARATOR + SEPARATOR) << "\n";
    }
	os << SEPARATOR << ")\n\n";
	return os.str();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction 
///	Routine: ARM_GramFunction
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::ARM_GramFunction( const ARM_GramFunction& Other)
:	itsDeclaration( Other.itsDeclaration ), itsFctor(Other.itsFctor)
{}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction 
///	Routine: ~ARM_GramFunction
///	Returns: nothing
///	Action : Destructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::~ARM_GramFunction()
{
	/// we do not delete itsDeclaration
	/// as it is given by a static table that
	/// should be deleted afterwards!
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction 
///	Routine: operator= 
///	Returns: *this with new value
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_GramFunction& ARM_GramFunction::operator=( const ARM_GramFunction &rhs )
{
	if( this != &rhs )
	{
		itsDeclaration	= rhs.itsDeclaration;
		itsFctor		= rhs.itsFctor;
	}
	return *this;
}


/////////////////////////////////////////////////////////////////
///	Routine: Clone
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_GramFunction::Clone() const
{
	return new ARM_GramFunction( * this );
}







/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction 
///	Routine: ARM_GramFunction
///	Returns: nothing
///	Action : Constructor
/////////////////////////////////////////////////////////////////

ARM_GramFunction::ARM_GramFunction(	const string& FuncName, const string& payModelName )
: itsFctor( ARM_GramFctorPtr( NULL ) )
{
	/// case insensitive search!
	MapStringToVoid* FTable			= ARM_GramFunction::TheFuncTable;
	MapStringToVoid::iterator iter	= FTable->find( stringGetUpper( FuncName ) );

	/// have we found it?
	if( iter == FTable->end() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			string( "ARM_GramFunction: cannot construct ARM_GramFunction: " ) + FuncName + "\n" );

	/// take it
	itsDeclaration	= reinterpret_cast< Declaration *>( (*iter).second );

	/// clone it to avoid problem!
	itsFctor		= itsDeclaration->Fctor()->Clone();

    /// Set payment model use for discounting in multi-currencies pricing
    itsFctor->SetPayModelName(payModelName);
}



	
/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction
///	Routine: ValidateArg
///	Returns: void
///	Action : throws an exception if a type is NULL
/////////////////////////////////////////////////////////////////

void ARM_GramFunction::ValidateArgStr( const string& t, int argNb, const string& typeName, 
	const string& FuncName )
{
	#ifdef __GP_STRICT_VALIDATION
		if( t == "" ) 
		{
			CC_Ostringstream os;
			os << "ARM_GramFunction:\nScrewed up argument number n° " << argNb << " for type " 
				<< typeName << "in function " << FuncName << "\n"; 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	#else
		/// nothing
	#endif
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramFunction
///	Routine: Reset
///	Returns: void
///	Action : basically reset the pre-computed arguments by taking the
///				original gramfunction.
/////////////////////////////////////////////////////////////////

void ARM_GramFunction::Reset( bool totalReset )
{
    /// Save payment model
    const string payModelName = itsFctor->GetPayModelName();

	/// clone it to avoid problem with change of the object!
	if ( totalReset ) 
	{
		itsFctor		= itsDeclaration->Fctor()->Clone();
	}
	else
	{
		itsFctor->SetAlreadyComputed(0);
	}

    /// Set payment model use for discounting in multi-currencies pricing
    itsFctor->SetPayModelName(payModelName);
}


/////////////////////////////////////////////////////////////////
///	Routine: CreateTheFuncTable
///	Returns: MapStringToVoid
///	Action : creates the static map table for all language functions
/////////////////////////////////////////////////////////////////

ARM_GramFunction::MapStringToVoid* ARM_GramFunction::CreateTheFuncTable()
{
	MapStringToVoid* theFuncTable = new MapStringToVoid;

	/// for easier typing
	typedef ARM_GramFunctionArgStruct GFAStruct; 
	typedef ARM_GramFunctionArgVector GFAVector;

	const int ARGSIZEMAX = 40;
	size_t i;

	for(i=0; i<FuncTableSize; ++i )
	{
		/// get variables
		string FuncName				= FunctionsTable[i].itsName;
		GramFuncArgType ReturnType	= FunctionsTable[i].itsReturnType;
		string Description			= FunctionsTable[i].itsDescription;
		const GFAStruct* ga			= FunctionsTable[i].itsArgs;
		ARM_GramFctor* Fctor		= FunctionsTable[i].itsFctor;
		string Category				= FunctionsTable[i].itsCategory;
		GramNodeBuilderFunc 
			gramNodeBuilderFunc		= FunctionsTable[i].itsGramNodeBuilder;

#ifdef __GP_STRICT_VALIDATION
		if( !gramNodeBuilderFunc )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"any function should have its gram node builder function\n" );
#endif

		int ArgSize = 0;
		while( ga[ArgSize].itsType != GFAT_TERMINATOR && ArgSize < ARGSIZEMAX )
			++ArgSize;

#ifdef __GP_STRICT_VALIDATION
		///  check for incoherence
		if( ArgSize == ARGSIZEMAX )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"ARM_GramFunction: you have probably a screwed table definition... please check it\n" );
#endif

		GFAVector Args(ArgSize);
		for( int j=0; j<ArgSize; ++j )
		{

#ifdef __GP_STRICT_VALIDATION
			/// validation checking
			ARM_GramFunction::ValidateArgStr( ga[j].itsName,			j, "Name",			FuncName );
			ARM_GramFunction::ValidateArgStr( ga[j].itsDescription,	j, "Description",	FuncName );
#endif				

			GramFuncArgType Type = ga[j].itsType;
			bool Required = ga[j].itsIsRequired;
			bool IsVaArg = ga[j].itsIsVaArg;
			string Name = ga[j].itsName;
			string Description = ga[j].itsDescription;
			string DefaultValue;
			
			/// if it is not required
			/// we should have a default!
			if( !Required )
			{
#ifdef __GP_STRICT_VALIDATION
				if( ga[j].itsDefaultValue == NULL )
				{
#if defined( WIN32 )
					/// if you end up here, it is because you have not set up 
					/// the default value for your argument
					__asm int 3
#else
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"no default but the table says that it is not a required argument, please advise!" );
#endif
				}
#endif
				/// currently default value can be only double or string
				/// case of double
				if( Type == GFAT_DOUBLE_TYPE || Type == GFAT_DATE_OR_DOUBLE_TYPE ||
                    Type == GFAT_VECTOR_TYPE || Type == GFAT_CURVE_TYPE || Type == GFAT_VECTOR_OR_CURVE_TYPE)
				{
					CC_Ostringstream os;
					os << ga[j].itsDefaultValue->itsDouble;
					DefaultValue = os.str();
				}
				else if( Type == GFAT_STRING_TYPE || Type == GFAT_MULTITOKENSTRING_TYPE ||
                    Type == GFAT_MATURITY_TYPE || Type == GFAT_STRING_OR_CURVE_TYPE )
					DefaultValue = ga[j].itsDefaultValue->itsChar? ga[j].itsDefaultValue->itsChar :"";
                else if( Type == GFAT_DATESTRIP_TYPE )
                    /// No default value, set 0
                    DefaultValue = string("0");
				else
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"defaulted argument can only be string or double" );
			}

			Args[j]= ARM_GramFunctionArg( ga[j].itsType, ga[j].itsIsRequired, ga[j].itsIsVaArg,  ga[j].itsName,
				ga[j].itsDescription, ga[j].itsIsVisible, DefaultValue, ga[j].itsAddRefNode);
		}

		theFuncTable->insert( pair< const string, void* >(
			(const string) StrUpper( FuncName ),
			(void*) new Declaration( FuncName, ReturnType, Description, Args, Category, Fctor, gramNodeBuilderFunc) ) );
	}

	return theFuncTable;
}


/////////////////////////////////////////////////////////////////
///	Routine: ReleaseTheFuncTable
///	Returns: MapStringToVoid
///	Action : deletes hte corresponding declaration.... functors have
///				not to be deleted as they are right now static variables
/////////////////////////////////////////////////////////////////

void ARM_GramFunction::ReleaseTheFuncTable()
{
	if (ARM_GramFunction::TheFuncTable)
	{
		MapStringToVoid::iterator iter = ARM_GramFunction::TheFuncTable->begin();
		
		while( iter != ARM_GramFunction::TheFuncTable->end() )
		{
			ARM_GramFunction::Declaration* theDeclaration = reinterpret_cast<ARM_GramFunction::Declaration*>( (*iter).second );
			delete theDeclaration;
			++iter;
		}

		delete ARM_GramFunction::TheFuncTable;
		ARM_GramFunction::TheFuncTable = NULL;
	}
}

/////////////////////////////////////////////////////////////////
///	Routine: CreateTheFuncTable
///	Returns: MapStringToVoid
///	Action : Create The FuncTable ... this is a static variable
///				and is created only once!
/////////////////////////////////////////////////////////////////

void ARM_GramFunction::InitTheFuncTable()
{
	ARM_GramFunction::TheFuncTable = ARM_GramFunction::CreateTheFuncTable();
}

/////////////////////////////////////////////////////////////////
/// declaration of the static variable
/////////////////////////////////////////////////////////////////

ARM_GramFunction::MapStringToVoid* ARM_GramFunction::TheFuncTable = NULL;



CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


