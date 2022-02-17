/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctiontable.h,v $
 * Revision 1.1  2003/10/08 16:48:26  ebenhamou
 * Initial revision
 *
 *
 */



/*! \file gramfunctiontable.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GRAMFUNCTIONTABLE_H
#define _INGPINFRA_GRAMFUNCTIONTABLE_H

#include "gpbase/port.h"
#include "gramargs.h"
#include "gramfunctor.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ExpNodeFunc;
class ARM_GramFunction;

extern const char* DEFAULT_CCY_CHAR;
extern double DEFAULT_CCY_DOUBLE;

///////////////////////////////////////////////////////////////////////////
/// \struct CFTableFunctionArg 
/// declaration of dafault value of an argument.
/// this object can be a string or a double
///////////////////////////////////////////////////////////////////////////
struct ARM_GramFunctionArgDef
{
    const char* itsChar;
    double      itsDouble;
};


/// these are default arguments
extern const ARM_GramFunctionArgDef DefaultPerCcyCharStruct;
extern const ARM_GramFunctionArgDef DefaultPerCcyDbleStruct;
extern const ARM_GramFunctionArgDef DefaultDbleZeroStruct;
extern const ARM_GramFunctionArgDef DefaultDbleOneStruct;
extern const ARM_GramFunctionArgDef DefaultResetTimingStruct;
extern const ARM_GramFunctionArgDef DefaultPaymentTimingStruct;

///////////////////////////////////////////////////////////////////////////
/// \struct CFTableFunctionArg 
/// declaration of an argument to a CFTable function.
/// the list of the function argument is defined as an array of arguments
/// Everything is public for easy C-style struct initialization with {}
///////////////////////////////////////////////////////////////////////////

struct ARM_GramFunctionArgStruct
{
	GramFuncArgType	itsType;

	bool itsIsRequired;

	bool itsIsVaArg;

	const char* itsName;

	const char* itsDescription;

    bool itsIsVisible;

	const ARM_GramFunctionArgDef* itsDefaultValue;

	bool itsAddRefNode;
};



///////////////////////////////////////////////////////////////////////////
/// \struct	  ARM_GramFunctionStruct
/// declaration of CF Table function. The description of its argument refers
/// to an array of CFTableFunctionArg
/// Everything is public for easy C-style struct initialization with {}
/// because we make the difference between MC and Backward induction
/// we have two functions mapping!
///////////////////////////////////////////////////////////////////////////

struct ARM_GramFunctionStruct
{
	typedef ARM_ExpNodeFunc* (*GramNodeBuilderFunc)( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate );

	const char* itsName;

	GramFuncArgType itsReturnType;

	const char* itsDescription;

	const ARM_GramFunctionArgStruct* itsArgs;

	ARM_GramFctor* itsFctor;

	const char* itsCategory;

	GramNodeBuilderFunc itsGramNodeBuilder;
};


///////////////////////////////////////////////////////////////////////////
/// table for global declaration of all generic pricer functions
///////////////////////////////////////////////////////////////////////////

extern const ARM_GramFunctionStruct FunctionsTable[];
extern const size_t FuncTableSize;


/// this is just for the help on the language
extern const ARM_GramFunctionStruct BuiltInFunctionsTable[];
extern const size_t BuiltInFuncTableSize;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

