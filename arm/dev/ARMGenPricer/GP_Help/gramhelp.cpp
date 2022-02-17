/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramhelp.cpp,v $
 * Revision 1.1  2003/10/30 07:43:31  ebenhamou
 * Initial revision
 *
 */

/*! \file gramhelp.cpp
 *
 *  \brief gramhelp is to provide help on the grammar, 
 *		its function and its syntax
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gphelp/gramhelp.h"
#include "gphelp/helpertext.h"

#include "gpinfra/gramfunction.h"
#include "gpinfra/gramfunctiontable.h"
#include "gpinfra/gramfunctionarg.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"


CC_BEGIN_NAMESPACE( ARM )

void ARM_GramHelper::FunctionTableHelpString( ARM_GramFunctionStruct* funcTable, size_t funcTableSize  )
{
	/// for easier typing
	typedef ARM_GramFunctionArgStruct GFAStruct;
	typedef ARM_GramFunctionArgVector GFAVector;

	/// variables for construction of the gramFunction
	string FuncName;
	string Category;
	GramFuncArgType ReturnType;
	string Description;
	const GFAStruct* ga;
	const int ARGSIZEMAX = 40;

	for(size_t i=0; i<funcTableSize; ++i)
	{
		/// get variables
		FuncName	= funcTable[i].itsName;
		ReturnType	= funcTable[i].itsReturnType;
		Description	= funcTable[i].itsDescription;
		ga			= funcTable[i].itsArgs;
		Category	= funcTable[i].itsCategory;

		/// get the arguments!
		int ArgSize = 0;
		while( ga[ArgSize].itsType != GFAT_TERMINATOR && ArgSize < ARGSIZEMAX )
			++ArgSize;

		int j;
		GFAVector Args(ArgSize);
		for(j=0; j<ArgSize; ++j)
		{
			string DefaultValue;
			GramFuncArgType Type = ga[j].itsType;

			if( !ga[j].itsIsRequired )
			{
#ifdef __GP_STRICT_VALIDATION
				if( ga[j].itsDefaultValue == NULL )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"no default but the table says that it is not a required argument, please advise!" );
#endif
				/// currently default value can be only double, string or maturity
				/// case of double
				if( Type == GFAT_DOUBLE_TYPE  || Type == GFAT_DATE_OR_DOUBLE_TYPE)
				{
					CC_Ostringstream os;
					os << ga[j].itsDefaultValue->itsDouble;
					DefaultValue = os.str();
				}
				else if( Type == GFAT_STRING_TYPE || Type == GFAT_MULTITOKENSTRING_TYPE || Type == GFAT_MATURITY_TYPE || Type == GFAT_STRING_OR_CURVE_TYPE)
					DefaultValue = ga[j].itsDefaultValue->itsChar? ga[j].itsDefaultValue->itsChar :"";
				else if( (Type == GFAT_VECTOR_TYPE) || (Type == GFAT_VECTOR_OR_CURVE_TYPE) )
					DefaultValue = "Vector Computed";
				else if (Type == GFAT_DATESTRIP_TYPE)
					DefaultValue = "Other dates";
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, "defaulted argument can only be string, maturity or double" );
			}

			Args[j]= ARM_GramFunctionArg( ga[j].itsType, ga[j].itsIsRequired, ga[j].itsIsVaArg, ga[j].itsName,
				ga[j].itsDescription,ga[j].itsIsVisible, DefaultValue );
		}

		/// Create the function
		/// and convert it to string
		ARM_GramFunction::Declaration func( FuncName, ReturnType, Description, Args, Category, NULL, NULL ); 
		itsHelpText += func.toString();
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramHelper
///	Routine: ARM_GramHelper
///	Returns: built object
///	Action : Constructor
/////////////////////////////////////////////////////////////////
ARM_GramHelper::ARM_GramHelper( const string& FuncName )
: itsHelpText()
{
	SetName(ARM_GRAMHELPER);

	if( FuncName == "ALLFUNC" )
	{
		/// third get the HelpText
		for(int i=0; i<AdditionalHelpTextSize; ++i)
			itsHelpText += AdditionalHelpText[i] + "\n";
		itsHelpText += "\n\n\n";

		itsHelpText += "List Of All The Functionalities And Tokens\n\n";

		/// first get the FunctionsTable
		itsHelpText += "******************************************************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "***********      1) FUNCTIONS         ****************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "\n\n";

		FunctionTableHelpString( (ARM_GramFunctionStruct*) FunctionsTable, FuncTableSize );

		/// second get the built-in Functions
		itsHelpText += "******************************************************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "***********    2) BUILT IN TOKENS     ****************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "******************************************************\n";
		itsHelpText += "\n\n";

		FunctionTableHelpString( (ARM_GramFunctionStruct*) BuiltInFunctionsTable, BuiltInFuncTableSize );
	}
	else
		itsHelpText = ARM_GramFunction(FuncName, "").toString();
}
	

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GramHelper 
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_GramHelper::Clone() const
{
	return new ARM_GramHelper( * this );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

