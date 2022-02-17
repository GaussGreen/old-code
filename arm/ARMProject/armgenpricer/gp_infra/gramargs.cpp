/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramargs.cpp,v $
 * Revision 1.1  2003/10/08 16:40:30  ebenhamou
 * Initial revision
 *
 *
 */



/*! \file gramargs.cpp
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpinfra/gramargs.h"

#include <glob/expt.h>

#include <string>
CC_USING_NS( std, string )



CC_BEGIN_NAMESPACE( ARM )

/// conversion to string
/// used to print an argument
string GramFuncArgTypeToString( GramFuncArgType t )
{
	switch( t )
	{
	case GFAT_TERMINATOR:
		return "TERMINATOR";
	case GFAT_UNKNOWN_TYPE:
		return "UNKNOWN";
	case GFAT_DOUBLE_TYPE:
		return "DOUBLE";
	case GFAT_VECTOR_TYPE:
		return "VECTOR";
	case GFAT_MATRIX_TYPE:
		return "MATRIX";
	case GFAT_STRING_TYPE:
		return "STRING";
	case GFAT_DATE_TYPE:
		return "DATE";
	case GFAT_FUTUREREF_TYPE:
		return "FUTURE_REFERENCE";
	case GFAT_VECTOR_OR_FUTUREREF_TYPE:
		return "VECTOR_OR_FUTURE_REFERENCE";
	case GFAT_MODEL_TYPE:
		return "MODEL";
	case GFAT_MATURITY_TYPE:
		return "MATURITY";
	case GFAT_DATEORMATU_TYPE:
		return "DATE_OR_MATURITY";
	case GFAT_MULTITOKENSTRING_TYPE:
		return "MULTI_TOKENS_STRING";
	case GFAT_DATE_OR_DOUBLE_TYPE:
		return "DATE_OR_DOUBLE";
	case GFAT_DATE_OR_VECTOR_TYPE:
		return "DATE_OR_VECTOR";
	case GFAT_DATE_VECTOR_OR_DOUBLE_TYPE:
		return "DATE_VECTOR_OR_DOUBLE";
	case GFAT_STRINGVEC_TYPE:
		return "STRING_VECTOR";
	case GFAT_STRINGVECTRANS_TYPE:
		return "STRING_VECTOR";
	case GFAT_CURVE_TYPE:
		return "CURVE";
	case GFAT_VECTOR_OR_CURVE_TYPE:
		return "CURVE_OR_VECTOR";
	case GFAT_DATESTRIP_TYPE:
		return "DATESTRIP";
	case GFAT_STRING_OR_CURVE_TYPE:
		return "STRING_OR_CURVE";
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"This is an unknown type not even defined as unknown or error..Please advise" );
	}
}



CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/



