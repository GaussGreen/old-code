/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_gp_local_excelManip.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */


////////////////////////////////////////////////
/// This file is for manipulation of excel type
/// pricer... Everything specific to the generic
/// pricer should be included here and not in ARM_local_interglob.cpp
/// to avoid exporting the generic pricer in the 
/// local dllarm project (activex project)
/// since the interglob.h file is included in the 
/// local dllarm project (activex project)
////////////////////////////////////////////////

#include <libCCxll\CCxll.h> /// for excel constants
#include "ARM_gp_local_excelManip.h"

/// remove namespace
using ARM::ARM_DOUBLE_TYPE;
using ARM::ARM_STRING_TYPE;
using ARM::ARM_BOOL_TYPE;
using ARM::ARM_MISSING_TYPE;
using ARM::ARM_MISSING_TYPE;
using ARM::ARM_INT_TYPE;
using ARM::ARM_ERR;
using ARM::ARM_DATE_TYPE;


////////////////////////////////////////////////////
//	Routine: ConvertXLtoARMType
//	Returns: ARM_GP_VALUE_TYPE
//	Action : convert an xl type to an ARM_GP_VALUE_TYPE
////////////////////////////////////////////////////

ARM_GP_VALUE_TYPE ConvertXLtoARMType( long type )
{
	switch( type )
	{
	case xltypeNum:			return ARM_DOUBLE_TYPE;
	case xltypeStr:			return ARM_STRING_TYPE;
	case xltypeBool:		return ARM_BOOL_TYPE;
	case xltypeMissing:		return ARM_MISSING_TYPE;
	case xltypeNil:			return ARM_MISSING_TYPE;
	case xltypeInt:			return ARM_INT_TYPE;
	default:				return ARM_ERR;
	}
}


////////////////////////////////////////////////////
//	Routine: ConvertARMTypeToXL
//	Returns: ARM_GP_VALUE_TYPE
//	Action : convert an ARM_GP_VALUE_TYPE to an xl type
////////////////////////////////////////////////////

long ConvertARMToXLType( ARM_GP_VALUE_TYPE type )
{
	switch( type )
	{
	case ARM_DOUBLE_TYPE :	return xltypeNum;
	case ARM_STRING_TYPE :	return xltypeStr;
	case ARM_BOOL_TYPE :	return xltypeBool;
	case ARM_INT_TYPE :		return xltypeInt;
	default:				return xltypeMissing;
	}
}


////////////////////////////////////////////////////
//	Routine: ConvertToGenPricerType
//	Returns: ARM_GP_VALUE_TYPE
//	Action : convert general ARM_GP_VALUE_TYPE to
//	simple ARM_GP_VALUE_TYPE supported by the generic
//	pricer
////////////////////////////////////////////////////

ARM_GP_VALUE_TYPE ConvertToGenPricerType( ARM_GP_VALUE_TYPE type, string& val )
{
	/// these are constant that should be moved to a common header
	/// a double is considered to be a date if between the 1Jan1990 (32874)
	/// and 1Jan 2150 = 91313

	const double XLDATEMIN		= 32874.00;
	const double XLDATEMAX		= 91313.00;
	const double JULIANDATEADD	= 2415019.0;

	switch( type )
	{
		case ARM_DOUBLE_TYPE: case ARM_INT_TYPE:
		{
			/// 
			if( atof( val.c_str() ) > XLDATEMIN && atof( val.c_str() ) < XLDATEMAX )
			{
				int result = atof( val.c_str() ) + JULIANDATEADD;
				char msg[20];
				/// put 3f toget some space!
				sprintf( msg, "%.1f", (double) result );
				val = string( msg );
				return ARM_DATE_TYPE;
			}
			else
				return ARM_DOUBLE_TYPE;
		}

		case ARM_STRING_TYPE:
		{
			/// blank are realling missing type
			if( val == "" || val == " ")
				return ARM_MISSING_TYPE;
			else
				return ARM_STRING_TYPE;
		}

		case ARM_BOOL_TYPE:
			return ARM_DOUBLE_TYPE;

		case ARM_MISSING_TYPE: case ARM_ERR: 
			return ARM_MISSING_TYPE;

		default:
			return ARM_STRING_TYPE;
	}
}

////////////////////////////////////////////////////
//	Routine: ConvertFromGenPricerType
//	Returns: ARM_GP_VALUE_TYPE
//	Action : convert general ARM_GP_VALUE_TYPE to
//	simple ARM_GP_VALUE_TYPE supported by the generic
//	pricer
////////////////////////////////////////////////////

ARM_GP_VALUE_TYPE ConvertFromGenPricerType( ARM_GP_VALUE_TYPE type)
{
	switch( type )
	{
		case ARM_DOUBLE_TYPE: case ARM_DATE_TYPE:
	        return ARM_DOUBLE_TYPE;

		default:
			return type;
	}
}

