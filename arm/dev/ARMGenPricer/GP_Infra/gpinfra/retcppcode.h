/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: retcppcode.h,v $
 *
 * Revision 1.1  2003/10/02 17:51:46  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file retcppcode.h
 *
 *  \brief c++ version of the retcode
 */
 
#ifndef _INGPINFRA_RETCPPODE_H
#define _INGPINFRA_RETCPPODE_H

#include "gpbase/valuetype.h"
#include <glob/expt.h>
#include "gpbase/port.h"
#include <string>


CC_BEGIN_NAMESPACE( ARM )
 
struct ARM_TYPE
{
	static CC_NS( std, string ) Name( ARM_GP_VALUE_TYPE t )
	{
		switch( t )
		{
		case ARM_UNKNOWN: 
			return "Unknown";
		case ARM_ERR :
			return "Error";
		case ARM_DOUBLE:
			return "Double";
		case ARM_INT:
			return "Integer";
		case ARM_STRING:
			return "String";
		case ARM_DATE_TYPE:
			return "Date";
		case ARM_TOBECOMPUTED:
			return "ToBeComputed";
		case ARM_BOOL_TYPE:
			return "BoolType";
		case ARM_MISSING_TYPE:
			return "Missing";
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"This is an unknown type not even defined as unknown or error..Please advise" );
		}
	}
};

CC_END_NAMESPACE()

#endif 

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
