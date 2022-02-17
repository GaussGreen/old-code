/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: valuetype.h,v $
 * Revision 1.1  2004/07/30 09:52:19  ebenhamou
 * Initial revision
 *
 */

#ifndef _INGPBASE_VALUETYPE_H
#define _INGPBASE_VALUETYPE_H

/// use our macro for namespace
#include "port.h"
#include "env.h"

CC_BEGIN_NAMESPACE( ARM )


typedef enum
{
    ARM_ERR = -1,
	ARM_UNKNOWN = 0,
    ARM_DOUBLE =1,		/// defined with value to make the two enums equal 
	ARM_DOUBLE_TYPE = 1,/// defined with value to make the two enums equal 
    ARM_INT = 2,		/// defined with value to make the two enums equal 
	ARM_INT_TYPE = 2,	/// defined with value to make the two enums equal 
    ARM_STRING = 3,		/// defined with value to make the two enums equal
	ARM_STRING_TYPE = 3,/// for type compliance                            
	ARM_DATE_TYPE,
	ARM_TOBECOMPUTED,
	ARM_BOOL_TYPE,
	ARM_MISSING_TYPE
		
} ARM_GP_VALUE_TYPE;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


