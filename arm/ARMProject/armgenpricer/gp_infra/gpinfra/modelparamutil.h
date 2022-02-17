/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: modelparamutil.h,v $
 * Revision 1.1  2003/10/08 16:45:26  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file modelparamutil.h
 *
 *  \brief 
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_MODELPARAMUTIL_H
#define _INGPINFRA_MODELPARAMUTIL_H

#include "gpbase/port.h"
#include <functional>
#include "modelparam.h"


CC_BEGIN_NAMESPACE( ARM )

/// the default param to compare with is the right param
/// therefore the argument is on the left side!


struct CompareModelParamWEnumUnaryLeftVersion : public CC_NS( std, unary_function )<ARM_ModelParam*,bool>
{
	CompareModelParamWEnumUnaryLeftVersion( ARM_ModelParam* param2 )
	:	itsParam2( param2 )
	{
		if( param2 == NULL )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model param2* is null!" );
	}
    
	bool operator()( ARM_ModelParam* param1 ) const
    {
#if defined(__GP_STRICT_VALIDATION )
		if( param1 == NULL )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model param1* is null!" );
#endif
        /// equality is only based on the type!
        return (param1->GetType() == itsParam2->GetType());
    }

	private:
		ARM_ModelParam* itsParam2;
};



struct FindModelParamWEnumUnaryVersion : public CC_NS( std, unary_function )<ARM_ModelParam*,bool>
{
	FindModelParamWEnumUnaryVersion( ARM_ModelParamType::ParamNb type )
	:	itsParamType( type ) 
	{}
    
	bool operator()( ARM_ModelParam* param1 ) const
    {
#if defined(__GP_STRICT_VALIDATION )
		if( param1 == NULL )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model param* is null!" );
#endif
        /// equality is only based on the type!
        return (param1->GetType() == itsParamType);
    }
private:
	ARM_ModelParamType::ParamNb itsParamType;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

