/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file ARM_ModelParamsSFRMFactory.cpp
 *  \brief
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 *
 */


/// this header comes first as it include some preprocessor constants
#include "gpmodels/ModelParamsSFRMFactory.h"
#include "gpmodels/ModelParamsSFRMDiag.h"
#include "gpmodels/ModelParamsSFRMRow.h"

#include "gpbase/singleton.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMFactory
///	Routine: CreateModelParamsSFRM 
///	Returns: 
///	Action : constructs model paramsSFRM of the correct type!
////////////////////////////////////////////////////
ARM_ModelParamsSFRM* ARM_ModelParamsSFRMFactoryImp::CreateModelParamsSFRM(
	const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb, size_t volType ) const
{
	switch( volType )
	{
	case K_DIAG:
		return new ARM_ModelParamsSFRMDiag( params, index, factorsNb );
	case K_ROW:
		return new ARM_ModelParamsSFRMRow( params, index, factorsNb );
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Unknown volType, permitted is ROW and DIAG");
	
	}
}




ARM_SingletonHolder<ARM_ModelParamsSFRMFactoryImp> ARM_ModelParamsSFRMFactory;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
