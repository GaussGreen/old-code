/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file AnalyticModelParams.cpp
 *
 *  \brief base class for the analytic model params
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"

/// gpbase
#include "gpbase/surface.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"

#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ARM_AnalyticModelParams
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_AnalyticModelParams::ARM_AnalyticModelParams( const ARM_ModelParamVector& modelParams )
:	ARM_ModelParams( modelParams )
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ARM_AnalyticModelParams
///	Returns: nothing
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_AnalyticModelParams::ARM_AnalyticModelParams( const ARM_AnalyticModelParams& rhs )
:	ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ARM_AnalyticModelParams
///	Returns: nothing
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_AnalyticModelParams& ARM_AnalyticModelParams::operator=( const ARM_AnalyticModelParams& rhs )
{
	if( this != &rhs )
		ARM_ModelParams::operator =(rhs);
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ValidateSurfaceModelParam
///	Returns: nothing
///	Action : validates the model params to check that it exists 
///				and that it is a surface model param
////////////////////////////////////////////////////
void ARM_AnalyticModelParams::ValidateSurfaceModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const
{
	if(!DoesModelParamExist(paramNb) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(paramNb) );
	if(!dynamic_cast<const ARM_SurfaceModelParam*>( &GetModelParam(paramNb) ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires surface model params" );
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ValidateCurveModelParam
///	Returns: nothing
///	Action : validates the model params to check that it exists 
///				and that it is a surface model param
////////////////////////////////////////////////////
void ARM_AnalyticModelParams::ValidateCurveModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const
{
	if(!DoesModelParamExist(paramNb) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(paramNb) );
	if(!dynamic_cast<const ARM_CurveModelParam*>( &GetModelParam(paramNb) ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires curve model params" );
}

////////////////////////////////////////////////////
///	Class  : ARM_AnalyticModelParams
///	Routine: ValidateSurfaceListModelParam
///	Returns: nothing
///	Action : validates the model params to check that it exists 
///				and that it is a surface list model param
////////////////////////////////////////////////////
void ARM_AnalyticModelParams::ValidateSurfaceListModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const
{
	if(!DoesModelParamExist(paramNb) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(paramNb) );
	if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(paramNb) ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + modelParamsName + " requires surface list model params" );
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

