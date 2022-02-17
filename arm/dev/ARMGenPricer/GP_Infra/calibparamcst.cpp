/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file calibparamcst.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#include "gpinfra/calibparamcst.h"

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////
///	Class		: ARM_CalibParamCst
///	Description	: Static variable of the struct ARM_CalibParamCst
///			that gives upper and lower bound
////////////////////////////////////////

const double ARM_CalibParamCst::CalibLowerBound	= -10.0;
const double ARM_CalibParamCst::CalibUpperBound	=  10.0;
const double ARM_CalibParamCst::CalibZeroBound  =  1.0e-10;



///////////////////////////////////////
///	Class   : ARM_CalibParamCst
///	Routines: ModelParamStdLowerBound
///	Returns :
///	Action  : provides the lower bound for a given model param
////////////////////////////////////////
double ARM_CalibParamCst::ModelParamStdLowerBound( ARM_ModelParamType::ParamNb paramNb )
{
	switch( paramNb )
	{
		case ARM_ModelParamType::Volatility:
        case ARM_ModelParamType::VolOfVol:
		case ARM_ModelParamType::Skew:
			return ARM_CalibParamCst::CalibZeroBound;
		default:
			return ARM_CalibParamCst::CalibLowerBound;
	}
}



///////////////////////////////////////
///	Class   : ARM_CalibParamCst
///	Routines: ModelParamStdUpperBound
///	Returns :
///	Action  : provides the upper bound for a given model param
////////////////////////////////////////
double ARM_CalibParamCst::ModelParamStdUpperBound( ARM_ModelParamType::ParamNb paramNb )
{
	return ARM_CalibParamCst::CalibUpperBound;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



