/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsQ1F_Fx.cpp
 *
 *  \brief Q model 1 factor FX version
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpvector.h"

/// gpmodels
#include "gpmodels/ModelParamsCEV_Fx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsCEV_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsCEV_Fx::ARM_ModelParamsCEV_Fx( const ARM_ModelParamVector& params, ARM_ZeroCurvePtr domCurve, ARM_ZeroCurvePtr fgnCurve, double spot  )
: ARM_ModelParams_Fx(domCurve, fgnCurve, spot),
ARM_ModelParamsHW1FStd(params)
{
	SetVolatilityType(ARM_ModelParamType::Volatility);

	if( params.size() != 3 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected 2 model parameters: Volatility, MeanReversion, Beta!" );

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsCEV_Fx
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_ModelParamsCEV_Fx::ValidateModelParams() const
{	
    /// checks that the cev model contains a volatility
	if( !DoesModelParamExist(ARM_ModelParamType::Volatility) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CEVModelParam: requires beta parameter!");

	/// checks that the cev model contains a beta parameter
	if( !DoesModelParamExist(ARM_ModelParamType::Beta) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CEVModelParam: requires a volatility!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsCEV_Fx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsCEV_Fx::ARM_ModelParamsCEV_Fx( const ARM_ModelParamsCEV_Fx& rhs )
: ARM_ModelParams_Fx(rhs),
ARM_ModelParamsHW1FStd(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsCEV_Fx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsCEV_Fx::~ARM_ModelParamsCEV_Fx()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsCEV_Fx
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////
string ARM_ModelParamsCEV_Fx::toString(const string& indent,const string& nextIndent) const
{
	string str;

	str += ARM_ModelParams_Fx::toString(indent,nextIndent);
	str += ARM_ModelParamsHW1FStd::toString(indent,nextIndent);

	return str;
}

CC_END_NAMESPACE()