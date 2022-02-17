/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsLN_Fx.cpp
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
#include "gpmodels/ModelParamsLN_Fx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsLN_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsLN_Fx::ARM_ModelParamsLN_Fx( const ARM_ModelParamVector& params, ARM_ZeroCurvePtr domCurve, ARM_ZeroCurvePtr fgnCurve, double spot  )
: ARM_ModelParams_Fx(domCurve, fgnCurve, spot),
ARM_ModelParamsHW1FStd(params),
itsFactorCount(1)
{
	SetVolatilityType(ARM_ModelParamType::Volatility);

	if( params.size() != 2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected 2 model parameters: Volatility, MeanReversion!" );

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsLN_Fx
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_ModelParamsLN_Fx::ValidateModelParams() const
{	
    /// checks that the cev model contains a volatility
	if( !DoesModelParamExist(ARM_ModelParamType::Volatility) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": LNModelParam: requires volatility!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsLN_Fx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsLN_Fx::ARM_ModelParamsLN_Fx( const ARM_ModelParamsLN_Fx& rhs )
: ARM_ModelParams_Fx(rhs),
ARM_ModelParamsHW1FStd(rhs),
itsFactorCount(rhs.itsFactorCount)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsLN_Fx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsLN_Fx::~ARM_ModelParamsLN_Fx()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsLN_Fx
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////
string ARM_ModelParamsLN_Fx::toString(const string& indent,const string& nextIndent) const
{
	string str;

	str += ARM_ModelParams_Fx::toString(indent,nextIndent);
	str += ARM_ModelParamsHW1FStd::toString(indent,nextIndent);

	return str;
}

CC_END_NAMESPACE()