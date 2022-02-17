/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsHeston_Fx.cpp
 *
 *  \brief Q model 1 factor FX version
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpvector.h"

/// gpmodels
#include "gpmodels/ModelParamsHeston_Fx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Fx::ARM_ModelParamsHeston_Fx( const ARM_ModelParamVector& params, ARM_ZeroCurvePtr domCurve, ARM_ZeroCurvePtr fgnCurve, double spot  )
: ARM_ModelParams_Fx(domCurve, fgnCurve, spot),
ARM_Heston_ModelParams(params)
{
	SetVolatilityType(ARM_ModelParamType::Sigma);

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHeston_Fx
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Heston model
////////////////////////////////////////////////////
void ARM_ModelParamsHeston_Fx::ValidateModelParams() const
{	
	ARM_Heston_ModelParams::ValidateModelParams();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Fx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Fx::ARM_ModelParamsHeston_Fx( const ARM_ModelParamsHeston_Fx& rhs )
: ARM_ModelParams_Fx(rhs),
ARM_Heston_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Fx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Fx::~ARM_ModelParamsHeston_Fx()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Fx
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////
string ARM_ModelParamsHeston_Fx::toString(const string& indent,const string& nextIndent) const
{
	string str;

	str += ARM_ModelParams_Fx::toString(indent,nextIndent);
	str += ARM_Heston_ModelParams::toString(indent,nextIndent);

	return str;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

