/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * \file Merton_ModelParams.cpp
 *
 *  \brief base class for Merton Model Params
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/Merton_ModelParams.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Merton_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_Merton_ModelParams::ARM_Merton_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_Merton_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_Merton_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 8 Model Param in total)
	static const string modelParamsName( "Merton Model Param" );

	/// -1) Volatility,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility	);

	/// -2) JumpProba,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::JumpProba	);

	/// -3) JumpSize,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::JumpSize	);

	/// -4) JumpVol,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::JumpSize	);


}




////////////////////////////////////////////////////
///	Class  : ARM_Merton_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Merton_ModelParams::ARM_Merton_ModelParams( const ARM_Merton_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Merton_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Merton_ModelParams::~ARM_Merton_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Merton_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Merton_ModelParams& ARM_Merton_ModelParams::operator=(const ARM_Merton_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Merton_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Merton_ModelParams::Clone() const
{
	return new ARM_Merton_ModelParams(*this);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

