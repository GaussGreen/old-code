/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file CEV_ModelParams.cpp
 *
 *  \brief base class for CEV Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/CEV_ModelParams.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_CEV_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_CEV_ModelParams::ARM_CEV_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_CEV_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_CEV_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 3 Model Param in total)
	static const string modelParamsName( "CEV Model Param" );

	/// -1) Volatility,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility	);

	///	-2) Beta,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Beta);

	///	-3) Drift,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Drift	);

}




////////////////////////////////////////////////////
///	Class  : ARM_CEV_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_CEV_ModelParams::ARM_CEV_ModelParams( const ARM_CEV_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CEV_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_CEV_ModelParams::~ARM_CEV_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_CEV_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_CEV_ModelParams& ARM_CEV_ModelParams::operator=(const ARM_CEV_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_CEV_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CEV_ModelParams::Clone() const
{
	return new ARM_CEV_ModelParams(*this);
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

