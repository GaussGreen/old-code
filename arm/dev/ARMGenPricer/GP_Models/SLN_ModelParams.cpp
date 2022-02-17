/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *	\file SLN_ModelParams.cpp
 *
 *  \brief base class for Shifted Log Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/SLN_ModelParams.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SLN_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_SLN_ModelParams::ARM_SLN_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_SLN_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_SLN_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 2 Model Param in total)
	static const string modelParamsName( "Shifted LogNormal Model Param" );

	/// -1) Volatility,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility	);

	///	-2) Shift,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Shift);

}




////////////////////////////////////////////////////
///	Class  : ARM_SLN_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SLN_ModelParams::ARM_SLN_ModelParams( const ARM_SLN_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_SLN_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SLN_ModelParams::~ARM_SLN_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_SLN_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SLN_ModelParams& ARM_SLN_ModelParams::operator=(const ARM_SLN_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_SLN_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_SLN_ModelParams::Clone() const
{
	return new ARM_SLN_ModelParams(*this);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

