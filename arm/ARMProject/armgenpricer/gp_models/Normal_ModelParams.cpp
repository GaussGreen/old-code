/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file Normal_ModelParams.cpp
 *
 *  \brief base class for NormalModel Model Params
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/Normal_ModelParams.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Normal_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_Normal_ModelParams::ARM_Normal_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_Normal_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_Normal_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 1 Model Param in total)
	static const string modelParamsName( "Normal Model Param" );

	/// -1) StandardDeviation,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility	);


}




////////////////////////////////////////////////////
///	Class  : ARM_Normal_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Normal_ModelParams::ARM_Normal_ModelParams( const ARM_Normal_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Normal_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Normal_ModelParams::~ARM_Normal_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Normal_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Normal_ModelParams& ARM_Normal_ModelParams::operator=(const ARM_Normal_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Normal_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Normal_ModelParams::Clone() const
{
	return new ARM_Normal_ModelParams(*this);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

