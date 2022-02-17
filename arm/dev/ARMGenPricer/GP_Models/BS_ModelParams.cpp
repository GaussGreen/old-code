/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *
 */

/*! \file BS_ModelParams.cpp
 *
 *  \brief base class for Black Scholes model params
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/BS_ModelParams.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"
CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_BS_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_BS_ModelParams::ARM_BS_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_BS_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_BS_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 1 Model Param in total)
	static const string modelParamsName( "BS Model Param" );

	/// -1) Volatility,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility	);

}




////////////////////////////////////////////////////
///	Class  : ARM_BS_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BS_ModelParams::ARM_BS_ModelParams( const ARM_BS_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_BS_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_BS_ModelParams::~ARM_BS_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_BS_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_BS_ModelParams& ARM_BS_ModelParams::operator=(const ARM_BS_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_BS_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_BS_ModelParams::Clone() const
{
	return new ARM_BS_ModelParams(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_BS_ModelParams
///	Routines: StateLocalVariance
///	Returns :
///	Action  : computes the state local variance
////////////////////////////////////////////////////
double ARM_BS_ModelParams::StateLocalVariance(double a,double b) const
{
	return 0.0;
}

////////////////////////////////////////////////////
///	Class   : ARM_BS_ModelParams
///	Routines: StateLocalDrift
///	Returns :
///	Action  : computes the state local drift
////////////////////////////////////////////////////
double ARM_BS_ModelParams::StateLocalDrift(double a,double b) const
{
	return 0.0;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

