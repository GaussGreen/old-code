/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_ModelParams.cpp
 *
 *  \brief file for the model params of Heston
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/ShiftedHeston_ModelParams.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ShiftedHeston_ModelParams::ARM_ShiftedHeston_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_ShiftedHeston_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 6 Model Param in total)
	static const string modelParamsName( "Shifted Heston Model Param" );

	/// -1) InitialVol,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::InitialVol	);

	///	-2) LongTermVol,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::LongTermVol);

	///	-3) VolOfVol,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::VolOfVol	);

	///	-4) VolMeanReversion,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::VolMeanReversion );

	/// -5) Correlation,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::Correlation);

	/// -6) Shift,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::Shift);
}




////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ShiftedHeston_ModelParams::ARM_ShiftedHeston_ModelParams( const ARM_ShiftedHeston_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ShiftedHeston_ModelParams::~ARM_ShiftedHeston_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ShiftedHeston_ModelParams& ARM_ShiftedHeston_ModelParams::operator=(const ARM_ShiftedHeston_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ShiftedHeston_ModelParams::Clone() const
{
	return new ARM_ShiftedHeston_ModelParams(*this);
}




CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

