/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *
 */

/*! \file SABR_CurveModelParams.cpp
 *
 *  \brief base class for SABR Model Params
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/SABR_CurveModelParams.h"
/// gpclosedforms
#include "gpclosedforms/extended_sabr_interface.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"
CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_SABR_CurveModelParams::ARM_SABR_CurveModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params)
{
	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_CurveModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_SABR_CurveModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 4 Model Param in total)
	static const string modelParamsName( "SABR Model Param" );

	/// -1) Alpha,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::Alpha );

	///	-2) Beta,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::Beta );

	///	-3) Correlation,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::Correlation );

	///	-4) VolOfVol,
	ARM_AnalyticModelParams::ValidateCurveModelParam( modelParamsName, ARM_ModelParamType::VolOfVol );

}

////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SABR_CurveModelParams::ARM_SABR_CurveModelParams( const ARM_SABR_CurveModelParams& rhs )
: ARM_AnalyticModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SABR_CurveModelParams::~ARM_SABR_CurveModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SABR_CurveModelParams& ARM_SABR_CurveModelParams::operator=(const ARM_SABR_CurveModelParams& rhs)
{
	if(this != &rhs)
		ARM_AnalyticModelParams::operator=(rhs);
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_CurveModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_SABR_CurveModelParams::Clone() const
{
	return new ARM_SABR_CurveModelParams(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_CurveModelParams
///	Routines: ComputeImpliedVolatility
///	Returns : double
///	Action  : compute implied Volatility
////////////////////////////////////////////////////
double ARM_SABR_CurveModelParams::ComputeImpliedVolatility(double underlying, 
       double strike, 
       double time,
       double tenor,
       int type,
       int IntegrationStep)
{
	double alpha			= GetModelParam(ARM_ModelParamType::Alpha).GetValue(time,tenor);
	double beta				= GetModelParam(ARM_ModelParamType::Beta).GetValue(time,tenor);
	double rho				= GetModelParam(ARM_ModelParamType::Correlation).GetValue(time,tenor);
	double nu				= GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(time,tenor);
    double expiry			= time/K_YEAR_LEN;

    double impliedVol = Export_SABR_ImplicitVol(underlying,strike,expiry,alpha,beta,rho,nu,type,IntegrationStep);

    return impliedVol;

}

////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_SABR_CurveModelParams::StateLocalVariance(double a,double b) const
{
	return 0;
}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_CurveModelParams
///	Routine: StateLocalDrift
///	Returns: value of the variance
///	Action : Drift in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_SABR_CurveModelParams::StateLocalDrift(double a,double b) const
{
	return 0;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

