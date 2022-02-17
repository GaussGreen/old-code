/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_ModelParams.cpp
 *
 *  \brief file for the model params of Heston
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date September 2004
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/Heston_ModelParams.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Heston_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_Heston_ModelParams::ARM_Heston_ModelParams( const ARM_ModelParamVector& params )
:	ARM_ModelParamsHW1FStd(params)
{
	if( params.size() < 7 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected at least 7 model parameters: Volatility, InitialVol, LongTermVol, VolOfVol, VolMeanReversion, Correlation and Betaparam! MeanReversion, Sigma and Alpha are optional.");

	if (!DoesModelParamExist(ARM_ModelParamType::MeanReversion))
	{
		ARM_CurveModelParam* meanreversion = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,0.0);
		SetModelParam(meanreversion);
	}
	if (!DoesModelParamExist(ARM_ModelParamType::Sigma))
	{
		ARM_CurveModelParam* sigma = new ARM_CurveModelParam(ARM_ModelParamType::Sigma,0.0);
		SetModelParam(sigma);
	}
	if (!DoesModelParamExist(ARM_ModelParamType::Alpha))
	{
		ARM_CurveModelParam* alpha = new ARM_CurveModelParam(ARM_ModelParamType::Alpha,0.0);
		SetModelParam(alpha);
	}	

	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_Heston_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_Heston_ModelParams::ValidateModelParams() const
{
	if( !DoesModelParamExist(ARM_ModelParamType::Volatility) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires volatility parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::MeanReversion) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires meanreversion parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::InitialVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires InitialVol parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::LongTermVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a long term vol!");

	if( !DoesModelParamExist(ARM_ModelParamType::VolOfVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a vol of vol!");

	if( !DoesModelParamExist(ARM_ModelParamType::VolMeanReversion) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a vol MR!");

	if( !DoesModelParamExist(ARM_ModelParamType::Correlation) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a correlation!");

	if( !DoesModelParamExist(ARM_ModelParamType::Beta) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires beta parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::Sigma) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires alpha parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::Alpha) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires alpha parameter!");

}




////////////////////////////////////////////////////
///	Class  : ARM_Heston_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Heston_ModelParams::ARM_Heston_ModelParams( const ARM_Heston_ModelParams& rhs )
: ARM_ModelParamsHW1FStd(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Heston_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Heston_ModelParams::~ARM_Heston_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Heston_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Heston_ModelParams& ARM_Heston_ModelParams::operator=(const ARM_Heston_ModelParams& rhs)
{
	if(this != &rhs)
		ARM_ModelParamsHW1FStd::operator=(rhs);
		//ARM_ModelParams::operator=(rhs);
	return *this;
}
double ARM_Heston_ModelParams::MC_scheme() const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::QParameter)).GetCurve()->Interpolate(0);
}
double ARM_Heston_ModelParams::InitialVol() const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::InitialVol)).GetCurve()->Interpolate(0);
}
double ARM_Heston_ModelParams::Scaling(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::LongTermVol(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::LongTermVol)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::VolOfVol(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::VolOfVol)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::VolMeanReversion(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::VolMeanReversion)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::Correlation(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::Beta(double t) const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->Interpolate(t);
}
double ARM_Heston_ModelParams::Alpha() const
{
	return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Alpha)).GetCurve()->Interpolate(0);
}

ARM_CurveModelParam ARM_Heston_ModelParams::GetVolAlpha() const
{
	double alpha = Alpha();

	std::vector<double> abs = ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Sigma)).GetCurve()->GetAbscisses();
	std::vector<double> ord = ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Sigma)).GetCurve()->GetOrdinates();

	ord *= alpha;

	ARM_CurveModelParam curveModelParam(ARM_ModelParamType::Sigma,&ord,&abs);

	return curveModelParam;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

