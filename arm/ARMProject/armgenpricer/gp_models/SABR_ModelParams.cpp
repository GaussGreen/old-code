/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *
 */

/*! \file SABR_ModelParams.cpp
 *
 *  \brief base class for SABR Model Params
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date February 2005
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/SABR_ModelParams.h"

/// gpbase
#include "gpbase/gpmatrix.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"

/// gpclosedforms
#include "gpclosedforms/extendedsabrformula.h"
/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/surfacemodelparam.h"

/// gpmodels
#include "gpmodels/AnalyticModelParams.h"
CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SABR_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_SABR_ModelParams::ARM_SABR_ModelParams( const ARM_ModelParamVector& params )
:	ARM_AnalyticModelParams(params),
	itsIsSigmaParam(false)
{
	ValidateModelParams();
	if(DoesModelParamExist(ARM_ModelParamType::Volatility )) itsIsSigmaParam = true;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_SABR_ModelParams::ValidateModelParams() const
{	
    /// checks that the model has the following model param :( 4 Model Param in total)
	static const string modelParamsName( "SABR Model Param" );

	/// -1) Alpha or Sigma
	if(!DoesModelParamExist(ARM_ModelParamType::Alpha ) && !DoesModelParamExist(ARM_ModelParamType::Volatility ))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + "SABR Model Param requires " 
		+ ARM_ModelParamType::GetTypeString(ARM_ModelParamType::Alpha) + " or " + ARM_ModelParamType::GetTypeString(ARM_ModelParamType::Volatility));

	if(DoesModelParamExist(ARM_ModelParamType::Alpha ) )
		ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Alpha );

	if(DoesModelParamExist(ARM_ModelParamType::Volatility ) )
		ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Volatility );
		
	///	-2) Beta,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Beta );

	///	-3) Correlation,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::Correlation );

	///	-4) VolOfVol,
	ARM_AnalyticModelParams::ValidateSurfaceModelParam( modelParamsName, ARM_ModelParamType::VolOfVol );

}

////////////////////////////////////////////////////
///	Class  : ARM_SABR_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SABR_ModelParams::ARM_SABR_ModelParams( const ARM_SABR_ModelParams& rhs )
: ARM_AnalyticModelParams(rhs),
  itsIsSigmaParam(rhs.itsIsSigmaParam)
{}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_SABR_ModelParams::Clone() const
{
	return new ARM_SABR_ModelParams(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_ModelParams
///	Routines: ComputeImpliedVolatility
///	Returns : double
///	Action  : compute implied Volatility
////////////////////////////////////////////////////
double ARM_SABR_ModelParams::ImpliedVol(double underlying, 
       double strike, 
       double time,
       double tenor,
       int type,
       int IntegrationStep)
{
	double beta				= GetModelParam(ARM_ModelParamType::Beta).GetValue(time,tenor);
	double rho				= GetModelParam(ARM_ModelParamType::Correlation).GetValue(time,tenor);
	double nu				= GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(time,tenor);
    double expiry			= time/K_YEAR_LEN;

	double alpha;			
	if(itsIsSigmaParam)
	{
		double sigma = GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);		
		alpha = SABR_ComputeAlphaFromSigmaATM(underlying,underlying,expiry,sigma,beta,rho,nu,type);
	}
	else
	{
		alpha = GetModelParam(ARM_ModelParamType::Alpha).GetValue(time,tenor) ;
	}

    double impliedVol = SABR_ComputeImpliedVol(underlying, strike, expiry, alpha, beta,  rho,  nu,type);

    return impliedVol;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_ModelParams
///	Routines: PartialDerivatives
///	Returns : double
///	Action  : to compute paratial derivative
////////////////////////////////////////////////////
double ARM_SABR_ModelParams::PartialDerivative(double underlying, 
       double strike, 
       double time,
       double tenor,
       int type,
       int IntegrationStep,
	   int modelParamType)
{
	double value;
	double alpha			= GetModelParam(ARM_ModelParamType::Alpha).GetValue(time,tenor);
	double beta				= GetModelParam(ARM_ModelParamType::Beta).GetValue(time,tenor);
	double rho				= GetModelParam(ARM_ModelParamType::Correlation).GetValue(time,tenor);
	double nu				= GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(time,tenor);
    double expiry			= time/K_YEAR_LEN;

   switch( modelParamType )
	{
    case ARM_ModelParamType::Alpha:
        {
            value = SABR_ComputePartialDerivative( underlying,
                    strike,
					time/K_YEAR_LEN,
					alpha,
					beta,
					rho,
					nu,
                    "Alpha",
                     type );
            break;
        }

    case ARM_ModelParamType::Beta:
        {
            value = SABR_ComputePartialDerivative( underlying,
                    strike,
					time/K_YEAR_LEN,
					alpha,
					beta,
					rho,
					nu,
                    "Beta",
                     type );
            break;
        }

    case ARM_ModelParamType::Correlation:
        {
            value = SABR_ComputePartialDerivative( underlying,
                    strike,
					time/K_YEAR_LEN,
					alpha,
					beta,
					rho,
					nu,
                    "Correlation",
                     type );
            break;
        }

    case ARM_ModelParamType::VolOfVol:
        {
            value = SABR_ComputePartialDerivative( underlying,
                    strike,
					time/K_YEAR_LEN,
					alpha,
					beta,
					rho,
					nu,
                    "VolOfVol",
                     type );
            break;
        }

    default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": " + ARM_ModelParamType::GetTypeString(ARM_ParamType(modelParamType)) + " can't be calibrated" );

    }
    return value;
}
////////////////////////////////////////////////////
///	Class   : ARM_SABR_ModelParams
///	Routine : PostProcessing
///	Returns : double
///	Action  : update the final surafec after calibration
////////////////////////////////////////////////////
void ARM_SABR_ModelParams::PostProcessing(const ARM_ModelFitter& modelFitter,
	 ARM_PricingModel* model ,
	 int factorNb)
{
    for(int i=0; i<size(); ++i)
    {
        if (ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(GetModelParams()[i]))
        {
            ARM_Surface* surface = surfaceCalibParam->GetSurface();
            for(size_t j=0; j<surface->GetX3().rows()  ;++j)
                for(size_t k=0; k<surface->GetX3().cols()  ;++k)
                {
                    double val = surface->Interpolate(surface->GetX1()[j],surface->GetX2()[k]);
                    surface->insertAtPoint(j, k, val );
                }
        }
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_SABR_ModelParams
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_SABR_ModelParams::StateLocalVariance(double a,double b) const
{
	return 0;
}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_ModelParams
///	Routine: StateLocalDrift
///	Returns: value of the variance
///	Action : Drift in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_SABR_ModelParams::StateLocalDrift(double a,double b) const
{
	return 0;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

