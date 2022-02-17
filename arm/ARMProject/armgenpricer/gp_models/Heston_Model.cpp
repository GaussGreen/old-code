/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Model.cpp
 *
 *  \brief base class for Heston Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

#include "gpmodels/Heston_Model.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"

/// gpclosedforms
#include "gpclosedforms/heston_interface.h"

/// gpmodels
#include "gpmodels/Heston_ModelParams.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Heston_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Heston_Model
////////////////////////////////////////////////////

ARM_Heston_Model::ARM_Heston_Model(const ARM_ZeroCurvePtr& zc,
								   const ARM_Heston_ModelParams& params,
								   size_t IntegrationStep  )
:	ARM_AnalyticIRModel( zc, params ), 
	itsIntegrationStep(IntegrationStep)
{}


////////////////////////////////////////////////////
///	Class   : ARM_Heston_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_Heston_Model::ARM_Heston_Model(const ARM_Heston_Model& rhs)
:	ARM_AnalyticIRModel( rhs ), 
	itsIntegrationStep( rhs.itsIntegrationStep)
{}

////////////////////////////////////////////////////
///	Class  : ARM_Heston_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Heston_Model::VanillaCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const std::vector<double>& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr liborRate	= Libor(curveName, evalTime,fwdStartTime, 
		fwdEndTime, period, fwdResetTime, payTime, states );
	ARM_VectorPtr DF		= DiscountFactor( curveName, evalTime, payTime, states );
	
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);
	double V0				= GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(time,tenor);
	double longtermV		= GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(time,tenor);
	double speed			= GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(time,tenor);
	double volvol			= GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(time,tenor);
	double rho				= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(time,tenor);
	double lambda			= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpSize).GetValue(time,tenor);
	double muJ				= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpProba).GetValue(time,tenor);
	double sigmaJ			= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpVol).GetValue(time,tenor);

	double optionValue = payNotional* (*DF)[0] * period * Export_GHeston_VanillaOption((*liborRate)[0], 
						strikesPerState[0], V0,	time/K_YEAR_LEN, longtermV, speed, volvol, rho, lambda, 
						muJ, sigmaJ, capFloor, itsIntegrationStep );

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



////////////////////////////////////////////////////
///	Class   : ARM_Heston_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_Heston_Model::VanillaSwaption(
	const string& curveName,
	double evalTime,
	double swapResetTime,
	const std::vector<double>& fixNotional,
	const std::vector<double>& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const ARM_GP_Matrix& strikesPerState,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
	/// TO BE UPDATED
	/// Check that the notional is constant
	double swapNotional = fixNotional[0];
	if (!(isConstantNotional&&isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );


	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!
	std::vector<double> dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;

	std::vector<double> margin = std::vector<double>(1,0.0);
	ARM_VectorPtr swapRate = SwapRate(
		curveName, 
		evalTime,
		floatStartTime, 
		floatEndTime, 
		fixPayTimes,
		fixPayPeriods,
		dummyFwdStartTimes,
        dummyFwdEndTimes,
        dummyFwdPayTimes,
        dummyFloatPayTimes,
        dummyFloatPayPeriods,
		margin,	/// margin
		true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
		ARM_PricingStatesPtr(NULL) );

	ARM_VectorPtr annuity = Annuity(
		curveName, 
        evalTime,
		fixPayTimes,
        fixPayPeriods,
		ARM_PricingStatesPtr(NULL) );

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double V0			= GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(time,tenor);
	double longtermV	= GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(time,tenor);
	double speed		= GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(time,tenor);
	double volvol		= GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(time,tenor);
	double rho			= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(time,tenor);
	double lambda		= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpSize).GetValue(time,tenor);
	double muJ			= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpProba).GetValue(time,tenor);
	double sigmaJ		= GetModelParams()->GetModelParam(ARM_ModelParamType::JumpVol).GetValue(time,tenor);

	double optionValue = swapNotional * (*annuity)[0] * Export_GHeston_VanillaOption(
		(*swapRate)[0],strikesPerState(0,0), V0,
		time/K_YEAR_LEN, longtermV, speed, volvol, rho, lambda, muJ, sigmaJ, callPut, itsIntegrationStep );

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



    
////////////////////////////////////////////////////
///	Class   : ARM_Heston_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Heston_Model::Clone() const
{
	return new ARM_Heston_Model(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_Heston_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Heston_Model::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Generalized Heston Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}



////////////////////////////////////////////////////
///	Class   : ARM_Heston_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Heston_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_Heston_ModelParams* Heston_ModelParams = dynamic_cast<const ARM_Heston_ModelParams*>(&params);
	if( !Heston_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_Heston_ModelParams" );
	return true;
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

