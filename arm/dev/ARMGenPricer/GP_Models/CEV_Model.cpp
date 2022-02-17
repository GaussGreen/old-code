/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CEV_Model.cpp
 *
 *  \brief base class for CEV Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

#include "gpmodels/CEV_Model.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"

/// gpclosedforms
#include "gpclosedforms/cev_interface.h"

/// gpmodels
#include "gpmodels/CEV_ModelParams.h"

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_CEV_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_CEV_Model
////////////////////////////////////////////////////

ARM_CEV_Model::ARM_CEV_Model(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, size_t IntegrationStep  )
:	ARM_AnalyticIRModel( zc, params ), itsIntegrationStep(IntegrationStep)
{}


////////////////////////////////////////////////////
///	Class   : ARM_CEV_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_CEV_Model::ARM_CEV_Model(const ARM_CEV_Model& rhs)
:	ARM_AnalyticIRModel( rhs ), itsIntegrationStep( rhs.itsIntegrationStep)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CEV_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_CEV_Model& ARM_CEV_Model::operator=(const ARM_CEV_Model& rhs)
{
	if( this != &rhs )
	{
		ARM_AnalyticIRModel::operator =(rhs);
		itsIntegrationStep = rhs.itsIntegrationStep;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_CEV_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_CEV_Model::~ARM_CEV_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_CEV_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_CEV_Model::VanillaCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const ARM_GP_Vector& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr liborRate	= Libor(curveName, evalTime,fwdStartTime, 
		fwdEndTime, period, fwdResetTime, payTime, states );
	ARM_VectorPtr DF		= DiscountFactor( curveName, evalTime, payTime, states );
	
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);
	double Volatility		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);
	double Beta				= GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(time,tenor);
	double Drift			= GetModelParams()->GetModelParam(ARM_ModelParamType::Drift).GetValue(time,tenor);

	double optionValue = payNotional * period * (*DF)[0] * Export_CEV_VanillaOption(
		(*liborRate)[0],strikesPerState[0],time/K_YEAR_LEN,Drift,Volatility,Beta,capFloor,itsIntegrationStep);

	return ARM_VectorPtr( new ARM_GP_Vector(1,optionValue));
}



////////////////////////////////////////////////////
///	Class   : ARM_CEV_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_CEV_Model::VanillaSwaption(
	const string& curveName,
	double evalTime,
	double swapResetTime,
	const ARM_GP_Vector& fixNotional,
	const ARM_GP_Vector& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const ARM_GP_Vector& floatResetTimes,
	const ARM_GP_Vector& floatStartTimes,
	const ARM_GP_Vector& floatEndTimes,
	const ARM_GP_Vector& floatIntTerms,
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
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
	ARM_GP_Vector dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;
	ARM_GP_Vector margin = ARM_GP_Vector(1,0.0);


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
	double Volatility		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);
	double Beta				= GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(time,tenor);
	double Drift			= GetModelParams()->GetModelParam(ARM_ModelParamType::Drift).GetValue(time,tenor);

	double optionValue = swapNotional * (*annuity)[0] * Export_CEV_VanillaOption(
		(*swapRate)[0],strikesPerState(0,0),time/K_YEAR_LEN,Drift,Volatility,Beta, callPut,itsIntegrationStep);

	return ARM_VectorPtr( new ARM_GP_Vector(1,optionValue));
}

    
////////////////////////////////////////////////////
///	Class   : ARM_CEV_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CEV_Model::Clone() const
{
	return new ARM_CEV_Model(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_CEV_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_CEV_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "CEV Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}



////////////////////////////////////////////////////
///	Class   : ARM_CEV_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_CEV_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_CEV_ModelParams* CEV_ModelParams = dynamic_cast<const ARM_CEV_ModelParams*>(&params);
	if( !CEV_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_CEV_ModelParams" );
	return true;
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

