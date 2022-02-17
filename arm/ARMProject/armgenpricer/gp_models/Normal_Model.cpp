/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Normal_Model.cpp
 *
 *  \brief base class for NormalModel Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

#include "gpmodels/Normal_Model.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"

/// gpclosedforms
#include "gpclosedforms/vanille_normal_interface.h"

/// gpmodels
#include "gpmodels/Normal_ModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Normal_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Normal_Model
////////////////////////////////////////////////////

ARM_Normal_Model::ARM_Normal_Model(const ARM_ZeroCurvePtr& zc, const ARM_Normal_ModelParams& params  )
:	ARM_AnalyticIRModel( zc, params )
{}


////////////////////////////////////////////////////
///	Class   : ARM_Normal_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_Normal_Model::ARM_Normal_Model(const ARM_Normal_Model& rhs)
:	ARM_AnalyticIRModel( rhs )
{}


////////////////////////////////////////////////////
///	Class  : ARM_Normal_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_Normal_Model& ARM_Normal_Model::operator=(const ARM_Normal_Model& rhs)
{
	if( this != &rhs )
	{
		ARM_AnalyticIRModel::operator =(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Normal_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Normal_Model::~ARM_Normal_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Normal_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Normal_Model::VanillaCaplet(
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
	double StandardDev		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);

	double optionValue = payNotional * period * (*DF)[0] * Export_normal_VanillaOption(
		(*liborRate)[0],strikesPerState[0],StandardDev,time/K_YEAR_LEN, capFloor);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}

////////////////////////////////////////////////////
///	Class   : ARM_Normal_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_Normal_Model::VanillaSwaption(
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
	double StandardDev		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);

	double optionValue = swapNotional * (*annuity)[0] * Export_normal_VanillaOption(
		(*swapRate)[0],strikesPerState(0,0),StandardDev,time/K_YEAR_LEN, callPut);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



    
////////////////////////////////////////////////////
///	Class   : ARM_Normal_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Normal_Model::Clone() const
{
	return new ARM_Normal_Model(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_Normal_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Normal_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Normal Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_Normal_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Normal_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_Normal_ModelParams* Normal_ModelParams = dynamic_cast<const ARM_Normal_ModelParams*>(&params);
	if( !Normal_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_Normal_ModelParams" );
	return true;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

