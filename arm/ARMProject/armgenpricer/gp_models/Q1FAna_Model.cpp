/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1FAna_Model.cpp
 *
 *  \brief Q model analytic part
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#include "gpmodels/Q1FAna_Model.h"

/// gpbase
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/curvemodelparam.h"

/// gpmodels
#include "gpmodels/QModelAnalytics.h"
#include "gpmodels/Q1FAna_ModelParams.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Q1FAna_Model::ARM_Q1FAna_Model(const ARM_ZeroCurvePtr& zc, const ARM_Q1FAna_ModelParams& params) 
:	ARM_AnalyticIRModel(zc,params)
{	CC_ARM_SETNAME(ARM_Q1F_MODEL);}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Q1FAna_Model::ARM_Q1FAna_Model(const ARM_Q1FAna_Model& rhs)
:	ARM_AnalyticIRModel(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Q1FAna_Model::~ARM_Q1FAna_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Q1FAna_Model& ARM_Q1FAna_Model::operator=(const ARM_Q1FAna_Model& rhs)
{
	if(this != &rhs)
		ARM_AnalyticIRModel::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Q1FAna_Model::Clone() const
{
	return new ARM_Q1FAna_Model(*this);
}




////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF (currently only forward value!)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Q1FAna_Model::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    // CurveName is not used (only for multi-currencies model)
    ARM_ZeroCurvePtr ZcCurve=GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);

	if( evalTime <= K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,zcT) );
    }
	else
	{
		double zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
		double zc = zcT/zct;
		size_t nbStates= states->size();
		ARM_VectorPtr values(new std::vector<double>(nbStates));
		for(size_t i=0;i<nbStates;++i)
	        (*values)[i] = zc;
		return values;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Q1FAna_Model::VanillaCaplet(
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
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
		return new std::vector<double>(1,0.0);
		
	ARM_VectorPtr zcPay,libor;
	size_t i,nbStates;

	zcPay	= GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

    if(-2 <= fwdEndTime - payTime && fwdEndTime - payTime <= 2)
    {
        /// No convexity : libor = (DFStart/DFEnd - 1.0)/ period
		libor	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		ARM_VectorPtr zcFwdEnd;

		if( fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			zcFwdEnd = zcPay;
        else
			zcFwdEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

		for(i=0;i<zcPay->size(); ++i)
			(*libor)[i] = ((*libor)[i]/(*zcFwdEnd)[i]-1.0)/fwdPeriod;
    }
    else
        /// Compute an adjusted forward rate
        libor = Libor(curveName,evalTime,fwdStartTime,fwdEndTime,fwdPeriod,fwdResetTime,payTime,states);

	nbStates = zcPay->size();

#if defined( __GP_STRICT_VALIDATION )
	if( strikesPerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

	ARM_VectorPtr caplet(new std::vector<double>(nbStates));

	/// do we just need to compute the intrinsic value!
	if( fwdResetTime <= K_NEW_DOUBLE_TOL || evalTime >= fwdResetTime-K_NEW_DOUBLE_TOL )
	{
		for(i=0;i<nbStates;i++)
			(*caplet)[i]=payNotional*period*CC_Max<double>(capFloor*(*libor)[i]-capFloor*strikesPerState[i],0) *(*zcPay)[i];
	}
	else
	{
		double time			= (fwdResetTime-evalTime)/K_YEAR_LEN;

		/// only in analytical case!
		if(nbStates==1)
			((ARM_Q1FAna_ModelParams*) GetModelParams())->CalibrateModelParams( (*libor)[0], time, capFloor );

		double qSigma		= ( (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetOrdinate(0);
		double qParameter	= ( (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		for(i=0;i<nbStates;i++)
			(*caplet)[i]=payNotional*period*QModelAnalytics::BSQFunction((*libor)[i], strikesPerState[i], qSigma, time, qParameter, (*zcPay)[i], capFloor, (*libor)[i] );

	}

	return caplet;
}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_Model
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////

ARM_VectorPtr ARM_Q1FAna_Model::VanillaSwaption(
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


	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
		return new std::vector<double>(1,0.0);

	/// get all the data
	std::vector<double> dummyVector;
	std::vector<double> margin = std::vector<double>(1,0.0);
	bool isDoubleNotional	= true;
	ARM_VectorPtr fixAnnuity= Annuity(curveName,evalTime,fixPayTimes,fixPayPeriods,states);
	ARM_VectorPtr fixAnnuityCloned( (std::vector<double>&) fixAnnuity->Clone() );
	ARM_VectorPtr swapfwd	= SwapRateInPlaceWithComputedAnnuity( curveName, evalTime, floatStartTime, floatEndTime, 
		fixPayTimes,fixPayPeriods, dummyVector, dummyVector, dummyVector, dummyVector, dummyVector, margin,
		isDoubleNotional, fixAnnuityCloned, states );

	size_t nbStates	 = states->size();
	ARM_VectorPtr swaption(new std::vector<double>(nbStates));
	size_t i;

	double time			= (swapResetTime-evalTime)/K_YEAR_LEN;
	/// only in analytical case!
	if(nbStates==1)
		((ARM_Q1FAna_ModelParams*) GetModelParams())->CalibrateModelParams( (*swapfwd)[0], time, callPut);

	double qSigma = ( (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetOrdinate(0);
	if( qSigma <= K_NEW_DOUBLE_TOL )
	{
		for(i=0;i<nbStates;i++)
			(*swaption)[i]=swapNotional*(*fixAnnuity)[i]*
				CC_Max<double>(callPut*((*swapfwd)[i]- strikesPerState(i,0)),0.0);
	}
	else
	{
		double qParameter	= ( (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		for(i=0;i<nbStates;i++)
			(*swaption)[i] = swapNotional* 
				QModelAnalytics::BSQFunction( (*swapfwd)[i], strikesPerState(i,0), qSigma, time, qParameter, (*fixAnnuity)[i], callPut, (*swapfwd)[i] );


	}

	return swaption;    
}


////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Q1FAna_Model::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F Q Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Q1FAna_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_Q1FAna_ModelParams* Q1FAna_ModelParams = dynamic_cast<const ARM_Q1FAna_ModelParams*>(&params);
	if( !Q1FAna_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_Q1FAna_ModelParams" );
	return true;
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

