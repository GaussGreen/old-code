/*!
 *
 * Copyright (c) IXIS CIB Paris 2005
 *	\file ForwardMarginIR.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2005
 */

/// gpinfra
#include "gpinfra/zccrvfunctor.h"

/// gpmodels
#include "gpmodels/ForwardMarginIR.h"

#define CHANGE_FUNCTORS() \
	ARM_ZeroCurveFunctor* oldDiscountFunctor=GetRefModel()->GetDiscountFunctor(); \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetDiscountFunctor( GetDiscountFunctor() ); \
	ARM_ZeroCurveFunctor* oldFixingFunctor=GetRefModel()->GetFixingFunctor(); \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetFixingFunctor( GetFixingFunctor() );

#define CAST_TO_IRMODEL() \
	ARM_PricingFunctionIR* IRMODEL = dynamic_cast<ARM_PricingFunctionIR*>( const_cast<ARM_PricingModel *>( GetRefModel() ));  \
	if( !IRMODEL ) \
		ARM_THROW(  ERR_INVALID_ARGUMENT, "could not find a reference model for " + modelName + "!" );

#define RESTORE_FUNCTORS() \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetDiscountFunctor(oldDiscountFunctor); \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetFixingFunctor(oldFixingFunctor);



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: Libor
///	Returns: a vector of libor rates
///	Action : Compute libor rates for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::Libor( 
	const string& modelName, 
	double evalTime,
	double fwdStartTime,
	double fwdEndTime,
	double period,
	double fwdResetTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->Libor( modelName, evalTime, fwdStartTime, fwdEndTime, 
		period, fwdResetTime, payTime, states );
	RESTORE_FUNCTORS();
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: SwapRate
///	Returns: a vector of swap rates
///	Action : Compute swap rates for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::SwapRate(
	const string& modelName, 
	double evalTime,
	double floatStartTime, 
	double floatEndTime, 
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const std::vector<double>& fwdStartTimes,
    const std::vector<double>& fwdEndTimes,
    const std::vector<double>& fwdPayPeriods,
	const std::vector<double>& floatPayTimes,
    const std::vector<double>& floatPayPeriods,
    const std::vector<double>& margin,
    bool isDbleNotional,
	const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->SwapRate( modelName, evalTime, floatStartTime, floatEndTime, fixPayTimes,
		fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods,
		margin, isDbleNotional, states );
	RESTORE_FUNCTORS();
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: NPVSwap
///	Returns: a vector of swap NPV
///	Action : Compute the swap NPV for the input model
////////////////////////////////////////////////////
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::NPVSwap(
	const string& modelName, 
	double evalTime,
	double floatStartTime,
	double floatEndTime, 
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const std::vector<double>& fwdStartTimes, 
	const std::vector<double>& fwdEndTimes, 
	const std::vector<double>& fwdPayPeriods, 
	const std::vector<double>& floatPayTimes, 
	const std::vector<double>& floatPayPeriods, 
	const std::vector<double>& margin,
	bool isDbleNotional,
	const std::vector<double>& FixNotional,
	const std::vector<double>& FloatNotional,
	const ARM_GP_Matrix& strikesPerState,
	int payRec,
	const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->NPVSwap( modelName, evalTime, floatStartTime, floatEndTime,fixPayTimes,
		fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, 
		floatPayTimes, floatPayPeriods, margin, isDbleNotional, FixNotional,FloatNotional,
		strikesPerState, payRec, states );
	RESTORE_FUNCTORS();
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwalLeg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_ForwardMarginIR::NPVSwapLeg(
		const string& modelName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& payTimes, 
		const std::vector<double>& payPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_GP_MatrixPtr result= DefaultNPVSwapLeg(modelName, evalTime,fwdStartTimes,fwdEndTimes, fwdPayPeriods, 
		 payTimes, payPeriods, margin, notional,  states); 

	RESTORE_FUNCTORS();

	return result;
}
////////////////////////////////////////////////////
///	Class   : ARM_ForwardMarginIR
///	Routines: NPVFixLeg
///	Returns : vector ptr
///	Action  : Fix leg computing
ARM_GP_MatrixPtr ARM_ForwardMarginIR::NPVFixLeg(
		const string& modelName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_GP_MatrixPtr result= IRMODEL->NPVFixLeg(modelName,evalTime,fixPayTimes,
		fixPayPeriods,FixNotional, strikesPerState,payRec,states);
	
	RESTORE_FUNCTORS();
	return result;
}
////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: VanillaCaplet
///	Returns: a vector of caplet prices
///	Action : Compute caplet price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::VanillaCaplet(
	const string& modelName, 
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
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaCaplet( modelName, evalTime, payTime, period, payNotional,
		fwdResetTime, fwdStartTime, fwdEndTime, fwdPeriod, strikesPerState,
		capFloor, states);
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: VanillaDigital
///	Returns: a vector of digital prices
///	Action : Compute digital price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::VanillaDigital(
	const string& modelName, 
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
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaDigital( modelName, evalTime, payTime, period,
		payNotional, fwdResetTime, fwdStartTime, fwdEndTime, fwdPeriod,
		strikesPerState, capFloor, states );
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: VanillaCorridorlet
///	Returns: a vector of VanillaCorridorlet prices
///	Action : Compute corridor leg price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::VanillaCorridorlet(
	const string& modelName, 
	double evalTime,
    double payTime,
    double resetTime,
    double startTime,
    double endTime,
    int indexPaymentType,
    double fwdPaymentPeriod,
    const std::vector<double>& refIdxResetTimes,
    const std::vector<double>& refIdxStartTimes,
    const std::vector<double>& refIdxEndTimes,
    const std::vector<double>& refFwdPeriods,
    const std::vector<double>& refIndexWeight,
    double couponMargin,
    const vector<const std::vector<double>*> downBarrierPerState,
    const vector<const std::vector<double>*> upBarrierPerState,
    double payNotional,
    int capFloor,
    const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaCorridorlet( modelName,  evalTime, payTime, resetTime,
		startTime, endTime, indexPaymentType, fwdPaymentPeriod, refIdxResetTimes,	
		refIdxStartTimes, refIdxEndTimes, refFwdPeriods, refIndexWeight,
		couponMargin, downBarrierPerState, upBarrierPerState,
		payNotional, capFloor, states );
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: VanillaSwaption
///	Returns: a vector of swaption prices
///	Action : Compute swaption price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginIR::VanillaSwaption(
	const string& modelName,
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
	const std::vector<double>& fixTimes,
    const std::vector<double>& fixPayPeriods,
    const ARM_GP_Matrix& strikesPerState,
    int callPut,
    const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();

	ARM_VectorPtr result = IRMODEL->VanillaSwaption( modelName, evalTime, swapResetTime, fixNotional,floatNotional,
		floatStartTime, floatEndTime, floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,fixTimes, fixPayPeriods, strikesPerState,
		callPut, states,isConstantNotional,isConstantSpread,isConstantStrike);
	
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class   : ARM_ForwardMarginIR
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardMarginIR::LocalDiscounts( size_t timeIdx, 
		double dt, 
		const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	const ARM_PricingModel* REFMODEL = GetRefModel();
	if( !REFMODEL )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "could not find a reference model!" );

	ARM_VectorPtr result = REFMODEL->LocalDiscounts( timeIdx, dt, states);

	RESTORE_FUNCTORS();
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginIR
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ForwardMarginIR::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_ForwardMarginIR Model\n";
    os << indent << "-------------------------\n";

    os << "Model name : " << GetModelName() << "\n";

	os << ARM_ForwardMargin::toString(indent,nextIndent);
	return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

