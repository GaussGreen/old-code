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

#include "gpmodels/ForwardMarginBasis.h"

/// gpinfra
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingstates.h"



#define CHANGE_FUNCTORS() \
	ARM_ZeroCurveFunctor* oldDiscountFunctor=GetRefModel()->GetDiscountFunctor(); \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetDiscountFunctor( GetDiscountFunctor() ); 

#define CAST_TO_IRMODEL() \
	ARM_PricingFunctionIR* IRMODEL = dynamic_cast<ARM_PricingFunctionIR*>( const_cast<ARM_PricingModel *>( GetRefModel() ));  \
	if( !IRMODEL ) \
		ARM_THROW(  ERR_INVALID_ARGUMENT, "could not find a reference model for " + modelName + "!" );

#define RESTORE_FUNCTORS() \
	const_cast<ARM_PricingModel*>( GetRefModel() )->SetDiscountFunctor(oldDiscountFunctor);



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginBasis
///	Routine: Libor
///	Returns: a vector of libor rates
///	Action : Compute libor rates for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::Libor( 
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: SwapRate
///	Returns: a vector of swap rates
///	Action : Compute swap rates for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::SwapRate(
	const string& modelName, 
	double evalTime,
	double floatStartTime, 
	double floatEndTime, 
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
	const ARM_GP_Vector& fwdStartTimes,
    const ARM_GP_Vector& fwdEndTimes,
    const ARM_GP_Vector& fwdPayPeriods,
	const ARM_GP_Vector& floatPayTimes,
    const ARM_GP_Vector& floatPayPeriods,
    const ARM_GP_Vector& margin,
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: NPVSwap
///	Returns: a vector of swap NPV
///	Action : Compute the swap NPV for the input model
////////////////////////////////////////////////////
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::NPVSwap(
	const string& modelName, 
	double evalTime,
	double floatStartTime,
	double floatEndTime, 
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
	const ARM_GP_Vector& fwdStartTimes, 
	const ARM_GP_Vector& fwdEndTimes, 
	const ARM_GP_Vector& fwdPayPeriods, 
	const ARM_GP_Vector& floatPayTimes, 
	const ARM_GP_Vector& floatPayPeriods, 
	const ARM_GP_Vector& margin,
	bool isDbleNotional,
	const ARM_GP_Vector& FixNotional,
	const ARM_GP_Vector& FloatNotional,
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwalLeg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_ForwardMarginBasis::NPVSwapLeg(
		const string& modelName, 
		double evalTime,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& payTimes, 
		const ARM_GP_Vector& payPeriods, 
		const ARM_GP_Vector& margin, 
		const ARM_GP_Vector& notional, 
		const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_GP_MatrixPtr result= IRMODEL->NPVSwapLeg(modelName, evalTime,fwdStartTimes,fwdEndTimes, fwdPayPeriods, 
		 payTimes, payPeriods, margin, notional,  states); 
	RESTORE_FUNCTORS();

	return result;
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMarginBasis
///	Routines: NPVFixLeg
///	Returns : matrix ptr
///	Action  : Fix leg computing
ARM_GP_MatrixPtr ARM_ForwardMarginBasis::NPVFixLeg(
		const string& modelName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: VanillaCaplet
///	Returns: a vector of caplet prices
///	Action : Compute caplet price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::VanillaCaplet(
	const string& modelName, 
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
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaCaplet( modelName, evalTime, payTime, period, payNotional,
		fwdResetTime, fwdStartTime, fwdEndTime, fwdPeriod, strikesPerState,
		capFloor, states);
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginBasis
///	Routine: VanillaDigital
///	Returns: a vector of digital prices
///	Action : Compute digital price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::VanillaDigital(
	const string& modelName, 
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
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaDigital( modelName, evalTime, payTime, period,
		payNotional, fwdResetTime, fwdStartTime, fwdEndTime, fwdPeriod,
		strikesPerState, capFloor, states );
	RESTORE_FUNCTORS();
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMarginBasis
///	Routine: VanillaCorridorlet
///	Returns: a vector of VanillaCorridorlet prices
///	Action : Compute corridor leg price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::VanillaCorridorlet(
	const string& modelName, 
	double evalTime,
    double payTime,
    double resetTime,
    double startTime,
    double endTime,
    int indexPaymentType,
    double fwdPaymentPeriod,
    const ARM_GP_Vector& refIdxResetTimes,
    const ARM_GP_Vector& refIdxStartTimes,
    const ARM_GP_Vector& refIdxEndTimes,
    const ARM_GP_Vector& refFwdPeriods,
    const ARM_GP_Vector& refIndexWeight,
    double couponMargin,
    const vector<const ARM_GP_Vector*> downBarrierPerState,
    const vector<const ARM_GP_Vector*> upBarrierPerState,
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: VanillaSwaption
///	Returns: a vector of swaption prices
///	Action : Compute swaption price for the input model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMarginBasis::VanillaSwaption(
	const string& modelName,
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
	const ARM_GP_Vector& fixTimes,
    const ARM_GP_Vector& fixPayPeriods,
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
///	Class   : ARM_ForwardMarginBasis
///	Routine: VanillaSpreadOptionLet
///	Returns : void
///	Action  : Computes the VanillaSpreadOptionLet
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_ForwardMarginBasis::VanillaSpreadOptionLet(const string& modelName,
														double evalTime,
														int callPut,
														double startTime,
														double endTime,
														double resetTime,
														double payTime,
														double payPeriod,
														double notional,
														double coeffLong,
														double coeffShort,
														const ARM_GP_Vector& strikes,
														double swapLongFloatStartTime,
														double swapLongFloatEndTime,
														const ARM_GP_Vector& swapLongFixPayTimes,
														const ARM_GP_Vector& swapLongFixPayPeriods,
														double swapShortFloatStartTime,
														double swapShortFloatEndTime,
														const ARM_GP_Vector& swapShortFixPayTimes,
														const ARM_GP_Vector& swapShortFixPayPeriods,
														const ARM_PricingStatesPtr& states) const
{
	CHANGE_FUNCTORS();
	CAST_TO_IRMODEL();
	ARM_VectorPtr result = IRMODEL->VanillaSpreadOptionLet( modelName, evalTime, callPut, startTime, endTime,resetTime,
		payTime, payPeriod, notional,coeffLong,coeffShort, strikes,swapLongFloatStartTime, swapLongFloatEndTime,
		swapLongFixPayTimes, swapLongFixPayPeriods, swapShortFloatStartTime, swapShortFloatEndTime, swapShortFixPayTimes,
		swapShortFixPayPeriods,states);
	
	RESTORE_FUNCTORS();
	return result;

}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMarginBasis
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardMarginBasis::LocalDiscounts( size_t timeIdx, 
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
///	Class  : ARM_ForwardMarginBasis
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ForwardMarginBasis::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_ForwardMarginBasis Model\n";
    os << indent << "----------------------------\n";

    os << "Model name : " << GetModelName() << "\n";

	os << ARM_ForwardMargin::toString(indent,nextIndent);

	return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

