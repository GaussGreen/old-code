/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingfunctionir.cpp
 *  \brief
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/curve.h"

#include "gpinfra/pricingfunctionir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/irrate.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: NPVSwapScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
/// 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::NPVSwapScalar(
		const string& curveName, 
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
		double margin,
		bool isDbleNotional,
		double Notional,
		double strike,
		int payRec,
		const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
	size_t nbFix = fixPayTimes.size();
	ARM_GP_Matrix strikeMatrix (statesSize,nbFix,strike);
	ARM_GP_Vector FixNotionalVector(fixPayPeriods.size(),Notional);
	ARM_GP_Vector FloatNotionalVector(floatPayPeriods.size(),Notional);
	ARM_GP_Vector MarginVector(floatPayPeriods.size(),margin);
	return NPVSwap(
		curveName, 
		evalTime,
		floatStartTime,
		floatEndTime, 
		fixPayTimes,
		fixPayPeriods,
		fwdStartTimes, 
		fwdEndTimes, 
		fwdPayPeriods, 
		floatPayTimes, 
		floatPayPeriods, 
		MarginVector,
		isDbleNotional,
		FixNotionalVector,
		FloatNotionalVector,
		strikeMatrix,
		payRec,
		states);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCapletScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCapletScalar(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        double strike,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
	ARM_GP_Vector strikeVector( statesSize, strike );
	
	return VanillaCaplet(curveName, evalTime, payTime, period, payNotional, fwdResetTime, fwdStartTime,
						fwdEndTime, fwdPeriod, strikeVector, capFloor, states );

}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Default for VanillaDigital is to have
///          numerical derivatives from a vanilla caplet
///          not the best but provides something in any case!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaDigital(
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
	const double epsilon  = 0.0001;
	ARM_GP_Vector upStrike	  = strikesPerState;
	upStrike			 += epsilon;
	ARM_VectorPtr upPrice = VanillaCaplet( curveName, 
		evalTime,
		payTime, 
		period,
        payNotional,
		fwdResetTime, 
		fwdStartTime,
        fwdEndTime,
        fwdPeriod,
        upStrike,
        capFloor,
		states);

	ARM_GP_Vector downStrike   = strikesPerState;
	downStrike			   -= epsilon;
	ARM_VectorPtr downPrice	=  VanillaCaplet( curveName, 
		evalTime,
		payTime, 
		period,
        payNotional,
		fwdResetTime, 
		fwdStartTime,
        fwdEndTime,
        fwdPeriod,
        downStrike,
        capFloor,
		states);


	/// does the minus in place
	ARM_GP_Vector::iterator 
		lhsBegin= upPrice->begin(),
		rhsBegin= downPrice->begin(),
		lhsEnd	= upPrice->end();

	const double invTwoEpsilon  = 5000.;
	for( ; lhsBegin!= lhsEnd; ++rhsBegin,++lhsBegin )
		*lhsBegin = ((*rhsBegin)-(*lhsBegin))*invTwoEpsilon;
	return upPrice;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaDigitalScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaDigitalScalar(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
        double fwdPeriod,
        double strike,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
	ARM_GP_Vector strikeVector( statesSize, strike );
	return VanillaDigital(curveName, evalTime,payTime, period, payNotional, fwdResetTime, 
		fwdStartTime, fwdEndTime, fwdPeriod, strikeVector, capFloor, states );
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCorridorletScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCorridorletScalar(
		const string& curveName, 
		double evalTime,
        double payTime,
        double ResetTime,
        double StartTime,
        double EndTime,
        int    indexPaymentType, 
        double fwdPaymentPeriod,
        const ARM_GP_Vector& RefIdxResettimes,
        const ARM_GP_Vector& RefIdxStarttimes,
        const ARM_GP_Vector& RefIdxEndtimes,
        const ARM_GP_Vector& RefFwdPeriods,
        const ARM_GP_Vector& RefIndexWeight,
        double CouponMargin,
        const ARM_GP_Vector& DownBarrier,
        const ARM_GP_Vector& UpBarrier,
        double payNotional,
        int capFloor,
        const ARM_PricingStatesPtr& states) const
{

    size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
	vector<const ARM_GP_Vector* > DownBarrierVector( statesSize, &( const_cast<ARM_GP_Vector&>(DownBarrier) ));
    vector<const ARM_GP_Vector*> UpBarrierVector( statesSize, &( const_cast<ARM_GP_Vector&>(UpBarrier) ));

    return VanillaCorridorlet(curveName, evalTime,payTime,ResetTime,StartTime,EndTime,indexPaymentType,fwdPaymentPeriod,
        RefIdxResettimes,RefIdxStarttimes,RefIdxEndtimes,RefFwdPeriods,RefIndexWeight,CouponMargin,
        DownBarrierVector,UpBarrierVector,payNotional,capFloor,states);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCorridor
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCorridor(
	const   string& curveName, 
	double  evalTime,
    const ARM_GP_Vector&  payTimes,
    const ARM_GP_Vector&  resetTimes,
    const ARM_GP_Vector&  startTimes,
    const ARM_GP_Vector&  endTimes,
    int     indexPaymentType,
    const ARM_GP_Vector&  fwdPaymentPeriods,
    const ARM_VectorVector& refIdxResettimes,
    const ARM_VectorVector& refIdxStarttimes,
    const ARM_VectorVector& refIdxEndtimes,
    const ARM_VectorVector& refFwdPeriods,
    const ARM_VectorVector& refIndexWeights,
    const ARM_GP_Vector&  couponMargins,
    const ARM_VectorVector& downBarriers,
    const ARM_VectorVector& upBarriers,
    const ARM_GP_Vector&  payNotionals,
    int     capFloor,
    const   ARM_PricingStatesPtr& states) const
{
	int size = states->size();
	ARM_GP_Vector* null = new ARM_GP_Vector(size,0.0); 
	ARM_VectorPtr result(null);
	ARM_GP_Vector::iterator
	resultEnd	= result->end(), corridorletBegin, resultBegin;

	size_t nbFlows = startTimes.size();

	for (size_t i = 0; i < nbFlows; ++i)
	{
		ARM_GP_VectorPtr corridorlet_i = VanillaCorridorletScalar(
			curveName,
			evalTime,
			payTimes[i],
			resetTimes[i],
			startTimes[i],
			endTimes[i],
			indexPaymentType,
			fwdPaymentPeriods[i],
			*refIdxResettimes[i],
			*refIdxStarttimes[i],
			*refIdxEndtimes[i],
			*refFwdPeriods[i],
			*refIndexWeights[i],
			couponMargins[i],
			*downBarriers[i],
			*upBarriers[i],
			payNotionals[i],
			capFloor,
			states);

		corridorletBegin = corridorlet_i->begin();
		resultBegin= result->begin();

		for( ; resultBegin!= resultEnd; ++corridorletBegin, ++resultBegin )
		*resultBegin = (*resultBegin)+(*corridorletBegin);
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCorridor
///	Returns: a vector
///	Action : indirects to VanillaCMSCorridor
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCMSCorridor(
	const string& curveName,
	double evalTime,
	const ARM_GP_Vector& payTimes,
	const ARM_GP_Vector& resetTimes,
	const ARM_GP_Vector& startTimes,
	const ARM_GP_Vector& endTimes,
	const ARM_VectorVector& refIdxResettimes,
	const ARM_VectorVector& refIndexWeights,
	const ARM_VectorVector& coeffs1,
	const ARM_SwapRatePtrVectorVector& firstIndexes,
	const ARM_VectorVector& coeffs2,
	const ARM_SwapRatePtrVectorVector& secondIndexes,
	const ARM_IntVector& payIndexType,		/// K_FIXED, K_CMS, K_LIBOR
	const ARM_GP_Vector& coupon,			/// in case of a fixed payment (K_FIXED)
	const ARM_SwapRatePtrVector& payRate,	/// rate description (if K_CMS or K_LIBOR)
	const ARM_GP_Vector& payIndexLeverage,
	const ARM_VectorVector& downBarriers,
    const ARM_VectorVector& upBarriers,
    const ARM_GP_Vector&  payNotionals,
    int     rcvPay,
	const ARM_SwapRatePtrVectorVector& thirdIndexes, // Single rate index for double condition
	const ARM_VectorVector& downBarriers3,
	const ARM_VectorVector& upBarriers3,
    const   ARM_PricingStatesPtr& states) const
{
	int size = states->size();
	ARM_GP_Vector* null = new ARM_GP_Vector(size,0.0); 
	ARM_VectorPtr result(null);
	ARM_GP_Vector::iterator
	resultEnd	= result->end(), corridorletBegin, resultBegin;

	size_t nbFlows = startTimes.size();

	for (size_t i = 0; i < nbFlows; ++i)
	{
		ARM_GP_VectorPtr corridorlet_i = VanillaCMSCorridorlet(
			curveName,
			evalTime,
			payTimes[i],
			resetTimes[i],
			startTimes[i],
			endTimes[i],
			*refIdxResettimes[i],
			*refIndexWeights[i],
			*coeffs1[i],
			firstIndexes[i],
			*coeffs2[i],
			secondIndexes[i],
			payIndexType[i],
			coupon[i],
			*payRate[i],
			payIndexLeverage[i],
			*downBarriers[i],
			*upBarriers[i],
			payNotionals[i],
			rcvPay,
			thirdIndexes[i],
			*downBarriers3[i],
			*upBarriers3[i],
			states);

		corridorletBegin = corridorlet_i->begin();
		resultBegin= result->begin();

		for( ; resultBegin!= resultEnd; ++corridorletBegin, ++resultBegin )
		*resultBegin = (*resultBegin)+(*corridorletBegin);
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCorridorlet
///	Returns: a vector
///	Action : Closed form function for CMS spread corridorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCMSCorridorlet(
	const string& curveName,
	double evalTime,
	double payTime,
	double resetTime,
	double startTime,
	double endTime,
	const ARM_GP_Vector& refIdxResettimes,
	const ARM_GP_Vector& refIndexWeights,
	const ARM_GP_Vector& coeff1,
	const ARM_SwapRatePtrVector& firstIndex,
	const ARM_GP_Vector& coeff2,
	const ARM_SwapRatePtrVector& secondIndex,
	int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
	double	coupon,					/// in case of a fixed payment (K_FIXED)
	const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
	double  payIndexLeverage,
	const ARM_GP_Vector& downBarriers,
    const ARM_GP_Vector& upBarriers,
    double  payNotional,
    int     rcvPay,
	const ARM_SwapRatePtrVector& thirdIndex, // 3rd index for double condition
	const ARM_GP_Vector& downBarriers3,
    const ARM_GP_Vector& upBarriers3,
    const   ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented yet!" );
	return ARM_VectorPtr(NULL); 
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCap
///	Returns: a vector
///	Action : sums its caplets
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCap(
	const string& curveName,
	double evalTime,
	const ARM_GP_Vector& payTimes,
	const ARM_GP_Vector& periods,
	double notional,
	const ARM_GP_Vector& fwdResetTimes,
	const ARM_GP_Vector& fwdStartTimes,
	const ARM_GP_Vector& fwdEndTimes,
	const ARM_GP_Vector& fwdPeriods,
	const ARM_GP_Vector& strikesPerState,
	int   capFloor,    
	const ARM_PricingStatesPtr& states) const
{
	
	int size = states->size();
	ARM_GP_Vector* null = new ARM_GP_Vector(size,0.0); 
	ARM_VectorPtr result(null);
	ARM_GP_Vector::iterator
	resultEnd	= result->end(), capletBegin, resultBegin;
	
	for(int i=0; i<fwdPeriods.size(); ++i )
	{
		ARM_VectorPtr caplet_i = VanillaCaplet( curveName, 
		evalTime,
		payTimes[i], 
		periods[i],
		notional,
		fwdResetTimes[i], 
		fwdStartTimes[i],
		fwdEndTimes[i],
		fwdPeriods[i],
		strikesPerState,
		capFloor,
		states);
		
		capletBegin = caplet_i->begin();
		resultBegin= result->begin();

		for( ; resultBegin!= resultEnd; ++capletBegin, ++resultBegin )
		*resultBegin = (*resultBegin)+(*capletBegin);
	}

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCapScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaCapScalar(
	const string& curveName,
	double evalTime,
	const ARM_GP_Vector& payTimes,
	const ARM_GP_Vector& periods,
	double notional,
	const ARM_GP_Vector& fwdResetTimes,
	const ARM_GP_Vector& fwdStartTimes,
	const ARM_GP_Vector& fwdEndTimes,
	const ARM_GP_Vector& fwdPeriods,
	double strike,
	int capFloor,    
	const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
	ARM_GP_Vector strikeVector( statesSize, strike );
	return VanillaCap(curveName, evalTime, payTimes, periods, notional, fwdResetTimes, fwdStartTimes,
		fwdEndTimes, fwdPeriods,strikeVector, capFloor, states);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSwaptionScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaSwaptionScalar(
	const string& curveName,
	double evalTime,
	double swapResetTime,
    const ARM_GP_Vector& swapFixNotional,
	const ARM_GP_Vector& swapFloatNotional,
	double floatStartTime,
	double floatEndTime,
	const ARM_GP_Vector& floatResetTimes,
	const ARM_GP_Vector& floatStartTimes,
	const ARM_GP_Vector& floatEndTimes,
	const ARM_GP_Vector& floatIntTerms,
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
    const ARM_GP_Vector& strikes,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
    ARM_GP_Matrix strikeMatrix;
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;
    for(int i=0; i<statesSize;++i)
        strikeMatrix.push_backRow(strikes);
    
	return VanillaSwaption( curveName, evalTime, swapResetTime,
		    swapFixNotional, swapFloatNotional, floatStartTime, floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms, fixPayTimes, 
			fixPayPeriods, strikeMatrix, callPut, states, isConstantNotional,isConstantSpread,isConstantStrike);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSpreadOptionScalar
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaSpreadOptionScalar(
	const string& curveName,
	double evalTime,
	int callPut,
	double startTime,
    double endTime,
	const ARM_GP_Vector& resetTimes,
	const ARM_GP_Vector& payTimes,
	const ARM_GP_Vector& payPeriods,
	const ARM_GP_Vector& notional,
	const ARM_GP_Vector& coeffLong,
	const ARM_GP_Vector& coeffShort,
	const ARM_GP_Vector& strikes,
	const ARM_GP_Vector& swapLongFloatStartTime,
	const ARM_GP_Vector& swapLongFloatEndTime,
	const ARM_VectorVector& swapLongFixPayTimes,
	const ARM_VectorVector& swapLongFixPayPeriods,
	const ARM_GP_Vector& swapShortFloatStartTime,
	const ARM_GP_Vector& swapShortFloatEndTime,
	const ARM_VectorVector& swapShortFixPayTimes,
	const ARM_VectorVector& swapShortFixPayPeriods,
	const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states != ARM_PricingStatesPtr(NULL)? states->size():1;

	if(statesSize!=1)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": only step-up deterministic strikes supported" );

	ARM_GP_Matrix strikeMatrix(1,strikes.size(),strikes);

	return VanillaSpreadOption(curveName,
								evalTime,
								callPut,
								startTime,
								endTime,
								resetTimes,
								payTimes,
								payPeriods,
								notional,
								coeffLong,
								coeffShort,
								strikeMatrix,
								swapLongFloatStartTime,
								swapLongFloatEndTime,
								swapLongFixPayTimes,
								swapLongFixPayPeriods,
								swapShortFloatStartTime,
								swapShortFloatEndTime,
								swapShortFixPayTimes,
								swapShortFixPayPeriods,
								0.0,
								states);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSpreadOption
///	Returns: a vector
///	Action : indirects to the vectorial version
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::VanillaSpreadOption(
	const string& curveName,
	double evalTime,
	int callPut,
	double startTime,
    double endTime,
	const ARM_GP_Vector& resetTimes,
	const ARM_GP_Vector& payTimes,
	const ARM_GP_Vector& payPeriods,
	const ARM_GP_Vector& notional,
	const ARM_GP_Vector& coeffLong,
	const ARM_GP_Vector& coeffShort,
	const ARM_GP_Matrix& strikes,
	const ARM_GP_Vector& swapLongFloatStartTime,
	const ARM_GP_Vector& swapLongFloatEndTime,
	const ARM_VectorVector& swapLongFixPayTimes,
	const ARM_VectorVector& swapLongFixPayPeriods,
	const ARM_GP_Vector& swapShortFloatStartTime,
	const ARM_GP_Vector& swapShortFloatEndTime,
	const ARM_VectorVector& swapShortFixPayTimes,
	const ARM_VectorVector& swapShortFixPayPeriods,
	double leveragePrev,
	const ARM_PricingStatesPtr& states) const
{
	size_t nbStates = states->size(),nbFlows = payPeriods.size();

	if(nbFlows != strikes.cols() || strikes.rows() < 1)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": wrong flow size or null state size for strikes" );

	ARM_VectorPtr result(new ARM_GP_Vector(nbStates,0.0));
	ARM_GP_Vector::iterator
	resultEnd	= result->end(), capletBegin, resultBegin;


	ARM_VectorPtr tmpVectorPtr;
	
	ARM_GP_Vector strikeVector(nbStates);
	for(int i=0; i<nbFlows; ++i )
	{
		for(int stateIdx=0; stateIdx<nbStates; ++stateIdx )
			strikeVector[stateIdx] = strikes(stateIdx<strikes.rows() ? stateIdx : strikes.rows()-1,i);

		ARM_VectorPtr spreadOptionLet_i = VanillaSpreadOptionLet( curveName, 
																evalTime,
																callPut,
																startTime,
																endTime,
																resetTimes[i],
																payTimes[i],
																payPeriods[i],
																notional[i],
																coeffLong[i],
																coeffShort[i],
																strikeVector,
																swapLongFloatStartTime[i],
																swapLongFloatEndTime[i],
																*swapLongFixPayTimes[i],
																*swapLongFixPayPeriods[i],
																swapShortFloatStartTime[i],
																swapShortFloatEndTime[i],
																*swapShortFixPayTimes[i],
																*swapShortFixPayPeriods[i],
																states);

		if (fabs(leveragePrev) > K_NEW_DOUBLE_TOL)
		{
			for (int j = 0; j < i; ++j)
			{
				tmpVectorPtr = VanillaSpreadOptionLet( curveName, 
																	evalTime,
																	callPut,
																	startTime,
																	endTime,
																	resetTimes[j],
																	payTimes[i],
																	payPeriods[j],
																	notional[j],
																	coeffLong[j],
																	coeffShort[j],
																	strikeVector,
																	swapLongFloatStartTime[j],
																	swapLongFloatEndTime[j],
																	*swapLongFixPayTimes[j],
																	*swapLongFixPayPeriods[j],
																	swapShortFloatStartTime[j],
																	swapShortFloatEndTime[j],
																	*swapShortFixPayTimes[j],
																	*swapShortFixPayPeriods[j],
																	states);

				(*spreadOptionLet_i) += (*tmpVectorPtr)*leveragePrev;
			}
		}
		
		capletBegin = spreadOptionLet_i->begin();
		resultBegin= result->begin();

		for( ; resultBegin!= resultEnd; ++capletBegin, ++resultBegin )
		*resultBegin = (*resultBegin)+(*capletBegin);
	}

	return result;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingFunctionIR
///	Routines: SwapRate
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::SwapRate(
		const string& curveName, 
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
	ARM_VectorPtr fixLegAnnuity = Annuity( curveName, evalTime, fixPayTimes, fixPayPeriods, states);
    return SwapRateInPlaceWithComputedAnnuity(	curveName,  evalTime,
		floatStartTime, floatEndTime, fixPayTimes,
		fixPayPeriods, fwdStartTimes, fwdEndTimes,
        fwdPayPeriods, floatPayTimes, floatPayPeriods,
        margin, isDbleNotional, fixLegAnnuity, states );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingFunctionIR
///	Routines: Spread
///	Returns : a vector of spread values
///	Action  : Default Spread computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::Spread(
		const string& curveName, 
		double evalTime,
		double coeff1,
		double floatStartTime1, 
		double floatEndTime1, 
		const ARM_GP_Vector& fixPayTimes1,
		const ARM_GP_Vector& fixPayPeriods1,
		const ARM_GP_Vector& fwdStartTimes1,
        const ARM_GP_Vector& fwdEndTimes1,
        const ARM_GP_Vector& fwdPayPeriods1,
		const ARM_GP_Vector& floatPayTimes1,
        const ARM_GP_Vector& floatPayPeriods1,
        const ARM_GP_Vector& margin1,
		double coeff2,
		double floatStartTime2, 
		double floatEndTime2, 
		const ARM_GP_Vector& fixPayTimes2,
		const ARM_GP_Vector& fixPayPeriods2,
		const ARM_GP_Vector& fwdStartTimes2,
        const ARM_GP_Vector& fwdEndTimes2,
        const ARM_GP_Vector& fwdPayPeriods2,
		const ARM_GP_Vector& floatPayTimes2,
        const ARM_GP_Vector& floatPayPeriods2,
        const ARM_GP_Vector& margin2,
		const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr fixLegAnnuity1 = Annuity( curveName, evalTime, fixPayTimes1, fixPayPeriods1, states);

	ARM_VectorPtr fixLegAnnuity2 = Annuity( curveName, evalTime, fixPayTimes2, fixPayPeriods2, states);

	ARM_VectorPtr swapRate1 = SwapRateInPlaceWithComputedAnnuity(	curveName,  evalTime,
		floatStartTime1, floatEndTime1, fixPayTimes1,
		fixPayPeriods1, fwdStartTimes1, fwdEndTimes1,
        fwdPayPeriods1, floatPayTimes1, floatPayPeriods1,
        margin1, true, fixLegAnnuity1, states );

	ARM_VectorPtr swapRate2 = SwapRateInPlaceWithComputedAnnuity(	curveName,  evalTime,
		floatStartTime2, floatEndTime2, fixPayTimes2,
		fixPayPeriods2, fwdStartTimes2, fwdEndTimes2,
        fwdPayPeriods2, floatPayTimes2, floatPayPeriods2,
        margin2, true, fixLegAnnuity2, states );

	ARM_GP_Vector::iterator iter1 = swapRate1->begin();
	ARM_GP_Vector::iterator iter2 = swapRate2->begin();

	for(; iter1!=swapRate1->end() ; ++iter1, ++iter2)
		(*iter1) = (*iter1)*coeff1 - (*iter2)*coeff2;
	
    return swapRate1;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaCorridorlet
///	Returns: a vector of Corridor let(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          Corridor caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////

ARM_VectorPtr ARM_PricingFunctionIR::VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     fwdPaymentType, 
        double  fwdPaymentPeriod,
        const ARM_GP_Vector& RefIdxResettimes,
        const ARM_GP_Vector& RefIdxStarttimes,
        const ARM_GP_Vector& RefIdxEndtimes,
        const ARM_GP_Vector& RefFwdPeriods,
        const ARM_GP_Vector& RefIndexWeight,
        double  couponMargin,
        const vector<const ARM_GP_Vector*> downBarrierPerState,
        const vector<const ARM_GP_Vector*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const 
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented yet!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSumOption
///	Returns: a vector of sum option values
///	Action : Closed form formula for standard
///          Sum Option cap/floor 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_PricingFunctionIR::VanillaSumOption(
		const string& curveName,
		double evalTime,
		int capFloor,
		const ARM_GP_Vector& coeffs,
		const ARM_GP_Vector& fwdResetTimes,
		const ARM_GP_Vector& fwdStartTimes,
		const ARM_GP_Vector& fwdEndTimes,
		double payTime,
		const ARM_GP_Vector& fwdPeriods,
		const ARM_GP_Vector& strikesPerState,
		double volatilityRatio,
		double* sumFwd,
		double* sumVol,
		const ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented yet!" );
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSmiledSwaption
///	Returns: a vector of sum option values
///	Action : Closed form formula for Smiled Swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_PricingFunctionIR::VanillaSmiledSwaption(
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
		const ARM_GP_Vector& fixTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_GP_Vector& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		const ARM_GP_Vector& data,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const 
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unimplemented yet!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: MaxRate
///	Returns: Max/min values of the rate
///	Action : Compute the Max/min rate or option on
///			 Max-Min given the start end the end of the
///			 lookback period and values of the rate
///			 at the beginning. Here default implementation
///			 assuming no volatility
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::MaxRate(
		const string& curveName, 
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
        double firstResetTime,
        double firstStartTime,
		const ARM_VectorPtr& firstRate,
		int MaxOrMin,
		int ResetFreq,
		const ARM_VectorPtr& strikes,
		int CapOrFloor,
		double RhoMinMax,
		bool IsAccrued,
		double MinAccrued,
		double MaxAccrued,
        const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr result = SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, 
								fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, 
								true, states);

	if(firstRate->size() == 0)
	{
		return result;
	}

	if(firstRate->size() == 1)
	{
		for(int i = 0; i < result->size(); i++) (*result)[i] = MaxOrMin * ((*result)[i] - (*firstRate)[0]) > 0 ? (*result)[i] : (*firstRate)[0];

		return result;
	}

	if(firstRate->size() == result->size())
	{
		for(int i = 0; i < result->size(); i++) (*result)[i] = MaxOrMin * ((*result)[i] - (*firstRate)[i]) > 0 ? (*result)[i] : (*firstRate)[i];

		return result;
	}

	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": first rate and swap rate have not the same size" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: DoubleDigital
///	Returns: double digital condition values
///	Action : Implicitly option expiry = eval date
///			 Computes the double digital price on
///			 both rates given their values
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionIR::DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const ARM_GP_Vector& firstStrikeDown,
        const ARM_GP_Vector& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const ARM_GP_Vector& secondStrikeDown,
        const ARM_GP_Vector& secondStrikeUp,
		double secondStrikeSpread,
        const ARM_PricingStatesPtr& states) const
 {
	ARM_THROW(ERR_INVALID_ARGUMENT,"Double Digital not implemented for this model");
}

ARM_VectorPtr ARM_PricingFunctionIR::ImpliedVol(
		const string& curveName, 
		double evalTime,
		double payTime,
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        const ARM_GP_Vector& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const
{
	ARM_THROW(ERR_INVALID_ARGUMENT,"Implied Vol not implemented for this model");
}

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
