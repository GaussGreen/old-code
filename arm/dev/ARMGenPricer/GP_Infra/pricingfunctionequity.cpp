/*!
 *
 * Copyright (c) CDC IXIS CM June 2004 Paris
 *
 *	\file pricingfunctionequity.cpp
 *  \brief
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#include "gpinfra/pricingmodelequity.h"


#include "gpbase/datestrip.h"

#include "gpinfra/pricingfunctionir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingcontext.h"

/// kernel
#include "crv/zerocurv.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: GetSettlementCalendar
///	Returns: string
///	Action : default implementation (static routine)
////////////////////////////////////////////////////
string ARM_PricingFunctionEquity::GetSettlementCalendar(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve)
{
	char setCal[7];
	strcpy(setCal, domCurve->GetCurrencyUnit()->GetCcyName());
	if(forCurve)
		strcat(setCal, forCurve->GetCurrencyUnit()->GetCcyName());
	return string(setCal);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: GetSettlementGap
///	Returns: double
///	Action : default implementation (static routine)
////////////////////////////////////////////////////
double ARM_PricingFunctionEquity::GetSettlementGap(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve)
{
	double domGap = domCurve->GetCurrencyUnit()->GetSpotDays();
	double forGap=domGap;
	if(forCurve)
		forGap = forCurve->GetCurrencyUnit()->GetSpotDays();
	return domGap > forGap ? domGap : forGap;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: CallScalar
///	Returns: ARM_VectorPtr
///	Action : compute a call or put with a scalar strike
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::CallScalar(
	const string& curveName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	double strike,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	size_t stateSize = states ==ARM_PricingStatesPtr(NULL)? 1 : states->size();
	ARM_GP_Vector strikePerState( stateSize, strike);
	return CallVectorial(curveName,evalTime,expiryTime,settlementTime,strikePerState,callPut,payTime,states,context);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: DigitalVectorial
///	Returns: ARM_VectorPtr
///	Action : compute a digital call or put with a vectorial strike
///          Default implementation with numerical derivate
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::DigitalVectorial(
	const string& curveName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikesPerState,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double eps,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	const double epsilon =1e-10;
	const double invTwoEpsilon = 1.0/(2.0*epsilon);
	
	/// computes in place the digital as (call(K-epsilon)-Call(K+epsilon))/(2.*epsilnon)
	ARM_GP_Vector upStrike   = strikesPerState;
	upStrike			    += epsilon;
	ARM_GP_Vector downStrike = strikesPerState;
	downStrike			 -= epsilon;

    ARM_VectorPtr callup,digital;
    if(context)
    {
        ARM_PricingDigitalContext* digitalContext = context->ToDigitalContext();
        ARM_PricingCallContext callContext;
	    digital = CallVectorial(curveName,evalTime,expiryTime,settlementTime,downStrike,callPut,payTime,states,&callContext);
        digitalContext->InsertDatas(callContext.GetVol(),callContext.GetShift(),callContext.GetForward());
	    callup  = CallVectorial(curveName,evalTime,expiryTime,settlementTime,upStrike,callPut,payTime,states,&callContext);
        digitalContext->InsertDatas(callContext.GetVol(),callContext.GetShift(),callContext.GetForward());
    }
    else
    {
	    digital = CallVectorial(curveName,evalTime,expiryTime,settlementTime,downStrike,callPut,payTime,states);
	    callup  = CallVectorial(curveName,evalTime,expiryTime,settlementTime,upStrike,callPut,payTime,states);
    }

	size_t stateSize = states == ARM_PricingStatesPtr(NULL)? 1 : states->size();
	for( size_t i=0; i<stateSize; ++i )
		(*digital)[i] =( (*digital)[i]-(*callup )[i] ) * invTwoEpsilon;
	return digital;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: DigitalScalar
///	Returns: ARM_VectorPtr
///	Action : compute a digital call or put with a scalar strike
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::DigitalScalar(
	const string& curveName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	double strike,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double eps,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	size_t stateSize = states ==ARM_PricingStatesPtr(NULL)? 1 : states->size();
	ARM_GP_Vector strikePerState( stateSize, strike);
	return DigitalVectorial(curveName,evalTime,expiryTime,settlementTime,strikePerState,notional,callPut,payTime,digitType,eps,states,context);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: DigitalScalar
///	Returns: ARM_VectorPtr
///	Action : compute a digital call or put with a scalar strike
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::CallStripVectorial(
		const string& modelName,
        double evalTime,
	    const ARM_GP_Vector& expiryTime,
	    const ARM_GP_Vector& settlementTime,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
	    const ARM_GP_Vector& payTime,
	    const ARM_GP_Vector& nominal,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{
	size_t i,nbStates = states ==ARM_PricingStatesPtr(NULL)? 1 : states->size();
    size_t j,nbCalls = expiryTime.size(); 
	ARM_GP_Vector strikes(nbStates);

    ARM_VectorPtr callStripValues(new ARM_GP_Vector(nbStates,0.0));
    ARM_VectorPtr callValues;
    double nomi;
    for(j=0;j<nbCalls;++j)
    {
        for(i=0;i<nbStates;++i)
            strikes[i] = strikesPerState(i,j);

	    callValues = CallVectorial(modelName,evalTime,expiryTime[j],settlementTime[j],strikes,callPut,payTime[j],states,context);

        nomi = nominal[j];
        for(i=0;i<nbStates;++i)
            (*callStripValues)[i] += nomi * (*callValues)[i];
    }

    return callStripValues;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: DigitalScalar
///	Returns: ARM_VectorPtr
///	Action : compute a digital call or put with a scalar strike
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::CallStripScalar(
		const string& modelName,
        double evalTime,
	    const ARM_GP_Vector& expiryTime,
	    const ARM_GP_Vector& settlementTime,
		const ARM_GP_Vector& strike,
        int callPut,
	    const ARM_GP_Vector& payTime,
	    const ARM_GP_Vector& nominal,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{
	size_t i,nbStates = states ==ARM_PricingStatesPtr(NULL)? 1 : states->size();
    size_t j,nbCalls = expiryTime.size(); 
	ARM_GP_Matrix strikesPerState( nbStates, nbCalls);
    double k;
    for(j=0;j<nbCalls;++j)
    {
        k=strike[j];
        for(i=0;i<nbStates;++i)
            strikesPerState(i,j) = k;
    }

	return CallStripVectorial(modelName,evalTime,expiryTime,settlementTime,strikesPerState,callPut,payTime,nominal,states,context);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: HybridCallVectorial
///	Returns: price an hybrid call that is a spread
///			 between a forward equity or FX strip
///			 and an IR swap
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::HybridCallVectorial(
	const string& modelName,
	double evalTime,
	double expiryTime,
	int callPut,
	const ARM_GP_Vector& strikesPerState,

	/// Strip of forwards FX (or equity)
	const ARM_GP_Vector& fxExpiryTimes,
	const ARM_GP_Vector& fxSettlementTimes,
	const ARM_GP_Vector& fxPayTimes,
	const ARM_GP_Vector& fxNotionals,

	/// IR Swap
	double swapResetTime,
	const ARM_GP_Vector& fixNotionals,
	const ARM_GP_Vector& floatNotionals,
	double floatStartTime,
	double floatEndTime,
	const ARM_GP_Vector& floatResetTimes,
	const ARM_GP_Vector& floatStartTimes,
	const ARM_GP_Vector& floatEndTimes,
	const ARM_GP_Vector& floatIntTerms,
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,

	const ARM_PricingStatesPtr& states,
	ARM_PricingContext* context) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ": Hybrid call not available with this model" ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: HybridCallScalar
///	Returns: vector ptr
///	Action : scalar version for hybrid call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::HybridCallScalar(
		const string& modelName,
		double evalTime,
		double expiryTime,
		int callPut,

		/// Strip of forwards FX (or equity)
		const ARM_GP_Vector& fxExpiryTimes,
		const ARM_GP_Vector& fxSettlementTimes,
		const ARM_GP_Vector& fxPayTimes,
		const ARM_GP_Vector& fxNotionals,

		/// Payer Swap
		double swapResetTime,
		const ARM_GP_Vector& fixNotionals,
		const ARM_GP_Vector& floatNotionals,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double strike,

		const ARM_PricingStatesPtr& states,
		ARM_PricingContext* context) const
{
	size_t nbStates = (states == ARM_PricingStatesPtr(NULL) ? 1 : states->size());
	ARM_GP_Vector strikesPerState(nbStates,strike);

	return HybridCallVectorial(modelName,evalTime,expiryTime,callPut,strikesPerState,
		fxExpiryTimes,fxSettlementTimes,fxPayTimes,fxNotionals,
		swapResetTime,fixNotionals,floatNotionals,floatStartTime,floatEndTime,
		floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
		fixPayTimes,fixPayPeriods,
		states,context);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: RangeAccrualVectorial
///	Returns: vector ptr
///	Action : compute a range accrual dual index : IR index and FX index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingFunctionEquity::RangeAccrualVectorial(
		const string& curveName,
		double evalTime,
		double startTime,
		double endTime,
		double payTime,
		const  ARM_GP_Vector& fixingTimes,
		int    payIndexType, 
        double payIndexTerm,
		const  string& fxModelName,
		int    irIndexType, 
		const  ARM_GP_Vector& irIndexResetTimes,
		const  ARM_GP_Vector& irIndexStartTimes,
		const  ARM_GP_Vector& irIndexEndTimes,
		const  ARM_GP_Vector& irIndexTerms,
		const  ARM_GP_Vector& fxDownBarriers,
		const  ARM_GP_Vector& fxUpBarriers,
		const  ARM_GP_Vector& irDownBarriers,
		const  ARM_GP_Vector& irUpBarriers,
		const  ARM_GP_Vector& notionals,
		const  ARM_PricingStatesPtr& states,
		ARM_Vector* eachFixingPrices,
        ARM_PricingContext* context) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ": Range Accrual not available with this model" ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: ValidateIRModel (static)
///	Returns: void
///	Action : check if an input model is an IR model
////////////////////////////////////////////////////
void ARM_PricingFunctionEquity::ValidateIRModel( const ARM_PricingModelPtr& irModel )
{
	if( !dynamic_cast<ARM_PricingFunctionIR*>( &*irModel ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ": this is not an interet rate model!" ); 
}

/***
////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionEquity
///	Routine: SetIRModel
///	Returns: void
///	Action : method to set an interest rate model
////////////////////////////////////////////////////
void ARM_PricingFunctionEquity::SetIRModel( const ARM_PricingModelPtr& irModel )
{ 
	ARM_THROW( ERR_INVALID_ARGUMENT, ": unimplemented method 'SetIRModel'" ); 
}
***/

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
