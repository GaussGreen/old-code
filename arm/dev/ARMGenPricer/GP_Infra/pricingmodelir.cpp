/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodelir.cpp
 *  \brief interest rate version of a pricing model!
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


///gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/surfacemodelparam.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelIR
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelIR::ARM_PricingModelIR( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params, ARM_DensityFunctor* densityFct)
:	ARM_PricingModel(zc, params, densityFct)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelIR
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelIR::ARM_PricingModelIR(const ARM_PricingModelIR& rhs)
:	ARM_PricingModel(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelIR
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PricingModelIR::~ARM_PricingModelIR()
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelIR
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PricingModelIR& ARM_PricingModelIR::operator=(const ARM_PricingModelIR& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModel::operator=(rhs);
        // Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Default Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{    
    /// Get Libor Zc through the fixing functor
#if defined( __GP_STRICT_VALIDATION)
	if( !GetFixingFunctor() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": FixingFunctor is null!" );
#endif

    return DefaultLibor(curveName,evalTime,fwdStartTime,fwdEndTime,period,fwdResetTime,payTime,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: LiborWithControlVariate
///	Returns : a vector of libor values
///	Action  : Default Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::LiborWithControlVariate( 
	const string& curveName,
    double evalTime,
	double fwdStartTime,
    double fwdEndTime,
	double period,
    double fwdResetTime,    // for convexity adjustment...
    double payTime,         //... in derived classes
    const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr libor = Libor(	curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
	return libor;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: Annuity
///	Returns : a vector of annuity
///	Action  : Default Annuity Computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::Annuity(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_PricingStatesPtr& states) const
{
     return DefaultAnnuity(curveName,evalTime,fixPayTimes,fixPayPeriods,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: Annuity
///	Returns : a vector of annuity with Nominal
///	Action  : Compute the Annuity with Nominal for Variabl Notional Swap
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::AnnuityWithNominal(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fixNominal,
		const ARM_PricingStatesPtr& states) const
{
     return DefaultAnnuityWithNominal(curveName,evalTime,fixPayTimes,fixPayPeriods,fixNominal,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: SwapRateInPlaceWithComputedAnnuity
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
///				WARNING: need to clone the annuity as the computation
///				is done in place... I REPEAT
///				the annuity argument will be modified by the computation as it is in place
///				so to avoid side effect on the annuity (or if you want to keep the value of the annuity)
///				one needs to first clone the value of the annuity
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::SwapRateInPlaceWithComputedAnnuity(
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
		const ARM_VectorPtr& FixedComputedAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
		const ARM_PricingStatesPtr& states ) const
{
	ARM_GP_VectorPtr FloatComputedAnnuity = ARM_GP_VectorPtr( NULL );
    if( isDbleNotional && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
		FloatComputedAnnuity = Annuity( curveName, evalTime, floatPayTimes, floatPayPeriods, states);

	return DefaultSwapRateInPlaceWithComputedAnnuity( curveName, evalTime, floatStartTime, floatEndTime, 
		fixPayTimes, fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods,
		margin, isDbleNotional, FixedComputedAnnuity, FloatComputedAnnuity, states);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: SwapRateInPlaceWithComputedAnnuity
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
///				WARNING: need to clone the annuity as the computation
///				is done in place... I REPEAT
///				the annuity argument will be modified by the computation as it is in place
///				so to avoid side effect on the annuity (or if you want to keep the value of the annuity)
///				one needs to first clone the value of the annuity
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::SwapRateInPlaceWithComputedAnnuityAndNominal(
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
		const ARM_VectorPtr& FixedComputedAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
		const ARM_GP_Vector& floatNotional,
		const ARM_PricingStatesPtr& states ) const
{
	ARM_GP_VectorPtr FloatComputedAnnuity = ARM_GP_VectorPtr( NULL );
	return DefaultSwapRateInPlaceWithComputedAnnuityAndNominal( curveName, evalTime, floatStartTime, floatEndTime, 
		fixPayTimes, fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods,
		margin, FixedComputedAnnuity, floatNotional,states);
}




////////////////////////////////////////////////////
///	Class  : ARM_PricingModelIR
///	Routine: NPVSwap
///	Returns: a vector of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwap computation
///           using double notional method
/// 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelIR::NPVSwap(
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
		const ARM_GP_Vector& FixNotional, 
		const ARM_GP_Vector& FloatNotional, 
		const ARM_GP_Matrix& strikesPerState,
		int payRec,
		const ARM_PricingStatesPtr& states) const
{
	/// Compute the fixed leg price
	//// If Strikes  are not constant we calculate the fixed Leg without Annuity
	//// ARM_GP_VectorPtr fixAnnuity = Annuity( curveName, evalTime, fixPayTimes, fixPayPeriods, states);
	
	///We calculate the margin at the same time as the float leg 
	ARM_GP_VectorPtr floatAnnuity = ARM_GP_VectorPtr (NULL);
    if( isDbleNotional && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) && margin[0] != 0.0)
		floatAnnuity = Annuity( curveName, evalTime, floatPayTimes, floatPayPeriods, states);
    
	return DefaultNPVSwapWithComputedAnnuity( curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes,
		fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, 
		margin, isDbleNotional, FixNotional,FloatNotional, strikesPerState, payRec, floatAnnuity,states);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwalLeg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PricingModelIR::NPVSwapLeg(
		const string& curveName, 
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
	return DefaultNPVSwapLeg(curveName, evalTime,fwdStartTimes,fwdEndTimes, fwdPayPeriods, 
		 payTimes, payPeriods, margin, notional,  states); 
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: NPVFixLeg
///	Returns : m,atrix ptr
///	Action  : To calculate a fix leg
ARM_GP_MatrixPtr ARM_PricingModelIR::NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const
{
	return DefaultNPVFixLeg(curveName,evalTime,fixPayTimes,fixPayPeriods,
		FixNotional, strikesPerState,payRec,states);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
double ARM_PricingModelIR::ImpliedVol(const ARM_VanillaArg& arg) const
{
    return DefaultImpliedVol(arg);
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_PricingModelIR::GetType() const
{
	return MT_INTEREST_RATE_MODEL;
}



CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
