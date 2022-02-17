/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingfunctionlir.h
 *
 *  \brief
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPINFRA_PRICINGFUNCTIONIR_H
#define _INGPINFRA_PRICINGFUNCTIONIR_H

#include "gpbase/env.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported

/// Forward declarations
class ARM_ModelParam;
struct ARM_VanillaArg;
struct ARM_SwapRate;


///////////////////////////////////////////////////////
/// \class ARM_PricingFunctionIR
/// \brief
/// This abstract class is the interface for IR keyword
/// functions of the GP
///////////////////////////////////////////////////////
class ARM_PricingFunctionIR
{
public:
	/// -----------------------------------------------------
    /// ----------- abstract function -----------------
	/// -----------------------------------------------------
	/// Libor function
	virtual ARM_VectorPtr Libor( 
		const string& curveName, double evalTime,
		double fwdStartTime, double fwdEndTime,
		double period, double fwdResetTime, double payTime,
        const ARM_PricingStatesPtr& states) const = 0;

	/// Annuity function
	virtual ARM_VectorPtr Annuity(
		const string& curveName, 
        double evalTime,
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuity(
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
		const ARM_VectorPtr& PreComputedAnnuity,
		const ARM_PricingStatesPtr& states) const = 0;

	/// Swap function (vectorial fixed rate version)
	virtual ARM_VectorPtr NPVSwap(
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
		const ARM_GP_Matrix&  strikesPerState,
        int payRec,
        const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing the floating leg with different fixing & dicount factors
	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& PayTimes, 
		const ARM_GP_Vector& PayPeriods, 
		const ARM_GP_Vector& margin, 
		const ARM_GP_Vector& notional, 
		const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing the fixed leg with different fixing & dicount factors
	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const =0;


    /// Implied volatility function to calculate 
    /// sqrt[1/resetTime * integral(evalTime-> resettime
    /// sigme*sigme) du]
    virtual double ImpliedVol(const ARM_VanillaArg& arg) const = 0;

    /// Vanilla caplet function (vectorial strike version)
    virtual ARM_VectorPtr VanillaCaplet(
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
        const ARM_PricingStatesPtr& states) const = 0; 

    /// Vanilla swaption function (vectorial strike version)
    virtual ARM_VectorPtr VanillaSwaption(
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
        const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const = 0;

	virtual ARM_VectorPtr  VanillaSpreadOption(
		const string& curveName,
		double evalTime,
		int callPut,
		double startTime,
		double endTime,
		const ARM_GP_Vector& resetDates,
		const ARM_GP_Vector& payDates,
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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
			const string& curveName,
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
			const ARM_PricingStatesPtr& states) const = 0;

	/// -----------------------------------------------------
	/// --------- function with default behavior
	/// -----------------------------------------------------
    /// Swap rate function
	virtual ARM_VectorPtr SwapRate(
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
        const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Spread(
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
        const ARM_PricingStatesPtr& states) const;

	/// Swap function (scalar fixed rate version)
    ARM_VectorPtr NPVSwapScalar(
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
		double strikes,
        int payRec,
        const ARM_PricingStatesPtr& states) const;

    /// Vanilla caplet function (scalar strike version)
    ARM_VectorPtr VanillaCapletScalar(
		const string& curveName, 
		double evalTime,
		double payTime,
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        double strike,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;  


    /// Default vanilla digital caplet function (vectorial strike version)
    virtual ARM_VectorPtr VanillaDigital(
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
        const ARM_PricingStatesPtr& states) const;  
    
    /// Vanilla digital caplet function (scalar strike version)
    ARM_VectorPtr VanillaDigitalScalar(
		const string& curveName, 
		double evalTime,
		double payTime,
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        double strike,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;  

    /// Vanilla corridorlet function (vectorial barriers version)
    virtual ARM_VectorPtr VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     indexPaymentType,
        double  fwdPaymentPeriod,
        const ARM_GP_Vector& refIdxResettimes,
        const ARM_GP_Vector& refIdxStarttimes,
        const ARM_GP_Vector& refIdxEndtimes,
        const ARM_GP_Vector& refFwdPeriods,
        const ARM_GP_Vector& refIndexWeight,
        double  couponMargin,
        const vector<const ARM_GP_Vector*> downBarrierPerState,
        const vector<const ARM_GP_Vector*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;   


    /// Vanilla corridorlet function (scalar barriers version)
    ARM_VectorPtr VanillaCorridorletScalar(
		const string& curveName, 
		double evalTime,
        double payTime,
        double ResetTime,
        double StartTime,
        double EndTime,
        int    indexPaymentType, 
        double fwdPaymentPeriod,
        const ARM_GP_Vector& refIdxResettimes,
        const ARM_GP_Vector& refIdxStarttimes,
        const ARM_GP_Vector& refIdxEndtimes,
        const ARM_GP_Vector& refFwdPeriods,
        const ARM_GP_Vector& refIndexWeight,
        double couponMargin,
        const ARM_GP_Vector& downBarrier,
        const ARM_GP_Vector& upBarrier,
        double payNotional,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;

	/// Vanilla corridor
    ARM_VectorPtr VanillaCorridor(
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
        const   ARM_PricingStatesPtr& states) const;

	// Vanilla Corridor function
	ARM_VectorPtr VanillaCMSCorridor(
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
		const ARM_SwapRatePtrVectorVector& thirdIndexes, // 3rd index for double condition
		const ARM_VectorVector& downBarriers3,
		const ARM_VectorVector& upBarriers3,
        const   ARM_PricingStatesPtr& states) const;

	// Vanilla Corridorlet function
	virtual ARM_VectorPtr VanillaCMSCorridorlet(
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
        const   ARM_PricingStatesPtr& states) const;

    /// Vanilla cap function (vectorial strike version)
	ARM_VectorPtr VanillaCap(
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
		const ARM_PricingStatesPtr& states) const;

    /// Vanilla cap function (scalar strike version)
    ARM_VectorPtr VanillaCapScalar(
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
		const ARM_PricingStatesPtr& states) const;


    /// Vanilla swaption function (scalar strike version)
    ARM_VectorPtr VanillaSwaptionScalar(
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
        const ARM_GP_Vector& strikes,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	ARM_VectorPtr VanillaSpreadOptionScalar(
		const string& curveName,
		double evalTime,
		int callPut,
		double startTime,
		double endTime,
		const ARM_GP_Vector& resetDates,
		const ARM_GP_Vector& payDates,
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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSumOption(
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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSmiledSwaption(
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
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	virtual ARM_VectorPtr MaxRate(
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
        const ARM_PricingStatesPtr& states) const;

    /// Double digital (expiry = eval date)
    virtual ARM_VectorPtr DoubleDigital(
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
        const ARM_PricingStatesPtr& states) const;
    
	virtual ARM_VectorPtr ImpliedVol(
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
        const ARM_PricingStatesPtr& states) const;
	
	
	//// Flag. True if the model has a closed formula for swaaptions
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantNominal&&isConstantSpread&&isConstantstrike;}


	////////////////////////////////////////////////////////////////////////////////////////////////
	/// Datas that any irmodel must provide to work in a multiasset environment when used with infaltion
	////////////////////////////////////////////////////////////////////////////////////////////////

	/// This is supposed to provide Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
	/// Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
	virtual double IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const { 
			        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "IntegratedBondSquaredVol not implemented in this model. Have a nice day. " ); }

	/// Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
	virtual double IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const { 
			        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "IntegratedBondCovariance not implemented in this model. Have a nice day. " ); }

	/// Scalar product: returns Int_startTime^endTime Gamma(s,bondMaturity) * otherModelVolatility(s) ds
	virtual double VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const { 
			        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "VolatilityScalarProduct not implemented in this model. Have a nice day. " ); }
};



CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
