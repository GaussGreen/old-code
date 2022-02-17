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
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
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
		const ARM_VectorPtr& PreComputedAnnuity,
		const ARM_PricingStatesPtr& states) const = 0;

	/// Swap function (vectorial fixed rate version)
	virtual ARM_VectorPtr NPVSwap(
		const string& curveName, 
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
		const ARM_GP_Matrix&  strikesPerState,
        int payRec,
        const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing the floating leg with different fixing & dicount factors
	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& PayTimes, 
		const std::vector<double>& PayPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const = 0;

	/// function to avoid computing the fixed leg with different fixing & dicount factors
	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
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
        const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const = 0; 

    /// Vanilla swaption function (vectorial strike version)
    virtual ARM_VectorPtr VanillaSwaption(
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
		const std::vector<double>& fixTimes,
        const std::vector<double>& fixPayPeriods,
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
		const std::vector<double>& resetDates,
		const std::vector<double>& payDates,
		const std::vector<double>& payPeriods,
		const std::vector<double>& notional,
		const std::vector<double>& coeffLong,
		const std::vector<double>& coeffShort,
		const ARM_GP_Matrix& strikes,
		const std::vector<double>& swapLongFloatStartTime,
		const std::vector<double>& swapLongFloatEndTime,
		const ARM_VectorVector& swapLongFixPayTimes,
		const ARM_VectorVector& swapLongFixPayPeriods,
		const std::vector<double>& swapShortFloatStartTime,
		const std::vector<double>& swapShortFloatEndTime,
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
			const std::vector<double>& strikes,
			double swapLongFloatStartTime,
			double swapLongFloatEndTime,
			const std::vector<double>& swapLongFixPayTimes,
			const std::vector<double>& swapLongFixPayPeriods,
			double swapShortFloatStartTime,
			double swapShortFloatEndTime,
			const std::vector<double>& swapShortFixPayTimes,
			const std::vector<double>& swapShortFixPayPeriods,
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
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
        const std::vector<double>& margin,
        bool isDbleNotional,
        const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Spread(
		const string& curveName, 
        double evalTime,
		double coeff1,
		double floatStartTime1,
        double floatEndTime1, 
		const std::vector<double>& fixPayTimes1,
        const std::vector<double>& fixPayPeriods1,
		const std::vector<double>& fwdStartTimes1,
        const std::vector<double>& fwdEndTimes1,
        const std::vector<double>& fwdPayPeriods1, 
		const std::vector<double>& floatPayTimes1,
        const std::vector<double>& floatPayPeriods1,
        const std::vector<double>& margin1,
		double coeff2,
        double floatStartTime2,
        double floatEndTime2, 
		const std::vector<double>& fixPayTimes2,
        const std::vector<double>& fixPayPeriods2,
		const std::vector<double>& fwdStartTimes2,
        const std::vector<double>& fwdEndTimes2,
        const std::vector<double>& fwdPayPeriods2, 
		const std::vector<double>& floatPayTimes2,
        const std::vector<double>& floatPayPeriods2,
        const std::vector<double>& margin2,
        const ARM_PricingStatesPtr& states) const;

	/// Swap function (scalar fixed rate version)
    ARM_VectorPtr NPVSwapScalar(
		const string& curveName, 
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
        const std::vector<double>& strikesPerState,
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
        const std::vector<double>& refIdxResettimes,
        const std::vector<double>& refIdxStarttimes,
        const std::vector<double>& refIdxEndtimes,
        const std::vector<double>& refFwdPeriods,
        const std::vector<double>& refIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
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
        const std::vector<double>&  payTimes,
        const std::vector<double>&  resetTimes,
        const std::vector<double>&  startTimes,
        const std::vector<double>&  endTimes,
        int     indexPaymentType,
        const std::vector<double>&  fwdPaymentPeriods,
        const ARM_VectorVector& refIdxResettimes,
        const ARM_VectorVector& refIdxStarttimes,
        const ARM_VectorVector& refIdxEndtimes,
        const ARM_VectorVector& refFwdPeriods,
        const ARM_VectorVector& refIndexWeights,
        const std::vector<double>&  couponMargins,
        const ARM_VectorVector& downBarriers,
        const ARM_VectorVector& upBarriers,
        const std::vector<double>&  payNotionals,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;

	// Vanilla Corridor function
	ARM_VectorPtr VanillaCMSCorridor(
		const string& curveName,
		double evalTime,
		const std::vector<double>& payTimes,
		const std::vector<double>& resetTimes,
		const std::vector<double>& startTimes,
		const std::vector<double>& endTimes,
		const ARM_VectorVector& refIdxResettimes,
		const ARM_VectorVector& refIndexWeights,
		const ARM_VectorVector& coeffs1,
		const ARM_SwapRatePtrVectorVector& firstIndexes,
		const ARM_VectorVector& coeffs2,
		const ARM_SwapRatePtrVectorVector& secondIndexes,
		const ARM_IntVector& payIndexType,		/// K_FIXED, K_CMS, K_LIBOR
		const std::vector<double>& coupon,			/// in case of a fixed payment (K_FIXED)
		const ARM_SwapRatePtrVector& payRate,	/// rate description (if K_CMS or K_LIBOR)
		const std::vector<double>& payIndexLeverage,
		const ARM_VectorVector& downBarriers,
        const ARM_VectorVector& upBarriers,
        const std::vector<double>&  payNotionals,
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
		const std::vector<double>& refIdxResettimes,
		const std::vector<double>& refIndexWeights,
		const std::vector<double>& coeff1,
		const ARM_SwapRatePtrVector& firstIndex,
		const std::vector<double>& coeff2,
		const ARM_SwapRatePtrVector& secondIndex,
		int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
		double	coupon,					/// in case of a fixed payment (K_FIXED)
		const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
		double  payIndexLeverage,
		const std::vector<double>& downBarriers,
        const std::vector<double>& upBarriers,
        double  payNotional,
        int     rcvPay,
		const ARM_SwapRatePtrVector& thirdIndex, // 3rd index for double condition
		const std::vector<double>& downBarriers3,
		const std::vector<double>& upBarriers3,
        const   ARM_PricingStatesPtr& states) const;

    /// Vanilla cap function (vectorial strike version)
	ARM_VectorPtr VanillaCap(
		const string& curveName,
		double evalTime,
		const std::vector<double>& payTimes,
		const std::vector<double>& periods,
		double notional,
		const std::vector<double>& fwdResetTimes,
		const std::vector<double>& fwdStartTimes,
		const std::vector<double>& fwdEndTimes,
		const std::vector<double>& fwdPeriods,
		const std::vector<double>& strikesPerState,
		int   capFloor,    
		const ARM_PricingStatesPtr& states) const;

    /// Vanilla cap function (scalar strike version)
    ARM_VectorPtr VanillaCapScalar(
		const string& curveName,
		double evalTime,
		const std::vector<double>& payTimes,
		const std::vector<double>& periods,
		double notional,
		const std::vector<double>& fwdResetTimes,
		const std::vector<double>& fwdStartTimes,
		const std::vector<double>& fwdEndTimes,
		const std::vector<double>& fwdPeriods,
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
		const std::vector<double>& coeffs,
		const std::vector<double>& fwdResetTimes,
		const std::vector<double>& fwdStartTimes,
		const std::vector<double>& fwdEndTimes,
		double payTime,
		const std::vector<double>& fwdPeriods,
		const std::vector<double>& strikesPerState,
		double volatilityRatio,
		double* sumFwd,
		double* sumVol,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSmiledSwaption(
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
		const std::vector<double>& fixTimes,
        const std::vector<double>& fixPayPeriods,
        const std::vector<double>& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		const std::vector<double>& data,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	virtual ARM_VectorPtr MaxRate(
		const string& curveName, 
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
        const std::vector<double>& firstStrikeDown,
        const std::vector<double>& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const std::vector<double>& secondStrikeDown,
        const std::vector<double>& secondStrikeUp,
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
        const std::vector<double>& strikesPerState,
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
