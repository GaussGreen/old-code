/*!
 *
 * Copyright (c) IXIS CIB Paris 2005
 *
 *	\file ForwardMarginIR.h
 *
 *  \brief
 *
 *	\author  E. Benhamou, J.M. Prié
 *	\version 1.0
 *	\date February 2005
 */

#ifndef _INGPMODELS_FORWARDMARGINIR_H
#define _INGPMODELS_FORWARDMARGINIR_H


/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "ForwardMargin.h"

CC_BEGIN_NAMESPACE( ARM )


///-----------------------------------------------------------------------------
/// \class ARM_ForwardMarginIR
/// \brief
///  Model for a deterministic margin between two yield curves.
///-----------------------------------------------------------------------------

class ARM_ForwardMarginIR : public ARM_ForwardMargin
{
public:
	/// ----------------  constructor
	ARM_ForwardMarginIR( const ARM_ZeroCurvePtr& shiftZcCurve,
        ARM_PricingModel* refModel = NULL,
        bool refModelDump = false)
		:	ARM_ForwardMargin( shiftZcCurve, refModel, refModelDump) {}
	
	ARM_ForwardMarginIR( const ARM_ForwardMarginIR& rhs ) : ARM_ForwardMargin( rhs ) {}
	virtual ~ARM_ForwardMarginIR() {}
	ARM_ForwardMarginIR& operator=( const ARM_ForwardMarginIR& rhs )
	{
		if( this != &rhs )
			ARM_ForwardMargin::operator=(rhs);	
		return *this;
	}


	/// ------------------ pricing function
	virtual ARM_VectorPtr Libor( 
		const string& modelName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr SwapRate(
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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr NPVSwap(
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
		const ARM_PricingStatesPtr& states) const;

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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaCaplet(
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
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaDigital(
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
        const ARM_PricingStatesPtr& states) const;



	virtual ARM_VectorPtr VanillaCorridorlet(
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
        const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSwaption(
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
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	virtual ARM_VectorPtr LocalDiscounts( size_t timeIdx, 
		double dt, 
		const ARM_PricingStatesPtr& states) const;

	/// -------- standard ARM root object support
	virtual ARM_Object* Clone() const { return new ARM_ForwardMarginIR(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
