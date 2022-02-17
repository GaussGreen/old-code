/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file pricingfunctioninflation.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

#ifndef _INGPINFRA_PRICINGFUNCTIONINFLATION_H
#define _INGPINFRA_PRICINGFUNCTIONINFLATION_H
 

#include "gpbase/env.h"
#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )
 /// macro for namespace ... define namespace only if supported


///////////////////////////////////////////////////////
/// \class ARM_PricingFunctionInflation
/// \brief
/// This abstract class is the interface for inflation 
/// functions of the GP
///////////////////////////////////////////////////////
class ARM_PricingFuncInflation
{
public:
	/// CPI Spot
	virtual ARM_GP_VectorPtr CPISpot( 
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, string DCFLag, long DailyInterp,
		string ResetLag,
		const ARM_PricingStatesPtr& states) const = 0;

	/// CPI Forward
	virtual ARM_GP_VectorPtr CPIForward(
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, 
		double FixingTime, 
		const ARM_PricingStatesPtr& states) const = 0;

	/// Convexity Adjustment
	virtual ARM_GP_VectorPtr ConvexityAdjustment(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const = 0;

	/// forward Ratio
	virtual ARM_GP_VectorPtr ForwardCPIRatio(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double CPITime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const = 0;

	/// YoY/OAT Swaps and SwapRates
 	virtual ARM_GP_VectorPtr YoYSwapRate(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const = 0;

 	virtual ARM_GP_VectorPtr YoYSwap(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const = 0;

	virtual ARM_GP_VectorPtr OATSwapRate(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsCoupon,
		const ARM_PricingStatesPtr& states) const = 0;

	virtual ARM_GP_VectorPtr OATSwap(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsCoupon,
		const ARM_PricingStatesPtr& states) const = 0;


	virtual ARM_GP_VectorPtr YoYCapFloor(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		int CapFloor,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const = 0;

	virtual ARM_GP_VectorPtr OATCapFloor(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		int CapFloor,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const = 0;

	virtual ARM_GP_VectorPtr ZCCap( ) const = 0;
};


CC_END_NAMESPACE()

#endif

