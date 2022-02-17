/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file HybridIRFX.cpp
 *
 *  \brief 1 IR or FX model  + 1 IR or FX model + 1 Payment model
 *
 *	\author  K.BELKHEIR
 * 	\version 1.0
 *	\date March 2007
 */



#ifndef _INGPMODELS_HybridIRFX_H
#define _INGPMODELS_HybridIRFX_H

//gpbase
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/numericconstant.h"

//gpmodels
#include "MultiAssets.h"
#include "gpmodels/EqFxBase.h"

//gpinfra
#include "gpinfra/modelnamemap.h"

//gpcalculator
//#include "gpcalculators\forexvanilla.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_GaussReplic2D;

////////////////////////////////////////////////////////////////////////////////
//
// This class defines a multi asset with 2 models:
// _ a 1 EqFXBase and  model associated to a FX1 rate
// _ a 1 EqFXBase model associated to a FX2 rate
// The two FX must have a common currency
//
//
////////////////////////////////////////////////////////////////////////////////

class ARM_HybridIRFX :	public ARM_MultiAssetsModel
{
public:
	//standard constructor
	ARM_HybridIRFX(const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix  = ARM_CurveMatrix());

	//copy constructor
	ARM_HybridIRFX(const ARM_HybridIRFX& rhs);
	ASSIGN_OPERATOR(ARM_HybridIRFX)
	virtual ~ARM_HybridIRFX()
	{
	}

	//Initialize the model
    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	///Financial vectorial functions

	// Digital function (vectorial strike version)
	virtual ARM_VectorPtr DigitalVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		double notional,
		int callPut,
	    double payTime,
		ARM_DigitType digitType,
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// CallSpread function (vectorial strike version)
	virtual ARM_VectorPtr CallSpreadVectorial(
		const string& Model1Name,
		const string& Model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const ARM_GP_Vector& strikePerState,
		double alpha,
		double beta,	
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// CallQuotient function (vectorial strike version)
	virtual ARM_VectorPtr CallQuotientVectorial(
		const string& Model1Name,
		const string& Model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const ARM_GP_Vector& strikePerState,
		double alpha,
		double beta,
		double strike2,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// RangeAccrual function (vectorial strike version)
	virtual ARM_VectorPtr RangeAccrualVectorial(
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
		ARM_Vector* eachFixingPrices=NULL,
        ARM_PricingContext* context=NULL) const;

	// SnowBall function
	ARM_VectorPtr SnowBallVectorial(
		const string& modelName,
		double evalTime,
		double payTime,
		const  ARM_GP_Vector& settelmentTimes,
		const  ARM_GP_Vector& fixingTimes,
		int	   callPut,
		const  ARM_GP_Vector& strikes,
		const  ARM_GP_Vector& coeffs,
		const  ARM_GP_Vector& leverages,
		const  ARM_GP_Vector& couponMin,
		const  ARM_GP_Vector& couponMax,
		const  ARM_PricingStatesPtr& states,
		bool isPerformOrNo = false,
		ARM_PricingContext* context=NULL) const;

	//Price of the call quanto
	double CallPrice( const ARM_EqFxBase& payoffModel, 
		const ARM_EqFxBase& quantoModel,
		ARM_GaussReplic2D* gaussReplic,
		double  fwd_payoff,
		double  fwd_quanto,
		double	expiryTime,
		const ARM_GP_Matrix& IntegParameters) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_HybridIRFX(*this);}
	virtual string ExportShortName() const { return "LIRFX";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
