/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingfunctionequity.h
 *
 *  \brief
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */

#ifndef _INGPINFRA_PRICINGFUNCTIONEQUITY_H
#define _INGPINFRA_PRICINGFUNCTIONEQUITY_H
 

#include "gpbase/env.h"
#include "typedef.h"
#include <glob/expt.h>

class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )
 /// macro for namespace ... define namespace only if supported

 struct ARM_PricingContext;
struct ARM_VanillaArg;

///////////////////////////////////////////////////////
/// \class ARM_PricingModelEquity
/// \brief
/// This abstract class is the interface for equity keyword
/// functions of the GP
///////////////////////////////////////////////////////
class ARM_PricingFunctionEquity
{
//private:
//	ARM_PricingModelPtr itsIRModel;

public:
	ARM_PricingFunctionEquity() {}; //:itsIRModel(NULL) {};
	ARM_PricingFunctionEquity(const ARM_PricingFunctionEquity& rhs)
	{};
	//: itsIRModel(rhs.itsIRModel){};
//	ARM_PricingFunctionEquity& operator=( const ARM_PricingFunctionEquity& rhs)
//	{
//		if( this!= & rhs )
			//itsIRModel =rhs.itsIRModel;
//		return *this;
//	}

	virtual ~ARM_PricingFunctionEquity() {};


	/// ----------- default functions
    /// Call function (scalar strike version)
	virtual ARM_VectorPtr CallScalar(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		double strike,
        int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// Digital call function (scalar strike version)
	virtual ARM_VectorPtr DigitalScalar(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		double strike,
		double notional,
        int callPut,
	    double payTime,
		ARM_DigitType digitType,
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// Default Digital call function (vectorial strike version)
	virtual ARM_VectorPtr DigitalVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikesPerState,
		double notional,
        int callPut,
	    double payTime,
		ARM_DigitType digitType, 
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// ----------- pure virtual functions
	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const = 0;

	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikesPerState,
        int callPut,
	    double payTime,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const = 0;

	/// CallStrip function (vectorial strike version)
	virtual ARM_VectorPtr CallStripVectorial(
		const string& modelName,
        double evalTime,
	    const ARM_GP_Vector& expiryTime,
	    const ARM_GP_Vector& settlementTime,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
	    const ARM_GP_Vector& payTime,
	    const ARM_GP_Vector& nominal,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	virtual ARM_VectorPtr CallStripScalar(
		const string& modelName,
        double evalTime,
	    const ARM_GP_Vector& expiryTime,
	    const ARM_GP_Vector& settlementTime,
		const ARM_GP_Vector& strike,
        int callPut,
	    const ARM_GP_Vector& payTime,
	    const ARM_GP_Vector& nominal,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	virtual ARM_VectorPtr HybridCallVectorial(
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
		ARM_PricingContext* context=NULL) const;

	virtual ARM_VectorPtr HybridCallScalar(
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
		ARM_PricingContext* contex=NULL) const;

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

	 virtual double ImpliedVol(const ARM_VanillaArg& arg) const = 0;


	/// for multi asset type equity model!
//	virtual void SetIRModel( const ARM_PricingModelPtr& irModel );
//	inline ARM_PricingModelPtr GetIRModel() const { return itsIRModel; }
	static void ValidateIRModel( const ARM_PricingModelPtr& irModel );

	/// convention support
	virtual string GetSettlementCalendar(const string& modelName="") const	= 0;
	virtual double GetSettlementGap(const string& modelName="") const		= 0;

	// default static versions
	static string GetSettlementCalendar(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve);
	static double GetSettlementGap(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve);
};


CC_END_NAMESPACE()

#endif

