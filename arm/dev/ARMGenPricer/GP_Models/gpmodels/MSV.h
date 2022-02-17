/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV.h
 *
 *  \brief 
 *
 *	\author  A  Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_MSV_H
#define _INGPMODELS_MSV_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "gpnumlib/odefunctions.h"
#include "gpmodels/SVModels.h"


#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_MarkovSV
// \brief
//  Interface class for Hull & White pricing models
//-----------------------------------------------------------------------------

/// Internal structures To Store Swaption Pricing Data

class ARM_MarkovSV : public ARM_SVModels
{
public:
	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

    static const double VOL_NORMAL_MAX;

	ARM_MarkovSV( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL );
	ARM_MarkovSV(const ARM_MarkovSV& rhs);
	virtual ~ARM_MarkovSV();

    ARM_MarkovSV& operator = (const ARM_MarkovSV& rhs);

    /// Closed form formulas
	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const;

    virtual ARM_VectorPtr VanillaCaplet(
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
		const ARM_PricingStatesPtr& states) const;

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
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;
    
    virtual ARM_VectorPtr VanillaDigital(
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
        const ARM_PricingStatesPtr& states) const;  


	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;



	virtual ARM_VectorPtr  VanillaSpreadOptionLet(const string& curveName,
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
													const ARM_PricingStatesPtr& states) const;
	
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return false;}

	/// no post init hence nothing
	virtual void PostInit(){};

	virtual ARM_GP_Vector* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );

	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	/// Volatilities and correlation time steps for PDEs
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;

	/// MSV is MEanReverting
	virtual bool IsMeanRevertingCompatible() const { return true;}

	virtual bool ValidateModelParams(const ARM_ModelParams& params) const = 0;

	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_MarkovSV : abstract class !"); }

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
