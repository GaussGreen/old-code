/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file IRFwdMod.h
 *
 *  \brief prototype model for the generic pricer
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_IRFWDMOD_H
#define _INGPINFRA_IRFWDMOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/pricingmodelir.h"

CC_BEGIN_NAMESPACE( ARM )


/// \class ARM_IrFwdMod
/// \brief simple interest rate forward curve model
/// mainly developped for testing the parser!

class  ARM_IrFwdMod: public ARM_PricingModelIR
{
public:
	ARM_IrFwdMod(const ARM_ZeroCurvePtr& zeroCurve);
	ARM_IrFwdMod(const ARM_IrFwdMod& rhs);
	ARM_IrFwdMod& operator=(const ARM_IrFwdMod& rhs);
	virtual ~ARM_IrFwdMod();

    /// Schedule initialisation, pre-computed datas
    /// and numerical method initialisation
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
        const ARM_PricingStatesPtr& states) const;


    // Default VanillaCaplet not implemented
    // (a Libor volatility function should be defined)
    virtual ARM_VectorPtr VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime,			/// not used for convexity correction
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;

	// Default VanillaSwaption not implemented
    // (a Swap volatility function should be defined)
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
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;

	/// induction is at the model stage ..
	/// this is exceptional as for other more complicated model
	/// induct calls the numerical method
    // Backward/forward induction to get prices at toTime time
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    /// Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;
	
	/// fonction to compute local drifts
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states ) const;

    /// Give the time to reach a given variance
    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

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
													const std::vector<double>& strikes,
													double swapLongFloatStartTime,
													double swapLongFloatEndTime,
													const std::vector<double>& swapLongFixPayTimes,
													const std::vector<double>& swapLongFixPayPeriods,
													double swapShortFloatStartTime,
													double swapShortFloatEndTime,
													const std::vector<double>& swapShortFixPayTimes,
													const std::vector<double>& swapShortFixPayPeriods,
													const ARM_PricingStatesPtr& states) const;

	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual std::vector<double>* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (25/05/2007 11:28:32):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const; 
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const; 

     /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalibSecIndex(size_t index, ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod){};
	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const {return true;}


	/// -------------------- standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LIRFM"; }

	/// model specific answer for generic numerical method!
	virtual bool SupportBackwardInduction() const { return true; }
	virtual bool SupportForwardInduction()  const { return true; }
	virtual bool SupportAnalyticMarginal()  const {	return true; }
	virtual int GetType() const;
	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	virtual size_t FactorCount() const { return 1; }
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

