/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QNF.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_QNF_H
#define _INGPMODELS_QNF_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "QBaseIR.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelParamsQNF;

//-----------------------------------------------------------------------------
// \class ARM_QModelNF
// \brief
//  Class for the Q Model analytics
//-----------------------------------------------------------------------------

class ARM_QModelNF : public ARM_QModelBaseIR
{
public:
	ARM_QModelNF(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQNF& params, bool degenerateInHW = false );
	ARM_QModelNF(const ARM_QModelNF& rhs);
    ARM_QModelNF& operator = (const ARM_QModelNF& rhs);
	virtual ~ARM_QModelNF();

    /// Reconstruction formula
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
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
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

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
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike= true) const;
    
   /// Default initialisation of the model
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
    virtual std::vector<double>* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos ) {return NULL;}
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	/// function for the generic tree
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// relative and absolute drifts
    virtual void EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    /// Give the local discount value between
    /// startTime and endTime for each states at startTime
    virtual ARM_VectorPtr LocalDiscounts(size_t timeIdx, 
		double dt, const ARM_PricingStatesPtr& states) const;

    // Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localvariances ) const;

	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& variances ) const;

    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0;}

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

    virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LQNFM";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
