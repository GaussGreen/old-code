/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ForwardMargin.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */

#ifndef _INGPMODELS_FORWARDMARGIN_H
#define _INGPMODELS_FORWARDMARGIN_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpinfra/pricingmodelir.h"
#include "gpinfra/typedef.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;

//-----------------------------------------------------------------------------
// \class ARM_ForwardMargin
// \brief
//  Model for a deterministic margin between two yield curves.
//  The model assumes that the multiplicative margin has no volatility  :
//  margin(t,T) = margin(0,T)/margin(0,t) = its forward value with
//  margin(0,u) = Bshift(0,u)/Bref(0,u)
//-----------------------------------------------------------------------------

class ARM_ForwardMargin :   public ARM_PricingModelIR
{
private :
    bool itsRefModelDump;
    ARM_PricingModel* itsRefModel;
    void CheckRefModel() const
    {
	    if( !itsRefModel )
		    ARM_THROW( ERR_INVALID_ARGUMENT, " refModel is NULL!" );
    }



public:
	ARM_ForwardMargin( 
		const ARM_ZeroCurvePtr& shiftZcCurve,
        ARM_PricingModel* refModel = NULL,
        bool refModelDump = false);

	ARM_ForwardMargin(const ARM_ForwardMargin& rhs);
	virtual ~ARM_ForwardMargin();
    ARM_ForwardMargin& operator = (const ARM_ForwardMargin& rhs);

    bool GetRefModelDump() const { return itsRefModelDump; }
    void SetRefModelDump(bool flag) { itsRefModelDump=flag; }

	virtual const ARM_PricingModel* GetRefModel() const { return itsRefModel; }
	virtual ARM_PricingModel* GetRefModel() { return itsRefModel; }

    void SetRefModel( ARM_PricingModel* refModel ) { itsRefModel=refModel; }

    /// Deterministic margin ratio computation
    inline double ComputeMarginRatio(double evalTime,double maturityTime) const;

    /// Closed form formulas
    virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName,
        double evalTime, 
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
    /// Always true because no model params
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const {return true;}

	/// typing of the model
	virtual int GetType() const;
	virtual size_t FactorCount() const{ return 0.0; } /// no model parameter in this non stochastic model!

	virtual void ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const = 0;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LFWMA";}

	/// function to communicate to the reference model
	virtual void SetNumeraire(const ARM_NumerairePtr& numerairePtr);
	virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);

	/// unsupported function
    /// --- calibration calibration part
	virtual void IntegratedLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const;
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter);
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
    virtual const ARM_ModelParams* const GetModelParams() const;
    virtual ARM_ModelParams* GetModelParams();

	virtual ARM_VectorPtr IntegratedRiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded=false ) const
    { CheckRefModel(); return itsRefModel->IntegratedRiskNeutralDrift(timeIdx,numMethodStates,isDriftAdded ); }
	virtual ARM_VectorPtr RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded=false ) const
    { CheckRefModel(); return itsRefModel->RiskNeutralDrift(timeIdx,numMethodStates,isDriftAdded ); }

	/// --- numerical functions for numerical method
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	virtual std::vector<double>* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (24/05/2007 10:41:59):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	virtual bool SupportBackwardInduction() const;
	virtual bool SupportForwardInduction() const;
	virtual bool SupportAnalyticMarginal() const;
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,		ARM_MatrixVector& localVariances ) const;
	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,	ARM_MatrixVector& localVariances ) const;
	virtual bool NeedsToCholeskyDecomposeFactors( ) const;
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const;
	
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

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual ARM_VectorPtr LocalDiscounts( size_t timeIdx, double dt, const ARM_PricingStatesPtr& states) const;

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );
};


double ARM_ForwardMargin::ComputeMarginRatio(double evalTime,double maturityTime) const
{
    /*double refDFT    = itsRefModel->GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
    double refDFt    = itsRefModel->GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
    double shiftedDFT= GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
    double shiftedDFt= GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
*/
    return 0;//(shiftedDFT*refDFt)/(shiftedDFt*refDFT);
}


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
