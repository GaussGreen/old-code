/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfFwdMod.h
 *
 *  \brief prototype model for the generic pricer
 *
 *	\author  N. Belgrade A. Schauly
 *	\version 1.0
 *	\date March 2005
 */

#ifndef _INGPMODELS_INFFWDMOD_H
#define _INGPMODELS_INFFWDMOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/pricingmodelinflation.h"
#include "gpinfra/pricingmodelir.h"

CC_BEGIN_NAMESPACE( ARM )

/// FwdDeclarations: 
class ARM_MultiAssetsModel;


/// \class ARM_InfFwdMod
/// \brief simple interest rate forward curve model
/// mainly developped for testing the parser!

class  ARM_InfFwdMod: public ARM_PricingModelInflation
{
public:

	/// ================================================================= ///
    /// ====================== Constructor/Destructor =================== ///
	/// ================================================================= ///

	ARM_InfFwdMod(const ARM_InfCurvPtr& infc);
	ARM_InfFwdMod(const ARM_InfFwdMod& rhs);
	ARM_InfFwdMod& operator=(const ARM_InfFwdMod& rhs);
	virtual ~ARM_InfFwdMod();


	/// ================================================================= ///
	/// ====================Pricing Model Support ======================= ///
	/// ================================================================= ///

    /// Schedule initialisation, pre-computed datas
    /// and numerical method initialisation
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	/// induction is at the model stage ..
	/// this is exceptional as for other more complicated model
	/// induct calls the numerical method
    // Backward/forward induction to get prices at toTime time
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    /// Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;
	
	/// fonction to compute local drifts
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states ) const;

    /// Give the time to reach a given variance
    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual ARM_GP_Vector* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (24/05/2007 10:46:11):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const; 
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const; 

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

	/// ================================================================= ///
    /// ====================== Calibration purpose ====================== ///
	/// ================================================================= ///

    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalibSecIndex(size_t index, ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod){};

	/// method to advise the break point times

    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const {return true;}

	virtual size_t FactorCount() const { return 1;}

	/// ===================================================================== ///
	/// ==================== standard ARM Object support ==================== ///
	/// ===================================================================== ///

	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LINFM"; }

	/// model specific answer for generic numerical method!

	virtual bool SupportBackwardInduction() const { return true; }
	virtual bool SupportForwardInduction()  const { return true; }
	virtual bool SupportAnalyticMarginal()  const {	return true; }
	virtual int GetType() const;
	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }


	/// ===================================================================== ///
    /// ================== elementary inflation pricing functions =========== ///
	/// ===================================================================== ///

	/// CPI Spot
	virtual ARM_GP_VectorPtr CPISpot( 
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, string DCFLag, long DailyInterp,
		string ResetLag,
		const ARM_PricingStatesPtr& states) const;

	/// CPI Forward
	virtual ARM_GP_VectorPtr CPIForward(
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, 
		double FixingTime, 
		const ARM_PricingStatesPtr& states) const;

	/// Convexity Adjustment
	virtual ARM_GP_VectorPtr ConvexityAdjustment(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

	/// forward Ratio
	virtual ARM_GP_VectorPtr ForwardCPIRatio(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double CPITime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

	/// YoY Cap
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
		const ARM_PricingStatesPtr& states) const;

	/// OAT Cap
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
		const ARM_PricingStatesPtr& states) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

