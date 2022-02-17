/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfFwdMod.h
 *
 *  \brief prototype model for the generic pricer
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */

#ifndef _INGPMODELS_INFEQMOD_H
#define _INGPMODELS_INFEQMOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/pricingmodelinflation.h"
#include "gpinfra/pricingmodelir.h"
#include "gpmodels/CPIManager.h"

CC_BEGIN_NAMESPACE( ARM )



/// \class ARM_InflationEquityMod
/// \brief simple interest rate forward curve model
/// mainly developped for testing the parser!

class  ARM_InflationEquityMod: public ARM_PricingModelInflation
{
private: 
	mutable double itsInitValue;
	double itsPublicationLag;
	size_t itsOtherModelRank;

	double integratedVolForBSYoY( double FromTime, double ToTime, double tenor ) const;
	ARM_GP_VectorPtr YoYCaplet( const string& irCurveName, const string& infCurveName, double EvalTime, double tenor, double CPITime, double maturityTime, double Strike, int CapFloor, const ARM_PricingStatesPtr& states ) const;
	CPIManagerPtr itsCPIManager;


	/// CPI Spot
	void ComputeCPISpotAndStore( 
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, string DCFLag, long DailyInterp,
		string ResetLag,
		const ARM_PricingStatesPtr& states) const;

public:

	/// ================================================================= ///
    /// ====================== Constructor/Destructor =================== ///
	/// ================================================================= ///

	ARM_InflationEquityMod(const ARM_InfCurvPtr& infc, const ARM_PricingModelIRPtr& irMod = ARM_PricingModelIRPtr(NULL), double publicationLag = 61, const ARM_ModelParams* params = NULL);
	ARM_InflationEquityMod(const ARM_InflationEquityMod& rhs);
	ARM_InflationEquityMod& operator=(const ARM_InflationEquityMod& rhs);
	virtual ~ARM_InflationEquityMod();


	/// ================================================================= ///
	/// ====================Pricing Model Support ======================= ///
	/// ================================================================= ///

    /// Schedule initialisation, pre-computed datas
    /// and numerical method initialisation
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    /// Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateLocalCovariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map );

	/// Inflation Equity is a kind of MeanReverting
	virtual bool IsMeanRevertingCompatible() const { return true;}
	
	/// fonction to compute local drifts
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    /// Give the time to reach a given variance
    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual ARM_GP_Vector* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (24/05/2007 10:46:52):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const; 
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

	virtual ARM_GP_Vector* PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos);

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

	/// ================================================================= ///
    /// ====================== Calibration purpose ====================== ///
	/// ================================================================= ///

    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index, ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
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
	virtual string ExportShortName() const { return "LIFEQ"; }

	/// model specific answer for generic numerical method!

	virtual bool SupportBackwardInduction() const { return true; }
	virtual bool SupportForwardInduction()  const { return true; }
	virtual bool SupportAnalyticMarginal()  const {	return true; }
	virtual int GetType() const;
	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(1, true); };

	/// ===================================================================== ///
	/// ========================= Internal Purpose ========================== ///
	/// ===================================================================== ///

	double IntegratedDividend( double FromTime, double ToTime ) const;

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

	/// OAT Cap Floors
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

