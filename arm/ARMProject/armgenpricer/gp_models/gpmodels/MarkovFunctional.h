/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarkovFunctional.h
 *
 *  \brief 
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_MF_H
#define _INGPMODELS_MF_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "MFDFMap.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

class ARM_DensityFunctor;


class ARM_MarkovFunctional : public ARM_PricingModelIR
{
private: 
	/// Map that contains discount factors
	mutable ARM_DiscountFactorMap itsDiscountFactorMap;

	/// To know if we are pricing with the same numerical method we used for calibration
	bool itsCurrentModelIsCalibrationModel;

	/// Functor For Calibration
	ARM_NumericalModelFitter * itsNumericalModelFitter;

	/// Are we calibrating or not ?
	bool itsCalibrationStatus;

	/// nb of stdev used for calibration
	double itsCalibrationNbStdev;

	/// previous Induct time
	double itsPrevToTime;
	
	/// DFs (payoffs) for previous induct time
	ARM_MemPool_Matrix itsPrevPayoffs;


	/// Does Backward calibration
	void BackwardCalibration( const ARM_PricingStatesPtr& states, double toTime );
	/// Builds Discount factors during backward calibatoin
	void BuildDFs( const ARM_PricingStatesPtr& states, double toTime, size_t lineIdx );
	/// Correction of the proba changes
	void CorrectProbaChanges( const ARM_PricingStatesPtr& states, double toTime, size_t lineIdx, bool correctNumeraireOnly = true) const;

	void buildNormalProbabilities( ARM_GP_VectorPtr& probas, const ARM_PricingStatesPtr& states, double toResetTime ) const;

public:
	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

    static const double VOL_NORMAL_MAX;

	ARM_MarkovFunctional( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL );
	ARM_MarkovFunctional(const ARM_MarkovFunctional& rhs);
	virtual ~ARM_MarkovFunctional();

    ARM_MarkovFunctional& operator = (const ARM_MarkovFunctional& rhs);

	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

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
		bool isConstantStrike = true) const;

	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return true;}
    
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
		const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;  


	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(FactorCount(), true); };

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
	
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);

    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return true;}
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalCorrels( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map );

	virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual void IntegratedLocalDrifts(	const std::vector<double>& timeSteps,	ARM_GP_MatrixPtr& relativeDrifts,	ARM_GP_MatrixPtr& absoluteDrifts ) const;

	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	/// no post init hence nothing
	virtual void PostInit(){};

	virtual std::vector<double>* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );

		
	/// Volatilities and correlation time steps for PDEs
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;
    virtual void EulerLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const;
	void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, ARM_GP_MatrixPtr& vols,ARM_GP_MatrixPtr& d1Vols,ARM_GP_MatrixPtr& correls, bool linearVol = false) const;

	/// MF is MEanReverting
	virtual bool IsMeanRevertingCompatible() const { return true;}

	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;
    virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);

	/// Markov Functional Calibration
	virtual void setNumericalModelFitter( ARM_NumericalModelFitter * numericalModelFitter ) ; 
	virtual void setCalibrationStatus( bool calibrationStatus ) { itsCalibrationStatus=calibrationStatus; }

	virtual size_t FactorCount() const { return 1; }


	/// Markov Functional Numeraire = always Terminal date
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;
	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_MarkovFunctional(*this); };
	virtual string ExportShortName() const { return "LMF1M";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
