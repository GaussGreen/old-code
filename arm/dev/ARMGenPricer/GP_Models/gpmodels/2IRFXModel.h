/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file 2IRFXModel.h
 *
 *  \brief
 *
 *  \brief 2 interest + fx multi assets model
 *
 *	\author  E. Benhamou, JM Prié
 *	\version 1.0
 *	\date January 2005
 */



#ifndef _INGPMODELS_2IRFXMODEL_H
#define _INGPMODELS_2IRFXMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_ModelParamsHW1FStd;
class ARM_ModelParamsQ1F;
class ARM_CurveModelParam;

class ARM_2IRFXModel :	public ARM_MultiAssetsMeanReverting
{
private:

    bool itsIsFlooredFxLocalModel;
    bool itsIsCappedFxLocalModel;
    bool itsIsRedemptionFxLocalModel;

    ARM_GP_MatrixPtr itsDriftCorrections;

    void AddIntegratedLocalCorrections( const ARM_GP_Vector& timeSteps,ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// to initialise the sub models
	void Validate();

	typedef void (ARM_2IRFXModel::*InitCalibrationFunc)(
		double evalTime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    void InitDomCalibration(double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const;
    void InitForCalibration(double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const;
	void InitSwapCalibration(
		double evalTime,
		double floatStartTime,
		double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& floatPayTimes, 
		const ARM_GP_Vector& floatPayPeriods, 
		const ARM_GP_Vector& margin,
		bool isDbleNotional,
		const ARM_PricingStatesPtr& states,
		const InitCalibrationFunc& initFunc) const;

    void AdjustVolToSchedule(const ARM_GP_VectorPtr& timeSteps, ARM_PricingModelPtr& domModel, ARM_PricingModelPtr& forModel, ARM_PricingModel& fxModel);

    bool IsBasisRefModel() const { return GetRefModel() == &*((*(GetModelMap()))[DomBasisModel]->Model()); }

    void ComputeIntegratedFwdFxVCV(double step, double nextStep,
            const ARM_ModelParamsHW1FStd* const domModelParams,
            const ARM_ModelParamsHW1FStd* const forModelParams,
            const ARM_ModelParamsHW1FStd* const fxModelParams,
            const ARM_CurveModelParam& fxVol,
            const ARM_GP_Matrix& correlMatrix,
            ARM_GP_Matrix& variances) const;

    void ComputeIntegratedSpotFxVCV(double step, double nextStep,
            const ARM_ModelParamsHW1FStd* const domModelParams,
            const ARM_ModelParamsHW1FStd* const forModelParams,
            const ARM_ModelParamsHW1FStd* const fxModelParams,
            const ARM_GP_Matrix& correlMatrix,
            ARM_GP_Matrix& variances) const;

	void ComputeIntegratedSpotFxVCVInQModel(double step, double nextStep,
            const ARM_ModelParamsHW1FStd* const domModelParams,
            const ARM_ModelParamsHW1FStd* const forModelParams,
            const ARM_ModelParamsQ1F* const fxModelParams,
            const ARM_GP_Matrix& correlMatrix,
            ARM_GP_Matrix& variances) const;

public:
	/// Specialised init method for tree because of 1D tree calibrations
	virtual ARM_PricingStatesPtr BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);


	/// Temporary specialised init method for MC waiting for actual use of sampler
	virtual ARM_PricingStatesPtr ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);
	
public:
    static const size_t NbLocalModels;

    enum modelsAlias
    {
        DomModel	=0,         /// the domestic stochastic IR model
        ForModel,			    /// the foreign stochastic IR model
        FxModel,			    /// the stochastic FX model
		DomBasisModel,		    /// pure basis model
		ForBasisModel,		    /// pure basis model
        FlooredFxLocalModel,    /// optional local volatility model for floored Fx coupon
        CappedFxLocalModel,     /// optional local volatility model for capped Fx coupon
        RedemptionFxLocalModel, /// optional local volatility model for redemption Fx coupon
        NbModels
    };

	ARM_2IRFXModel(
		const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlCurveMatrix );

	ARM_2IRFXModel(const ARM_2IRFXModel& rhs);
	ASSIGN_OPERATOR(ARM_2IRFXModel)
	virtual ~ARM_2IRFXModel(){}
	
    const ARM_PricingModelPtr& GetModel(modelsAlias modelIdx) const { return (* GetModelMap())[modelIdx]->Model(); }

    virtual ARM_BoolVector NeedMCIntegProcess() const;

    /// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { return false; }

    virtual bool SupportAnalyticMarginal() const;

    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	void ComputeMeanRevertingVariances( const ARM_GP_Vector& timeSteps,
	    ARM_MatrixVector& localVariances,
	    ARM_MatrixVector& variances ) const;

    void NumMethodStateLocalGlobalVariances( const ARM_GP_Vector& timeSteps,
	    ARM_MatrixVector& localVariances,
	    ARM_MatrixVector& variances ) const;

	void NumMethodStateGlobalVariances( 
		const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& variances ) const;

    /// No more used because variance initialisation is done through the numerical method sampler
    void ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps) {}
    void NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps ) {}

    void AdjustVol();

    virtual ARM_VectorPtr LocalDiscounts(size_t timeIdx, 
		double dt, const ARM_PricingStatesPtr& states) const;

    /// Give the local payoff for numerical method calibration
    virtual ARM_VectorPtr LocalPayoffs(size_t timeIdx, 
		double dt, const ARM_PricingStatesPtr& states) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
	virtual void PdeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex, double lambda) const;
    virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

    virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;

	/// function to propagate the markovian drift from one to another model
    virtual ARM_GP_MatrixPtr IntegratedMarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, const ARM_GP_VectorPtr& driftCorrection ) const;
    virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates) const;

    /// To compute drifts of diffused processes from 0 to each time steps 
    virtual void IntegratedGlobalDrifts(const ARM_GP_Vector& timeSteps,ARM_GP_MatrixPtr& drifts);

	/// function to compute the relative and absolute drifts with a Euler scheme
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// function to compute the relative and absolute drifts with an integrated scheme
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	//function to compute (t,T) -> q(t,T) and (t,T) -> sigma(t,T) of fwd FX
	virtual ARM_VectorPtr ComputeFwdFXModelParam(
	    double evalTime,
	    double settlementTime,
	    ARM_PricingContext* context=NULL) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_2IRFXModel(*this);}
	virtual string ExportShortName() const { return "L2IFX";}

    virtual ARM_VectorPtr DiscountFactor( 
		const string& modelName,
        double evalTime, 
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    virtual ARM_VectorPtr Libor( 
		const string& modelName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr NPVSwap(
	const string& modelName, 
	double evalTime,
	double floatStartTime,
	double floatEndTime, 
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
	const ARM_GP_Vector& fwdStartTimes, 
	const ARM_GP_Vector& fwdEndTimes, 
	const ARM_GP_Vector& fwdPayPeriods, 
	const ARM_GP_Vector& floatPayTimes, 
	const ARM_GP_Vector& floatPayPeriods, 
	const ARM_GP_Vector& margin,
	bool isDbleNotional,
	const ARM_GP_Vector& FixNotional,
	const ARM_GP_Vector& FloatNotional,
	const ARM_GP_Matrix& strikesPerState,
	int payRec,
	const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr NPVBasisSwap( 
		const string& domCurveName,
		const string& forCurveName, 
		const string& fxCurveName, 
		double evalTime,
		double startTime,
		double endTime,
		int		payRec,
		const ARM_GP_Vector& domResetTimes,	    
		const ARM_GP_Vector& domFwdStartTimes,
		const ARM_GP_Vector& domFwdEndTimes,
		const ARM_GP_Vector& domFlowStartTimes,			
		const ARM_GP_Vector& domFlowEndTimes,	
		const ARM_GP_Vector& domFwdPayPeriods,	
		const ARM_GP_Vector& domPayTimes,
		const ARM_GP_Vector& domPayPeriods,
		const ARM_GP_Vector& domMarginVector,
		const ARM_GP_Vector& domNotionalVector,
		bool                 isDomFlottant,
		const ARM_GP_Vector& forResetTimes,       
		const ARM_GP_Vector& forFwdStartTimes,
		const ARM_GP_Vector& forFwdEndTimes,
		const ARM_GP_Vector& forFlowStartTimes,   		
		const ARM_GP_Vector& forFlowEndTimes,	    
		const ARM_GP_Vector& forFwdPayPeriods,	
		const ARM_GP_Vector& forPayTimes,
		const ARM_GP_Vector& forPayPeriods,
		const ARM_GP_Vector& forMarginVector,
		const ARM_GP_Vector& forNotionalVector,
		bool                 isForFlottant,
		const string&        exNotionalType,   
		const ARM_GP_Vector& fxResetTimes,
		const ARM_GP_Vector& fxSettlTimes,  
		const ARM_GP_Vector& fxPayTimes,
		const ARM_GP_Matrix& domStrikePerState,
		const ARM_GP_Matrix& forStrikePerState,
		const ARM_PricingStatesPtr& states ) const;

    virtual ARM_VectorPtr Forward(
	    const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

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

	virtual double UnderlyingCovariance (  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const;

	void ARM_2IRFXModel::UpdatePDE3DCoeffs(
		size_t timeIdx,
		const ARM_PricingStatesPtr& states,
		ARM_GP_VectorPtr& qxx,
		ARM_GP_VectorPtr& qyy,
		ARM_GP_VectorPtr& qzz,
		ARM_GP_VectorPtr& qxy,
		ARM_GP_VectorPtr& qyz,
		ARM_GP_VectorPtr& qzx,
		ARM_GP_VectorPtr& px,
		ARM_GP_VectorPtr& py,
		ARM_GP_VectorPtr& pz,
		ARM_GP_VectorPtr& o,
		double lambda,
		bool IsInit
		) const;

	/// Just to allow correct branching but Q1F pricing function doesn't work for
	/// variable notional swaption... (to do)
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

	// Compute the correlation used for the 1IRFX
	void ComputeVolatilitiesAndCorrelMatrix(
			const ARM_GP_Vector& resetTimes,
			const ARM_GP_Vector& settlementTimes,
			ARM_GP_Matrix& volatilities,
			ARM_MatrixVector& correlMatrixVector,
			ARM_GP_Vector& totalVolatilities,
			bool withDomesticIR) const;

	// Compute the volatilities used for the 1IRFX
	void ComputeVolatilities(
			const ARM_GP_Vector& resetTimes,
			const ARM_GP_Vector& settlementTimes,
			ARM_GP_Matrix& volatilities,
			ARM_GP_Vector& totalVar) const;


	virtual string toString(const string& indent="",const string& nextIndent="") const;
	
	double ComputeTimeLag(double expiryTime, double payTime, bool isDomestic = true) const;
	double ComputeAbsTimeLag(double expiryTime, double payTime, bool isDomestic = true) const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
