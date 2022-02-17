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

    void AddIntegratedLocalCorrections( const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& absoluteDrifts) const;

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
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& floatPayTimes, 
		const std::vector<double>& floatPayPeriods, 
		const std::vector<double>& margin,
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

	void ComputeMeanRevertingVariances( const std::vector<double>& timeSteps,
	    ARM_MatrixVector& localVariances,
	    ARM_MatrixVector& variances ) const;

    void NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	    ARM_MatrixVector& localVariances,
	    ARM_MatrixVector& variances ) const;

	void NumMethodStateGlobalVariances( 
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& variances ) const;

    /// No more used because variance initialisation is done through the numerical method sampler
    void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps) {}
    void NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps ) {}

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
    virtual void IntegratedGlobalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& drifts);

	/// function to compute the relative and absolute drifts with a Euler scheme
    virtual void EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// function to compute the relative and absolute drifts with an integrated scheme
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
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
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const std::vector<double>& fwdStartTimes, 
	const std::vector<double>& fwdEndTimes, 
	const std::vector<double>& fwdPayPeriods, 
	const std::vector<double>& floatPayTimes, 
	const std::vector<double>& floatPayPeriods, 
	const std::vector<double>& margin,
	bool isDbleNotional,
	const std::vector<double>& FixNotional,
	const std::vector<double>& FloatNotional,
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
		const std::vector<double>& domResetTimes,	    
		const std::vector<double>& domFwdStartTimes,
		const std::vector<double>& domFwdEndTimes,
		const std::vector<double>& domFlowStartTimes,			
		const std::vector<double>& domFlowEndTimes,	
		const std::vector<double>& domFwdPayPeriods,	
		const std::vector<double>& domPayTimes,
		const std::vector<double>& domPayPeriods,
		const std::vector<double>& domMarginVector,
		const std::vector<double>& domNotionalVector,
		bool                 isDomFlottant,
		const std::vector<double>& forResetTimes,       
		const std::vector<double>& forFwdStartTimes,
		const std::vector<double>& forFwdEndTimes,
		const std::vector<double>& forFlowStartTimes,   		
		const std::vector<double>& forFlowEndTimes,	    
		const std::vector<double>& forFwdPayPeriods,	
		const std::vector<double>& forPayTimes,
		const std::vector<double>& forPayPeriods,
		const std::vector<double>& forMarginVector,
		const std::vector<double>& forNotionalVector,
		bool                 isForFlottant,
		const string&        exNotionalType,   
		const std::vector<double>& fxResetTimes,
		const std::vector<double>& fxSettlTimes,  
		const std::vector<double>& fxPayTimes,
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
	    const std::vector<double>& strikePerState,
	    int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	virtual ARM_VectorPtr HybridCallVectorial(
		const string& modelName,
		double evalTime,
		double expiryTime,
		int callPut,
		const std::vector<double>& strikesPerState,

		/// Strip of forwards FX (or equity)
		const std::vector<double>& fxExpiryTimes,
		const std::vector<double>& fxSettlementTimes,
		const std::vector<double>& fxPayTimes,
		const std::vector<double>& fxNotionals,

		/// IR Swap
		double swapResetTime,
		const std::vector<double>& fixNotionals,
		const std::vector<double>& floatNotionals,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,

		const ARM_PricingStatesPtr& states,
		ARM_PricingContext* context=NULL) const;

	// RangeAccrual function (vectorial strike version)
	virtual ARM_VectorPtr RangeAccrualVectorial(
		const string& curveName,
		double evalTime,
		double startTime,
		double endTime,
		double payTime,
		const  std::vector<double>& fixingTimes,
		int    payIndexType, 
        double payIndexTerm,
		const  string& fxModelName,
		int    irIndexType, 
		const  std::vector<double>& irIndexResetTimes,
		const  std::vector<double>& irIndexStartTimes,
		const  std::vector<double>& irIndexEndTimes,
		const  std::vector<double>& irIndexTerms,
		const  std::vector<double>& fxDownBarriers,
		const  std::vector<double>& fxUpBarriers,
		const  std::vector<double>& irDownBarriers,
		const  std::vector<double>& irUpBarriers,
		const  std::vector<double>& notionals,
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
			const std::vector<double>& resetTimes,
			const std::vector<double>& settlementTimes,
			ARM_GP_Matrix& volatilities,
			ARM_MatrixVector& correlMatrixVector,
			std::vector<double>& totalVolatilities,
			bool withDomesticIR) const;

	// Compute the volatilities used for the 1IRFX
	void ComputeVolatilities(
			const std::vector<double>& resetTimes,
			const std::vector<double>& settlementTimes,
			ARM_GP_Matrix& volatilities,
			std::vector<double>& totalVar) const;


	virtual string toString(const string& indent="",const string& nextIndent="") const;
	
	double ComputeTimeLag(double expiryTime, double payTime, bool isDomestic = true) const;
	double ComputeAbsTimeLag(double expiryTime, double payTime, bool isDomestic = true) const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
