/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file MultiAssets.h
 *
 *  \brief
 *
 *  \brief multi asset model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */



#ifndef _INGPMODELS_MULTIASSETS_H
#define _INGPMODELS_MULTIASSETS_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/curvematrix.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingmodelequity.h"
#include "gpinfra/pricingfunctioninflation.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelNameMap;
struct ARM_PricingContext;

class ARM_MultiAssetsModel :	public ARM_PricingModelIR,
								public ARM_PricingFunctionEquity,
								public ARM_PricingFuncInflation
{
private :
    ARM_ModelNameMap* itsModelMap;
	ARM_CurveMatrix* itsCorrelMatrix;
	ARM_PricingModel*	itsRefModel;

	typedef void (ARM_PricingModel::*ComputeDriftFunc)(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// common function for Euler and Integrated drift!
	void ComputeDriftCommon( 
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts,
		size_t rowSize,
		const ComputeDriftFunc& func ) const;

	void UpdateSubModelLinks();
	void SetModelParamsVec();
	void SetMultiFactorFlagOnModel();

protected:
	/// to give access only to derived classes!
    virtual ARM_PricingStatesPtr BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);
    virtual ARM_PricingStatesPtr ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);

public:

    ARM_MultiAssetsModel(
		const ARM_ModelNameMap*	modelNameMap	= NULL,
		const ARM_GP_Matrix* correlationMatrix	= NULL );

	ARM_MultiAssetsModel(
		const ARM_ModelNameMap*	modelNameMap,
		const ARM_CurveMatrix* correlationMatrix);

	ARM_MultiAssetsModel(const ARM_MultiAssetsModel& rhs);
	ASSIGN_OPERATOR(ARM_MultiAssetsModel)
	virtual ~ARM_MultiAssetsModel();

	void InitMultiAsset();

    /// -------------------------------------------
	/// pricing function
    virtual ARM_VectorPtr DiscountFactor( 
		const string& modelName,
        double evalTime, 
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    /// Forward is overiden to use the linked Forward Forex model
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



    /// Specialized IR model functions
    /// ------------------------------
    /// Libor is overriden to allow reference model payment lag adjustment 
    virtual ARM_VectorPtr Libor( 
		const string& modelName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const;

    /// SwapRate is overriden to allow different fixing & discount fctors
    virtual ARM_VectorPtr SwapRate(
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
		const ARM_PricingStatesPtr& states) const;

	/// Spread is overriden to allow different fixing & discount fctors
    virtual ARM_VectorPtr Spread(
		const string& curveName, 
        double evalTime,
		double coeff1,
		double floatStartTime1,
        double floatEndTime1, 
		const std::vector<double>& fixPayTimes1,
        const std::vector<double>& fixPayPeriods1,
		const std::vector<double>& fwdStartTimes1,
        const std::vector<double>& fwdEndTimes1,
        const std::vector<double>& fwdPayPeriods1, 
		const std::vector<double>& floatPayTimes1,
        const std::vector<double>& floatPayPeriods1,
        const std::vector<double>& margin1,
		double coeff2,
        double floatStartTime2,
        double floatEndTime2, 
		const std::vector<double>& fixPayTimes2,
        const std::vector<double>& fixPayPeriods2,
		const std::vector<double>& fwdStartTimes2,
        const std::vector<double>& fwdEndTimes2,
        const std::vector<double>& fwdPayPeriods2, 
		const std::vector<double>& floatPayTimes2,
        const std::vector<double>& floatPayPeriods2,
        const std::vector<double>& margin2,
        const ARM_PricingStatesPtr& states) const;

    /// NPVSwap is overriden to allow different fixing & discount fctors
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

	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& PayTimes, 
		const std::vector<double>& PayPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const ;

	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const;

///	Action  : Default Swap Rate computation with null strike and exchanging notional flow by flow
	virtual ARM_VectorPtr NPVBasisSwap( 
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

    /// Vanilla caplet/floorlet formula is overiden to allow different fixing & discount fctors
    virtual ARM_VectorPtr VanillaCaplet(
		const string& modelName, 
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

    /// Vanilla digital caplet/floorlet is overiden to allow different fixing & discount fctors
    virtual ARM_VectorPtr VanillaDigital(
		const string& modelName, 
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

    /// Vanilla corridorlet formula is overiden to allow different fixing & discount fctors
    virtual ARM_VectorPtr VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     indexPaymentType,
        double  fwdPaymentPeriod,
        const std::vector<double>& refIdxResetTimes,
        const std::vector<double>& refIdxStartTimes,
        const std::vector<double>& refIdxEndTimes,
        const std::vector<double>& refFwdPeriods,
        const std::vector<double>& refIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;

    /// Vanilla swaption formula is overiden to allow different fixing & discount fctors
    virtual ARM_VectorPtr VanillaSwaption(
		const string& modelName,
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

    /// Specialized Inflation model functions
    /// -------------------------------------
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

/// YoY/OAT Swaps and SwapRates
 	virtual ARM_GP_VectorPtr YoYSwapRate(const string& irCurveName, const string& infCurveName, 
														double evalTime, const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														const ARM_DateStripPtr& fixedDateStrip,
														double itsSpread,
														const ARM_PricingStatesPtr& states) const;

 	virtual ARM_GP_VectorPtr YoYSwap(const string& irCurveName, const string& infCurveName, 
														double evalTime, double Strike, double FloatMargin, 
														const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														const ARM_DateStripPtr& fixedDateStrip,
														double itsSpread,
														const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr OATSwapRate(const string& irCurveName, const string& infCurveName, 
														double evalTime, const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														const ARM_DateStripPtr& fixedDateStrip,
														double itsCoupon,
														const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr OATSwap(const string& irCurveName, const string& infCurveName, 
														double evalTime, double Strike, double FloatMargin, 
														const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														const ARM_DateStripPtr& fixedDateStrip,
														double itsCoupon,
														const ARM_PricingStatesPtr& states) const;


	virtual ARM_GP_VectorPtr YoYCapFloor( const string& irCurveName,	const string& infCurveName, 
														double evalTime,
														double Strike,
														double FloatMargin, 
														int CapFloor,
														const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														double itsSpread,
														const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr OATCapFloor( const string& irCurveName, const string& infCurveName, 
														double evalTime,
														double Strike,
														double FloatMargin, 
														int CapFloor,
														const ARM_DateStripPtr& numDateStrip,
														const ARM_DateStripPtr& denomDateStrip,
														double itsSpread,
														const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr ZCCap( ) const;

	virtual void ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;
	void InitModelNb();
    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	virtual void SetModelStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool shareStateLocalStdDevs = false );
	virtual void SetModelStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool shareStateLocalVars =false );
	virtual void SetNumMethodStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool isNumMethodStateLocalVarsShared = false);
	virtual void SetNumMethodStateLocalStdDevs( const ARM_MatrixVector& stateLocalVars,bool isNumMethodStateLocalVarsShared = false);
    
	virtual void SetNumeraire(const ARM_NumerairePtr& numerairePtr); 
    virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);

    virtual const ARM_PricingModel* GetRefModel() const;
    virtual ARM_PricingModel* GetRefModel();
	virtual void SetRefModel( ARM_PricingModel* refModel ) { itsRefModel=refModel; }

	/// nb of factors
	virtual size_t FactorCount() const;

	// nb of ModelStates (sum of all modelstates for all underlying models)
	virtual size_t ModelStatesSize() const;

    /// Multi-currency retriever
	virtual ARM_Currency* GetCurrency( const string& modelName  ) const;

    /// Virtual fcts are overiden to call the reference model...
    ///...for initialisation if multi-loop pricing is used
	virtual ARM_PricingStatesPtr ReInit();

    ///...for induction
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);

    
	/// Pure virtual fcts calling the reference model...
    ///...for calibration
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter);
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

    ///...for pricing process   
	virtual void PostInit();
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual std::vector<double>* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;
	virtual bool SupportBackwardInduction() const;
	virtual bool SupportForwardInduction() const;
	virtual bool SupportAnalyticMarginal() const;
	virtual bool NeedArrowDebreuPrices() const;

    ///...for pricing computation
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;
	
	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const;

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual ARM_VectorPtr LocalDiscounts( size_t timeIdx, double dt, 
		const ARM_PricingStatesPtr& states) const;

	/// function for the markovian drift (no need for integrated covariance)
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// relative and absolute drifts
    virtual void EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

	virtual ARM_BoolVector NeedMCIntegProcess() const;


    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const;

	virtual double UnderlyingCorrelation(  string	underlyingType,
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

	virtual double UnderlyingCovariance(   string	underlyingType,
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

	// Vanilla Corridorlet function
	virtual ARM_VectorPtr VanillaCMSCorridorlet(
		const string& curveName,
		double evalTime,
		double payTime,
		double resetTime,
		double startTime,
		double endTime,
		const std::vector<double>& refIdxResettimes,
		const std::vector<double>& refIndexWeights,
		const std::vector<double>& coeff1,
		const ARM_SwapRatePtrVector& firstIndex,
		const std::vector<double>& coeff2,
		const ARM_SwapRatePtrVector& secondIndex,
		int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
		double	coupon,					/// in case of a fixed payment (K_FIXED)
		const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
		double  payIndexLeverage,
		const std::vector<double>& downBarriers,
        const std::vector<double>& upBarriers,
        double  payNotional,
        int     rcvPay,
		const ARM_SwapRatePtrVector& thirdIndex, // 3rd index for double condition
		const std::vector<double>& downBarriers3,
		const std::vector<double>& upBarriers3,
        const   ARM_PricingStatesPtr& states) const;

	// RangeAccrual function : FX and Libor double condition.
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

	// Double digital condition
	virtual ARM_VectorPtr DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const std::vector<double>& firstStrikeDown,
        const std::vector<double>& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const std::vector<double>& secondStrikeDown,
        const std::vector<double>& secondStrikeUp,
		double secondStrikeSpread,
        const ARM_PricingStatesPtr& states) const;

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    /// Always validated because no specific model parameters
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const;
	virtual ARM_Date GetAsOfDate() const;

	/// accessors for the mean reverting multi asset
	ARM_ModelNameMap* GetModelMap() { return itsModelMap; }
	const ARM_ModelNameMap* const GetModelMap() const { return itsModelMap; }

	ARM_CurveMatrix* GetCorrelMatrix() { return itsCorrelMatrix; }
	const ARM_CurveMatrix* const GetCorrelMatrix() const { return itsCorrelMatrix; }
	const ARM_GP_MatrixPtr GetCorrelSubMatrix( const string& ModelName1, const string& ModelName2 ) const;

	/// Set ModelMap
    virtual void SetModelMap(ARM_ModelNameMap* RefModelMap);
	virtual void SetModelMapNoClone(ARM_ModelNameMap* RefModelMap);

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
	virtual string GetModelNamePerFactor(size_t ) const;
	virtual void SetPayModelName(const string& modelName );

	/// Volatilities and correlation time steps for PDEs
	ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;

	virtual double ImpliedVol(const ARM_VanillaArg& arg) const ;

	/// typing of the model
	virtual int GetType() const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_MultiAssetsModel(*this);}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LMAMO";}	
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
