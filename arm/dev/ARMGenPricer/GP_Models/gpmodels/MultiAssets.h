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
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// common function for Euler and Integrated drift!
	void ComputeDriftCommon( 
		const ARM_GP_Vector& timeSteps,
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
		const ARM_GP_Vector& strikePerState,
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
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes,
        const ARM_GP_Vector& fwdEndTimes,
        const ARM_GP_Vector& fwdPayPeriods,
		const ARM_GP_Vector& floatPayTimes,
        const ARM_GP_Vector& floatPayPeriods,
        const ARM_GP_Vector& margin,
        bool isDbleNotional,
		const ARM_PricingStatesPtr& states) const;

	/// Spread is overriden to allow different fixing & discount fctors
    virtual ARM_VectorPtr Spread(
		const string& curveName, 
        double evalTime,
		double coeff1,
		double floatStartTime1,
        double floatEndTime1, 
		const ARM_GP_Vector& fixPayTimes1,
        const ARM_GP_Vector& fixPayPeriods1,
		const ARM_GP_Vector& fwdStartTimes1,
        const ARM_GP_Vector& fwdEndTimes1,
        const ARM_GP_Vector& fwdPayPeriods1, 
		const ARM_GP_Vector& floatPayTimes1,
        const ARM_GP_Vector& floatPayPeriods1,
        const ARM_GP_Vector& margin1,
		double coeff2,
        double floatStartTime2,
        double floatEndTime2, 
		const ARM_GP_Vector& fixPayTimes2,
        const ARM_GP_Vector& fixPayPeriods2,
		const ARM_GP_Vector& fwdStartTimes2,
        const ARM_GP_Vector& fwdEndTimes2,
        const ARM_GP_Vector& fwdPayPeriods2, 
		const ARM_GP_Vector& floatPayTimes2,
        const ARM_GP_Vector& floatPayPeriods2,
        const ARM_GP_Vector& margin2,
        const ARM_PricingStatesPtr& states) const;

    /// NPVSwap is overriden to allow different fixing & discount fctors
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

	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& PayTimes, 
		const ARM_GP_Vector& PayPeriods, 
		const ARM_GP_Vector& margin, 
		const ARM_GP_Vector& notional, 
		const ARM_PricingStatesPtr& states) const ;

	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
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
        const ARM_GP_Vector& strikesPerState,
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
        const ARM_GP_Vector& strikesPerState,
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
        const ARM_GP_Vector& refIdxResetTimes,
        const ARM_GP_Vector& refIdxStartTimes,
        const ARM_GP_Vector& refIdxEndTimes,
        const ARM_GP_Vector& refFwdPeriods,
        const ARM_GP_Vector& refIndexWeight,
        double  couponMargin,
        const vector<const ARM_GP_Vector*> downBarrierPerState,
        const vector<const ARM_GP_Vector*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;

    /// Vanilla swaption formula is overiden to allow different fixing & discount fctors
    virtual ARM_VectorPtr VanillaSwaption(
		const string& modelName,
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
		const ARM_GP_Vector& fixTimes,
        const ARM_GP_Vector& fixPayPeriods,
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
	virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;
	virtual bool SupportBackwardInduction() const;
	virtual bool SupportForwardInduction() const;
	virtual bool SupportAnalyticMarginal() const;
	virtual bool NeedArrowDebreuPrices() const;

    ///...for pricing computation
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;
	
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual ARM_VectorPtr LocalDiscounts( size_t timeIdx, double dt, 
		const ARM_PricingStatesPtr& states) const;

	/// function for the markovian drift (no need for integrated covariance)
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// relative and absolute drifts
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
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

	// Vanilla Corridorlet function
	virtual ARM_VectorPtr VanillaCMSCorridorlet(
		const string& curveName,
		double evalTime,
		double payTime,
		double resetTime,
		double startTime,
		double endTime,
		const ARM_GP_Vector& refIdxResettimes,
		const ARM_GP_Vector& refIndexWeights,
		const ARM_GP_Vector& coeff1,
		const ARM_SwapRatePtrVector& firstIndex,
		const ARM_GP_Vector& coeff2,
		const ARM_SwapRatePtrVector& secondIndex,
		int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
		double	coupon,					/// in case of a fixed payment (K_FIXED)
		const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
		double  payIndexLeverage,
		const ARM_GP_Vector& downBarriers,
        const ARM_GP_Vector& upBarriers,
        double  payNotional,
        int     rcvPay,
		const ARM_SwapRatePtrVector& thirdIndex, // 3rd index for double condition
		const ARM_GP_Vector& downBarriers3,
		const ARM_GP_Vector& upBarriers3,
        const   ARM_PricingStatesPtr& states) const;

	// RangeAccrual function : FX and Libor double condition.
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

	// Double digital condition
	virtual ARM_VectorPtr DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const ARM_GP_Vector& firstStrikeDown,
        const ARM_GP_Vector& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const ARM_GP_Vector& secondStrikeDown,
        const ARM_GP_Vector& secondStrikeUp,
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
