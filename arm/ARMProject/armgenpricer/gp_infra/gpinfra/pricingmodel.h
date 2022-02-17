
#ifndef _INGPINFRA_PRICINGMODEL_H
#define _INGPINFRA_PRICINGMODEL_H

/// this header has to come first
#include "gpbase/removeidentifiedwarning.h"

/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "ccy/currency.h"
#include "typedef.h"

/// to avoid include every where about the numMethod and the numeraire... due to the smart pointor choice on these!
#include "numeraire.h"
#include "nummethod.h"
#include "dates.h"
#include <string>
CC_USING_NS(std,string)

/// forward declaration
class ARM_ZeroCurve;


CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported

/// forward declaration
class ARM_TriangularMatrix;
class ARM_CalibMethod;
class ARM_ModelParams;
class ARM_MultiAssetsModel;
class ARM_NumericalModelFitter;
class ARM_DensityFunctor;

/// struct forward declaration
struct ARM_DiscretisationScheme;
struct ARM_PricerInfo;
struct ARM_ZeroCurveFunctor;
struct ARM_VanillaArg;

///////////////////////////////////////////////////////
/// \class ARM_PricingModel
/// \brief
/// This abstract class is the standard
/// interface for pricing models of the generic pricer.
/// Derived from the ARM_Model
///////////////////////////////////////////////////////
class ARM_PricingModel : public ARM_RootObject
{
private:
	/// model name
	string itsModelName;
	/// Gensec payment model name
	string itsPayModelName;

	/// the discount curve
	ARM_ZeroCurvePtr itsZeroCurve;

    /// Model parameters (vol,MRS,calibrated drift...)
	ARM_ModelParams*	            itsParams;
	ARM_DensityFunctor*				itsDensityFunctor;

    /// The numeraire,
	// the numerical method used for
    /// diffusion purpose 
    ARM_NumerairePtr            itsNumeraire;
	ARM_NumMethodPtr	        itsNumMethod;
    ARM_ZeroCurveFunctor*       itsDiscountFunctor;
    ARM_ZeroCurveFunctor*       itsFixingFunctor;

    /// Market datas container & selector keys
    ARM_MarketData_ManagerRep*  itsMktDataManager;
    ARM_StringVector			itsMDMKeys;

	/// the state local var and std dev can be shared... this is only in the case of hybrid model
	/// to avoid any problem... these data member are declared private and the multi asset model is
	/// declared friend. A better solution would have been to use ARM_GP_TensorPtr but this would have 
	/// required to rewrite the numerical method, which we do not want yet!
	ARM_MatrixVector			itsModelStateLocalVars;
	ARM_MatrixVector			itsNumMethodStateLocalVars;
	bool						itsModelStateLocalVarsIsShared;
    ARM_MatrixVector			itsModelStateLocalStdDevs;
	ARM_MatrixVector			itsNumMethodStateLocalStdDevs;
	bool						itsModelStateLocalStdDevsIsShared;
	bool                        itsNumMethodStateLocalVarsIsShared;
	ARM_GP_MatrixPtr			itsNumMethodAbsoluteDrifts;
	ARM_GP_MatrixPtr			itsNumMethodRelativeDrifts;

	/// to tell whether this is from the multi factor
	bool						itsFromMultiFactor;

    void CleanUp();
    void CopyNoCleanUp(const ARM_PricingModel& rhs);

	size_t itsModelNb;
	size_t itsModelRank;

public:
	/// constructor
	ARM_PricingModel( const ARM_ZeroCurvePtr& zc=ARM_ZeroCurvePtr(NULL), const ARM_ModelParams* params=NULL, const ARM_DensityFunctor* densityFct=NULL );
    ARM_PricingModel( const ARM_ObjectVector& marketDatas, const ARM_StringVector& mdmKeys, const ARM_ModelParams* params=NULL, const ARM_DensityFunctor* densityFct=NULL );
	ARM_PricingModel( const ARM_PricingModel& rhs);
	virtual ~ARM_PricingModel();
    ARM_PricingModel& operator = (const ARM_PricingModel& rhs);

    /// Accessors
    virtual const ARM_ModelParams* const GetModelParams() const {return itsParams;}
    virtual ARM_ModelParams* GetModelParams() {return itsParams;}
    void SetModelParams(const ARM_ModelParams& params);
	/// pure virtual to force implementation
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const=0;

    inline const ARM_NumerairePtr const GetNumeraire() const {return itsNumeraire;}
    inline ARM_NumerairePtr GetNumeraire() {return itsNumeraire;}
    virtual void SetNumeraire(const ARM_NumerairePtr& numerairePtr) { itsNumeraire=numerairePtr; }

    inline const ARM_NumMethodPtr const GetNumMethod() const {return itsNumMethod;}
    inline ARM_NumMethodPtr GetNumMethod() {return itsNumMethod;}
    virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);

	/// accessor for zero curve
	inline ARM_ZeroCurvePtr GetZeroCurve() const { return itsZeroCurve; }
	virtual void SetZeroCurve( const ARM_ZeroCurvePtr& zc);
	/// Reset the DF Map
	virtual void ResetDFMap () const {}

    inline const ARM_MarketData_ManagerRep* const GetMktDataManager() const {return itsMktDataManager;}
    inline ARM_MarketData_ManagerRep* GetMktDataManager() {return itsMktDataManager;}
    inline const ARM_StringVector& GetKeys() const {return itsMDMKeys;}
    inline void SetKeys(const ARM_StringVector& mdmKeys) {itsMDMKeys=mdmKeys;}
    void SetMktDataManager(const ARM_MarketData_ManagerRep& mktDataManager,const ARM_StringVector& mdmKeys);

    /// ================== calibration for any model========================
    /// initialise the news parameters to calibrate
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter);
    virtual void PreProcessing(ARM_ModelFitter& modelFitter)=0;
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter)=0;
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)=0;
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter)=0;
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod)=0;
	virtual bool HasClosedFormsDerivatives( ARM_ModelParamType::ParamNb paramType, size_t factorNb ){ return false; }
	virtual double PartialDerivative( const ARM_ModelParam& modelParam, size_t number, size_t factorNb,
        const ARM_VanillaArg& arg, ARM_MktTargetType targetFuncType = ARM_CalibrationTarget::PriceTarget );

	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0 ) = 0;
    ARM_PricerInfo* CreatePricerInfo() const;

    /// ================== generic pricing ========================
	/// if a model has not set its numerical method, it is only the analytical part of
	/// the model that we can use!
	bool IsCurrentlyOnlyAnalyticalModel() const;

	/// Get the times used in pricing
	virtual std::vector<double>* PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos);

    /// Schedule initialisation, pre-computed datas
    /// and numerical method initialisation
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)=0;

	/// the ReInit method is for multi-loop method to do a second intialisation lighter 
	/// than the first one ... the default is to reinit the numerical method only
	virtual ARM_PricingStatesPtr ReInit();

	/// The ReInit method is useful for BackwardForward pricing (direction change)
	virtual ARM_PricingStatesPtr ReInitLoop();

	/// Finalize
	void Finalize();

	/// nb of factors
	virtual size_t FactorCount() const;

	/// Sizeof ModelStates
	virtual size_t ModelStatesSize() const { return FactorCount();}

	/// Does model advise to use other payoffs ?
	bool GetOtherPayoffsFlag() const;

	/// method to intialize after the numerical method has been intialised
	/// using the numerical discretisation ... postinit is mainly responsible
	/// for caching model dependent data
	virtual void PostInit() {};

	// With this function the model tells to the numerical method for each 
	// dimensions of pocess simulated if they are stored as the integrated
	// process or increments.
	virtual ARM_BoolVector NeedMCIntegProcess() const;

	///for var and stddev 
	inline const ARM_MatrixVector& GetModelStateLocalVars() const	 { return itsModelStateLocalVars;}
	inline const ARM_MatrixVector& GetModelStateLocalStdDevs() const { return itsModelStateLocalStdDevs;}
	virtual void SetModelStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool shareStateLocalStdDevs = false );
	virtual void SetModelStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool shareStateLocalVars =false );
	inline const ARM_MatrixVector& GetNumMethodStateLocalVars() const	 { return itsNumMethodStateLocalVars;}
	inline const ARM_MatrixVector& GetNumMethodStateLocalStdDevs() const { return itsNumMethodStateLocalStdDevs;}
	virtual void SetNumMethodStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool isNumMethodStateLocalVarsShared = false);
	virtual void SetNumMethodStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool isNumMethodStateLocalVarsShared = false);
	virtual ARM_GP_MatrixPtr GetNumMethodAbsoluteDrifts() const { return itsNumMethodAbsoluteDrifts; };
	virtual void SetNumMethodAbsoluteDrifts( const ARM_GP_MatrixPtr& drifts ) { itsNumMethodAbsoluteDrifts = drifts; };
	virtual ARM_GP_MatrixPtr GetNumMethodRelativeDrifts() const { return itsNumMethodRelativeDrifts; };
	virtual void SetNumMethodRelativeDrifts( const ARM_GP_MatrixPtr& drifts ) { itsNumMethodRelativeDrifts = drifts; };
	virtual bool IsMeanRevertingCompatible() const { return false;}

    /// Backward/forward induction to get prices at toTime date
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);

	/// for looping method 
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const=0;

	/// Default implementation of functions to process paid payoffs (useful for change of measure  ... etc)
	virtual void ProcessPaidPayoffs(const string& payCurveName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payCurveName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;
	virtual void ProcessUnPaidPayoffs(const string& payCurveName, ARM_GP_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const;

	/// function for discretisation non const because model may cache things
	/// compute model times is responsible for giving specific model times (case of model
	///	requiring specific dates like BGM)
	virtual std::vector<double>* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos )=0;
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const=0;
	/// model Fix Time Step gives number of fix point for model supporting non analytical Marginal!
	virtual int ModelFixTimeStep( int fixTimeStep ) const;

	// Every model must answer if it supports
    /// bacward/forward induction
	virtual bool SupportBackwardInduction() const=0;
	virtual bool SupportForwardInduction() const=0;	

    virtual bool NeedLocalDiscount() const;
    virtual bool NeedArrowDebreuPrices() const { return false; }

	/// Every model must answer if it supports or not analytic Marginal
	/// allow to avoid discretizing in certain numerical scheme (case of HW model
	/// but not SFRM that needs intermediate points to diffuse)
	virtual bool SupportAnalyticMarginal() const=0;	

	// Do we need to add the model time in discretisation scheme of the numerical
	// method ?
	virtual bool NeedModelTime() const {return false;}
	// Do we need to evaluate the states at his timeIndex? 
	virtual bool NeedStatesEval(int timeIndex) const {return false;}

    /// Every model should implement a discount factor function
    /// to support its discount & fixing functors
	/// NDC : We MUST uncomment the following line 
//protected:_
    virtual ARM_VectorPtr DiscountFactor( 
											const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states
										) const = 0;
//public:
	virtual ARM_VectorPtr DiscountFactor( 
											const string& curveName,
											double evalTime, 
											const std::vector<double>&  maturityTime,
											const ARM_PricingStatesPtr& states
										) const;

	/// Functions to compute stochastic part of the integrated or instantaneous
    /// zero-coupon risk neutral drift
	virtual ARM_VectorPtr IntegratedRiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded=false ) const;
	virtual ARM_VectorPtr RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded=false ) const;


    /// Give the local discount value between startTime
    /// and endTime for each states at startTime
    virtual ARM_VectorPtr LocalDiscounts(size_t timeIdx, 
		double dt, const ARM_PricingStatesPtr& states) const;

    /// Give the local payoff for numerical method calibration
    virtual ARM_VectorPtr LocalPayoffs(size_t timeIdx, 
		double dt, const ARM_PricingStatesPtr& states) const
    { return LocalDiscounts(timeIdx,dt,states); }


    /// Give global drifts w.r.t. a given schedule
    /// Default implementation is to return NULL (meanning centred processes)
    virtual void IntegratedGlobalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& drifts) const { drifts = ARM_GP_MatrixPtr( NULL ); }

    /// Give local drifts w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    virtual void EulerLocalDrifts(
		const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// compute the markovian drift
	virtual ARM_GP_MatrixPtr IntegratedMarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, const ARM_GP_VectorPtr& driftCorrection) const;
	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;


	/// computes the vols, the vols derivatives and the correlation
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;
	
	/// Volatlities And Correlation TimeSteps: for PDE discretization
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Volatilities And Correlation TimeSteps not implemented"); }

	/// Coefficients of a PDE3D
	virtual void UpdatePDE3DCoeffs(
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

	virtual bool NeedLocalDiscounts() const;

	/// function to compute variances with full matrices
	/// the design is for j variables (Xj) with k factors at time step i (Ti)
	/// the ARM_MatrixVector (which is vector< ARM_Matrix* > contains for 
	/// localVar( Xj between Ti-1 and Ti ) = Sum( k=0 ... nbCol) (*localVariances[i])(j,k)
	/// we have both localVariances and StDev to avoid taking sqrt each time!
	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const = 0;

	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Unimplemented Method NumMethodStateGlobalVariances"); }

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const = 0;

	virtual void NumMethodStateLocalCovariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map ) {}

	/// function to cache in stdDev (for performance reason!)
	virtual void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);
	virtual void NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);
	void NumMethodsDrifts( const std::vector<double>& timeSteps);


	/// to compute the std Dev cholesky or elementwise of a matrix!
	void ComputestdDevMatrixVector( 
		ARM_MatrixVector& varCovarMatrixVector,
		ARM_MatrixVector& stdDevMatrixVector,
		bool needToCholeskyDecompose,
		bool skipFirst = false ) const;

	/// in order to get independent variables does the Cholesky decomposition if necessary
	virtual bool NeedsToCholeskyDecomposeFactors( ) const = 0;

	/// fonction to compute variances for model with the same nb of model
	///	variables and factors ... hence the used of triangular matrices
	/// the design is for j variables (Xj) with k factors at time step i (Ti)
	/// the ARM_TriangularMatrixVector (which is vector< ARM_TriangularMatrixVector* > contains for 
	/// localVar( Xj between Ti-1 and Ti ) = Sum( k=j ... nbCol) (*localVariances[i])(j,k)
	/// we have both localVariances and StDev to avoid taking sqrt each time!
    virtual void NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
        ARM_MatrixVector& localVariances,
        ARM_MatrixVector& variances ) const;

	virtual void NumMethodStateLocalGlobalVariances( const ARM_GP_Vector& timeSteps,
        ARM_MatrixVector& localVariances,
        ARM_MatrixVector& variances ) const
	{
		NumMethodStateLocalGlobalVariances(timeSteps.GetValues(), localVariances, variances);
	};

	/// version that does the variances and stdDev
    virtual void NumMethodStateLocalGlobalVariancesAndStdDev(
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances,
		ARM_MatrixVector& globlaVariances,
		ARM_MatrixVector& localStDev,
        ARM_MatrixVector& globalStdDev ) const;

    /// Give the time to reach a given variance
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const=0;

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

	/// function given a MonteCarlo
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const=0;

	/// function given a Tree
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const=0;

	/// function given a Pde
	virtual void PdeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex, double lambda) const {};

    /// ================== convention functions (ccies & dates) =============
	/// Calendar management
    /// A multi-currency model will override this function
    virtual ARM_Currency* GetCurrency( const string& modelName  ) const { return nullptr;/*GetZeroCurve()->GetCurrencyUnit();*/ }

    /// --- for basis swap discounting purpose
    void SetDiscountFunctor( ARM_ZeroCurveFunctor* newDiscountFunctor ){ itsDiscountFunctor = newDiscountFunctor; }
    ARM_ZeroCurveFunctor* GetDiscountFunctor() const{ return itsDiscountFunctor; }
    
    void SetFixingFunctor( ARM_ZeroCurveFunctor* newFixingFunctor ){ itsFixingFunctor = newFixingFunctor; }
    ARM_ZeroCurveFunctor* GetFixingFunctor() const{ return itsFixingFunctor; }

	/// conversion from date to time and vice versa
    /// For multi-currency model, there is only one asOfDate
	double GetTimeFromDate( const ARM_Date& d ) const { return d.GetJulian()-GetAsOfDate().GetJulian(); }
	ARM_Date GetDateFromTime( double d ) const { return ARM_Date( d+GetAsOfDate().GetJulian() ); }
	virtual ARM_Date GetAsOfDate() const { return ARM_Date()/*GetZeroCurve()->GetAsOfDate()*/; }

    /// ================== Standard ARM object support ========================
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_CLASS_NAME GetRootName() { return ARM_PRICINGMODEL; }

    /// ================== default functions for model ir =============
    /// ================== elementary pricing functions =============
	/// Defauft Libor provided but may be redefined by
    /// a Libor Market Model
	ARM_VectorPtr DefaultLibor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double fwdResetTime,    // for convexity adjustment...
        double payTime,         //... in derived classes
        const ARM_PricingStatesPtr& states) const;

	/// Defauft annuity provided but may be redefined
	ARM_VectorPtr DefaultAnnuity(
		const string& curveName, 
        double evalTime,
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const;

	/// Defauft annuity provided but may be redefined (With Nominal for Variable Notional Swaps)
	ARM_VectorPtr DefaultAnnuityWithNominal(
		const string& curveName, 
        double evalTime,
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fixNominal,
        const ARM_PricingStatesPtr& states) const;
	
	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	ARM_VectorPtr DefaultSwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
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
		const ARM_VectorPtr& FixedComputedAnnuity,
		const ARM_VectorPtr& FloatComputedAnnuity,
		const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr DefaultSwapRateInPlaceWithComputedAnnuityAndNominal(
		const string& curveName, 
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
		const ARM_VectorPtr& FixedComputedAnnuity,
		const std::vector<double>& FloatNotional,
		const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr DefaultNPVSwapWithComputedAnnuity(
		const string& curveName, 
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
		const ARM_VectorPtr& floatAnnuity,
		const ARM_PricingStatesPtr& states) const;

	ARM_GP_MatrixPtr DefaultNPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& PayTimes, 
		const std::vector<double>& PayPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const;

	ARM_GP_MatrixPtr DefaultNPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const;


    /// By default implied volatility
    double DefaultImpliedVol(const ARM_VanillaArg& arg) const;

	/// function for multi-assets
	inline size_t GetModelNb() const { return itsModelNb; }
	inline void SetModelNb(size_t modelNb) { itsModelNb=modelNb; }
	inline size_t GetModelRank() const { return itsModelRank; }
	inline void SetModelRank(size_t modelRank) { itsModelRank = modelRank; }
	size_t OffsetTimeIndexWithModelNb( size_t timeIdx, size_t modelNb ) const;
	double GetLocalMatrixElemWithModelNb( const ARM_MatrixVector& matrix, size_t timeIdx, size_t modelNb, size_t i, size_t j ) const;
	virtual void SetModelName(const string& modelName ) { itsModelName=modelName; }
	inline string GetModelName() const { return itsModelName; }
	virtual void SetPayModelName(const string& modelName ) { itsPayModelName=modelName; }
	inline string GetPayModelName() const { return itsPayModelName; }
	virtual void UpdateLinks( const ARM_MultiAssetsModel& map ) {}; /// default is to do nothing!
	virtual string GetModelNamePerFactor(size_t ) const { return itsModelName; }
	inline bool GetFromMultiFactor() const { return itsFromMultiFactor; }
	inline void SetFromMultiFactor(bool value) { itsFromMultiFactor=value; }

	/// Markov Functional Calibration
	virtual void setNumericalModelFitter( ARM_NumericalModelFitter * ) { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": setNumericalModelFitter not implemented!"); }
	virtual void setFwdResetStartEndDates( const ARM_GP_VectorPtr& ResetDates, const ARM_GP_VectorPtr& StartDates, const ARM_GP_VectorPtr& EndDates ) { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": setDensityResetDates not implemented!"); }
	virtual void setStorageDates( const ARM_GP_VectorPtr& StorageDates ) { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": setStorageDates not implemented!"); }
	virtual void setCalibrationStatus( bool calibrationStatus ) { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": setCalibrationStatus not implemented!"); }
	virtual ARM_DensityFunctor* GetDensityFunctor() const {return itsDensityFunctor;};//  Accessor
	virtual void UpdateDensityFunctor(double fwd, double expiryTime, double tenor=0.0) { ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": UpdateDensityFunctor is not implemented!"); };//  Update the density functor at expiryTime and tenor


	/// tells the type of the model
	virtual int GetType() const = 0;
	virtual const ARM_PricingModel* GetRefModel() const { return this; }
	virtual ARM_PricingModel* GetRefModel() { return this; }
};


CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
