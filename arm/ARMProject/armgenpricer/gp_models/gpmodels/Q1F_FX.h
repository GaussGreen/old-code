/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1F_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_Q1F_FX_H
#define _INGPMODELS_Q1F_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParamsQ1F_Fx.h"
#include "gpmodels/typedef.h"
#include "gpmodels/Local_Functional.h"
#include "gpmodels/typedef.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )

class ARM_QModel1F;
struct ARM_PricingContext;

class ARM_QModel1F_Fx :  public ARM_EqFxBase
{
private:
	bool itsIntegratedVersion;
	bool itsIsForexDiffusion;
	ARM_GP_VectorPtr itsPrecomputedFwds;
	ARM_PricingModelPtr itsIRDomModel;
	ARM_PricingModelPtr itsIRForModel;

    bool itsIsImpliedVolCalc;
    mutable std::vector<double> itsFwdFxVol;

	static ARM_QModel1F* GetStochasticModel( const ARM_PricingModelPtr& model );

	ARM_ModelParamsQ1F_Fx* GetQ1FFXModelParams();

    void MarkovianDriftData( size_t timeIdx,
        double& qFx, 
		double& fx0, 
		double* qLinearCoef=NULL,
		double* fx0Deriv=NULL, 
		double* dt=NULL ) const;

    static double ComputeForwardFxQ(double t,
		double T,
		const std::vector<double>& times,
		const std::vector<double>& sigmas, 
		const std::vector<double>& qs, 
		const std::vector<double>& dqs,
		double& qDrift);

	ARM_LocalFunctionalPtr itsLocalFunctional;

public:
    enum modelsAlias /// for correlation matrix
    {
        DomModel	=0,     /// the domestic stochastic IR model
        ForModel,			/// the foreign stochastic IR model
        FxModel,			/// the stochastic FX model (myself !)
    };

	ARM_QModel1F_Fx(const ARM_ZeroCurvePtr& zc, 
		ARM_ModelParamsQ1F_Fx* modelParam, 
		const ARM_CurveMatrix& correlMatrix, 
		bool IntegratedVersion = true,
		bool isSpotDiffusion = false);

	ARM_QModel1F_Fx( const ARM_QModel1F_Fx& rhs );
	ASSIGN_OPERATOR(ARM_QModel1F_Fx)
	virtual ~ARM_QModel1F_Fx();

    /// To say if forward model(s) is(are) set to generate stochastic values of forward to
    /// use an integrated diffusion version
    virtual bool IsFwdModels() const { return itsIRDomModel != ARM_PricingModelPtr(NULL ) && itsIRForModel != ARM_PricingModelPtr(NULL ); }

	/// accessors
	ARM_GP_VectorPtr GetPrecomputedFwds() const { return itsPrecomputedFwds; }
	bool GetIntegratedVersion() const { return itsIntegratedVersion; }
    void SetIntegratedVersion(bool integratedVersion) { itsIntegratedVersion=integratedVersion; }
	bool GetIsForexDiffusion() const { return itsIsForexDiffusion; }
    void SetIsForexDiffusion(bool isForexDiffusion) { itsIsForexDiffusion=isForexDiffusion; }

	/// monte carlo part
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual void SetIRForeignModel( const ARM_PricingModelPtr& irModel ); 
	virtual void SetIRDomesticModel( const ARM_PricingModelPtr& irModel );

    bool GetIsImpliedVolCalc() const { return itsIsImpliedVolCalc; }
    void SetIsImpliedVolCalc(bool isImpliedVolCalc) { itsIsImpliedVolCalc=isImpliedVolCalc; }
    const std::vector<double>& GetFwdFxVol() const { return itsFwdFxVol; }

	void SetLocalFunctional(const ARM_LocalFunctionalPtr& localFunctional) { itsLocalFunctional = localFunctional; }

	virtual double MappingFunction( double x, double x0, double q0 ) const;
	virtual double InverseMappingFunction( double fx, double fx0, double q0, bool& status ) const;

	/// By default the domestic model is used to compute a discount factor
	virtual ARM_VectorPtr DiscountFactor( 
		const string& modelName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
    { return itsIRDomModel->DiscountFactor(modelName,evalTime,maturityTime,states); }

    /// Forward function
	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

	/// Call function (vectorial strike version)
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

	/// Computation of Q and sigma of Fwd FW from Q and sigma of spot FX
	virtual ARM_VectorPtr ComputeFwdFXModelParam(
		double evalTime,
	    double expiryTime,
	    double settlementTime,
		ARM_PricingContext* context=NULL) const; 

	/// Computation of integral between 0 and evalTime of (Gamma(s,Tsettl) - Gamma(s,Texp))² ds)
	virtual ARM_VectorPtr ComputeFwdZCModelParam(
		double t1,
		double t2,
	    double T1,
	    double T2,
		ARM_PricingContext* context=NULL) const; 

	/// Default Digital call function (vectorial strike version)
	virtual ARM_VectorPtr DigitalVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const std::vector<double>& strikePerState,
		double notional,
		int callPut,
	    double payTime,
		ARM_DigitType digitType,
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// Hybriud Call function (vectorial strike version)
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
	
	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	/// three functions to compute the Markovian drift
    virtual ARM_GP_MatrixPtr IntegratedMarkovianDrift(size_t timeIdx, 
		const ARM_GP_MatrixPtr& numMethodStates,
		const ARM_GP_VectorPtr& driftCorrection ) const;
	virtual ARM_GP_MatrixPtr MarkovianDrift(size_t timeIdx, 
		const ARM_GP_MatrixPtr& numMethodStates ) const;
	// Just used by the model
	virtual void MarkovianDriftPDE(size_t timeIdx, 
		const ARM_GP_MatrixPtr& numMethodStates, ARM_GP_VectorPtr result ) const;

	virtual ARM_VectorPtr LocalDiscounts(
		size_t timeIdx, 
		double dt, 
		const ARM_PricingStatesPtr& states) const;

	virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;


	virtual bool SupportAnalyticMarginal() const;

	void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);
	virtual ARM_PricingStatesPtr Init(const string& payModelName, 
		const ARM_TimeInfoPtrVector& timeInfos);
	virtual void PostInit();

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

	/// function for the generic tree
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	void EulerLocalDrifts( 
		const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    /// Calibration purpose
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

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

	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	double VarianceFwdFx(double a, 
		double b, 
		double settlementTime,
        const ARM_QModel1F* domRefModel,
        const ARM_QModel1F* forRefModel,
        const ARM_ModelParamsQ1F* fxParams,
        bool isFwdZcFwdFxCov, 
		double& fwdZcFwdFxCov) const;

	/// ARM Support
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_QModel1F_Fx( *this ); }
	virtual string ExportShortName() const { return "LQ1FX";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
