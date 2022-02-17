/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file EqFxBase.h
 *
 *  \brief 
 *
 *	\author  R.Guillemot
 *	\version 1.0
 *	\date May 2006
 */


#ifndef _INGPMODELS_EQFXBASE_H
#define _INGPMODELS_EQFXBASE_H


/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/curvetypedef.h"
#include "gpbase/typedef.h"
#include "gpbase/curvematrix.h"

#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingmodelequity.h"

#include "gpcalib/densityfunctors.h"
#include "gpcalib/typedef.h"

#include "gpclosedforms/gaussian_integrals.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_EqFxBase :	public ARM_PricingModel,
						public ARM_PricingFunctionEquity
{
public:
	enum CallType
	{
		ClosedFormula =0,
		Quanto,
		InvClosedFormula,
		InvQuanto
	};
	typedef ARM_EqFxBase::CallType  ARM_CallType; 

protected:
	string itsSettlementCalendar;
	double itsSettlementGap;
	ARM_CurveMatrix itsCorrelMatrix;
	ARM_CallType itsCallType;
	ARM_PricingModelIR* itsConvAdjustModel;
	ARM_Curve itsRho; // for payment lag convexity adjustment

public:
	ARM_EqFxBase( const ARM_ZeroCurvePtr& zc, ARM_ModelParams* modelParam,		
		const ARM_CurveMatrix& correlMatrix = ARM_CurveMatrix(ARM_GP_Matrix(NULL)),
		ARM_DensityFunctor* densityFct = NULL,
		ARM_PricingModelIR* convAdjustModel = NULL);
	ARM_EqFxBase( const ARM_EqFxBase& rhs );
	ASSIGN_OPERATOR(ARM_EqFxBase)
	virtual ~ARM_EqFxBase() {}; 

	// Initialize the model
	void Init();

	/// convention support : calendar + gap support
	virtual string ComputeSettlementCalendar(const string& modelName="") const;
	virtual double ComputeSettlementGap(const string& modelName="") const;

	string GetSettlementCalendar(const string& modelName="") const { return itsSettlementCalendar; }
	double GetSettlementGap(const string& modelName="") const { return itsSettlementGap; }

	virtual double ComputeFwdAtTime( double evalTime  ) const;

	/// virtual pure Call price without discounting: E[(S-K)^{+}]
	virtual double CallPrice(double Fwd,double Strike,double Expiry, int callPut) const {return 0;};

	/// NormalDistribution function: ie S(fwd,T)=Normaldistribution(X) with X~N(0,1)
	virtual ARM_GP_VectorPtr NormalDistribution(GaussLegendre_Coefficients glc, double forward, double expiryTime) const;

	//  Normal(invNormal)Mapping of the forward at expiryTime in a Black world of volatility vol 
	virtual ARM_VectorPtr Black_Forward_Mapping(
		double	fwd,
		double	vol,
		double	expiryTime,
		GaussLegendre_Coefficients glc,
		bool IsInv) const;

	//  Normal(invNormal)Mapping of the forward at expiryTime in a smiled world  
	virtual ARM_GP_Matrix Forward_Mapping(
		double	fwd,
		double	expiryTime,
		const ARM_GP_Vector& glparams) const;

	// Discount Factor
	virtual ARM_VectorPtr DiscountFactor( 
		const string& modelName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

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
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const;

	/// Default Digital call function (vectorial strike version)
	virtual ARM_VectorPtr DigitalVectorial(
		const string& modelName,
		double evalTime,
		double expiryTime,
		double settlementTime,
		const ARM_GP_Vector& strikePerState,
		double notional,
		int callPut,
		double payTime,
		ARM_DigitType digitType,
		double epsilon,
		const ARM_PricingStatesPtr& states,
		ARM_PricingContext* context) const;

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(const string& modelName,
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

	// RangeAccrual function (vectorial strike version)
	virtual ARM_VectorPtr RangeAccrualVectorial(
		const string& model1Name,
		double evalTime,
		double startTime,
		double endTime,
		const ARM_GP_Vector& fixingTimes,
		double payTime,
		const ARM_GP_Vector& downBarrierVect,
		const ARM_GP_Vector& upBarrierVect,
		const ARM_GP_Vector& notionalVect,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// ----------------- function for numerical method ----------------
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	/// Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

    virtual void NumMethodStateGlobalVariances(
        const ARM_GP_Vector& timeSteps,
        ARM_MatrixVector& variances) const;

    virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos ) {return NULL;}
// FIXMEFRED: mig.vc8 (22/05/2007 18:06:27):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0;}

	/// function for the generic tree
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// multi-factor support
    void SetCorrelMatrix( ARM_CurveMatrix& correlMatrix ) { itsCorrelMatrix = correlMatrix; }
	const ARM_CurveMatrix& GetCorrelMatrix() const { return itsCorrelMatrix; }

	/// relative and absolute drifts
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;
	virtual ARM_VectorPtr LocalDiscounts(
		size_t timeIdx, 
		double dt, 
		const ARM_PricingStatesPtr& states) const;

	/// --------- numerical method bool flags + validation
    virtual bool SupportBackwardInduction() const {	return true; }
    virtual bool SupportForwardInduction()  const {	return true; }
	virtual bool SupportAnalyticMarginal()  const {	return true;}
    virtual bool NeedArrowDebreuPrices() const { return true; }

    /// Default initialisation of the model
	virtual void TreeStatesToModelStates( ARM_PricingStatesPtr& states, int timeIndex ) const;
	virtual void MCModelStatesFromToNextTime( ARM_PricingStatesPtr& states,int timeIndex ) const;
    virtual bool NeedsToCholeskyDecomposeFactors() const {return false;}

    /// Calibration purpose
    virtual void Re_InitialiseCalibParams( ARM_ModelFitter& modelFitter ){};
    virtual void PreProcessing( ARM_ModelFitter& modelFitter ) {}
    virtual void PostProcessing( const ARM_ModelFitter& modelFitter ) {}
    virtual void AdviseCurrentCalibSecIndex( size_t index,ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib( ARM_ModelFitter& modelFitter ){};
	virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	 /// Convexity adjustment purpose
	virtual double PaymentLagAdjst(
		double expiryTime,
		double settlementTime,
		double payTime) const; 

	virtual double ImpliedVol(const ARM_VanillaArg& arg) const ;

	inline int GetType() const { return MT_FX_MODEL; };

	//accessors
	inline ARM_CallType GetCallType() const {return itsCallType;};
	inline void SetCallType( const ARM_CallType& calltype ) { itsCallType = calltype; }

	inline ARM_PricingModelIR* GetConvAdjustModel() const {return itsConvAdjustModel;};
	inline void SetConvAdjustModel( ARM_PricingModelIR* convAdjustModel ) { itsConvAdjustModel = convAdjustModel; };

	inline ARM_Curve GetRho() const {return itsRho;};
	inline void SetRho( const ARM_Curve& rho ) { itsRho = rho; };

	/// Standard ARM object support
	virtual ARM_Object* Clone()  const { return new ARM_EqFxBase(*this); };
	virtual string ExportShortName() const { return "LFX";}

};

CC_END_NAMESPACE()

#endif
