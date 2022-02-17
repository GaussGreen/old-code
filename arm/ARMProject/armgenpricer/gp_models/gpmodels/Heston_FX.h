/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */


#ifndef _INGPMODELS_HESTON_FX_H
#define _INGPMODELS_HESTON_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpinfra/pricingmodelir.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsHeston_Fx.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/Local_Functional.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_HestonModel_Fx :  public ARM_EqFxBase
{
public:
	enum MCScheme
	{
		Euler,
		Andreasen,
		Unknown
	};
private:
	ARM_PricingModelPtr itsIRDomModel;
	ARM_PricingModelPtr itsIRForModel;
	ARM_GP_Matrix itsCorrelMatrix;
	MCScheme itsMCScheme;

	ARM_LocalFunctionalPtr itsLocalFunctional;

	static ARM_QModel1F* GetStochasticModel( const ARM_PricingModelPtr& model );

	string ComputeSettlementCalendar(const string& modelName="") const;
	double ComputeSettlementGap(const string& modelName="") const;

	ARM_ModelParamsHeston_Fx* GetHestonFXModelParams();

public:
	enum modelsAlias /// for correlation matrix
    {
        DomModel	=0,     /// the domestic stochastic IR model
        ForModel,			/// the foreign stochastic IR model
        FxModel,			/// the stochastic FX model (myself !)
    };

	enum modelsAliasHeston /// for the model states
	{
		SpotHeston = 1,			/// the spot heston
		VarHeston				/// the var heston
	};

	ARM_HestonModel_Fx(
		const ARM_ZeroCurvePtr& zc,
		ARM_ModelParamsHeston_Fx* modelParam,
		MCScheme itsMCScheme);
	
	ARM_HestonModel_Fx( const ARM_HestonModel_Fx& rhs )
	:	ARM_EqFxBase(rhs),
		itsIRDomModel(NULL),
		itsIRForModel(NULL),
		itsMCScheme(rhs.itsMCScheme),
		itsLocalFunctional(rhs.itsLocalFunctional)
	{}
	ASSIGN_OPERATOR(ARM_HestonModel_Fx)
	virtual ~ARM_HestonModel_Fx(){};

	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

	// Call function (vectorial strike version)
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

	//  Update the density functor at expiryTime
	virtual void UpdateDensityFunctor(double fwd, double expiryTime);

	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual bool NeedsToCholeskyDecomposeFactors() const {return false;}

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const {};

	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const {};

	virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void SetIRForeignModel( const ARM_PricingModelPtr& irModel ); 
	virtual void SetIRDomesticModel( const ARM_PricingModelPtr& irModel );

	virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	/// multi-factor support
    virtual void SetCorrelMatrix( ARM_GP_Matrix& correlMatrix ) { itsCorrelMatrix = correlMatrix; }

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

	double VarianceFwdFx(double a, 
		double b, 
		double settlementTime,
        const ARM_QModel1F* domRefModel,
        const ARM_QModel1F* forRefModel,
        const ARM_Heston_ModelParams* fxParams,
        bool isFwdZcFwdFxCov, 
		double& fwdZcFwdFxCov) const;
	

	double Proba(
		double settlementTime,
		double strike) const;

	void Probas(
		double settlementTime,
		const std::vector<double>& strikes,
		std::vector<double>& probasVec) const;

	double VolATM(double settlementTime) const;

	void ARM_HestonModel_Fx::HestonParameter(
		const std::vector<double>& resetTimes,
		std::vector<double>& fwds,
		std::vector<double>& volATMs,
		std::vector<double>& times,
		std::vector<double>& levels,
		std::vector<double>& kappas,
		std::vector<double>& v0s,
		std::vector<double>& thetas,
		std::vector<double>& rhos,
		std::vector<double>& nus,
		std::vector<double>& shifts,
		std::vector<double>& sigmas
		) const;

	double Quantile(
		double settlementTime,
		double proba) const;

	void SetLocalFunctional(const ARM_LocalFunctionalPtr& localFunctional) { itsLocalFunctional = localFunctional; }

	/// ARM Support
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_HestonModel_Fx( *this ); }
	virtual string ExportShortName() const { return "LHFXM";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
