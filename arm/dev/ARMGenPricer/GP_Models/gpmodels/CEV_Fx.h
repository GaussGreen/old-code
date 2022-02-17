/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CEV_Fx.h
 *
 *  \brief 
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2006
 */

#ifndef _INGPMODELS_CEV_FX_H
#define _INGPMODELS_CEV_FX_H

// gpbase
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

//gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingmodelequity.h"

//gpmodel
#include "gpmodels/EqFXBase.h"
#include "gpmodels/CEV_Model.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsCEV_Fx.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )

class ARM_QModel1F;

class ARM_CEVModel_Fx :  public ARM_EqFxBase
{
private:
	ARM_PricingModelPtr itsIRDomModel;
	ARM_PricingModelPtr itsIRForModel;

	static ARM_QModel1F* GetStochasticModel( const ARM_PricingModelPtr& model );

	string ComputeSettlementCalendar(const string& modelName="") const;
	double ComputeSettlementGap(const string& modelName="") const;

	ARM_ModelParamsCEV_Fx* GetCEVFXModelParams();

    double VarianceFwdFx(double a, 
		double b, 
		double settlementTime,
        const ARM_QModel1F* domRefModel,
        const ARM_QModel1F* forRefModel,
        const ARM_ModelParamsCEV_Fx* fxParams,
        bool isFwdZcFwdFxCov, 
		double& fwdZcFwdFxCov) const;

public:
    enum modelsAlias /// for correlation matrix
    {
        DomModel	=0,     /// the domestic stochastic IR model
        ForModel,			/// the foreign stochastic IR model
        FxModel,			/// the stochastic FX model (myself !)
    };

	ARM_CEVModel_Fx(
		const ARM_ZeroCurvePtr& zc, 
		ARM_ModelParamsCEV_Fx* modelParam,
		const ARM_CurveMatrix& correlMatrix);
	ARM_CEVModel_Fx( const ARM_CEVModel_Fx& rhs );
	ASSIGN_OPERATOR(ARM_CEVModel_Fx)
	virtual ~ARM_CEVModel_Fx();

    /// To say if forward model(s) is(are) set to generate stochastic values of forward to
    /// use an integrated diffusion version
    virtual bool IsFwdModels() const { return itsIRDomModel != ARM_PricingModelPtr(NULL ) && itsIRForModel != ARM_PricingModelPtr(NULL ); }

	/// ----------------- function for numerical method ----------------
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual bool NeedsToCholeskyDecomposeFactors() const {return false;}

	virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const {};

	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const {};

	virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const { return itsSettlementCalendar; }
	virtual double GetSettlementGap(const string& modelName="") const { return itsSettlementGap; }

	virtual void SetIRForeignModel( const ARM_PricingModelPtr& irModel ); 
	virtual void SetIRDomesticModel( const ARM_PricingModelPtr& irModel );

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
		const ARM_GP_Vector& strikePerState,
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
		const ARM_GP_Vector& strikePerState,
		double notional,
		int callPut,
	    double payTime,
		ARM_DigitType digitType,
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const;

	virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

	/// ARM Support
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_CEVModel_Fx( *this ); }
	virtual string ExportShortName() const { return "LQ1FX";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/