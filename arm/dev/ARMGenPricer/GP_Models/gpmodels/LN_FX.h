/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file LN_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_LN_FX_H
#define _INGPMODELS_LN_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpinfra/pricingmodelir.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParamsLN_Fx.h"
#include "gpmodels/typedef.h"


/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_LN_Fx :  public ARM_EqFxBase
{
private:
	ARM_PricingModelPtr itsIRDomModel;
	ARM_PricingModelPtr itsIRForModel;

public:
	enum modelsAlias /// for correlation matrix
    {
        DomModel	=0,     /// the domestic stochastic IR model
        ForModel,			/// the foreign stochastic IR model
        FxModel,			/// the stochastic FX model (myself !)
    };

	ARM_LN_Fx(
		const ARM_ZeroCurvePtr& zc, 
		ARM_ModelParamsLN_Fx* modelParam,
		const ARM_CurveMatrix& correlMatrix);
	ARM_LN_Fx( const ARM_LN_Fx& rhs );
	ASSIGN_OPERATOR(ARM_LN_Fx)
	virtual ~ARM_LN_Fx(){};

	/// To say if forward model(s) is(are) set to generate stochastic values of forward to
    /// use an integrated diffusion version
    virtual bool IsFwdModels() const { return itsIRDomModel != ARM_PricingModelPtr(NULL ) && itsIRForModel != ARM_PricingModelPtr(NULL ); }

		/// multi Factor support
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );

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
        ARM_PricingContext* context=NULL) const;

	/// ----------------- function for numerical method ----------------
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
    
	/// Give local drifts and variances w.r.t. a given schedule
	virtual ARM_GP_MatrixPtr MarkovianDrift(size_t timeIdx, 
		const ARM_GP_MatrixPtr& numMethodStates ) const;
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual ARM_VectorPtr LocalDiscounts(
		size_t timeIdx, 
		double dt, 
		const ARM_PricingStatesPtr& states) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// function for the monte carlo
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void TreeStatesToModelStates( ARM_PricingStatesPtr& states, int timeIndex ) const;

    virtual void SetFactorCount(size_t fc);

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_LN_Fx( *this ); }
	virtual string ExportShortName() const { return "LLNFX";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
