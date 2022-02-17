/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SABR_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */


#ifndef _INGPMODELS_SABR_FX_H
#define _INGPMODELS_SABR_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpinfra/pricingmodelir.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/SABR_ModelParams.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_SABRModel_Fx :  public ARM_EqFxBase
{
private:
	ARM_PricingModelIR* itsIRModel;
	ARM_PricingModelPtr itsIRModelPtr;

	double itsIRSpotCorrel;
	double itsIRVolCorrel;
	double itsSpotVolCorrel;

	size_t itsIntegrationStep;
	int itsTypeOfModel;

public:
	typedef ARM_ModelParams_Fx_T<ARM_SABR_ModelParams> ARM_ModelParamsSABR_Fx;

	ARM_SABRModel_Fx(const ARM_ZeroCurvePtr& zc,
		ARM_ModelParamsSABR_Fx* modelParam);
	
	ARM_SABRModel_Fx( const ARM_SABRModel_Fx& rhs )
	:	ARM_EqFxBase(rhs),
		itsIntegrationStep(rhs.itsIntegrationStep),
	    itsTypeOfModel(rhs.itsTypeOfModel)
	{}

	ASSIGN_OPERATOR(ARM_SABRModel_Fx)
	virtual ~ARM_SABRModel_Fx(){};

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

	//  Update the density functor at expiryTime
	virtual void UpdateDensityFunctor(double fwd, double expiryTime);

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;
	
	/// to compute implied volatility
	virtual double ImpliedVol(const ARM_VanillaArg& arg) const;

	 /// To calculate Implied Volatility partial derivatives
    /// in ralation to each ModelParam
	virtual bool HasClosedFormsDerivatives( ARM_ModelParamType::ParamNb paramType, size_t factorNb ){ return true; }
	virtual double PartialDerivative( const ARM_ModelParam& modelParam, size_t number, 
        size_t factorNb,const ARM_VanillaArg& arg,
        ARM_MktTargetType targetFuncType = ARM_CalibrationTarget::PriceTarget);
    

	/// function for the monte carlo
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_SABRModel_Fx( *this ); }
	virtual string ExportShortName() const { return "LSFXM";}

	virtual size_t FactorCount() const { return 2; }

	/// Accessors
	inline ARM_PricingModelIR* getIRModel() const { return itsIRModel; }
	inline void setIRModel( ARM_PricingModelIR* IRModel ) { itsIRModel = IRModel;}
	inline void setIRModel( const ARM_PricingModelPtr& IRModel )
		{ itsIRModelPtr = IRModel; itsIRModel = dynamic_cast<ARM_PricingModelIR*> (&*IRModel); if (!itsIRModel) ARM_THROW( ERR_INVALID_ARGUMENT, " SABR_Eq Model is supposed to work with IRModels. " ); }

	/// To keep links with its IRModel
	//virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
