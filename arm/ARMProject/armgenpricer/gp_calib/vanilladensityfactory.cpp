/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file vanilladensityfactory.h
 *  \brief factory class for vanilla density
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date March 2007
 */


#include "gpcalib/vanilladensityfactory.h"
#include "gpcalib/densityfunctors.h"
#include "gpcalib/vanillasecuritydensity.h"

#include "gpinfra/pricingmodel.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Mixture_Fx.h"
#include "gpbase/singleton.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_VanillaDensityFactorImp
///	Routine: CreateVanillaDenstityFactory
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_VanillaSecurityDensity* ARM_VanillaDensityFactorImp::CreateVanillaDenstityFactory(
		const ARM_PricingModel&  model,
		const ARM_Date&  expiryDate,
		bool isDirect)
{
	ARM_VanillaSecurityDensity* vanillaDensity = NULL;
	double expiryTime = model.GetTimeFromDate(expiryDate);
	ARM_EqFxBase* eqfxmodel = dynamic_cast<ARM_EqFxBase*> (const_cast<ARM_PricingModel*>(&model));
	if(eqfxmodel)
	{
		ARM_DensityFunctor* density = NULL;
		ARM_ModelParams_Fx*  modelParams = dynamic_cast<ARM_ModelParams_Fx*>(eqfxmodel->GetModelParams());
		ARM_ZeroCurvePtr domZc = modelParams->GetDomCurve();
		ARM_ZeroCurvePtr forZc = modelParams->GetForCurve();
		double spot = modelParams->GetSpot();

		double fwd = eqfxmodel->ComputeFwdAtTime(expiryTime) ;
		ARM_MixtureModel_Fx* mixtureModel = dynamic_cast<ARM_MixtureModel_Fx*>(eqfxmodel);
		if(mixtureModel)
		{
			double volATM = model.GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
			double decVol = model.GetModelParams()->GetModelParam(ARM_ModelParamType::Smile).GetValue(expiryTime);
			double alpha = model.GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(expiryTime);
			double lambda = model.GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(expiryTime);

			density = new ARM_MixtureDensityFunctor(fwd, expiryTime/K_YEAR_LEN, volATM, decVol, alpha, lambda,isDirect);
		}
		vanillaDensity = new ARM_VanillaSecurityDensityFX( expiryDate.GetJulian(),ARM_DensityFunctorPtr(density),domZc,forZc,spot);
	}

	return vanillaDensity;
}

ARM_SingletonHolder<ARM_VanillaDensityFactorImp> ARM_VanillaDensityFactor;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
