/*!
 *
 * Copyright (c) CDC IXIS CM July 2005 Paris
 *
 *	\file LN Model Factory.cpp
 *
 *  \brief LN Model Factory
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


/// gpmodel
#include "gpmodels/EqFx_ModelFactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/typedef.h"

/// BS Model 
#include "gpmodels/BS_ModelParams.h"
#include "gpmodels/LN_Eq.h"
#include "gpmodels/LN_Fx.h"

/// Q1F Model 
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/Q1F_Eq.h"
#include "gpmodels/Q1F_Fx.h"

// CEV Model
#include "gpmodels/CEV_Fx.h"
#include "gpmodels/CEV_ModelParams.h"

/// Heston Model
#include "gpmodels/Heston_ModelParams.h"
#include "gpmodels/Heston_Fx.h"
#include "gpmodels/Heston_Eq.h"

/// SABR Model
#include "gpmodels/SABR_ModelParams.h"
#include "gpmodels/SABR_Fx.h"
#include "gpmodels/SABR_Eq.h"

/// Mixture Model
#include "gpmodels/Mixture_Fx.h"


/// gpbase
#include "gpbase/singleton.h"
/// gpinfra
#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_EqFx_ModelFactoryImp
///	Routine: CreateModel
///	Returns: ARM_PricingModel*
///	Action : creates the corresponding pricing Model
////////////////////////////////////////////////////

ARM_PricingModel* ARM_EqFx_ModelFactoryImp::CreateModel( 
	const ARM_ZeroCurvePtr& domesticZc, 
	const ARM_ModelParamVector& params,
	double spot, // 1foriegn Ccy = X domestic Ccy
	const ARM_ZeroCurvePtr& foreignZc, 
	const ARM_CurveMatrix& correlMatrix,
	ModelType modelType,
	int mcScheme)
{
	switch( modelType )
	{
	case BS_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsLN_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_LN_Fx( domesticZc, &modelParams,correlMatrix);
			}
			else
			{
				ARM_ModelParamsBS_Eq modelParams( params, domesticZc, spot );
				return new 	ARM_LN_Eq( domesticZc, &modelParams);
			}
		}
		break;
	case Q1F_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsQ1F_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_QModel1F_Fx( domesticZc, &modelParams,correlMatrix);
			}
			else
			{
				ARM_ModelParamsQ1F_Eq modelParams( params, domesticZc, spot );
				return new 	ARM_QModel1F_Eq( domesticZc, &modelParams);
			}
		}
		break;
	case HESTON_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsHeston_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_HestonModel_Fx( domesticZc, &modelParams,(ARM_HestonModel_Fx::MCScheme)mcScheme);
			}
			else
			{
				ARM_ModelParamsHeston_Eq modelParams( params, domesticZc, spot );
				return new 	ARM_HestonModel_Eq( domesticZc, &modelParams);
			}
		}
		break;
	case SABR_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsSABR_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_SABRModel_Fx( domesticZc, &modelParams);
			}
			else
			{
				ARM_ModelParamsSABR_Eq modelParams( params, domesticZc, spot );
				return new 	ARM_SABR_Eq( domesticZc, &modelParams);
			}
		}
		break;

	case CEV_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsCEV_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_CEVModel_Fx( domesticZc, &modelParams,correlMatrix);
			}
			else
			{				
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": there is no CEV equity model" );
			}
		}
		break;
	case Mixture_Model:
		{
			if( foreignZc != ARM_ZeroCurvePtr(NULL) )
			{
				ARM_ModelParamsMixture_Fx modelParams( params, domesticZc, foreignZc, spot );
				return new 	ARM_MixtureModel_Fx( domesticZc, &modelParams);
			}
			else
			{				
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": there is no CEV equity model" );
			}
		}
		break;

	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unknown type" );
	}
}


ARM_SingletonHolder<ARM_EqFx_ModelFactoryImp> ARM_EqFx_ModelFactory;

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

