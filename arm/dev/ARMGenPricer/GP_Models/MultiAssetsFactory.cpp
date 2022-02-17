/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file MultiAssetsFactory.cpp
 *
 *  \brief multi asset factory
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/MultiAssetsFactory.h"

#include "gpbase/singleton.h"
#include "gpmodels/argconvdefault.h"
#include "gpmodels/typedef.h"
#include "gpmodels/MultiAssets.h"
#include "gpmodels/MultiAssetsMeanReverting.h"
#include "gpmodels/2IRFXModel.h"
#include "gpmodels/1IRFXModel.h"
#include "gpmodels/HybridIRFX.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsFactoryImp
///	Routine : CreateMultiAssetsModel
///	Returns : ARM_MultiAssetsModel
///	Action  : creates the corresponding multi assets model
////////////////////////////////////////////////////

ARM_MultiAssetsModel* ARM_MultiAssetsFactoryImp::CreateMultiAssetsModel(const ARM_ModelNameMap& modelNameMap, 
																		const ARM_CurveMatrix& correlationMatrix, 
																		const string& NameStr)
{
    ARM_MultiAssetsType type = (ARM_MultiAssetsType)ARM_ArgConv_MultiAssetsType.GetNumber(NameStr);

	switch(type)
	{
	case ARM_MultiAssetsModelType::twoirfx:
		return  new ARM_2IRFXModel(modelNameMap, correlationMatrix);
		break;
	case ARM_MultiAssetsModelType::oneirfx:
		return  new ARM_1IRFXModel(modelNameMap, correlationMatrix);
		break;
	case ARM_MultiAssetsModelType::irfx:		
		return new ARM_HybridIRFX(modelNameMap, correlationMatrix);
		break;
	default:
		{
			if( ARM_MultiAssetsMeanReverting::MeanRevertingCompatible( modelNameMap ) )
				return new ARM_MultiAssetsMeanReverting(&modelNameMap, &correlationMatrix);
			else
				return new ARM_MultiAssetsModel(&modelNameMap, &correlationMatrix);
		}
	}
}

ARM_SingletonHolder<ARM_MultiAssetsFactoryImp> ARM_MultiAssetsFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

