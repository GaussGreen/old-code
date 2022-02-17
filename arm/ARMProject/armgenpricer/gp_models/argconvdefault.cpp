/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date October 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/argconvdefault.h"
#include "gpmodels/enummodel.h"
#include "gpinfra/pricingstates.h"
#include "gpmodels/marketirmodel.h"
#include "gpmodels/smiledfrm.h"
#include "gpmodels/smiledfrmfactory.h"
#include "gpmodels/hwsv.h"
#include "gpmodels/Heston_Fx.h"



CC_BEGIN_NAMESPACE( ARM )

ARGConvTable PricingModelTable[] =
{	
    /// Type Name       /// number
    "HWM1F",      ARM_PricingModelType::HWM1F,
    "HWM2F",      ARM_PricingModelType::HWM2F,
    "QGM1F",      ARM_PricingModelType::QGM1F,
    "SFRM1F",     ARM_PricingModelType::SFRM1F,
    "SFRM2F",     ARM_PricingModelType::SFRM2F,
    "QM",         ARM_PricingModelType::QM,
    "SBGM",       ARM_PricingModelType::SBGM,
    "HK",	      ARM_PricingModelType::HK,
    "Unknown",    ARM_PricingModelType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_PricingModelType( PricingModelTable, "PricingModel" );
const ARM_ArgConvReverse ARM_ArgConvReverse_PricingModelType( PricingModelTable, "PricingModel" );



ARGConvTable VnsPricingMethodTable[] =
{	
    /// Type Name       /// number
    "MONEYNESS",   ARM_MarketIRModel::MONEYNESS,
    "ATM",         ARM_MarketIRModel::ATM,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_VnsPricingMethod( VnsPricingMethodTable, "VnsPricingMethod" );
const ARM_ArgConvReverse ARM_ArgConvReverse_VnsPricingMethod( VnsPricingMethodTable, "VnsPricingMethod" );

ARGConvTable MMCalibProxyTable[] =
{	
    /// Type Name       /// number
    "ATM",			ARM_ModelParamsSmiled::ATM,
    "ATM0",			ARM_ModelParamsSmiled::AtmBlack,
    "MOMENT",		ARM_ModelParamsSmiled::MomentMatching,
	"LOCAL",		ARM_ModelParamsSmiled::LocalVolatility,
	"LOCAL0",		ARM_ModelParamsSmiled::LocalVolatilityBlack,
	"LOCAL+",		ARM_ModelParamsSmiled::LocalVolatilityWithRescaling,
	"LOCAL0+",		ARM_ModelParamsSmiled::LocalVolatilityBlackWithRescaling,
	"EFF",			ARM_ModelParamsSmiled::EffectiveSkew,
	"EFF+",			ARM_ModelParamsSmiled::EffectiveSkewWithRescaling,
	"GAUSS0",		ARM_ModelParamsSmiled::GaussBasketAtm,
	"GAUSS",		ARM_ModelParamsSmiled::GaussBasketMoneyness,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MMCalibProxy( MMCalibProxyTable, "MMCalibProxy" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MMCalibProxy( MMCalibProxyTable, "MMCalibProxy" );

ARGConvTable MMCorrelTypeTable[] =
{	
    /// Type Name       /// number
    "BETA",			ARM_ModelParamsSmiled::Beta,
    "THETA",		ARM_ModelParamsSmiled::Theta,
	"CORRELMATRIX",	ARM_ModelParamsSmiled::CorrelMatrix,
	"FWD",			ARM_ModelParamsSmiled::Fwd,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MMCorrelType( MMCorrelTypeTable, "MMCorrelType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MMCorrelType( MMCorrelTypeTable, "MMCorrelType" );


ARGConvTable MMCalibPatternTable[] =
{	
    /// Type Name       /// number
    "LIBOR",			ARM_SmiledFRMfactoryImp::LIBOR,
    "CMS_OLD",			ARM_SmiledFRMfactoryImp::CMS_OLD,
	"CMS",				ARM_SmiledFRMfactoryImp::CMS,
	"VMS_OLD",			ARM_SmiledFRMfactoryImp::VMS_OLD,
	"VMS",				ARM_SmiledFRMfactoryImp::VMS,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MMCalibPattern( MMCalibPatternTable, "MMCalibPattern" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MMCalibPattern( MMCalibPatternTable, "MMCalibPattern" );

ARGConvTable HWSVFormulaTable[] =
{	
    /// Type Name       /// number
    "Heston",			HWSVNumericals::Heston,
    "Lewis",			HWSVNumericals::Lewis,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_HWSVFormula( HWSVFormulaTable, "HWSVFormula" );
const ARM_ArgConvReverse ARM_ArgConvReverse_HWSVFormula( HWSVFormulaTable, "HWSVFormula" );


ARGConvTable MultiAssetsModelTable[] =
{	
    /// Type Name       /// number
    "2IRFX",			ARM_MultiAssetsModelType::twoirfx,
    "1IRFX",			ARM_MultiAssetsModelType::oneirfx,
	"IRFX",				ARM_MultiAssetsModelType::irfx,
	"UNKNOWN",          ARM_MultiAssetsModelType::unknown,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MultiAssetsType( MultiAssetsModelTable, "MultiAssetsModel" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MultiAssetsType( MultiAssetsModelTable, "MultiAssetsModel" );


ARGConvTable HestonMCSchemeTable[] =
{	
    /// Type Name       /// number
    "Euler",			ARM_HestonModel_Fx::Euler,
    "Andreasen",		ARM_HestonModel_Fx::Andreasen,
	"UNKNOWN",          ARM_HestonModel_Fx::Unknown,
    
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_HestonMCScheme( HestonMCSchemeTable, "HestonMCScheme" );
const ARM_ArgConvReverse ARM_ArgConvReverse_HestonMCScheme( HestonMCSchemeTable, "HestonMCScheme" );

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

