/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file vanillapricer.h
 *	\brief pricer for vanilla instruments
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

/// gpcalib
#include "gpcalib/vanillapricer.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/vanillaarg.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"

/// gpbase
#include "gpbase/autocleaner.h"

/// gpmodels
#include "gpmodels/MultiAssets.h"


///kernel
#include <inst/swaption.h>
#include <inst/armdigital.h>
#include <inst/corridorleg.h>
#include <inst/option.h>
#include <inst/forex.h>
#include <inst/sumopt.h>
#include <inst/stripoption.h>
#include <inst/stripdigitaloption.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Struct   : ARM_VanillaPricer
///	Name     : Price (static function)
///	Action   : price instruments
/// Arguments: ARM_Security* sec, const ARM_PricingModel& mod 
////////////////////////////////////////////////////
string ARM_VanillaPricer::GetDefaultModelName( ARM_Security* Security, ARM_PricingModel* mod )
{
	ARM_MultiAssetsModel* multiAssetMod = dynamic_cast<ARM_MultiAssetsModel*>(mod);
	if( multiAssetMod && multiAssetMod->GetRefModel() )
		return multiAssetMod->GetRefModel()->GetModelName();

	/// test the type of the security!
	/// support only 
    if(		dynamic_cast<ARM_CapFloor*>(Security)
		||  dynamic_cast<ARM_Swaption*>(Security)
		||	dynamic_cast<ARM_Digital*>(Security)
		||	dynamic_cast<ARM_CorridorLeg*>(Security)
		||	dynamic_cast<ARM_SwapLeg*>(Security)
		||	dynamic_cast<ARM_SumOpt*>(Security))
    {
		return Security->GetCurrencyUnit()->GetCcyName();
	}
	else if( ARM_Option* option = dynamic_cast<ARM_Option*>(Security))
	{
		ARM_Security* underlying = option->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			string forCcy = forex->GetMainCurrency()->GetCcyName();
			string domCcy = forex->GetMoneyCurrency()->GetCcyName();
			return  forCcy + domCcy;
		}
		else
		{
			return "NO Equity Name";
		}
	}
	else if( ARM_StripOption* option = dynamic_cast<ARM_StripOption*>(Security) )
	{
		ARM_Security* underlying = option->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			string forCcy = forex->GetMainCurrency()->GetCcyName();
			string domCcy = forex->GetMoneyCurrency()->GetCcyName();
			return  forCcy + domCcy;
		}
		else
		{
			return "NO Equity Name";
		}
	}
	else if( ARM_Portfolio* portfolio = dynamic_cast<ARM_Portfolio*>(Security) )
	{
		return GetDefaultModelName(portfolio->GetAsset(0), mod);
	}
	else
    {
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"security not supported!");
	}
}


////////////////////////////////////////////////////
///	Struct   : ARM_VanillaPricer
///	Name     : Price (static function)
///	Action   : price instruments
/// Arguments: ARM_Security* sec, const ARM_PricingModel& mod 
////////////////////////////////////////////////////
double ARM_VanillaPricer::Price( ARM_Security* sec, ARM_PricingModel* pmod )
{
	/// create this corresponding vanilla Arg
	string modelName = ARM_VanillaPricer::GetDefaultModelName( sec, pmod );
	ARM_VanillaArg* vanillaArg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject( (ARM_Security*)sec, pmod->GetAsOfDate().GetJulian(),modelName);
	vanillaArg->SetCurveName( modelName );

	/// put an autocleaner on vanilla arg to ensure proper deletion after use!
	/// very simple smart pointor!
	CC_NS(ARM,ARM_AutoCleaner)<ARM_VanillaArg> Hold(vanillaArg);
	double price = vanillaArg->Price(pmod);

	ARM_ConverterFromKernel::TransferToSecurity(vanillaArg,sec);

	return price;
}



////////////////////////////////////////////////////
///	Struct   : ARM_VanillaPricer
///	Name     : Price (static function)
///	Action   : price instruments
/// Arguments: const ARM_VanillaArg& sec, const ARM_PricingModel& mod 
////////////////////////////////////////////////////
double ARM_VanillaPricer::Price( const ARM_VanillaArg& sec, ARM_PricingModel* pmod )
{
	return sec.Price(pmod);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

