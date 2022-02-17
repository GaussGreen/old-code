/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file pricerbuilder.cpp
 *  \brief file for pricer builder
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#include "xxxproject/pricerbuilder.h"
#include "xxxproject/vanillapricer.h"

#include <inst/capfloor.h>
#include <inst/armdigital.h>
#include <inst/swaption.h>
#include <inst/corridorleg.h>
#include <inst/option.h>
#include <inst/forex.h>
#include <inst/spreadoption.h>

#include <gpinflation/infcorridorLeg.h>

#include <gpcalculators/gencalculator.h>

CC_BEGIN_NAMESPACE( ARM )

ARM_Pricer* ARM_PricerBuilder::BuildPricer(ARM_Object* security)
{
	if (ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(security) )
		return new ARM_CAP_Pricer(capFloor);

	else if (ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(security) )
		return new ARM_OSW_Pricer(swaption);
	
//	else if (ARM_Option* fxOption = dynamic_cast<ARM_Option*>(security) )
//		return new ARM_FX_Pricer(fxOption);

	else if (ARM_InfCorridorLeg* infOption = dynamic_cast<ARM_InfCorridorLeg*>(security) )
		return new ARM_INF_Pricer(infOption);

//	else if (ARM_GenCalculator* genCalculator = dynamic_cast<ARM_GenCalculator* >(security) )
//		return new ARM_CalculatorPricer(genCalculator);

	return NULL;
}


CC_END_NAMESPACE()





/*


  ARM_Pricer* ARM_PricerBuilder::BuildPricer(ARM_Security* security)
{
	if (ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(security))
	{
		if (ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>(security))
		{
			return new ARM_SpreadOptionPricer(spreadOption);
		}
		else
		{
			return new ARM_CAP_Pricer(capFloor);
	//	}
	}
	if (ARM_Digital* digital = dynamic_cast<ARM_Digital*>(security))
	{
		return new ARM_DigitalPricer(digital);
	}
	else if (ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(security))
	{
		return new ARM_OSW_Pricer(swaption);
	}
	else if (ARM_CorridorLeg* corridorLeg = dynamic_cast<ARM_CorridorLeg*>(security))
	{
		return new ARM_CorridorPricer(corridorLeg);
	}	else if (ARM_Option* fxOption = dynamic_cast<ARM_Option*>(security))
	{
		ARM_Forex* forex = dynamic_cast<ARM_Forex*>(fxOption->GetUnderlying());

		if (!forex)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Just FX Options are managed by XXXProject." );
		}

		return new ARM_FX_Pricer(fxOption);
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "This security is not managed in the XXXproject." );
	}

	return NULL;
}

*/