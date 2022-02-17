/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillairfxswaption.cpp
 *
 *  \brief vanilla IR/FX  swaption
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date March 2006
 */

#include "gpcalib/vanillairfxswaption.h"
#include "gpcalib/vanillafxoption.h"
#include "gpcalib/vanillaswaption.h"

#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

ARM_VanillaIrFxSwaption::~ARM_VanillaIrFxSwaption()
{
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIrFxSwaption
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action : No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaIrFxSwaption::ImpliedVol(ARM_PricingModel* model) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : No formula is valid to calculate Implied Volatilty");
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIrFxSwaption
///	Routine: Price
///	Returns: 
///	Action : price an fx option with a model checking that it
///				is derived from an fx model
////////////////////////////////////////////////////

double ARM_VanillaIrFxSwaption::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionEquity* irFxModel = dynamic_cast<ARM_PricingFunctionEquity*>(model);
    ARM_VectorPtr Price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    

    if(irFxModel)
	{
        Price = irFxModel->HybridCallScalar(
			GetCurveName(),
			GetEvalTime(),
			GetExpiry(),
			GetCallPut(),
			itsFxResetTimes,
			itsFxSettlementTimes,
			itsFxPayTimes,
			itsFxNominals,
			itsIrSwaption->GetResetTime(),
			*(itsIrSwaption->GetFixNotional()),
			*(itsIrSwaption->GetFloatNotional()),
			itsIrSwaption->GetStartTime(),
			itsIrSwaption->GetEndTime(),
			*(itsIrSwaption->GetFloatResetTimes()),
			*(itsIrSwaption->GetFloatStartTimes()),
			*(itsIrSwaption->GetFloatEndTimes()),
			*(itsIrSwaption->GetFloatIntTerms()),
			*(itsIrSwaption->GetFixPayTimes()),
			*(itsIrSwaption->GetFixPayPeriods()),
			itsStrike,
			dumStates);
	}
    else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : input model is not a Market IR/FX model !");
  
    return (*Price)[0];
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaFxOption
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaIrFxSwaption::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaIrFxSwaption"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

