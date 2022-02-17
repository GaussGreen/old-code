/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillaequityoption.cpp
 *
 *  \brief vanilla equity option
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#include "gpcalib/vanillaequityoption.h"

#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaEqOption
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaEqOption::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaEqOption
///	Routine: Price
///	Returns: 
///	Action : price an equity option with a model checking that it
///				is derived from an equity model
////////////////////////////////////////////////////

double ARM_VanillaEqOption::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionEquity* EqModel = dynamic_cast<ARM_PricingFunctionEquity*>(model);
    ARM_VectorPtr Price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    

    if(EqModel)
        Price = EqModel->CallScalar(
		GetCurveName(),
		GetEvalTime(),
		GetExpiry(),
		GetFwdTime(),
		GetStrike(),
		GetCallPut(),
        GetPayTime(),
		dumStates);
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Model is not an equity model: not derived from ARM_PricingFunctionEquity. So cannot price Equity option");
  
    return (*Price)[0];
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaEqOption
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaEqOption::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaEqOption"; 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

