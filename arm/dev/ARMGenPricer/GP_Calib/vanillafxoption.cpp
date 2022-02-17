/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillafxoption.cpp
 *
 *  \brief vanilla fx  option
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#include "gpcalib/vanillafxoption.h"

#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpmodels/ModelParams_EqFxBase.h"



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaFxOption
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaFxOption::ImpliedVol(ARM_PricingModel* model) const
{
  
	/// create a dumState to avoid dummy price..
	double impliedVol;
	ARM_PricingFunctionEquity* fxModel = dynamic_cast<ARM_PricingFunctionEquity*>(model);
	
    if(fxModel)
    {
		impliedVol = fxModel->ImpliedVol(*this);
    }
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an forex model: not derived from ARM_PricingFunctionEquity. So cannot price fx option, please advise");

    return impliedVol;

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaFxOption
///	Routine: Price
///	Returns: 
///	Action : price an fx option with a model checking that it
///				is derived from an fx model
////////////////////////////////////////////////////

double ARM_VanillaFxOption::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionEquity* FxModel = dynamic_cast<ARM_PricingFunctionEquity*>(model);
    ARM_PricingModel* IRModel = dynamic_cast<ARM_PricingModel*>(model);
	const ARM_ModelParams_Fx* fxModelParams = dynamic_cast<const ARM_ModelParams_Fx* const>(model->GetModelParams());
    double Price;

    if (FxModel)
	{
       if ( GetExpiry() <= GetEvalTime() ) // option expired            
       {
          ARM_ZeroCurvePtr ZcCurve = IRModel->GetZeroCurve();
		  double zcT = ZcCurve->DiscountPrice(GetPayTime()/K_YEAR_LEN);   
          Price = MAX(GetCallPut()*(GetFixedFx()- GetStrike()), 0.0)*zcT;
       }
       else
       {
          ARM_VectorPtr PriceVec;
          ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    

          PriceVec = FxModel->CallScalar(GetCurveName(),
		                              GetEvalTime(),
		                              GetExpiry(),
		                              GetFwdTime(),
		                              GetStrike(),
		                              GetCallPut(),
                                      GetPayTime(),
		                              dumStates);
          Price = (*PriceVec)[0];
       }

       if (( itsAccountingPricingFlag == ARM_DISC_ACCOUNTING_METH ) &&  ( GetPayTime() == GetEvalTime() ) )
       {
           Price = 0.0;
       }

	   if (itsPayFgnCcy) /// Wrong wrong should be fixed
		  Price /= fxModelParams->GetSpot();
	}
    else
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Model is not an FX Model: not derived from ARM_PricingFunctionEquity. So cannot price Fx option, please advise");

    return Price;

}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaFxOption
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaFxOption::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaFxOption"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

