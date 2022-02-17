/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillasumopt.cpp
 *
 *  \brief vanilla sum option
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2005
 */


#include "gpcalib/vanillasumopt.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"


/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_VanillaSumOptArg
///	Routine: ARM_VanillaSumOptArg
///	Returns: 
///	Action : copy constructor
////////////////////////////////////////////////////

ARM_VanillaSumOptArg::ARM_VanillaSumOptArg(const ARM_VanillaSumOptArg& rhs)
:	ARM_VanillaArg(rhs),
	itsCoeffs ( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsCoeffs->Clone())) ),
	itsFwdResetTimes( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdResetTimes->Clone())) ),  
	itsFwdStartTimes( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdStartTimes->Clone())) ),
	itsFwdEndTimes( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdEndTimes->Clone())) ),
	itsPayTime ( rhs.itsPayTime ),
	itsStrike ( rhs.itsStrike ),
	itsVolatilityRatio(rhs.itsVolatilityRatio)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSumOptArg
///	Routine: ARM_VanillaSumOptArg
///	Returns: 
///	Action : Affectation operator
////////////////////////////////////////////////////

ARM_VanillaSumOptArg& ARM_VanillaSumOptArg::operator=(const ARM_VanillaSumOptArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		
		itsCoeffs = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsCoeffs->Clone()));
		itsFwdResetTimes = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdResetTimes->Clone()));  
		itsFwdStartTimes = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdStartTimes->Clone()));
		itsFwdEndTimes = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsFwdEndTimes->Clone()));
		itsPayTime = rhs.itsPayTime;
		itsStrike = rhs.itsStrike;
		itsVolatilityRatio = rhs.itsVolatilityRatio;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSumOptArg
///	Routine: ARM_VanillaSumOptArg
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_VanillaSumOptArg::~ARM_VanillaSumOptArg()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSumOptArg
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Make a deep copy of the object
////////////////////////////////////////////////////

ARM_Object* ARM_VanillaSumOptArg::Clone() const
{
	return new ARM_VanillaSumOptArg(*this);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwapArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaSumOptArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSumOptArg
///	Routine: Price
///	Returns: 
///	Action : price a sum option with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaSumOptArg::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
    ARM_VectorPtr Price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
	ARM_GP_Vector strikeVec( 1,itsStrike);
    if(IRModel)
        Price = IRModel->VanillaSumOption(
		GetCurveName(),
		GetEvalTime(),
		GetCallPut(),
		*itsCoeffs,
		*itsFwdResetTimes,
		*itsFwdStartTimes,
		*itsFwdEndTimes,
		itsPayTime,
		*itsFwdPeriods,
		strikeVec,
		itsVolatilityRatio,
		const_cast<double*>(&itsSumFwd),
		const_cast<double*>(&itsSumVol),
		dumStates);
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price sum option, please advise");
  
    return (*Price)[0];
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSumOptArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaSumOptArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_SumOption"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

