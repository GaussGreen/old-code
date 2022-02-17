/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillasecuritydensity.cpp
 *
 *  \brief vanilla security descr + density for HK calib
 *
 *	\author  A. Chaix
 *	\version 1.0
 *	\date November 2005
 */

/// First To Include
#include "gpbase/removeidentifiedwarning.h"
/// gpcalib
#include "gpcalib/vanillaswaptionsmile.h"

/// gpbase
#include "gpbase/datestrip.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/zccrvfunctor.h"

/// kernel
#include <crv/zerocurv.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSmiledSwaption
///	Routine: constructor (default + contextual)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSmiledSwaption::ARM_VanillaSmiledSwaption(const ARM_VanillaSwaptionArg& refSwaption,
							const ARM_GP_Vector& target)
:	ARM_VanillaSwaptionArg (refSwaption),
	itsTargetParams	(target)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSmiledSwaption
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSmiledSwaption::ARM_VanillaSmiledSwaption( const ARM_VanillaSmiledSwaption& rhs)
:	ARM_VanillaSwaptionArg(rhs)
{
	itsTargetParams = rhs.itsTargetParams;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSmiledSwaption
///	Routine: destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSmiledSwaption::~ARM_VanillaSmiledSwaption()
{
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: Price
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaSmiledSwaption::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
    ARM_VectorPtr Price;
	double price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
    if(IRModel)
	{
		Price = IRModel->VanillaSmiledSwaption(
			GetCurveName(),
			GetEvalTime(),
			itsResetTime,
			*itsFixNominal,
			*itsFloatNominal,
			itsStartTime,
			itsEndTime,
			*itsFloatResetTimes,
			*itsFloatStartTimes,
			*itsFloatEndTimes,
			*itsFloatIntTerms,
			*itsFixPayTimes,
			*itsFixPayPeriods,
			*itsStrikes,
			GetCallPut(),
			dumStates,
			itsTargetParams,
			itsIsConstantNotional,
			itsIsConstantSpread,
			itsIsConstantStrike);
			price = (*Price)[0];
	}
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price swaption, please advise");
  
    return price;
}



////////////////////////////////////////////////////
///	Class  : ARM_VanillaSmiledSwaption
///	Routine: toString
///	Returns: string
///	Action : toString
////////////////////////////////////////////////////
string ARM_VanillaSmiledSwaption::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_VanillaSmiledSwaption\n";
    os << indent << "---------------------------\n";

	ARM_VanillaSwaptionArg::toString(indent,nextIndent);
	os << CC_NS(std,endl);
	
	return os.str();	

}


CC_END_NAMESPACE()


