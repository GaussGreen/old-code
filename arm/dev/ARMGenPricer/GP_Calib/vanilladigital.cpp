/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanilladigital.cpp
 *
 *  \brief vanilla digital
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/vanilladigital.h"

/// gpbase
#include "gpbase/gpvector.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArg
///	Routine: PriceOplet
///	Returns: double
///	Action : price a caplet
////////////////////////////////////////////////////

double ARM_VanillaDigitalArg::PriceOplet(
	ARM_PricingFunctionIR* IRModel,
	const string& curveName, 
	double evalTime,
	double payTime,
	double period,
	double payNotional,
	double fwdResetTime,	/// used for volatility computation
	double fwdStartTime,
	double fwdEndTime,
	double fwdPeriod,
	double strike,
	int capFloor,
	const ARM_PricingStatesPtr& states ) const
{
	static ARM_GP_VectorPtr result;

	result = IRModel->VanillaDigitalScalar(
		curveName,
		evalTime,
		payTime,
		period,
		payNotional,
		fwdResetTime,
		fwdStartTime,
		fwdEndTime,
		fwdPeriod,
		strike,
		capFloor,
		states);
	return (*result)[0];
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaDigitalArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaDigitalArg"; 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

