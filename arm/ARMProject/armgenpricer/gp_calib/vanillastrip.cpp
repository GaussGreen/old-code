/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillacap.cpp
 *
 *  \brief vanilla cap is a vanilla interest rate cap
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date November 2003
 */

#include "gpcalib/vanillastrip.h"

/// gpbase
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaStripArg
///	Routine: operator=
///	Returns: 
///	Action :  Copy constructor
////////////////////////////////////////////////////

ARM_VanillaStripArg::ARM_VanillaStripArg(const ARM_VanillaStripArg& rhs)
:	
	ARM_VanillaArg(rhs),
	itsVanillas(rhs.itsVanillas),
	itsCoeffs(rhs.itsCoeffs),
    itsAccountingPricingFlag(rhs.itsAccountingPricingFlag)
{
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaStripArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaStripArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaStripArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaStripArg::Price(ARM_PricingModel* model) const
{
    double price = 0.0;

	size_t i;

	for (i = 0; i < itsVanillas.size(); ++i)
	{
		price += itsCoeffs[i]*itsVanillas[i]->Price(model);
	}

	return price;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaStripArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaStripArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaStripArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

