/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file PricingModelType.h
 *
 *  \brief
 *
 *  \brief type for the pricing model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPINFRA_PRICINGMODELTYPE_H
#define _INGPINFRA_PRICINGMODELTYPE_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

/// type of the model on the bit of the integer
///
///   bit	modeltype
///		1	interest rate
///		2	equity
///		3	fx
///		4	commodity
///		5	credit
///		6	emerging markets
///		7	non stochastic

const int MT_INTEREST_RATE_MODEL	= 0x0001;
const int MT_EQUITY_MODEL			= 0x0002;
const int MT_FX_MODEL				= 0x0004;
const int MT_COMMODITY_MODEL		= 0x0008;
const int MT_CREDIT_MODEL			= 0x0010;
const int MT_EMERGING_MKT_MODEL		= 0x0020;
const int MT_NON_STOCHASTIC_MODEL	= 0x0040;
const int MT_MULTIASSET_MODEL		= 0x0080;

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
