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


#ifndef _INGPCALIB_VANILLAPRICER_H
#define _INGPCALIB_VANILLAPRICER_H

/// gpinfra
#include "gpbase/port.h"
/// STL
#include <string>
CC_USING_NS(std,string)

/// forward declaration
class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;
struct ARM_VanillaArg;


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaPricer
///////////////////////////////////////////////////////////////
struct ARM_VanillaPricer
{
	static double Price( ARM_Security* sec, ARM_PricingModel* mod );
	static double Price( const ARM_VanillaArg& sec, ARM_PricingModel* mod );
	static string GetDefaultModelName( ARM_Security* Security, ARM_PricingModel* mod );

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
