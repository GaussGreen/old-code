/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file pricerfactory.h
 *	\brief pricer factory
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPCALCULATORS_PRICERFACTORY_H
#define _INGPCALCULATORS_PRICERFACTORY_H

/// gpinfra
#include "gpbase/port.h"
#include <utility> /// for pair

/// forward declaration
class ARM_Object;

CC_BEGIN_NAMESPACE( ARM )


struct ARM_PricerFactory
{
	static CC_NS(std,pair)<bool,double> Price( ARM_Object* sec, ARM_Object* mod	);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
