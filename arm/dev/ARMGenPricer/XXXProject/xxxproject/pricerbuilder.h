/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file pricerbuilder.h
 *  \brief file for pricer builder
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#ifndef INXXXPROJECT_PRICERBUILDER
#define INXXXPROJECT_PRICERBUILDER

#include "xxxproject/pricer.h"
#include <gpbase/port.h>

class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricerBuilder
{
public:
	static ARM_Pricer* BuildPricer(ARM_Object *);
};

CC_END_NAMESPACE()

#endif