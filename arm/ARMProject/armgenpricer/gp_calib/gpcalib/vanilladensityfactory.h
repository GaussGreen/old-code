/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file vanilladensityfactory.h
 *  \brief factory class for vanilla density
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date March 2007
 */


#ifndef _INGPCALIB_VANILLADENSITYFACTORY_H
#define _INGPCALIB_VANILLADENSITYFACTORY_H

#include "gpbase/port.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;

class ARM_VanillaSecurityDensity;
class ARM_PricingModel;

struct ARM_VanillaDensityFactorImp
{
	ARM_VanillaSecurityDensity* CreateVanillaDenstityFactory(const ARM_PricingModel&  model,
		const ARM_Date&  expiryDate,
		bool isDirect	= false);
private:
	/// to forbid client from using it except for the singleton holder
	ARM_VanillaDensityFactorImp() {};
	friend class ARM_SingletonHolder<ARM_VanillaDensityFactorImp>;
};

extern ARM_SingletonHolder<ARM_VanillaDensityFactorImp> ARM_VanillaDensityFactor;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

