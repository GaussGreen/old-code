/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file FXVanillaFactory.h
 *
 *  \brief
 *
 *  \brief fxvanilla factory
 *
 *	\author  K.Belkheir
 *	\version 1.0
 *	\date February 2007
 */


#ifndef _INGPCALCULATORS_FXVANILLAFACTORY_H
#define _INGPCALCULATORS_FXVANILLAFACTORY_H

#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_GaussReplic2D;

struct ARM_FXVanillaFactoryImp
{
	ARM_GaussReplic2D CreateFXVanillaAndGaussReplic(const string& payName,
			const string& fx1Name,
			int callPut,
			double strike,
			ARM_VanillaType vanillaType,
			const string& fx2Name = "",
			double rho = 0.0,
			double alpha = 1.0,
			double beta = 0.0,
			double strike2 = 1.0);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_FXVanillaFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_FXVanillaFactoryImp>;
};

extern ARM_SingletonHolder<ARM_FXVanillaFactoryImp> ARM_FXVanillaFactory;

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
