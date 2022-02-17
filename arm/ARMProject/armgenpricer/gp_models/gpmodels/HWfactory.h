/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file modelparamsfactory.h
 *  \brief factory class for model params
 *	\author  E. BEzzine
 *	\version 1.0
 *	\date June 2005
 */


#ifndef _INGPCALIB_MODELPARAMSFACTORY_H
#define _INGPCALIB_MODELPARAMSFACTORY_H

#include "gpbase/port.h"
#include "typedef.h"
#include "gpinfra/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;

class ARM_HullWhite;

struct ARM_HWModelFactoryImp
{
	ARM_HullWhite* CreateHWModel(const ARM_ZeroCurvePtr& zc,
        const ARM_ModelParamVector& params = ARM_ModelParamVector() );
private:
	/// to forbid client from using it except for the singleton holder
	ARM_HWModelFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_HWModelFactoryImp>;
};

extern ARM_SingletonHolder<ARM_HWModelFactoryImp> ARM_HWModelFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

