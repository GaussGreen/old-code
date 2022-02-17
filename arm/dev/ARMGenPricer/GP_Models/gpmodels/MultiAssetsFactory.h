/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file MultiAssetsFactory.h
 *
 *  \brief
 *
 *  \brief multi asset model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_MULTIASSETSFACTORY_H
#define _INGPMODELS_MULTIASSETSFACTORY_H

#include "gpbase/port.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_MultiAssetsModel;
class ARM_ModelNameMap;
class ARM_CurveMatrix;

struct ARM_MultiAssetsFactoryImp
{
	ARM_MultiAssetsModel* CreateMultiAssetsModel(const ARM_ModelNameMap& modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix, 
		const string& multiAssetsModelName = "UNKNOWN");

private:
	/// to forbid client from using it except for the singleton holder
	ARM_MultiAssetsFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_MultiAssetsFactoryImp>;
};

extern ARM_SingletonHolder<ARM_MultiAssetsFactoryImp> ARM_MultiAssetsFactory;

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
