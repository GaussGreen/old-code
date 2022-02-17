/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsSFRMFactory.h
 *
 *  \brief class to control the creation of model params of SFRM!
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPMODELS_MODELPARAMSSFRMFACTORY_H
#define _INGPMODELS_MODELPARAMSSFRMFACTORY_H

#include "gpbase/port.h"
#include "ModelParamsSFRM.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;

///-----------------------------------------------------------------------------
/// \class ARM_ModelParamsSFRMFactory
/// \brief Factory class to create model params for SFRM implemented 
/// as a singleton
///-----------------------------------------------------------------------------
struct ARM_ModelParamsSFRMFactoryImp
{
	ARM_ModelParamsSFRM* CreateModelParamsSFRM(
		const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb, size_t volType ) const;

private:
	/// to forbid client from using it except for the singleton holder
	ARM_ModelParamsSFRMFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_ModelParamsSFRMFactoryImp>;
};

extern ARM_SingletonHolder<ARM_ModelParamsSFRMFactoryImp> ARM_ModelParamsSFRMFactory;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

