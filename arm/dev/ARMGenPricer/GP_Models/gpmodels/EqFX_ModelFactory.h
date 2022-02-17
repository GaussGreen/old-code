/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file EqFx_ModelFactory.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_EQFX_MODELFACTORY_H
#define _INGPMODELS_EQFX_MODELFACTORY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/curvetypedef.h"
#include "gpbase/typedef.h"
#include "gpbase/curvematrix.h"

/// gpinfra
#include "gpinfra/modelparams.h"

/// forward declaration in the global namespace
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_PricingModel;


class ARM_EqFx_ModelFactoryImp
{
public:
	enum ModelType
	{ 
		BS_Model = 0,
		Q1F_Model,
		HESTON_Model,
		SABR_Model,
		CEV_Model,
		Mixture_Model

	};

    ARM_PricingModel* CreateModel(
		const ARM_ZeroCurvePtr& zc, 
		const ARM_ModelParamVector& params	= ARM_ModelParamVector(),
		double spot							= ARM_NULL_OBJECT,
		const ARM_ZeroCurvePtr& forCurve	= ARM_ZeroCurvePtr(NULL),
		const ARM_CurveMatrix& correlMatrix   = ARM_CurveMatrix(),
		ModelType modelType					= Q1F_Model,
		int mcScheme						= 0);

private:

	/// to forbid client from using it except for the singleton holder
	ARM_EqFx_ModelFactoryImp(){};
	friend class ARM_SingletonHolder<ARM_EqFx_ModelFactoryImp>;
};

extern ARM_SingletonHolder<ARM_EqFx_ModelFactoryImp> ARM_EqFx_ModelFactory;


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
