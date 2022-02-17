/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsHeston_Fx.cpp
 *
 *  \brief Q model 1 factor FX version
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpvector.h"

/// gpmodels
#include "gpmodels/ModelParamsHeston_EqFx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_EqFx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHeston_EqFx::ARM_ModelParamsHeston_EqFx( const ARM_ModelParamVector& params, double spot  )
: ARM_ModelParams_EqFxBase(spot),
ARM_Heston_ModelParams(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_EqFx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_EqFx::ARM_ModelParamsHeston_EqFx( const ARM_ModelParamsHeston_EqFx& rhs )
: ARM_ModelParams_EqFxBase(rhs),
ARM_Heston_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_EqFx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_EqFx::~ARM_ModelParamsHeston_EqFx()
{}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

