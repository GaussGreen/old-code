/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsSABR_Fx.cpp
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
#include "gpmodels/ModelParamsSABR_EqFx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSABR_EqFx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsSABR_EqFx::ARM_ModelParamsSABR_EqFx( const ARM_ModelParamVector& params, double spot  )
: ARM_ModelParams_EqFxBase(spot),
ARM_SABR_ModelParams(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSABR_EqFx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSABR_EqFx::ARM_ModelParamsSABR_EqFx( const ARM_ModelParamsSABR_EqFx& rhs )
: ARM_ModelParams_EqFxBase(rhs),
ARM_SABR_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSABR_EqFx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsSABR_EqFx::~ARM_ModelParamsSABR_EqFx()
{}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

