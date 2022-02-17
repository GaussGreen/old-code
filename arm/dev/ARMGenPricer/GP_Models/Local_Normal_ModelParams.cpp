/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file Local_Normal_ModelParams.cpp
 *
 *  \brief base class for NormalModel Model Params
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/local_normal_modelparams.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"

/// gpmodels
#include "gpinfra/surfacelistmodelparam.h"

/// gpbase
#include "gpbase/surface.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector & param types
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams::ARM_Local_Normal_ModelParams( const ARM_ModelParamVector& params,
														    const ARM_IntVector& paramTypes)
:	ARM_ModelParams(params)
{
	if(paramTypes.size()>0)
		ValidateModelParams(paramTypes);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_ModelParams
///	Routines: Validate
///	Returns : void
///	Action  : checks the validy ²
////////////////////////////////////////////////////
void ARM_Local_Normal_ModelParams::ValidateModelParams(const ARM_IntVector& paramTypes) const
{
	size_t nbParams = paramTypes.size();
	if(nbParams<=2)
	{
		/// Default validation ARM_Local_Normal_ModelParams is required to have
		/// o 1 surface list of type ForwardAdjustment
		/// o 1 surface list of type Volatility
    
		/// Adjustments
		if(!DoesModelParamExist(ARM_ModelParamType::ForwardAdjustment) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");

		if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::ForwardAdjustment) ) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");

		/// Volatility
		if(!DoesModelParamExist(ARM_ModelParamType::Volatility) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");

		if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::Volatility) ) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");
	}
	else
	{
		/// General validation
		for(size_t i=0;i<nbParams;++i)
		{
			if(!DoesModelParamExist(paramTypes[i]) )
				ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");

			if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(paramTypes[i]) ) )
				ARM_THROW( ERR_INVALID_ARGUMENT, "Invalid Local Normal Model Params");
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams::ARM_Local_Normal_ModelParams( const ARM_Local_Normal_ModelParams& rhs )
: ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams::~ARM_Local_Normal_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams& ARM_Local_Normal_ModelParams::operator=(const ARM_Local_Normal_ModelParams& rhs)
{
		if (&rhs != this)
	{ 
		this->~ARM_Local_Normal_ModelParams();
		new (this) ARM_Local_Normal_ModelParams (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Local_Normal_ModelParams::Clone() const
{
	return new ARM_Local_Normal_ModelParams(*this);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

