/*!
 *
 * Copyright (c) IXIS CIB July 2005 Paris
 *
 *	\file Local_SLN_ModelParams.cpp
 *
 *  \brief model params for local normal model
 *	\author  J-M Prié
 *	\version 1.0
 *	\date July 2005
 */

#include "gpmodels/Local_SLN_ModelParams.h"

/// gpbase (for ARM_CountedPtr deletion)
#include "gpbase/surface.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_Local_SLN_ModelParams::ARM_Local_SLN_ModelParams( const ARM_ModelParamVector& params )
:	ARM_ModelParams(params)
{
	ValidateModelParams();
}




////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_ModelParams
///	Routines: Validate
///	Returns : void
///	Action  : checks the validy
////////////////////////////////////////////////////
void ARM_Local_SLN_ModelParams::ValidateModelParams() const
{	
	// ARM_Local_SLN_ModelParams is required to have
	// o 1 surface list of type ForwardAdjustment (only first surface used at the moment)
	// o 1 surface list of type Volatility (only first surface used at the moment)
	// o 1 surface list of type Shift (only first surface used at the moment)
    
	/// Adjustments
	if(!DoesModelParamExist(ARM_ModelParamType::ForwardAdjustment) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Forward Adjustment surface is missing in Local_SLN_ModelParams");

	if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::ForwardAdjustment) ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Forward Adjustment surface is invalid in Local_SLN_ModelParams");

	/// Volatility or QVol
	if( DoesModelParamExist(ARM_ModelParamType::Volatility) )
    {
	    if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::Volatility) ) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, "Volatility surface is invalid in Local_SLN_ModelParams");
    }
    else if( DoesModelParamExist(ARM_ModelParamType::QVol) )
    {
	    if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::QVol) ) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, "QVol surface is invalid in Local_SLN_ModelParams");
    }
    else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Volatility or QVol surface is missing in Local_SLN_ModelParams");


	/// Shift or QParameter
	if( DoesModelParamExist(ARM_ModelParamType::Shift) )
    {
	    if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::Shift) ) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, "Shift surface is invalid in Local_SLN_ModelParams");
    }
    else if( DoesModelParamExist(ARM_ModelParamType::QParameter) )
    {
	    if(!dynamic_cast<const ARM_SurfaceListModelParam*>( &GetModelParam(ARM_ModelParamType::QParameter) ) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, "QParameter surface is invalid in Local_SLN_ModelParams");
    }
    else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Shift or QParameter surface is missing in Local_SLN_ModelParams");

}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Local_SLN_ModelParams::ARM_Local_SLN_ModelParams( const ARM_Local_SLN_ModelParams& rhs )
: ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Local_SLN_ModelParams::~ARM_Local_SLN_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Local_SLN_ModelParams& ARM_Local_SLN_ModelParams::operator=(const ARM_Local_SLN_ModelParams& rhs)
{
		if (&rhs != this)
	{ 
		this->~ARM_Local_SLN_ModelParams();
		new (this) ARM_Local_SLN_ModelParams (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Local_SLN_ModelParams::Clone() const
{
	return new ARM_Local_SLN_ModelParams(*this);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

