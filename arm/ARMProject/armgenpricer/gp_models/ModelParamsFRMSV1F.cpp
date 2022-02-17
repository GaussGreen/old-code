/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_ModelParamsMSV.cpp
 *
 *  \brief
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */


/// this header comes first as it includes some preprocessor constants!

#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"

#include "gpmodels/VolDiffusionParams.h"
#include "gpmodels/ModelParamsFRMSV1F.h"
#include "gpmodels/ModelParamsSFRM.h"
#include "gpbase/curve.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsFRMSV1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsFRMSV1F::ARM_ModelParamsFRMSV1F( const ARM_ModelParamsFRMSV1F& rhs )
: ARM_ModelParamsVec(rhs)
{
	itsSFRM_param_index = rhs.itsSFRM_param_index;
	itsSV_param_index = rhs.itsSV_param_index;	
//	GenerateSchedule();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsFRMSV1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsFRMSV1F::ARM_ModelParamsFRMSV1F(  const vector<ARM_ModelParams*>& paramsVec)
: ARM_ModelParamsVec(paramsVec),itsSFRM_param_index(0),itsSV_param_index(1)
{
	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsFRMSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsFRMSV1F::Clone() const
{
	return new ARM_ModelParamsFRMSV1F(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsFRMSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsFRMSV1F::~ARM_ModelParamsFRMSV1F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsFRMSV1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsFRMSV1F& ARM_ModelParamsFRMSV1F::operator=(const ARM_ModelParamsFRMSV1F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsVec::operator=(rhs);
		itsSFRM_param_index = rhs.itsSFRM_param_index;
		itsSV_param_index = rhs.itsSV_param_index;			/// Copy class attributes if any
	}
	return *this;
}

void ARM_ModelParamsFRMSV1F::ValidateModelParams()
{
	if (dynamic_cast<ARM_ModelParamsSFRM*> (itsParamsVector[0]))
	{
		if (!dynamic_cast<ARM_VolParams*> (itsParamsVector[1]))
				ARM_THROW( ERR_INVALID_ARGUMENT, "Stochastic Volatility Params are not valid" );
	}
	else
	{
		if (dynamic_cast<ARM_ModelParamsSFRM*> (itsParamsVector[1]))
		{
			itsSFRM_param_index = 1;
			if (dynamic_cast<ARM_VolParams*> (itsParamsVector[0]))
				itsSV_param_index = 0;
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "Stochastic Volatility Params are not valid" );
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "SFRM Params are not valid" );
	}
}



////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsFRMSV1F
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsFRMSV1F::ModelParamsTimeSteps() const
{
    //std::vector<double> sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)
	//return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
	return ARM_GP_VectorPtr( NULL );
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsFRMSV1F
///	Routines: 
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_ModelParamsSFRM* ARM_ModelParamsFRMSV1F::GetSFRMParams() const
{
	return ((ARM_ModelParamsSFRM*) itsParamsVector[itsSFRM_param_index]);
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsFRMSV1F
///	Routines: 
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_VolParams* ARM_ModelParamsFRMSV1F::GetVolParams() const
{
	return ((ARM_VolParams*) itsParamsVector[itsSV_param_index]);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

