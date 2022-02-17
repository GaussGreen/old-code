/*!
 *
 * Copyright (c) CDC IXIS CM May 2006 Paris
 *
 *	\file ModelParamsQGM2F.cpp
 *
 *  \brief 
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date May 2006
 */



/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/modelparamsqgm2f.h"
#include "gpmodels/ouprocess.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/bootstrap1d.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/gpmatrix.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM2FOneDim
///	Routine: Constructor, copy , egal, clone, PreProcessing
///	Returns: 
////////////////////////////////////////////////////
ARM_ModelParamsQGM2FOneDim::ARM_ModelParamsQGM2FOneDim( const ARM_ModelParamVector& params )
: ARM_ModelParams(params)
{}

ARM_ModelParamsQGM2FOneDim::ARM_ModelParamsQGM2FOneDim( const ARM_ModelParamsQGM2FOneDim& rhs )
: ARM_ModelParams(rhs)
{}

ARM_ModelParamsQGM2FOneDim::~ARM_ModelParamsQGM2FOneDim()
{}

ARM_ModelParamsQGM2FOneDim& ARM_ModelParamsQGM2FOneDim::operator=(const ARM_ModelParamsQGM2FOneDim& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}

ARM_Object* ARM_ModelParamsQGM2FOneDim::Clone() const
{
	return new ARM_ModelParamsQGM2FOneDim(*this);
}

void ARM_ModelParamsQGM2FOneDim::PreProcessing(ARM_ModelFitter& modelFitter,int factorNb)
{
    /// we put a typeid to implement in a sense a double dispatcher...
    /// we did not use a dynamic to avoid throwing exception of type std::bad_cast
    if( typeid(modelFitter) == typeid(ARM_Bootstrap1D) && 
        modelFitter.GetPortfolio()->GetAsset(0)->GetName() == ARM_SWAPTION &&
		(modelFitter.GetCalibParam())->GetType() == ARM_ModelParamType::Skew) 
    {
        modelFitter.SetCalibDirection(CalibDirection_Backward);
    }
}

///Fin ClassARM_ModelParamsQGM2FOneDim




////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM2F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsQGM2F::ARM_ModelParamsQGM2F( const vector<ARM_ModelParams*>& paramsVec  )
: ARM_ModelParamsVec(paramsVec)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM2F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsQGM2F::ARM_ModelParamsQGM2F( const ARM_ModelParamsQGM2F& rhs )
: ARM_ModelParamsVec(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM2F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsQGM2F::~ARM_ModelParamsQGM2F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM2F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsQGM2F& ARM_ModelParamsQGM2F::operator=(const ARM_ModelParamsQGM2F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsVec::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQGM2F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsQGM2F::Clone() const
{
	return new ARM_ModelParamsQGM2F(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

