/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsMSV1F.cpp
 *
 *  \brief
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/modelparamshwsv1f.h"

/// gpbase headers
#include "gpbase/gpmatrix.h"
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/comparisonfunctor.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/curvemodelparam.h"

/// gpmodel 
#include "gpmodels/TargetFuncHW.h"

/// gpcalib
#include "gpcalib/modelfitter.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"

/// kernel
#include <inst/portfolio.h>
#include <inst/irindex.h>


#include <algorithm>
CC_USING_NS( std, sort )

#include <functional>
CC_USING_NS( std, ptr_fun )

CC_BEGIN_NAMESPACE( ARM )

const double VOL_LIMIT      = 0.000001;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHWSV1F
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHWSV1F::ARM_ModelParamsHWSV1F( const ARM_ModelParamsHWSV1F& rhs )
: ARM_ModelParams(rhs)
{
	GenerateSchedule();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHWSV1F::ARM_ModelParamsHWSV1F( const ARM_ModelParamVector& params)
: ARM_ModelParams(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHWSV1F::~ARM_ModelParamsHWSV1F()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHWSV1F& ARM_ModelParamsHWSV1F::operator=(const ARM_ModelParamsHWSV1F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHWSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHWSV1F::Clone() const
{
	return new ARM_ModelParamsHWSV1F(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHWSV1F
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsHWSV1F::ModelParamsTimeSteps() const
{
	return ARM_GP_VectorPtr( NULL );
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ModelParamsHWSV1F::GenerateSchedule()
{
	/// Compute a merged schedule of model parameters
    const ARM_GP_Vector& schedAt  = GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedBt  = GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedCt = GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->GetAbscisses();
    ARM_GP_Vector tmpSched(schedAt.size()+schedBt.size());
	itsSchedule = ARM_GP_Vector(tmpSched.size()+schedCt.size());
    CC_NS(std,merge)(schedAt.begin(),schedAt.end(),schedBt.begin(),schedBt.end(),tmpSched.begin());
    CC_NS(std,merge)(tmpSched.begin(),tmpSched.end(),schedCt.begin(),schedCt.end(),itsSchedule.begin());
    ARM_GP_Vector::iterator last=CC_NS(std,unique)(itsSchedule.begin(),itsSchedule.end());
	itsSchedule.resize( last-itsSchedule.begin() );
	if(itsSchedule[0] > K_NEW_DOUBLE_TOL)
        itsSchedule.insert(itsSchedule.begin(),0.0);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

