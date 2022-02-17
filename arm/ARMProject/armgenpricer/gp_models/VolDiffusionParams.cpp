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
#include "gpmodels/VolDiffusionParams.h"

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
// ARM_VolParams
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_VolParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_VolParams::ARM_VolParams( const ARM_VolParams& rhs )
: ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_VolParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_VolParams::ARM_VolParams( const ARM_ModelParamVector& params)
: ARM_ModelParams(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_VolParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_VolParams::~ARM_VolParams()
{}

////////////////////////////////////////////////////
///	Class  : ARM_VolParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_VolParams& ARM_VolParams::operator=(const ARM_VolParams& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_VolParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_VolParams::Clone() const
{
	return new ARM_VolParams(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_VolParams
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_VolParams::ModelParamsTimeSteps() const
{
    //std::vector<double> sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)
	//return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
	return ARM_GP_VectorPtr( NULL );
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

