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
#include "gpmodels/modelparamsmsv1f.h"

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
// ARM_ModelParamsMSV1F
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsMSV1F::ARM_ModelParamsMSV1F( const ARM_ModelParamsMSV1F& rhs )
: ARM_ModelParamsMSV(rhs)
{
	itsIRIndex = rhs.itsIRIndex?(ARM_IRIndex*)rhs.itsIRIndex->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsMSV1F::ARM_ModelParamsMSV1F( const ARM_ModelParamVector& params, 	ARM_IRIndex* index )
: ARM_ModelParamsMSV(params)
{
	itsIRIndex = index? (ARM_IRIndex*) index->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsMSV1F::~ARM_ModelParamsMSV1F()
{
	delete itsIRIndex;
    itsIRIndex = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsMSV1F& ARM_ModelParamsMSV1F::operator=(const ARM_ModelParamsMSV1F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsMSV::operator=(rhs);
		/// Copy class attributes if any
		itsIRIndex		= rhs.itsIRIndex?(ARM_IRIndex*)rhs.itsIRIndex->Clone() : NULL;	
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsMSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsMSV1F::Clone() const
{
	return new ARM_ModelParamsMSV1F(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsMSV1F
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsMSV1F::ModelParamsTimeSteps() const
{
    //std::vector<double> sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)
	//return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
	return ARM_GP_VectorPtr( NULL );
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: BetatT
///	Returns: value fo beta(t,T)
///	Action : beta(t,T)=(1-exp(-MRS*(T-t))/MRS
///                   = Integ{t->T,exp(-MRS*(u-t))du}
////////////////////////////////////////////////////
double ARM_ModelParamsMSV1F::BetatT(double t,double T) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

    if(fabs(MRSValue)>K_NEW_DOUBLE_TOL)
        return (1.0-exp(-MRSValue*(T-t)/K_YEAR_LEN))/MRSValue;
    else
        return (T-t)/K_YEAR_LEN;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV1F
///	Routine: Deriv_BetatT
///	Returns: value fo beta(t,T)
///	Action : beta(t,T)=exp(-MRS*(T-t))
////////////////////////////////////////////////////
double ARM_ModelParamsMSV1F::Deriv_BetatT(double t,double T) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
        return (exp(-MRSValue*(T-t)/K_YEAR_LEN));
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

