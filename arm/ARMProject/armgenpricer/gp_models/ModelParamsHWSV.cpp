/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsHWSV.cpp
 *
 *  \brief
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date Septembre 2006
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/modelparamshwsv.h"

/// gpbase headers
#include "gpbase/gpmatrix.h"
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/comparisonfunctor.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/curvemodelparam.h"

#include <algorithm>
CC_USING_NS( std, sort )

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHWSV
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHWSV::ARM_ModelParamsHWSV( const ARM_ModelParamsHWSV& rhs )
: ARM_ModelParams(rhs)
{
	GenerateSchedule();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHWSV::ARM_ModelParamsHWSV( const ARM_ModelParamVector& params)
: ARM_ModelParams(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHWSV::~ARM_ModelParamsHWSV()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHWSV& ARM_ModelParamsHWSV::operator=(const ARM_ModelParamsHWSV& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHWSV
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHWSV::Clone() const
{
	return new ARM_ModelParamsHWSV(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHWSV
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsHWSV::ModelParamsTimeSteps() const
{
	return ARM_GP_VectorPtr( NULL );
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: FactorCount
///	Returns: number of diffusion factors
///	Action : 3 (1F) or 6 (2F)
////////////////////////////////////////////////////
size_t ARM_ModelParamsHWSV::FactorCount() const
{
	/// 1F : X, V
	/// 2F : X1, X2, V
	return	DoesModelParamExist( ARM_ModelParamType::VolatilityRatio ) &&
			DoesModelParamExist( ARM_ModelParamType::MeanReversionSpread ) ? 3 : 2;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHWSV
///	Routine: GenerateSchedule
///	Returns: nothing
///	Action : merge all model parameter schedules
////////////////////////////////////////////////////
void ARM_ModelParamsHWSV::GenerateSchedule()
{
	/// Compute a merged schedule of model parameters : vol, correl, VoV and if ncessary LTVar & VolRatio
    const std::vector<double>& schedAt	= GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->GetAbscisses();
    const std::vector<double>& schedBt	= GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve()->GetAbscisses();
    const std::vector<double>& schedCt	= GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->GetAbscisses();

    std::vector<double> tmpSched(schedAt.size()+schedBt.size());
    CC_NS(std,merge)(schedAt.begin(),schedAt.end(),schedBt.begin(),schedBt.end(),tmpSched.begin());

	std::vector<double> tmpSched2(tmpSched.size()+schedCt.size());
    CC_NS(std,merge)(tmpSched.begin(),tmpSched.end(),schedCt.begin(),schedCt.end(),tmpSched2.begin());

	std::vector<double> tmpSched3;
	if(DoesModelParamExist( ARM_ModelParamType::LongTermVol))
	{
		const std::vector<double>&  schedDt = GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->GetAbscisses();
		tmpSched3.resize(tmpSched2.size()+schedDt.size());
		CC_NS(std,merge)(tmpSched2.begin(),tmpSched2.end(),schedDt.begin(),schedDt.end(),tmpSched3.begin());
	}
	else
		tmpSched3 = tmpSched2;

	if(DoesModelParamExist( ARM_ModelParamType::VolatilityRatio ))
	{
		const std::vector<double>& schedDt  = GetModelParam( ARM_ModelParamType::VolatilityRatio).ToCurveModelParam().GetCurve()->GetAbscisses();
		itsSchedule = std::vector<double>(tmpSched3.size()+schedDt.size());
		CC_NS(std,merge)(tmpSched3.begin(),tmpSched3.end(),schedDt.begin(),schedDt.end(),itsSchedule.begin());
	}
	else
		itsSchedule = tmpSched3;

    std::vector<double>::iterator last=CC_NS(std,unique)(itsSchedule.begin(),itsSchedule.end());
	itsSchedule.resize( last-itsSchedule.begin() );
	if(itsSchedule[0] > K_NEW_DOUBLE_TOL)
        itsSchedule.insert(itsSchedule.begin(),0.0);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

