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
#include "gpmodels/ModelParamsMSV.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpbase/curve.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsMSV::ARM_ModelParamsMSV( const ARM_ModelParamsMSV& rhs )
: ARM_ModelParams(rhs)
{
	GenerateSchedule();
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsMSV::ARM_ModelParamsMSV( const ARM_ModelParamVector& params )
: ARM_ModelParams(params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsMSV::~ARM_ModelParamsMSV()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsMSV& ARM_ModelParamsMSV::operator=(const ARM_ModelParamsMSV& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsMSV
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ModelParamsMSV::GenerateSchedule()
{
	/// Compute a merged schedule of model parameters
    const std::vector<double>& schedAt  = GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->GetAbscisses();
    const std::vector<double>& schedBt  = GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->GetAbscisses();
    const std::vector<double>& schedCt = GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->GetAbscisses();
    std::vector<double> tmpSched(schedAt.size()+schedBt.size());
	itsSchedule = std::vector<double>(tmpSched.size()+schedCt.size());
    CC_NS(std,merge)(schedAt.begin(),schedAt.end(),schedBt.begin(),schedBt.end(),tmpSched.begin());
    CC_NS(std,merge)(tmpSched.begin(),tmpSched.end(),schedCt.begin(),schedCt.end(),itsSchedule.begin());
    std::vector<double>::iterator last=CC_NS(std,unique)(itsSchedule.begin(),itsSchedule.end());
	itsSchedule.resize( last-itsSchedule.begin() );
    if(itsSchedule[0] < K_NEW_DOUBLE_TOL && itsSchedule.size() > 1)
        itsSchedule.erase(itsSchedule.begin());
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

