/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsHW.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/ModelParamsHW.h"

CC_BEGIN_NAMESPACE( ARM )

const double ARM_ModelParamsHW::MrsMinValue = 1.0e-5;
const double ARM_ModelParamsHW::VOL_LIMIT   = 0.000001;


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW::ARM_ModelParamsHW( const ARM_ModelParamsHW& rhs )
: ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW::ARM_ModelParamsHW( const ARM_ModelParamVector& params )
: ARM_ModelParams(params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW::~ARM_ModelParamsHW()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW& ARM_ModelParamsHW::operator=(const ARM_ModelParamsHW& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

