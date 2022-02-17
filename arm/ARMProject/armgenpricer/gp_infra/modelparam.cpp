/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelparam.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// ARM Kernel
#include <glob/expt.h>			/// necessary for exception throwing

/// gpbase
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/modelparam.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParam
///	Routine: ARM_ModelParam
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParam::ARM_ModelParam( ParamType type, bool adviseBreakPointTimes )
:	ARM_RootObject(), itsType( type ), itsAdviseBreakPointTimes( adviseBreakPointTimes )
{
    CC_ARM_SETNAME(ARM_MODELPARAM);
};


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: ARM_ModelParam
///	Returns: void
///	Action : copy constructor
///////////////////////////////////////
ARM_ModelParam::ARM_ModelParam( const ARM_ModelParam& rhs )
: ARM_RootObject( rhs ), itsType( rhs.itsType ), itsAdviseBreakPointTimes( rhs.itsAdviseBreakPointTimes )
{}


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: assignment operator
///	Returns: ARM_ModelParam&
///	Action : 
///////////////////////////////////////
ARM_ModelParam& ARM_ModelParam::operator =( const ARM_ModelParam& rhs )
{
	if( this !=	 &rhs )
	{
		ARM_RootObject::operator=( rhs );
		itsType  = rhs.itsType;
		itsAdviseBreakPointTimes = rhs.itsAdviseBreakPointTimes;
	}
	return *this;
}


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: assignment operator
///	Returns: destructor
///	Action : 
///////////////////////////////////////
ARM_ModelParam::~ARM_ModelParam()
{}


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: ToCurveModelParam
///	Returns: ARM_CurveModelParam
///	Action : downcast to an ARM_CurveModelParam
///////////////////////////////////////
ARM_CurveModelParam& ARM_ModelParam::ToCurveModelParam()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to ARM_CurveModelParam" );
}

const ARM_CurveModelParam& ARM_ModelParam::ToCurveModelParam() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to const ARM_CurveModelParam" );
}


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: ToSurfaceModelParam
///	Returns: ARM_SurfaceModelParam
///	Action : downcast to an ARM_SurfaceModelParam
///////////////////////////////////////

ARM_SurfaceModelParam& ARM_ModelParam::ToSurfaceModelParam()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to ARM_SurfaceModelParam" );
}

const ARM_SurfaceModelParam& ARM_ModelParam::ToSurfaceModelParam() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to const ARM_SurfaceModelParam" );
}


////////////////////////////////////////
///	Class  : ARM_ModelParam 
///	Routine: ToSurfaceListModelParam
///	Returns: ARM_SurfaceListModelParam
///	Action : downcast to an ARM_SurfaceListModelParam
///////////////////////////////////////

ARM_SurfaceListModelParam& ARM_ModelParam::ToSurfaceListModelParam()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to ARM_SurfaceListModelParam" );
}

const ARM_SurfaceListModelParam& ARM_ModelParam::ToSurfaceListModelParam() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to const ARM_SurfaceListModelParam" );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

