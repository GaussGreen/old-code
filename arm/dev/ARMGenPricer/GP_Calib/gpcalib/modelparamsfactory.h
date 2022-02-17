/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file modelparamsfactory.h
 *  \brief factory class for model params
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPCALIB_MODELPARAMSFACTORY_H
#define _INGPCALIB_MODELPARAMSFACTORY_H

#include "gpbase/port.h"
#include "gpbase/surfacetypedef.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;

class ARM_ModelParam;
class ARM_SurfaceModelParam;


struct ARM_ModelParamFactoryImp
{
	ARM_ModelParam* CreateModelParam(
		int type,
		ARM_GP_Vector* values    = NULL,
		ARM_GP_Vector* breakPointTimes	= NULL,
		const string& name              = "",
		const string& interpolatorName  = "STEPUPRIGHT",
		ARM_GP_Vector* lowerbound		= NULL ,
		ARM_GP_Vector* upperbound		= NULL,
		bool adviseBreakPointTimes		= false,
		const string & currency			= "");

	ARM_SurfaceModelParam* CreateSurfaceModelParam(
		int modelParamType,
		ARM_Surface* surface,
		const string& modelParamName,
		double lowerBound,
		double upperBound,
		bool adviseBreakPointTimes		= false);

	ARM_ModelParam* CreateTrigoCorrelMatParam( 
		double theta,
		const ARM_DateStrip& datestrip,
		double asOf,
		const string& interpolatorName,
		ARM_GP_Vector* lowerBound   = NULL,
		ARM_GP_Vector* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );

	ARM_ModelParam* CreateTrigoCorrelMatParam( 
		double theta,
		const ARM_GP_Vector& resetDates,
		double asOf,
		const string& interpolatorName,
		ARM_GP_Vector* lowerBound   = NULL,
		ARM_GP_Vector* upperBound   = NULL,
		bool adviseBreakPointTimes	= false );


private:
	/// to forbid client from using it except for the singleton holder
	ARM_ModelParamFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_ModelParamFactoryImp>;
};

extern ARM_SingletonHolder<ARM_ModelParamFactoryImp> ARM_ModelParamFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

