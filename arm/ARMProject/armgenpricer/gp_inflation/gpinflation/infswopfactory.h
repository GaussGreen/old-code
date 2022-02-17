/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file discretisationscheme.cpp
 *  \brief computes the inflation swaption volatility from year to year volatility
 *	\author  E Benhamou, N. Belgrade
 *	\version 1.0
 *	\date September 2004
 *
 *	\version 2.0
 *	\date March 2005
 */
#ifndef _INGPINFLATION_INFSWOPTFACTORY_H
#define _INGPINFLATION_INFSWOPTFACTORY_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"

/// forward declaration
class ARM_VolCurve;
class ARM_VolCube;
class ARM_Date;


CC_BEGIN_NAMESPACE( ARM )

class ARM_InfBSModel;


class ARM_InfVolComputation_Factory
{
public:
	enum ComputationMethodType
	{
		INF_SWOPT_STD,
		INF_SWOPT_EQUAL_WEIGHT, 
		INF_SWOPT_DF_WEIGHT, 
		INF_SWOPT_DF_WEIGHT_SQUARE, 
	};

	enum choice
	{
		Yes,
		No,
	};
	
	static ARM_VolCurve* GenerateSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_InfBSModel*	infIRBSModel,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ComputationMethodType	method);

	static ARM_VolCurve* GenerateOATSwopVolCurve( 
		const ARM_Date&			asOfDate,
		ARM_InfBSModel*			infIRBSModel,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ComputationMethodType	method,
		double					coupon,
		long					choice);
};

//YK
//for vol cube
class ARM_InfVolCubeComputation_Factory
{
public:
	static ARM_VolCube* GenerateSwopVolCube( 
		const ARM_Date&			asOfDate,
		ARM_InfBSModel*			infIRBSModel,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ARM_GP_Vector*			smiledTenors,
		ARM_GP_Vector*			strikes,
		ARM_InfVolComputation_Factory::ComputationMethodType	method);
};
//YK

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
