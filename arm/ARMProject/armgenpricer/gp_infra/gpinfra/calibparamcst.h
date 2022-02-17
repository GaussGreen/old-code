/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file calibparamcst.h
 *  \brief file for the constant for all calib params
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#ifndef _INGPCALIB_CALIBPARAMCST_H
#define _INGPCALIB_CALIBPARAMCST_H

#include "gpbase/port.h"
#include "gpinfra/modelparamtype.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_CalibParamCst
{
	static double ModelParamStdLowerBound( ARM_ModelParamType::ParamNb paramNb );
	static double ModelParamStdUpperBound( ARM_ModelParamType::ParamNb paramNb );
	static const double CalibLowerBound;
	static const double CalibUpperBound;
	static const double CalibZeroBound;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

