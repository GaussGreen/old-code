/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 * $Log: solver.cpp,v $
 * Revision 1.1  2004/09/22 10:15:09  aschauly
 * Initial revision
 *
 *
 *
 */


/*! \file solver.cpp
 *
 *  \brief
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnumlib/optimizer.h"

CC_BEGIN_NAMESPACE( ARM )

const double OptimizerConstant::DefaultXTolerance		= 1e-8;
const double OptimizerConstant::DefaultFxTolerance		= 1e-8;
const size_t OptimizerConstant::DefaultMax_Iter			= 100;
const double OptimizerConstant::DefaultStepMax          =2.0;
const double OptimizerConstant::UpperInfiniteBound		= +1e+20;
const double OptimizerConstant::LowerInfiniteBound		= -1e+20;
const double OptimizerConstant::Zero					= 0;
const double OptimizerConstant::DoubleTolerance			= 1e-14;
const bool	 OptimizerConstant::DefaultLocalSearch		= false;
const bool   OptimizerConstant::DefaultPrint			= false;

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

