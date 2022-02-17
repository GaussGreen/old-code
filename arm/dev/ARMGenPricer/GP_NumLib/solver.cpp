/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: solver.cpp,v $
 * Revision 1.1  2004/09/22 10:15:09  ebenhamou
 * Initial revision
 *
 *
 *
 */


/*! \file solver.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#include "gpnumlib/solver.h"

CC_BEGIN_NAMESPACE( ARM )

const double SolverConstant::DefaultXTolerance		= 1e-12;
const double SolverConstant::DefaultFxTolerance     = 1.0e-14;
const double SolverConstant::DefaultGradTolerance   = 1.0e-10;
const double SolverConstant::DoubleTolerance		= 1e-12;
const size_t SolverConstant::DefaultMax_Iter		= 100;

const double SolverConstant::DefaultDichoXTolerance = 1e-4;
const size_t SolverConstant::DefaultDichoMax_Iter	= 10;

const double SolverConstant::UpperInfiniteBound		= +1e+20;
const double SolverConstant::LowerInfiniteBound		= -1e+20;
const double SolverConstant::DefaultStepMax			= 100.;
const bool   SolverConstant::DefaultPrint			= false;


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

