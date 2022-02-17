/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file nummultidimfunction.cpp
 *
 *  \brief template function to compute easily the gradient
 *		of multi dimensional function
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date Ocotber 2004
 */


#include "gpnumlib/nummultidimfunction.h"

/// kernel
#include "expt.h"
#include "gpbase/utilityport.h"

/// standard libraries
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : GradientFunction
///	Routine: operator()( const std::vector<double>& x ) const
///	Returns: ARM_GP_Matrix 
///	Action : computes the gradient of a multi-dimensional function
////////////////////////////////////////////////////




CC_END_NAMESPACE()

/*----------------------------------------------------------------------*/
/*---- End of file ----*/
