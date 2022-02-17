/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file smothingfunc.h
 *
 *  \brief file for the smoothing or the constraint function
 *		in the calibration problem
 *	\author  E.M Ezzine
 *	\version 1.0
 *	\date December 2005
 */


#ifndef _INGPCALIB_SMOOTHINGFUNC_H
#define _INGPCALIB_SMOOTHINGFUNC_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/functor.h"
#include "gpbase/port.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"

#include "typedef.h"
#include <functional>



CC_BEGIN_NAMESPACE( ARM )


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
