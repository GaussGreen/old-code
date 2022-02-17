/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_calibration.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"


#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"


#include "gpclosedforms/optimization1.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_TOLERANCE 1.0e-13
#define ARM_CF_MAXIT 2000
#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.05





CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
