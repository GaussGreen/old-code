/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file polynomial_rootfinding.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date fev 2007
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/polynomial_rootfinding.h"

#include "gpbase/numericconstant.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000



/////////////////////////////////////////////////////////////////////////
///
///    
///
/////////////////////////////////////////////////////////////////////////


CC_END_NAMESPACE()


 
#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/