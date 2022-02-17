/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file lambert_function.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date September 2005
 */
 
#ifndef _GP_CF_LAMBERT_FUNCTION_H
#define _GP_CF_LAMBERT_FUNCTION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "long_double.h"
#include <complex>

using std::complex;
//using std::complex<double>;

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////
///
///     solution in w of  w exp(w) == x
///
/////////////////////////////////////////////////////////////////////////

double lambertfunc (double z);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

