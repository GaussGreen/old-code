/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file incompletebeta.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_INCOMPLETEBETA_H
#define _GP_CF_INCOMPLETEBETA_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"

#include "long_double.h"
#include <complex>

using std::complex;
//using std::complex<double>;

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////
///
///    incomplete beta function (x,a,b) 
///     = 1/Beta(a,b)  Integral({t,0,x}, t^(a-1) (1-t)^(b-1) )
///
/////////////////////////////////////////////////////////////////////////


double IncompleteBeta(const double& a, const double& b, const double& x);
double IncompleteBeta(const double& a, const double& b, const double& x1,const double& x2);
/////////////////////////////////////////////////////////////////////////
///
///    incomplete beta inverse function (x,a,b) 
///     solution in y of IncompleteBeta(a,b,y)==x
///
/////////////////////////////////////////////////////////////////////////


double IncompleteBetaInverse(const double& a, const double& b, const double& y);
double IncompleteBetaInverse(const double& a, const double& b, const double& y0,const double& y);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

