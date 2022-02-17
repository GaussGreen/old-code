/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Complex Numbers
 * Revision 1.1  2003/10/08 16:45:06  atriki
 * Initial revision
 *
 *
 *
 */
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/cplx.h"

#include <complex>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: robust_quotient
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_Cplx::robust_quotient (double re1, double im1, double re2, double im2, /* result-> */ double& re_quot,  double& im_quot)
{
	std::complex<double> z1 (re1, im1);
	std::complex<double> z2 (re2, im2);
	std::complex<double> z = z1/z2;
	re_quot = z.real();
	im_quot = z.imag();
}

CC_END_NAMESPACE()
