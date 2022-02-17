#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

    /* Compute x^p, x complex and p real */
#ifdef ANSI
Mf_complex imsl_cf_power(Mf_complex x, Mfloat p_arg)
#else
Mf_complex imsl_cf_power(x, p_arg)
    Mf_complex	x;
    Mfloat	p_arg;
#endif
{
	Mf_complex    temp1, temp3;
	Mfloat	    temp2;
	Mfloat	    p = p_arg;

	temp2 = imsl_c_abs(x);
	if (temp2 == F_ZERO) return(x);
	temp1.re = F_ZERO;
	temp1.im = p*imsl_c_arg(x);
	temp1    = imsl_c_exp(temp1);
	temp3.re = pow(temp2,p)*temp1.re;
	temp3.im = pow(temp2,p)*temp1.im;
	return(temp3);
}
