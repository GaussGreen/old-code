#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

    /* Compute x^p, x complex and p complex */
#ifdef ANSI
Mf_complex imsl_cc_power(Mf_complex x, Mf_complex p)
#else
Mf_complex imsl_cc_power(x, p)
    Mf_complex	x;
    Mf_complex	p;
#endif
{
	Mf_complex    temp1;
	Mf_complex    ans;
	Mfloat	    temp2;
	Mfloat	    r;
	Mfloat	    absx;

	temp2 = imsl_c_abs(x);
	if (temp2 == F_ZERO) {
	    temp1.re = F_ZERO;
	    temp1.im = F_ZERO;
	    return(temp1);
	}
	temp1.re = F_ZERO;
	absx     = imsl_c_arg(x);
	temp1.im = p.im*log(temp2) + p.re*absx;
	temp1    = imsl_c_exp(temp1);

	r = fabs(pow(temp2,p.re)*exp(-p.im*absx));
	ans.re = r*temp1.re;
	ans.im = r*temp1.im;

	return(ans);
}
