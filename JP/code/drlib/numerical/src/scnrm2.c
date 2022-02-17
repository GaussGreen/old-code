#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* INCX is not used, but leave in calling sequence. */
#ifdef ANSI
Mfloat imsl_scnrm2(Mint *n, Mf_complex *x, Mint *incx)
#else
Mfloat imsl_scnrm2(n, x, incx)
Mint *n;
Mf_complex x[];
Mint *incx;
#endif
{
	Mfloat sum = 0;
	Mint i;

	for(i=0; i<=*n-1; ++i)
		sum = sum + x[i].re*x[i].re + x[i].im*x[i].im;

	return sqrt(sum);
}
