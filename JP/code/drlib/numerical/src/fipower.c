#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

    /* Compute x^k, x real and k an integer */
#ifdef ANSI
Mfloat imsl_fi_power(Mfloat x_arg, Mint k)
#else
Mfloat imsl_fi_power(x_arg, k)
    Mfloat	x_arg;
    Mint	k;
#endif
{
    Mint	n;
    Mfloat	y;
    Mfloat	z;
    Mfloat	x = x_arg;

    n = abs(k);
    if (k<0) x = F_ONE/x;
    switch (n) {
	case 0:
	    y = F_ONE;
	    break;
	case 1:
	    y = x;
	    break;
	case 2:
	    y = x*x;
	    break;
	case 3:
	    y = x*x;
	    y *= x;
	    break;
	case 4:
	    y = x*x;
	    y *= y;
	    break;
	case 5:
	    y = x*x;
	    y *= x*y;
	    break;
	case 6:
	    y = x*x;
	    y *= y*y;
	    break;
	default:
		/* see Knuth, seminumerical Algorithms, page 399-400 */
		/* Algorithm A */
	    y = F_ONE;
	    z = x;
	    for (;;) {
		if (n & 0x1 == 0x1) y *= z;
		n = n >> 1;
		if (n == 0) goto RETURN;
		z *= z;
	    }
    }
    RETURN:
	return y;
}
