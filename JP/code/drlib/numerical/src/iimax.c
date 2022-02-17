#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

    /* Compute x^k, x integer and k an integer */
#ifdef ANSI
Mint imsl_ii_power(Mint x, Mint k)
#else
Mint imsl_ii_power(x, k)
    Mint	x;
    Mint	k;
#endif
{
    Mint	y;
    Mint	z;

    /********************* if (k<0) ERROR *********************/
    switch (k) {
	case 0:
	    y = 1;
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
	    y = x*y*y;
	    break;
	case 6:
	    y = x*x;
	    y *= y*y;
	    break;
	default:
		/* see Knuth, seminumerical Algorithms, page 399-400 */
		/* Algorithm A */
	    y = 1;
	    z = x;
	    for (;;) {
		if (k & 0x1 == 0x1) y *= z;
		k = k >> 1;
		if (k == 0) goto RETURN;
		z *= z;
	    }
    }
    RETURN:
	return y;
}
