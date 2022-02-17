#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

    /* Compute x^p, x real and p real */
#ifdef ANSI
Mfloat imsl_ff_power(Mfloat x_arg, Mfloat p_arg)
#else
Mfloat imsl_ff_power(x_arg, p_arg)
    Mfloat	x_arg;
    Mfloat	p_arg;
#endif
{
#ifndef DOUBLE
    Mfloat	x = x_arg;
    Mfloat	p = p_arg;

#if 0
#ifdef sparc
#define SPOW
    union {double df;  float rf} uf;
    uf.df   = r_pow_(&x, &p);
    return  uf.rf;
#endif

#ifdef mc68020
#define SPOW
    float val;
    *(int *)(&val) = r_pow_(&x, &p);
    return val;
#endif

#endif	/* not DOUBLE */
#endif  /* 0 */

#ifndef SPOW
    return pow(x_arg, p_arg);
#endif
}
