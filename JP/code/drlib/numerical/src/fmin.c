#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_vmin(Mint ncount, va_list argptr);
#else
static VA_LIST_HACK l_vmin();
#endif

static Mfloat lv_min;



#ifdef ANSI
Mfloat imsl_f_min(Mfloat a, Mfloat b)
#else
Mfloat imsl_f_min(a, b)
    Mfloat  a;
    Mfloat  b;
#endif
{
    return (a < b) ? a : b;
}



#ifdef ANSI
Mfloat imsl_f_vmin(Mint ncount, ...)
#else
Mfloat imsl_f_vmin(ncount, va_alist)
    Mint        ncount;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, ncount);

    E1PSH("imsl_f_vmin", "imsl_d_vmin");

    lv_min = F_ZERO;
    IMSL_CALL(l_vmin(ncount, argptr));
    va_end(argptr);

    E1POP("imsl_f_vmin", "imsl_d_vmin");

    return lv_min;
}


#ifdef ANSI
static VA_LIST_HACK l_vmin(Mint ncount, va_list argptr)
#else
static VA_LIST_HACK l_vmin(ncount, argptr)
    Mint        ncount;
    va_list	argptr;
#endif
{
    Mint	    k;
    Mint	    arg_number = 1;
    Mfloat	    x;

    if (ncount < 1) {
	lv_min = F_ZERO;
	return (argptr);
    }

    arg_number++;
    lv_min = va_arg(argptr, Mdouble);

    for (k=1;  k<ncount;  k++) {
	arg_number++;
	x = va_arg(argptr, Mdouble);
	if (lv_min > x) lv_min = x;
    }
    return (argptr);
}
