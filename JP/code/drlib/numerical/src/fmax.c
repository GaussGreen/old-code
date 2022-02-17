#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_vmax(Mint ncount, va_list argptr);
#else
static VA_LIST_HACK l_vmax();
#endif

static Mfloat lv_max;



#ifdef ANSI
Mfloat imsl_f_max(Mfloat a, Mfloat b)
#else
Mfloat imsl_f_max(a, b)
    Mfloat  a;
    Mfloat  b;
#endif
{
    return (a > b) ? a : b;
}



#ifdef ANSI
Mfloat imsl_f_vmax(Mint ncount, ...)
#else
Mfloat imsl_f_vmax(ncount, va_alist)
    Mint        ncount;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, ncount);

    E1PSH("imsl_f_vmax", "imsl_d_vmax");

    lv_max = F_ZERO;
    IMSL_CALL(l_vmax(ncount, argptr));
    va_end(argptr);

    E1POP("imsl_f_vmax", "imsl_d_vmax");

    return lv_max;
}


#ifdef ANSI
static VA_LIST_HACK l_vmax(Mint ncount, va_list argptr)
#else
static VA_LIST_HACK l_vmax(ncount, argptr)
    Mint        ncount;
    va_list	argptr;
#endif
{
    Mint	    k;
    Mint	    arg_number = 1;
    Mfloat	    x;

    if (ncount < 1) {
	lv_max = F_ZERO;
	return (argptr);
    }

    arg_number++;
    lv_max = va_arg(argptr, Mdouble);

    for (k=1;  k<ncount;  k++) {
	arg_number++;
	x = va_arg(argptr, Mdouble);
	if (lv_max < x) lv_max = x;
    }
    return (argptr);
}
