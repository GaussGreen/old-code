#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_vmax(Mint ncount, va_list argptr);
#else
static VA_LIST_HACK l_vmax();
#endif

static Mint lv_max;



#ifdef ANSI
Mint imsl_i_max(Mint a, Mint b)
#else
Mint imsl_i_max(a, b)
    Mint  a;
    Mint  b;
#endif
{
    return (a > b) ? a : b;
}



#ifdef ANSI
Mint imsl_i_vmax(Mint ncount, ...)
#else
Mint imsl_i_vmax(ncount, va_alist)
    Mint        ncount;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, ncount);

    imsl_e1psh("imsl_i_vmax");

    lv_max = 0;
    IMSL_CALL(l_vmax(ncount, argptr));
    va_end(argptr);

    imsl_e1pop("imsl_i_vmax");

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
    Mint	    x;

    if (ncount < 1) {
	lv_max = 0;
	return (argptr);
    }

    arg_number++;
    lv_max = va_arg(argptr, Mint);

    for (k=1;  k<ncount;  k++) {
	arg_number++;
	x = va_arg(argptr, Mint);
	if (lv_max < x) lv_max = x;
    }
    return (argptr);
}
