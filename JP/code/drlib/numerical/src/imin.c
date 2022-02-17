#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_vmin(Mint ncount, va_list argptr);
#else
static VA_LIST_HACK l_vmin();
#endif

static Mint lv_min;



#ifdef ANSI
Mint imsl_i_min(Mint a, Mint b)
#else
Mint imsl_i_min(a, b)
    Mint  a;
    Mint  b;
#endif
{
    return (a < b) ? a : b;
}



#ifdef ANSI
Mint imsl_i_vmin(Mint ncount, ...)
#else
Mint imsl_i_vmin(ncount, va_alist)
    Mint        ncount;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, ncount);

    imsl_e1psh("imsl_i_vmin");

    lv_min = 0;
    IMSL_CALL(l_vmin(ncount, argptr));
    va_end(argptr);

    imsl_e1pop("imsl_i_vmin");

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
    Mint	    x;

    if (ncount < 1) {
	lv_min = 0;
	return (argptr);
    }

    arg_number++;
    lv_min = va_arg(argptr, Mint);

    for (k=1;  k<ncount;  k++) {
	arg_number++;
	x = va_arg(argptr, Mint);
	if (lv_min > x) lv_min = x;
    }
    return (argptr);
}
