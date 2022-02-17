#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_bessel_Kx(Mfloat xnu, Mf_complex z,
                Mint n, va_list argptr);
static VA_LIST_HACK l_bessel_Kx_adr(Mfloat *xnu, Mf_complex *z,
                Mint n, va_list argptr);
static void l_c2ks(Mfloat *xnu, Mf_complex *z, Mint *n,
		Mf_complex cbs[], Mf_complex fi[]);
#else
static VA_LIST_HACK l_bessel_Kx();
static VA_LIST_HACK l_bessel_Kx_adr();
static void l_c2ks();
#endif

static Mf_complex *lv_cbs = NULL;
#ifdef ANSI
Mf_complex *imsl_c_bessel_Kx(Mfloat xnu, Mf_complex z,
		Mint n, ...)
#else
 Mf_complex *imsl_c_bessel_Kx(xnu, z, n, va_alist)
     Mfloat 	xnu;
    Mf_complex 	z;
    Mint        n;
    va_dcl
#endif
/*                        User callable section */
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Kx", "imsl_z_bessel_Kx");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Kx(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Kx", "imsl_z_bessel_Kx");

    return(lv_cbs);
}
#ifdef ANSI
Mf_complex *imsl_c_bessel_Kx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, ...)
#else
 Mf_complex *imsl_c_bessel_Kx_adr(xnu, z, n, va_alist)
     Mfloat 	*xnu;
    Mf_complex 	*z;
    Mint        n;
    va_dcl
#endif
/*                        User callable section */
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Kx_adr", "imsl_z_bessel_Kx_adr");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Kx_adr(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Kx_adr", "imsl_z_bessel_Kx_adr");

    return(lv_cbs);
}
/*                          End User callable section */
/*                          This routine only has one optional argument.*/
#ifdef ANSI
static VA_LIST_HACK l_bessel_Kx(Mfloat xnu, Mf_complex z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Kx(xnu, z, n, argptr)
    Mfloat 	xnu;
    Mf_complex 	z;
    Mint        n;
    va_list     argptr;
#endif
{
    Mfloat float_xnu;
    Mint        code, user_cbs = 0, arg_number = 3;
    Mf_complex      *cwork      = NULL;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_cbs = va_arg(argptr, Mf_complex*);
        user_cbs = 1;
    }
/*                               Not a valid variable argument */
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

	if (n < 1) {
		imsl_e1sti(1, n);

		imsl_e1mes(5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");
/*		goto L_9000; */
                goto RETURN;
	}


    if (!user_cbs) {
        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));
        if (!lv_cbs){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    cwork       = (Mf_complex *) imsl_malloc(n*sizeof(*cwork));
/* make sure you free this space */
    if (cwork==NULL ) {
        imsl_e1stl(1, "n");
        imsl_e1sti(1, n);
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto FREE_SPACE; 
    }
    float_xnu = xnu;
    l_c2ks (&float_xnu, &z, &n, lv_cbs, cwork);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }
/*                         Only free space actually allocated. */
FREE_SPACE:
    if (cwork != NULL) imsl_free(cwork); 

RETURN:

    return(argptr);
}

/*                          This routine only has one optional argument.*/
#ifdef ANSI
static VA_LIST_HACK l_bessel_Kx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Kx_adr(xnu, z, n, argptr)
    Mfloat 	*xnu;
    Mf_complex 	*z;
    Mint        n;
    va_list     argptr;
#endif
{
    Mfloat float_xnu;
    Mint        code, user_cbs = 0, arg_number = 3;
    Mf_complex      *cwork      = NULL;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_cbs = va_arg(argptr, Mf_complex*);
        user_cbs = 1;
    }
/*                               Not a valid variable argument */
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

	if (n < 1) {
		imsl_e1sti(1, n);

		imsl_e1mes(5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");
/*		goto L_9000; */
                goto RETURN;
	}


    if (!user_cbs) {
        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));
        if (!lv_cbs){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    cwork       = (Mf_complex *) imsl_malloc(n*sizeof(*cwork));
/* make sure you free this space */
    if (cwork==NULL ) {
        imsl_e1stl(1, "n");
        imsl_e1sti(1, n);
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto FREE_SPACE; 
    }
    float_xnu = *xnu;
    l_c2ks (&float_xnu, z, &n, lv_cbs, cwork);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }
/*                         Only free space actually allocated. */
FREE_SPACE:
    if (cwork != NULL) imsl_free(cwork); 

RETURN:

    return(argptr);
}



/*Translated by FOR_C++, v0.1, on 06/10/91 at 14:15:22 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/10/91 at 14:15:20
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C2KS/DC2KS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    June 10, 1991

    Purpose:    Evaluate a sequence of Modified Bessel functions of
                the second kind with real order and d_complex arguments.

    Usage:      CALL C2KS (XNU, Z, N, CBS, FI)

    Arguments:
       XNU    - Real argument which is the lowest order desired.
                (Input)
                XNU must be greater than -1/2.
       Z      - Complex argument for which the sequence of Bessel
                functions is to be evaluated.  (Input)
       N      - Number of elements in the sequence.  (Input)
       CBS    - Vector of length N containing the values of the
                function through the series.  (Output)
       FI     - Complex vector of length N used as workspace.  (Output)

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2ks(Mfloat *xnu, Mf_complex *z, Mint *n,
		Mf_complex cbs[], Mf_complex fi[])
#else
static void l_c2ks(xnu, z, n, cbs, fi)
	Mfloat          *xnu;
	Mf_complex      *z;
	Mint            *n;
	Mf_complex       cbs[], fi[];
#endif
{
	Mint             mode;
	Mf_complex       fip[1], fkp[1];


	mode = 2;
	imsl_c3is(z, xnu, n, fi, cbs, fip, fkp, &mode);

	return;
}				/* end of function */
