#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_bessel_Jx(Mfloat xnu, Mf_complex z,
                Mint n, va_list argptr);
static VA_LIST_HACK l_bessel_Jx_adr(Mfloat *xnu, Mf_complex *z,
                Mint n, va_list argptr);
static void l_cbjs(Mfloat *xnu, Mf_complex *z, Mint *n, Mf_complex cbs[]);
#else
static VA_LIST_HACK l_bessel_Jx();
static VA_LIST_HACK l_bessel_Jx_adr();
static void l_cbjs();
#endif

static Mf_complex *lv_cbs = NULL;
#ifdef ANSI
Mf_complex *imsl_c_bessel_Jx(Mfloat xnu, Mf_complex z,
		Mint n, ...)
#else
Mf_complex *imsl_c_bessel_Jx(xnu, z, n, va_alist)
    Mfloat 	xnu;
    Mf_complex 	z;
    Mint        n;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Jx", "imsl_z_bessel_Jx");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Jx(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Jx", "imsl_z_bessel_Jx");

    return(lv_cbs);
}

#ifdef ANSI
Mf_complex *imsl_c_bessel_Jx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, ...)
#else
Mf_complex *imsl_c_bessel_Jx_adr(xnu, z, n, va_alist)
    Mfloat 	*xnu;
    Mf_complex 	*z;
    Mint        n;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Jx_adr", "imsl_z_bessel_Jx_adr");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Jx_adr(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Jx_adr", "imsl_z_bessel_Jx_adr");

    return(lv_cbs);
}

#ifdef ANSI
static VA_LIST_HACK l_bessel_Jx(Mfloat xnu, Mf_complex z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Jx(xnu, z, n, argptr)
    Mfloat 	xnu;
    Mf_complex 	z;
    Mint        n;
    va_list     argptr;
#endif
{
    Mfloat float_xnu;
    Mint        code, user_cbs = 0, ner, arg_number = 3;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_cbs = va_arg(argptr, Mf_complex*);
        user_cbs = 1;
    }
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

    if (n <= 0) {
        imsl_e1sti (1, n);

        imsl_e1mes (5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");
    }
    if (!user_cbs) {
        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));
        if (!lv_cbs){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    float_xnu = xnu;
    l_cbjs (&float_xnu, &z, &n, lv_cbs);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }

RETURN:

    return(argptr);
}


#ifdef ANSI
static VA_LIST_HACK l_bessel_Jx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Jx_adr(xnu, z, n, argptr)
    Mfloat 	*xnu;
    Mf_complex 	*z;
    Mint        n;
    va_list     argptr;
#endif
{
    Mfloat float_xnu;
    Mint        code, user_cbs = 0, ner, arg_number = 3;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_cbs = va_arg(argptr, Mf_complex*);
        user_cbs = 1;
    }
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

    if (n <= 0) {
        imsl_e1sti (1, n);

        imsl_e1mes (5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");
    }
    if (!user_cbs) {
        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));
        if (!lv_cbs){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    float_xnu = *xnu;
    l_cbjs (&float_xnu, z, &n, lv_cbs);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }

RETURN:

    return(argptr);
}




/*Translated by FOR_C++, v0.1, on 06/06/91 at 15:21:52 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/06/91 at 15:21:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CBJS/DCBJS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    June 7, 1991

    Purpose:    Evaluate a sequence of Bessel functions of the first
                kind with real order and d_complex arguments.

    Usage:      CALL CBJS (XNU, Z, N, CBS)

    Arguments:
       XNU    - Real argument which is the lowest order desired.
                (Input)
                XNU must be greater than -1/2.
       Z      - Complex argument for which the sequence of Bessel
                functions is to be evaluated.  (Input)
       N      - Number of elements in the sequence.  (Input)
       CBS    - Vector of length N containing the values of the
                function through the series.  (Output)
                CBS(I) contains the value of the Bessel function of
                order XNU+I-1 at Z for I = 1 to N.

    Remark:
       Informational errors
       Type Code
         3   1  One the continued fractions failed.
         4   2  Only the first several entries in CBS are valid.

    GAMS:       C10a4

    Chapter:    MATH/LIBRARY SPECIAL FUNCTIONS Bessel Functions

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_cbjs(Mfloat *xnu, Mf_complex *z, Mint *n, Mf_complex cbs[])
#else
static void l_cbjs(xnu, z, n, cbs)
	Mfloat          *xnu;
	Mf_complex      *z;
	Mint            *n;
	Mf_complex       cbs[];
#endif
{
	Mint             _l0, k, mode;
	Mfloat           t, x, y;
	Mf_complex       _cx0, fip[1], fk[1], fkp[1], zfact, zi;
	static Mfloat    pi2 = 1.57079632679489661923132169164e0;



	imsl_e1psh("CBJS ");
	/*
	 * Use Abramowitz & Stegun equation 9.6.3, imsl_page 375 to relate J
	 * to I.
	 */
	_cx0.re = 0.0;
	_cx0.im = 0.0;
	_l0 = 1;
	imsl_cset(n, &_cx0, cbs, &_l0);	
	x = imsl_fc_convert(*z);
	y = imsl_c_aimag(*z);
/* orig 	zi = imsl_c_ftocf(y, -x); */
        zi = imsl_cf_convert(y,-x);
	mode = 4;
	imsl_c3is(&zi, xnu, n, cbs, fk, fip, fkp, &mode);
	if (imsl_n1rty(1) >= 4)
		goto L_9000;

	if (x < 0.0 && y < 0.0) {
		for (k = 1; k <= *n; k++) {
			t = -3.0 * pi2 * (*xnu + k - 1);
/*orig			zfact = imsl_c_ftocf(cos(t), sin(t));*/
                        zfact = imsl_cf_convert(cos(t),sin(t));
			cbs[k - 1] = imsl_c_mul(zfact, cbs[k - 1]);
		}
	} else {
		for (k = 1; k <= *n; k++) {
			t = pi2 * (*xnu + k - 1);
/*orig			zfact = imsl_c_ftocf(cos(t), sin(t));*/
                        zfact = imsl_cf_convert(cos(t),sin(t));
			cbs[k - 1] = imsl_c_mul(zfact, cbs[k - 1]);
		}
	}

L_9000:
	imsl_e1pop("CBJS ");

	return;
}				/* end of function */
