#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_bessel_Yx(Mfloat xnu, Mf_complex z,
                Mint n, va_list argptr);
static VA_LIST_HACK l_bessel_Yx_adr(Mfloat *xnu, Mf_complex *z,
                Mint n, va_list argptr);
static void l_c2ys(Mfloat *xnu, Mf_complex *z, Mint *n,
		Mf_complex cbs[], Mf_complex fk[]);
#else
static VA_LIST_HACK l_bessel_Yx();
static VA_LIST_HACK l_bessel_Yx_adr();
static void l_c2ys();
#endif

static Mf_complex *lv_cbs = NULL;
#ifdef ANSI
Mf_complex *imsl_c_bessel_Yx(Mfloat xnu, Mf_complex z,
		Mint n, ...)
#else
Mf_complex *imsl_c_bessel_Yx(xnu, z, n, va_alist)
    Mfloat 	xnu;
    Mf_complex 	z;
    Mint        n;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Yx", "imsl_z_bessel_Yx");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Yx(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Yx", "imsl_z_bessel_Yx");

    return(lv_cbs);
}

#ifdef ANSI
Mf_complex *imsl_c_bessel_Yx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, ...)
#else
Mf_complex *imsl_c_bessel_Yx_adr(xnu, z, n, va_alist)
    Mfloat 	*xnu;
    Mf_complex 	*z;
    Mint        n;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n);

    E1PSH("imsl_c_bessel_Yx_adr", "imsl_z_bessel_Yx_adr");

    lv_cbs = NULL;
    IMSL_CALL(l_bessel_Yx_adr(xnu, z, n, argptr));
    va_end(argptr);

    E1POP("imsl_c_bessel_Yx_adr", "imsl_z_bessel_Yx_adr");

    return(lv_cbs);
}

#ifdef ANSI
static VA_LIST_HACK l_bessel_Yx(Mfloat xnu, Mf_complex z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Yx(xnu, z, n, argptr)
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
    l_c2ys (&float_xnu, &z, &n, lv_cbs, cwork);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }

FREE_SPACE:
    if (cwork != NULL) imsl_free(cwork); 

RETURN:

    return(argptr);
}

#ifdef ANSI
static VA_LIST_HACK l_bessel_Yx_adr(Mfloat *xnu, Mf_complex *z,
		Mint n, va_list argptr)
#else
static VA_LIST_HACK l_bessel_Yx_adr(xnu, z, n, argptr)
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
    l_c2ys (&float_xnu, z, &n, lv_cbs, cwork);
    if (imsl_n1rty(0)>3 && !user_cbs) {
        imsl_free(lv_cbs);
        lv_cbs = NULL;
    }

FREE_SPACE:
    if (cwork != NULL) imsl_free(cwork); 

RETURN:

    return(argptr);
}

/*Translated by FOR_C++, v0.1, on 06/10/91 at 12:33:22 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/10/91 at 12:33:20
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C2YS/DC2YS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    January 6, 1990

    Purpose:    Evaluate a sequence of Bessel functions of the second
                kind with real order and d_complex arguments.

    Usage:      CALL C2YS (XNU, Z, N, CBS, FK)

    Arguments:
       XNU    - Real argument which is the lowest order desired.
                (Input)
                XNU must be greater than -1/2.
       Z      - Complex argument for which the sequence of Bessel
                functions is to be evaluated.  (Input)
       N      - Number of elements in the sequence.  (Input)
       CBS    - Vector of length N containing the values of the
                function through the series.  (Output)
       FK     - Complex vector of length N used as workspace.  (Output)

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2ys(Mfloat *xnu, Mf_complex *z, Mint *n,
		Mf_complex cbs[], Mf_complex fk[])
#else
static void l_c2ys(xnu, z, n, cbs, fk)
	Mfloat          *xnu;
	Mf_complex      *z;
	Mint            *n;
	Mf_complex       cbs[], fk[];
#endif
{
	Mint       quad4;
	Mint             k, mode;
	Mfloat           t, x, y;
	Mf_complex       fact, facti, factk, fip[1], fkp[1], zi;
        Mf_complex       ctmp1, ctmp2, ctmp3;
        Mfloat           den;
	static Mfloat    pi2 = 1.57079632679489661923132169164e0;



	imsl_e1psh("C2YS ");
	/*
	 * Use Abramowitz & Stegun equation 9.6.5, imsl_page 375 to relate Y
	 * to I and K.
	 */
#if 0
	x = imsl_fc_convert(*z);
	y = imsl_c_aimag(*z);
	zi = imsl_cf_convert(y, -x);
#endif
	x = (*z).re;
	y = (*z).im;
        zi.re =  y;
        zi.im = -x;

	/*
	 * Use conjugate formula for points in the fourth quadrant
	 */
	quad4 = x < 0.0 && y < 0.0;
#if 0
	if (quad4)
		zi = imsl_cf_convert(-y, x);
#endif
	if (quad4) {
	        zi.re = -y;
                zi.im =  x;
	      }

	mode = 2;
	imsl_c3is(&zi, xnu, n, cbs, fk, fip, fkp, &mode);
	if (imsl_n1rty(1) >= 4)
		goto L_9000;

	for (k = 1; k <= *n; k++) {
		t = pi2 * (*xnu + k - 1);
#if 0
		fact = imsl_cf_convert(cos(t), sin(t));
		factk = imsl_c_neg(imsl_c_div(imsl_c_conjg(fact), imsl_cf_convert(pi2, 0.)));
		facti = imsl_cf_convert(-imsl_c_aimag(fact), imsl_fc_convert(fact));
		if (quad4) {
			facti = imsl_c_conjg(facti);
			factk = imsl_c_conjg(factk);
		}
		cbs[k - 1] = imsl_c_add(imsl_c_mul(facti, cbs[k - 1]), imsl_c_mul(factk, fk[k - 1]));
#endif
		fact.re = cos(t);
		fact.im = sin(t);
                ctmp1.re    =  pi2;
                ctmp1.im    =  0.0;
                ctmp2.re    =  fact.re;
                ctmp2.im    = -fact.im;
                if( ctmp1.re == F_ZERO && ctmp1.im == F_ZERO ){
                    imsl_ermes(IMSL_TERMINAL, IMSL_COMPLEX_DIVIDE_BY_ZERO);
                    return;
                }
        
                if (fabs(ctmp1.re) > fabs(ctmp1.im)){
                    t = ctmp1.im/ctmp1.re;
                    den = ctmp1.re + ctmp1.im*t;
                    ctmp3.re = (ctmp2.re + ctmp2.im*t)/den;
                    ctmp3.im = (ctmp2.im - ctmp2.re*t)/den;
                }
                else{
                    t = ctmp1.re/ctmp1.im;
                    den = ctmp1.im + ctmp1.re*t;
                    ctmp3.re = (ctmp2.im + ctmp2.re*t)/den;
                    ctmp3.im = (ctmp2.im*t - ctmp2.re)/den;
                }
                factk.re    = -ctmp3.re;
                factk.im    = -ctmp3.im;
                

		facti.re = -fact.im;
		facti.im =  fact.re;
		if (quad4) {
			facti.im = -facti.im;
			factk.im = -factk.im;
		}
                ctmp1.re      = factk.re*((fk[k-1]).re) - factk.im*((fk[k-1]).im);
                ctmp1.im      = factk.re*((fk[k-1]).im) + factk.im*((fk[k-1]).re);
                ctmp2.re      = facti.re*((cbs[k-1]).re) - facti.im*((cbs[k-1]).im);
                ctmp2.im      = facti.re*((cbs[k-1]).im) + facti.im*((cbs[k-1]).re);
                (cbs[k-1]).re = ctmp1.re + ctmp2.re;
                (cbs[k-1]).im = ctmp1.im + ctmp2.im;
	}

L_9000:
	imsl_e1pop("C2YS ");

	return;
}				/* end of function */
