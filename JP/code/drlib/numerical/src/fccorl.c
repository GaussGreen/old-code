#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#define ADR(t,x)    ( t = x, &t )
#ifndef COMPUTER_LINUX
#define nint(x) (x < 0 ? (int)((x)-.5) : (int)((x)+.5))
#endif
#ifdef ANSI
static VA_LIST_HACK l_convolution (Mint nx, Mf_complex *x, Mint ny,
                        Mf_complex *y, Mint *nz, va_list argptr);
static void l_c2onv (Mint *, Mint *, Mf_complex [], Mint *, Mf_complex [],
	Mint *, Mint *, Mf_complex [], Mf_complex [], Mf_complex [], 
        Mf_complex [],
	Mfloat []);
static void l_c2orl (Mint *, Mint *, Mf_complex [], Mf_complex [],
	Mint *, Mint *, Mf_complex [], Mf_complex [], Mf_complex [], 
        Mf_complex [],
	Mfloat []); 
static void l_f3tci(Mint *n, Mfloat wa[], Mfloat imsl_fac[]);
static void l_fftci(Mint *n, Mfloat wfftc[]);

#else
static VA_LIST_HACK l_convolution ();
static void l_c2onv ();
static void l_c2orl ();
static void l_fftci ();
static void l_f3tci ();
#endif


static Mf_complex *lv_z = NULL;
static Mf_complex *lv_temp_zhat = NULL;
static Mf_complex **lv_zhat = NULL;
static Mf_complex *lv_zhat_user = NULL;

/* copied declarations from ffftc.c */
extern  void imsl_f2tcf (Mint *n, Mf_complex seq[], Mf_complex coef[], Mfloat wfftc[], Mfloat cpy[]);
extern void imsl_f2tcb (Mint *n, Mf_complex coef[], Mf_complex seq[], Mfloat wfftc[], Mfloat cpy[]);


#ifdef ANSI
Mf_complex     *imsl_c_convolution (Mint nx, Mf_complex *x, Mint ny,
			Mf_complex *y, Mint *nz, ...)
#else
Mf_complex     *imsl_c_convolution (nx, x, ny, y, nz, va_alist)
    Mint        nx;
    Mf_complex     *x;
    Mint        ny;
    Mf_complex     *y;
    Mint       *nz;

va_dcl

#endif
/*                                   User Callable part */
{
    va_list     argptr;
    VA_START (argptr, nz);
    E1PSH ("imsl_c_convolution", "imsl_z_convolution");
    lv_z = NULL;
    lv_temp_zhat = NULL;
    lv_zhat = NULL;
    lv_zhat_user = NULL;
    IMSL_CALL (l_convolution (nx, x, ny, y, nz, argptr));
    va_end (argptr);
    E1POP ("imsl_c_convolution", "imsl_z_convolution");
    return lv_z;
}


#ifdef ANSI
static VA_LIST_HACK l_convolution (Mint nx, Mf_complex *x, Mint ny,
			Mf_complex *y, Mint *nz, va_list argptr)
#else
static VA_LIST_HACK l_convolution (nx, x, ny, y, nz, argptr)
    Mint        nx;
    Mf_complex     *x;
    Mint        ny;
    Mf_complex     *y;
    Mint       *nz;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 5;
    Mint        user_space = 0;
    Mint        return_zhat = 0;
    Mint        zhat_user = 0;
    Mint        periodic = 0;
    Mint        correlation = 0;
    Mint        first_call = 0;
    Mint        continue_call = 0;
    Mint        last_call = 0;
    Mint        ido = 0;
    Mint        ipad;
    Mint        krlse;
    Mint        n;
    Mint        nz2;
    Mint        valold;
    Mint        i;
    Mint        j;
    Mint        k;
    Mint        valnew;
    static 	Mfloat     *irwksp;
    static 	Mf_complex *icwksp;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_z = va_arg (argptr, Mf_complex *);
	    user_space = 1;
	    arg_number++;
	    break;
	case IMSL_PERIODIC:
	    periodic = 1;
	    break;
	case IMSL_CORRELATION:
	    correlation = 1;
	    break;
	case IMSL_Z_TRANS:
	    lv_zhat = va_arg (argptr, Mf_complex **);
	    return_zhat = 1;
	    arg_number++;
	    break;
	case IMSL_Z_TRANS_USER:
	    lv_zhat_user = va_arg (argptr, Mf_complex *);
	    zhat_user = 1;
	    return_zhat = 1;
	    arg_number++;
	    break;
	case IMSL_FIRST_CALL:
	    first_call = 1;
	    break;
	case IMSL_CONTINUE_CALL:
	    continue_call = 1;
	    first_call = 0;
	    break;
	case IMSL_LAST_CALL:
	    last_call = 1;
	    continue_call = 0;
	    break;
	case 0:
	    break;
	default:
	    /* Argument number %(I2) is an unknown */
	    /* optional argument %(I1). */
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }

    if (imsl_n1rty (0))
	goto RETURN;


    if (user_space && (lv_z == NULL)) {
	imsl_e1stl (1, "z");
	imsl_e1stl (2, "IMSL_RETURN_USER");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }

    if (zhat_user && (lv_zhat_user == NULL) ) {
	imsl_e1stl (1, "lv_zhat");
	imsl_e1stl (2, "IMSL_PARAMS");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }

    if (x == NULL) {
	imsl_e1stl (1, "x");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
	goto RETURN;
    }

    if (y == NULL) {
	imsl_e1stl (1, "y");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
	goto RETURN;
    }

    if (first_call)
	ido = 1;
    if (continue_call)
	ido = 2;
    if (last_call)
	ido = 3;

    if (&x[0] != &y[0]){
	if (periodic)
	    ipad = 0;
	else
	    ipad = 1;
    }
    if (&x[0] == &y[0]) {
	if (periodic)
	    ipad = 2;
	else
	    ipad = 3;
    }

    /* Check NX */
    if (nx <= 0) {
	imsl_e1sti (1, nx);
/*
	imsl_e1mes (5, 3, "The length of the vector X must be positive while NX = %(i1) is given.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_X_NOT_POSITIVE);
    }
    /* Check NY */
    if (ny <= 0) {
	imsl_e1sti (1, ny);
/*
	imsl_e1mes (5, 4, "The length of the vector Y must be positive while NY = %(i1) is given.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_Y_NOT_POSITIVE);
    }

    if (imsl_n1rty (0) != 0)
	goto L_9000;

if (!correlation) {


/*  cconv follows here.  This computation is needed to determine
    how much space to allocate before call c2onv */


    if (ipad == 2) ipad = 0;
    if (ipad == 3) ipad = 1;

        /* Compute correct NZ */
        if (ipad == 0) {
                n = imsl_i_max(nx, ny);
        } else {
                nz2 = nx + ny - 1;
                valold = 2.0 * (nz2);
                for (i = 0; i <= nint(log((float) (nz2)) / log(5.0e0)); i++) {
                        for (j = 0; j <= nint((log((float) (nz2) / (imsl_fi_power(5.0e0, i)))) /
                                              (log(3.0e0))); j++) {
                                for (k = nint((log((float) (nz2) / ((imsl_fi_power(5.0e0, i)) *
                                                                    (imsl_fi_power(3.0e0, j))))) / log(2.0e0)); k <= nint((log(((float) (nz2) /
                                                                                                                                (imsl_fi_power(5.0e0, i))) / (imsl_fi_power(3.0e0, j)))) / (log(2.0e0))); k++) {
                                        valnew = (imsl_fi_power(2.0e0, k)) * (imsl_fi_power(3.0e0, j)) * (imsl_fi_power(5.0e0, i));
                                        if ((valnew >= nz2) && (valnew < valold)) {
                                                valold = valnew;
                                                if (valold == (nz2))
                                                        goto L_40;
                                        }
                                }
                        }
                }
        }
L_40:
        if (ipad == 1)
            n = valold;
        *nz = n;
        if (imsl_n1rty(0) != 0)
                goto L_9000;
        /*
         * Allocate workspace based on the minimum value for NZ as computed
         * above and stored in N.  We use N instead of the input NZ since we
         * will be using the minimum value for NZ when padding the vectors X
         * and Y.
         */
        if ((ido == 0) || (ido == 1)) {
            icwksp = (Mf_complex*) imsl_malloc(2*n*sizeof(*icwksp));
            irwksp = (Mfloat*) imsl_malloc ((6 * n + 15) * sizeof (*irwksp));
	}

        if( !zhat_user && return_zhat)
	   *lv_zhat = (Mf_complex*) imsl_malloc (n * sizeof (**lv_zhat));

        if (!return_zhat)
            lv_temp_zhat = (Mf_complex*) imsl_malloc (n*sizeof(*lv_temp_zhat));

        if ( (!zhat_user && return_zhat && *lv_zhat == NULL) || 
                ( (ido==0 || ido==1) && (irwksp==NULL || icwksp==NULL)) || 
                (!return_zhat && lv_temp_zhat == NULL) ) {
	    /* Not enough memory, with %(L1) = %(I1). */
	    imsl_e1sti (1, n);
	    imsl_e1stl (1, "n");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}

	if (!user_space) {
	    lv_z = (Mf_complex *) imsl_malloc (n * sizeof (*lv_z));
	    if (lv_z == NULL) {
		/* Not enough memory, with %(L1) = %(I1). */
		imsl_e1sti (1, n);
		imsl_e1stl (1, "n");
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		goto FREE_SPACE;
	    }
	}

     /* Call computational routine */

    if (zhat_user) 
        l_c2onv(&ido, &nx, x, &ny, y, &ipad, nz, lv_z, lv_zhat_user, 
                icwksp, icwksp + n , irwksp);
    else if (return_zhat)
        l_c2onv(&ido, &nx, x, &ny, y, &ipad, nz, lv_z, *lv_zhat, 
                icwksp, icwksp + n , irwksp);
    else 
        l_c2onv(&ido, &nx, x, &ny, y, &ipad, nz, lv_z, lv_temp_zhat, 
                icwksp, icwksp + n , irwksp);

    } /* end if convolution */

    else /* correlation */ {



/*  ccorl follows here.  This computation is needed to determine
    how much space to allocate before call c2orl */

        /* Compute correct NZ */
/*    Put this in */
        n = nx;            
        nz2 = n;
        if ((ipad == 1) || (ipad == 3)) {
                nz2 = 2*n - 1;
                valold = 2.0 * (nz2);
                for (i = 0; i <= nint(log((float) (nz2)) / log(5.0e0)); i++) {
                        for (j = 0; j <= nint((log((float) (nz2) / (imsl_fi_power(5.0e0, i)))) /
                                              (log(3.0e0))); j++) {
                                for (k = nint((log((float) (nz2) / ((imsl_fi_power(5.0e0, i)) *
                                                                    (imsl_fi_power(3.0e0, j))))) / log(2.0e0)); k <= nint((log(((float) (nz2) /
                                                                                                                                (imsl_fi_power(5.0e0, i))) / (imsl_fi_power(3.0e0, j)))) / (log(2.0e0))); k++) {
                                        valnew = (imsl_fi_power(2.0e0, k)) * (imsl_fi_power(3.0e0, j)) * (imsl_fi_power(5.0e0, i));
                                        if ((valnew >= nz2) && (valnew < valold)) {
                                                valold = valnew;
                                                if (valold == (nz2))
                                                        goto L_42;
                                        }
                                }
                        }
                }
        }
L_42:
        if ((ipad == 1) || (ipad == 3))
                nz2 = valold;
        *nz = nz2;
         /*
         * Allocate workspace based on the minimum value fo NZ as computed
         * above and stored in NZ2.  We use NZ2 instead of the input NZ since
         * we will be using the minimum value for NZ when padding the vectors
         * X and Y.
         */
        if ((ido == 0) || (ido == 1)) {
        	icwksp = (Mf_complex *) imsl_malloc (2*nz2 * sizeof(*icwksp));
        	irwksp = (Mfloat *)imsl_malloc ((6*nz2+15) * sizeof(*irwksp));
	}

        if( !zhat_user && return_zhat)
           *lv_zhat = (Mf_complex*) imsl_malloc (nz2 * sizeof (**lv_zhat));

        if (!return_zhat)
            lv_temp_zhat = (Mf_complex*) imsl_malloc (nz2*sizeof(*lv_temp_zhat));

        if ( (!zhat_user && return_zhat && *lv_zhat == NULL) ||
                ( (ido==0 || ido==1) && (irwksp==NULL || icwksp==NULL)) ||
                (!return_zhat && lv_temp_zhat == NULL) ) {
            /* Not enough memory, with %(L1) = %(I1). */
            imsl_e1sti (1, n);
            imsl_e1stl (1, "n");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE;
        }

        if (!user_space) {
            lv_z = (Mf_complex *) imsl_malloc (nz2 * sizeof (*lv_z));
            if (lv_z == NULL) {
                /* Not enough memory, with %(L1) = %(I1). */
                imsl_e1sti (1, n);
                imsl_e1stl (1, "n");
                imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
                goto FREE_SPACE;
            }
        }

        /* Call computational routine */

    if (zhat_user) 
        l_c2orl(&ido, &n, x, y, &ipad, nz, lv_z, lv_zhat_user, 
                   icwksp , icwksp + nz2 ,irwksp);
    else if (return_zhat)
        l_c2orl(&ido, &n, x, y, &ipad, nz, lv_z, *lv_zhat,
                   icwksp , icwksp + nz2 ,irwksp);
    else
        l_c2orl(&ido, &n, x, y, &ipad, nz, lv_z, lv_temp_zhat,
                   icwksp , icwksp + nz2 ,irwksp);
    }

L_9000:
    if (ido == 0 || ido == 3) {
       if (irwksp != NULL) { 
		imsl_free (irwksp);
		irwksp = NULL;
	}
       if (icwksp != NULL) {
		imsl_free (icwksp);
		icwksp = NULL;
	}
    }
    if (lv_temp_zhat != NULL && !return_zhat)
       imsl_free (lv_temp_zhat);

FREE_SPACE:
    ;
RETURN:
    return (argptr);
}
/* Structured by FOR_STRUCT, v0.2, on 06/21/91 at 10:46:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C2ONV/DC2ONV (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    December 1, 1989

    Purpose:    Compute the convolution of two d_complex vectors.

    Usage:      CALL C2ONV (IDO, NX, X, NY, Y, IPAD, NZ, Z,
                            ZHAT, XWK, YWK, WK)

    Arguments:  (See CCONV)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI 
static void l_c2onv(Mint *ido, Mint *nx, Mf_complex x[], Mint *ny, 
                    Mf_complex y[], Mint *ipad, Mint *nz, Mf_complex z[], 
                    Mf_complex zhat[], Mf_complex xwk[], Mf_complex ywk[], 
                    Mfloat wk[])
#else
static void l_c2onv(ido, nx, x, ny, y, ipad, nz, z,
                    zhat, xwk, ywk, wk)
        Mint            *ido, *nx;
        Mf_complex       x[];
        Mint            *ny;
        Mf_complex       y[];
        Mint            *ipad, *nz;
        Mf_complex       z[], zhat[], xwk[], ywk[];
        Mfloat           wk[];
#endif 
{
        Mint             _l0,i,_l1, icpy, iwfftc, j, k, n, nz2, nzold;
        Mfloat           valnew, valold;
        Mf_complex       _cx0;

        imsl_e1psh("C2ONV ");
        /* Check IDO */
        if ((*ido < 0) || (*ido > 3)) {
                imsl_e1sti(1, *ido);

                imsl_e1mes(5, 1, "The value of IDO must be zero, one, two, or three but IDO = %(i1) was given.");
        }
        /* Check IPAD */
        if ((*ipad < 0) || (*ipad > 1)) {
                imsl_e1sti(1, *ipad);

                imsl_e1mes(5, 2, "The value of IPAD must be either zero or one, but IPAD = %(i1) was given.");
        }
        /* Check NX */
        if (*nx <= 0) {
                imsl_e1sti(1, *nx);

                imsl_e1mes(5, 3, "The length of the vector X must be positive while NX = %(i1) is given.");
        }
        /* Check NY */
        if (*ny <= 0) {
                imsl_e1sti(1, *ny);

                imsl_e1mes(5, 4, "The length of the vector Y must be positive while NY = %(i1) is given.");
        }
        if (imsl_n1rty(0) != 0)
                goto RETURN;
        /* Compute correct NZ */
        nzold = *nz;
        if (*ipad == 0) {
                n = imsl_i_max(*nx, *ny);
        } else {
                nz2 = *nx + *ny - 1;
                valold = 2.0 * (nz2);
                for (i = 0; i <= nint(log((float) (nz2)) / log(5.0e0)); i++) {
                        for (j = 0; j <= nint((log((float) (nz2) / (imsl_fi_power(5.0e0, i)))) /
                                              (log(3.0e0))); j++) {
                                for (k = nint((log((float) (nz2) / ((imsl_fi_power(5.0e0, i)) *
                                                                    (imsl_fi_power(3.0e0, j))))) / log(2.0e0)); k <= nint((log(((float) (nz2) /
                                                                                                                                (imsl_fi_power(5.0e0, i))) / (imsl_fi_power(3.0e0, j)))) / (log(2.0e0))); k++) {
valnew = (imsl_fi_power(2.0e0, k)) * (imsl_fi_power(3.0e0, j)) * (imsl_fi_power(5.0e0, i));
                                        if ((valnew >= nz2) && (valnew < valold)) {
                                                valold = valnew;
                                                if (valold == (nz2))
                                                        goto L_43;
                                        }
                                }
                        }
                }
        }
L_43:
        if (*ipad == 1)
                n = valold;
        /* Check input NZ. */
        if (*nz < n) {
                imsl_e1sti(1, n);
                imsl_e1sti(2, *nz);
                imsl_e1sti(3, *nx);
                imsl_e1sti(4, *ny);

                imsl_e1mes(4, 1, "The length of the vector Z must be at least %(i1) while NZ = %(i2) is given, where NZ is based on NX = %(i3) and NY = %(i4).");
        }
        *nz = n;
        if (imsl_n1rty(0) != 0)
                goto RETURN;
        /* Fill and pad work vectors. */
        imsl_ccopy(nx, x, ADR(_l0, 1), xwk, ADR(_l1, 1));
        if (*nz > *nx)
                imsl_cset(ADR(_l0, *nz - *nx), ADR(_cx0, imsl_cf_convert(0.0, 0.0)), &xwk[*nx],
                          ADR(_l1, 1));

        imsl_ccopy(ny, y, ADR(_l0, 1), ywk, ADR(_l1, 1));
        if (*nz > *ny)
                imsl_cset(ADR(_l0, *nz - *ny), ADR(_cx0, imsl_cf_convert(0.0, 0.0)), &ywk[*ny],
                          ADR(_l1, 1));
        /* Compute the forward transforms. */
        iwfftc = 1;
        icpy = iwfftc + 4 ** nz + 15;
        if ((*ido == 0) || (*ido == 1))
                l_fftci(nz, &wk[iwfftc - 1]);
        imsl_f2tcf(nz, xwk, xwk, &wk[iwfftc - 1], &wk[icpy - 1]);
        imsl_f2tcf(nz, ywk, ywk, &wk[iwfftc - 1], &wk[icpy - 1]);
        /*
         * Multiply elements of transformed vectors.
         */
        for (i = 1; i <= *nz; i++) {
                zhat[i - 1] = imsl_c_mul(xwk[i - 1], ywk[i - 1]);
        }
        /* Compute backward transform of ZWK. */
        imsl_f2tcb(nz, zhat, z, &wk[iwfftc - 1], &wk[icpy - 1]);
        /* Scale answer */
        for (i = 1; i <= *nz; i++) {
                z[i - 1] = imsl_c_div(z[i - 1], imsl_cf_convert((double) *nz, 0.));        }
        /* Fill rest of Z with zeros */
        if (nzold > *nz)
                imsl_cset(ADR(_l0, nzold - *nz), ADR(_cx0, imsl_cf_convert(0.0e0, 0.0e0)),&z[*nz], ADR(_l1, 1));
RETURN:
        imsl_e1pop("C2ONV ");
        return;
}                               /* end of function */
/*Translated by FOR_C++, v0.1, on 06/21/91 at 11:13:30 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/21/91 at 11:13:28
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  FFTCI/DFFTCI (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    June 10, 1991

    Purpose:    Compute parameters needed by FFTCF and FFTCB.

    Usage:      CALL FFTCI (N, WFFTC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WFFTC  - Array of length 4N+15 containing parameters needed by
                FFTCF and FFTCB.  (Output)

    Remark:
       Different WFFTC arrays are needed for different values of N.

    Keywords:   Complex exponential FFT; Initialization

    GAMS:       J1a2

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
*/
#ifdef ANSI
static void l_fftci(Mint *n, Mfloat wfftc[])
#else
static void l_fftci(n, wfftc)
        Mint            *n;
        Mfloat           wfftc[];
#endif
{

        /* CHECK ARGUMENT N */
        if (*n < 1) {
                imsl_e1psh("FFTCI ");
                imsl_e1sti(1, *n);

                imsl_e1mes(5, 1, "The length of the sequence N = %(i1).  It must be at least 1.");
                imsl_e1pop("FFTCI ");
                goto RETURN;
        }
        if (*n > 1) {
                l_f3tci(n, &wfftc[*n * 2], &wfftc[*n * 4]);
        }
RETURN:
        ;
        return;
}                               /* end of function */
/*Translated by FOR_C++, v0.1, on 06/21/91 at 11:17:06 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/21/91 at 11:17:04
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F3TCI/DF3TCI (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute parameters needed by FFTCI and FFTCB.

    Usage:      CALL F3TCI (N, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WA     - Vector of length 2*N.  (Output)
       FAC    - Vector of length 14.  (Output)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f3tci(Mint *n, Mfloat wa[], Mfloat imsl_fac[])
#else
static void l_f3tci(n, wa, imsl_fac)
        Mint            *n;
        Mfloat           wa[], imsl_fac[];
#endif
{
        Mint             i, i1, ido, idot, ii, ip, ipm, j, k1, l1, l2, ld,
                        nf, nl, nq, nr, ntry;
        Mfloat           arg, argh, argld, fi, tpi;
        static Mint      ntryh[4] = {3, 4, 2, 5};

        /*
         * FAC WILL BE A VECTOR OF LENGTH 14. FAC(1) WILL CONTAIN N, FAC(2)
         * THRU FAC(14) WILL HAVE THE PRIME FACTORS OF N, IN ASCENDING ORDER.
         * WA WILL BE A VECTOR OF LENGTH 2*N, CONTAINING THE REAL AND
         * IMAGINARY PARTS OF CEXP(2*PI*?/N)
         */
        nl = *n;
        nf = 0;
        j = 0;
L_10:
        j += 1;
        if (j - 4 > 0)
                goto L_30;
        goto L_20;
L_20:
        ntry = ntryh[j - 1];
        goto L_40;
L_30:
        ntry += 2;
L_40:
        nq = nl / ntry;
        nr = nl - ntry * nq;
        if (nr != 0)
                goto L_10;
        goto L_50;
L_50:
        nf += 1;
        imsl_fac[nf + 1] = (float) (ntry);
        nl = nq;
        if (ntry != 2)
                goto L_60;
        if (nf == 1)
                goto L_60;
        imsl_scopy(nf - 1, &imsl_fac[2], -1, &imsl_fac[3], -1);
        imsl_fac[2] = 2.0;
L_60:
        if (nl != 1)
                goto L_40;
        imsl_fac[0] = (float) (*n);
        imsl_fac[1] = (float) (nf);
        tpi = 2.0 * 3.1415926535897932384626433831;
        argh = tpi / (float) (*n);
        i = 2;
        l1 = 1;
        for (k1 = 1; k1 <= nf; k1++) {
                ip = nint(imsl_fac[k1 + 1]);
                ld = 0;
                l2 = l1 * ip;
                ido = *n / l2;
                idot = ido + ido + 2;
                ipm = ip - 1;
                for (j = 1; j <= ipm; j++) {
                        i1 = i;
                        wa[i - 2] = 1.0;
                        wa[i - 1] = 0.0;
                        ld += l1;
                        fi = 0.0;
                        argld = (float) (ld) * argh;
                        for (ii = 4; ii <= idot; ii += 2) {
i += 2;
                                fi += 1.0;
                                arg = fi * argld;
                                wa[i - 2] = cos(arg);
                                wa[i - 1] = sin(arg);
                        }
                        if (ip > 5) {
                                wa[i1 - 2] = wa[i - 2];
                                wa[i1 - 1] = wa[i - 1];
                        }
                }
                l1 = l2;
        }
        return;
}                               /* end of function */


/*Translated by FOR_C++, v0.1, on 06/21/91 at 10:45:44 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/21/91 at 10:45:42
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C2ORL/DC2ORL (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    December 1, 1989

    Purpose:    Compute the correlation of two d_complex vectors.

    Usage:      CALL C2ORL (IDO, N, X, Y, IPAD, NZ, Z, ZHAT, XWK, YWK,
                            WK)

    Arguments:  (See CCORL)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2orl (Mint *ido, Mint *n, Mf_complex x[], Mf_complex y[],
	Mint *ipad, Mint *nz, Mf_complex z[], Mf_complex zhat[], 
	Mf_complex xwk[], Mf_complex ywk[], Mfloat wk[])
#else 
static void l_c2orl(ido, n, x, y, ipad, nz, z, zhat, xwk, ywk,
	   wk)

	Mint            *ido, *n;
	Mf_complex       x[], y[];
	Mint            *ipad, *nz;
	Mf_complex       z[], zhat[], xwk[], ywk[];
	Mfloat           wk[];
#endif 
{
	Mint             _l0, _l1, i, icpy, iwfftc, j, k, nz2, nzold;
	Mfloat           valnew, valold;
	Mf_complex       _cx0;


	imsl_e1psh("C2ORL ");
	/* Check IDO */
	if ((*ido < 0) || (*ido > 3)) {
		imsl_e1sti(1, *ido);

		imsl_e1mes(5, 1, "The value of IDO must be zero, one, two, or three but IDO = %(i1) was given.");
	}
	/* Check IPAD */
	if ((*ipad < 0) || (*ipad > 3)) {
		imsl_e1sti(1, *ipad);

		imsl_e1mes(5, 2, "The value of IPAD must be either zero one, two, or three but IPAD = %(i1) was given.");
	}
	/* Check N */
	if (*n <= 0) {
		imsl_e1sti(1, *n);

		imsl_e1mes(5, 3, "The length of the input vectors must be positive while N = %(i1) is given.");
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Compute correct NZ */
	nzold = *nz;
	nz2 = *n;
	if ((*ipad == 1) || (*ipad == 3)) {
		nz2 = 2 ** n - 1;
		valold = 2.0 * (nz2);
		for (i = 0; i <= nint(log((float) (nz2)) / log(5.0e0)); i++) {
			for (j = 0; j <= nint((log((float) (nz2) / (imsl_fi_power(5.0e0, i)))) /
					      (log(3.0e0))); j++) {
				for (k = nint((log((float) (nz2) / ((imsl_fi_power(5.0e0, i)) *
								    (imsl_fi_power(3.0e0, j))))) / log(2.0e0)); k <= nint((log(((float) (nz2) /
																(imsl_fi_power(5.0e0, i))) / (imsl_fi_power(3.0e0, j)))) / (log(2.0e0))); k++) {
					valnew = (imsl_fi_power(2.0e0, k)) * (imsl_fi_power(3.0e0, j)) * (imsl_fi_power(5.0e0, i));
					if ((valnew >= nz2) && (valnew < valold)) {
						valold = valnew;
						if (valold == (nz2))
							goto L_40;
					}
				}
			}
		}
	}
L_40:
	if ((*ipad == 1) || (*ipad == 3))
		nz2 = valold;
	/* Check input NZ */
	if (*nz < nz2) {
		imsl_e1sti(1, nz2);
		imsl_e1sti(2, *nz);
		imsl_e1sti(3, *n);

		imsl_e1mes(4, 1, "The length of the vector Z must be at least %(i1) while NZ = %(i2) is given, where NZ is based on N = %(i3).");
	}
	*nz = nz2;
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * If ipad is zero or one, we have to transform both vectors.  If
	 * ipad is two or three we can compute the correlation using only X.
	 */
	if ((*ipad == 0) || (*ipad == 1)) {
		/* Fill and pad work vectors. */
		imsl_ccopy(n, x, ADR(_l0, 1), xwk, ADR(_l1, 1));
		imsl_ccopy(n, y, ADR(_l0, 1), ywk, ADR(_l1, 1));
		imsl_cset(ADR(_l0, *nz - *n), ADR(_cx0, imsl_cf_convert(0.0, 0.0)), &xwk[*n],
			  ADR(_l1, 1));
		imsl_cset(ADR(_l0, *nz - *n), ADR(_cx0, imsl_cf_convert(0.0, 0.0)), &ywk[*n],
			  ADR(_l1, 1));
		/* Compute the forward transforms. */
		iwfftc = 1;
		icpy = iwfftc + 4 ** nz + 15;
		if ((*ido == 0) || (*ido == 1))
			l_fftci(nz, &wk[iwfftc - 1]);
		imsl_f2tcf(nz, xwk, xwk, &wk[iwfftc - 1], &wk[icpy - 1]);
		imsl_f2tcf(nz, ywk, ywk, &wk[iwfftc - 1], &wk[icpy - 1]);
		/*
		 * Multiply the elements of the transformed vectors.
		 */
		for (i = 1; i <= *nz; i++) {
			zhat[i - 1] = imsl_c_mul(xwk[i - 1], imsl_c_conjg(ywk[i - 1]));
		}

		/*
		 * This else-clause is executed if IPAD is two or three.
		 */
	} else {
		/* Fill and pad work vector. */
		imsl_ccopy(n, x, ADR(_l0, 1), xwk, ADR(_l1, 1));
		imsl_cset(ADR(_l0, *nz - *n), ADR(_cx0, imsl_cf_convert(0.0, 0.0)), &xwk[*n],
			  ADR(_l1, 1));
		/* Compute the forward transform. */
		iwfftc = 1;
		icpy = iwfftc + 4 ** nz + 15;
		if ((*ido == 0) || (*ido == 1))
			l_fftci(nz, &wk[iwfftc - 1]);
		imsl_f2tcf(nz, xwk, xwk, &wk[iwfftc - 1], &wk[icpy - 1]);
		/*
		 * Multiply the elements of the transformed vector.
		 */
		for (i = 1; i <= *nz; i++) {
			zhat[i - 1] = imsl_c_mul(xwk[i - 1], imsl_c_conjg(xwk[i - 1]));
		}
	}
	/*
	 * Compute backward transform of ZHAT and store it in Z.
	 */
	imsl_f2tcb(nz, zhat, z, &wk[iwfftc - 1], &wk[icpy - 1]);

	/* Scale answer. */
	for (i = 1; i <= *nz; i++) {
		z[i - 1] = imsl_c_div(z[i - 1], imsl_cf_convert((float) (*nz), 0.));
	}
	/* Fill rest of Z with zeros */
	if (nzold > *nz)
		imsl_cset(ADR(_l0, nzold - *nz), ADR(_cx0, imsl_cf_convert(0.0e0, 0.0e0)),
			  &z[*nz], ADR(_l1, 1));

L_9000:
	imsl_e1pop("C2ORL ");
	return;
}				/* end of function */
