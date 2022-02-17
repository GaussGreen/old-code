#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_convolution (Mint, Mfloat *, Mint, Mfloat *,
		Mint *, va_list argptr);
static void l_r2onv (Mint *, Mint *, Mfloat [], Mint *, Mfloat [],
	Mint *, Mint *, Mfloat [], Mfloat [], Mfloat [], Mfloat [],
	Mfloat []);
static void l_r2orl (Mint *, Mint *, Mfloat [], Mfloat [],
		Mint *, Mint *, Mfloat [], Mfloat [],
		Mfloat [], Mfloat [], Mfloat []);
static void l_fftri (Mint *, Mfloat []);
static void l_f3tri (Mint *, Mfloat [], Mfloat []);
#else
static VA_LIST_HACK l_convolution ();
static void l_r2onv ();
static void l_r2orl ();
static void l_fftri ();
static void l_f3tri ();
#endif

static Mfloat *lv_z = NULL;
static Mfloat **lv_zhat = NULL;
static Mfloat *lv_zhat_user = NULL;

#ifdef ANSI
Mfloat     *imsl_f_convolution (Mint nx, Mfloat *x, Mint ny,
		Mfloat *y, Mint *nz, ...)
#else
Mfloat     *imsl_f_convolution (nx, x, ny, y, nz, va_alist)
    Mint        nx;
    Mfloat     *x;
    Mint        ny;
    Mfloat     *y;
    Mint       *nz;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, nz);
    E1PSH ("imsl_f_convolution", "imsl_d_convolution");
    lv_z = NULL;
    lv_zhat = NULL;
    lv_zhat_user = NULL;
    IMSL_CALL (l_convolution (nx, x, ny, y, nz, argptr));
    va_end (argptr);
    E1POP ("imsl_f_convolution", "imsl_d_convolution");
    return lv_z;
}


#ifdef ANSI
static VA_LIST_HACK l_convolution (Mint nx, Mfloat *x, Mint ny, Mfloat *y,
		Mint *nz, va_list argptr)
#else
static VA_LIST_HACK l_convolution (nx, x, ny, y, nz, argptr)
    Mint        nx;
    Mfloat     *x;
    Mint        ny;
    Mfloat     *y;
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
    Mint        n;
    Mint        nz2;
    Mint        valold;
    Mint        i;
    Mint        j;
    Mint        k;
    Mint        valnew;
    static Mfloat     *irwksp;
    Mfloat *lv_temp_zhat = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_z = va_arg (argptr, Mfloat *);
	    user_space = 1;
	    arg_number++;
	    break;
	case IMSL_PERIODIC:
	    periodic = 1;
	    break;
	case IMSL_CORRELATION:
	    correlation = 1;
	    arg_number++;
	    break;
	case IMSL_Z_TRANS:
	    lv_zhat = va_arg (argptr, Mfloat **);
	    return_zhat = 1;
	    arg_number++;
	    break;
	case IMSL_Z_TRANS_USER:
	    lv_zhat_user = va_arg (argptr, Mfloat *);
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

    if (zhat_user && (lv_zhat_user == NULL)) {
	imsl_e1stl (1, "zhat");
	imsl_e1stl (2, "IMSL_Z_TRANS_USER");
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

    if (&x[0] != &y[0]) {
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
	imsl_e1mes (5, 3, "The length of the vector X must be positive 
			while NX = %(i1) is given.");
*/
	imsl_ermes(IMSL_TERMINAL, IMSL_X_NOT_POSITIVE);
    }
    /* Check NY */
    if (ny <= 0) {
	imsl_e1sti (1, ny);
/*
	imsl_e1mes (5, 4, "The length of the vector Y must be positive
			 while NY = %(i1) is given.");
*/
	imsl_ermes(IMSL_TERMINAL, IMSL_Y_NOT_POSITIVE);
    }

    if (imsl_n1rty (0) != 0)
	goto L_9000;

if (!correlation) {


/*  rconv follows here.  This computation is needed to determine
    how much space to allocate before call r2onv */


    if (ipad == 2) ipad = 0;
    if (ipad == 3) ipad = 1;

    /* Compute correct NZ */
    if (ipad == 0) {
	n = imsl_i_max (nx, ny);
    }
    else {
	nz2 = nx + ny - 1;
	valold = 2.0 * (nz2);
	for (i = 0; i <= nint (log ((Mfloat) (nz2)) / log (5.0e0)); i++) {
	    for (j = 0; j <= nint ((log ((Mfloat) (nz2) / (imsl_fi_power (5.0e0, i)))) /
		    (log (3.0e0))); j++) {
		for (k = nint ((log ((Mfloat) (nz2) / ((imsl_fi_power (5.0e0, i)) *
				    (imsl_fi_power (3.0e0, j))))) / log (2.0e0)); k <= nint ((log (((Mfloat) (nz2) /
				    (imsl_fi_power (5.0e0, i))) / (imsl_fi_power (3.0e0, j)))) / (log (2.0e0))); k++) {
		    valnew = (imsl_fi_power (2.0e0, k)) * (imsl_fi_power (3.0e0, j)) * (imsl_fi_power (5.0e0, i));
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
    /*
     * Allocate workspace based on the minimum value for NZ as computed above
     * and stored in N.  We use N instead of the input NZ since we will be
     * using the minimum value for NZ when padding the vectors X and Y.
     */
	if (ido == 0 || ido == 1)
	    irwksp = (Mfloat *) imsl_malloc ((4 * n + 15) * sizeof (*irwksp));


	if (!zhat_user && return_zhat)
	    *lv_zhat = (Mfloat *) imsl_malloc (n * sizeof (**lv_zhat));

	if (!return_zhat)
	    lv_temp_zhat = (Mfloat *) imsl_malloc (n * sizeof(*lv_temp_zhat));

	if ( (!zhat_user && return_zhat && *lv_zhat == NULL) || 
		( (ido == 0 || ido == 1) && irwksp == NULL) || 
		(!return_zhat && lv_temp_zhat == NULL) ) {
	    /* Not enough memory, with %(L1) = %(I1). */
	    imsl_e1sti (1, n);
	    imsl_e1stl (1, "n");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}

	if (!user_space) {
	    lv_z = (Mfloat *) imsl_malloc (n * sizeof (*lv_z));
	    if (lv_z == NULL) {
		/* Not enough memory, with %(L1) = %(I1). */
		imsl_e1sti (1, n);
		imsl_e1stl (1, "n");
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		goto FREE_SPACE;
	    }
	}

    /* CALL COMPUTATIONAL ROUTINE */

    if (zhat_user) 
	l_r2onv (&ido, &nx, x, &ny, y, &ipad, nz, lv_z, lv_zhat_user,
				irwksp, (irwksp+n), (irwksp+n*2));
    else if (return_zhat)
        l_r2onv (&ido, &nx, x, &ny, y, &ipad, nz, lv_z, *lv_zhat, irwksp,
				(irwksp + n), (irwksp + n * 2));
    else 
        l_r2onv (&ido, &nx, x, &ny, y, &ipad, nz, lv_z, lv_temp_zhat, irwksp,
				(irwksp + n), (irwksp + n * 2));

    } /* end-if convolution */

    else /* correlation */ {



/*  rcorl follows here.  This computation is needed to determine
    how much space to allocate before call r2orl */


    /* Compute correct NZ */


    n = nx;
    nz2 = n;
    if ((ipad == 1) || (ipad == 3)) {
        nz2 = 2 * n - 1;
        valold = 2.0 * (nz2);
        for (i = 0; i <= nint (log ((Mfloat) (nz2)) / log (5.0e0)); i++) {
            for (j = 0; j <= nint ((log ((Mfloat) (nz2) / (imsl_fi_power (5.0e0, i)))) /
                    (log (3.0e0))); j++) {
                for (k = nint ((log ((Mfloat) (nz2) / ((imsl_fi_power (5.0e0, i)) *
                                    (imsl_fi_power (3.0e0, j))))) / log (2.0e0)); k <= nint ((log (((Mfloat) (nz2) /
                                    (imsl_fi_power (5.0e0, i))) / (imsl_fi_power (3.0e0, j)))) / (log (2.0e0))); k++) {
                    valnew = (imsl_fi_power (2.0e0, k)) * (imsl_fi_power (3.0e0, j)) * (imsl_fi_power (5.0e0, i));
                    if ((valnew >= nz2) && (valnew < valold)) {
                        valold = valnew;
                        if (valold == (nz2))
                            goto L_41;
                    }
                }
            }
        }
    }
L_41:
    if ((ipad == 1) || (ipad == 3))
        nz2 = valold;
    *nz = nz2;
    /*
     * Allocate workspace based on the minimum value for NZ as computed above
     * and stored in NZ2.  We use NZ2 instead of the input NZ since we will
     * be using the minimum value for NZ when padding the vectors X and Y.
     */

	if (ido == 0 || ido == 1)
	irwksp = (Mfloat *) imsl_malloc((4*nz2+15)*sizeof(*irwksp));

	if (!zhat_user && return_zhat)
            *lv_zhat = (Mfloat *) imsl_malloc (nz2 * sizeof (**lv_zhat));

        if (!return_zhat)
            lv_temp_zhat = (Mfloat *) imsl_malloc (nz2*sizeof(*lv_temp_zhat));

        if ( (!zhat_user && return_zhat && *lv_zhat == NULL) || 
		( (ido == 0 || ido == 1) && irwksp == NULL) ||
		(!return_zhat && lv_temp_zhat == NULL) ) {
            /* Not enough memory, with %(L1) = %(I1). */
            imsl_e1sti (1, n);
            imsl_e1stl (1, "n");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE;
        }
        if (!user_space) {
            lv_z = (Mfloat *) imsl_malloc (nz2 * sizeof (*lv_z));
            if (lv_z == NULL) {
                /* Not enough memory, with %(L1) = %(I1). */
                imsl_e1sti (1, n);
                imsl_e1stl (1, "n");
                imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
                goto FREE_SPACE;
            }
        }


    /* CALL COMPUTATIONAL ROUTINE */

    if (zhat_user) 
	l_r2orl (&ido, &n, x, y, &ipad, nz, lv_z, lv_zhat_user,
				irwksp, (irwksp+nz2), (irwksp+nz2*2));
    else if (return_zhat)
        l_r2orl (&ido, &n, x, y, &ipad, nz, lv_z, *lv_zhat, irwksp, 
				(irwksp + nz2), (irwksp + nz2 * 2) );
    else
        l_r2orl (&ido, &n, x, y, &ipad, nz, lv_z, lv_temp_zhat, irwksp, 
				(irwksp + nz2), (irwksp + nz2 * 2) );

    }


L_9000:

if ((ido == 0 || ido == 3) && irwksp != NULL) {
	imsl_free (irwksp);
	irwksp = NULL;
}
if (lv_temp_zhat != NULL)
	imsl_free(lv_temp_zhat);
/*
	if (!return_zhat && lv_zhat != NULL)
	    imsl_free(*lv_zhat);
*/
FREE_SPACE:
    ;
RETURN:
    return (argptr);
}





/*----------------------------------------------------------------------- */

/*  IMSL Name:  R2ONV/DR2ONV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 1, 1989

    Purpose:    Compute the convolution of two real vectors.

    Usage:      CALL R2ONV (IDO, NX, X, NY, Y, IPAD, NZ, Z, ZHAT,
                            XWK, YWK, WK)

    Arguments:  (See RCONV)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2onv (Mint *ido, Mint *nx, Mfloat x[], Mint *ny,
	Mfloat y[], Mint *ipad, Mint *nz, Mfloat z[], Mfloat zhat[],
	Mfloat xwk[], Mfloat ywk[], Mfloat wk[])
#else
static void l_r2onv (ido, nx, x, ny, y, ipad, nz, z, zhat, xwk,
                ywk, wk)
    Mint        *ido, *nx;
    Mfloat       x[];
    Mint        *ny;
    Mfloat       y[];
    Mint        *ipad, *nz;
    Mfloat       z[], zhat[], xwk[], ywk[], wk[];
#endif
{
    Mint         i, iwfftc, j, k, n, nz2, nzold;
    Mfloat       valnew, valold;


    imsl_e1psh ("R2ONV ");
    /* COMPUTE CORRECT NZ */
    nzold = *nz;
    if (*ipad == 0) {
	n = imsl_i_max (*nx, *ny);
    }
    else {
	nz2 = *nx + *ny - 1;
	valold = 2.0 * (nz2);
	for (i = 0; i <= nint (log ((Mfloat) (nz2)) / log (5.0e0)); i++) {
	    for (j = 0; j <= nint ((log ((Mfloat) (nz2) / (imsl_fi_power (5.0e0, i)))) /
		    (log (3.0e0))); j++) {
		for (k = nint ((log ((Mfloat) (nz2) / ((imsl_fi_power (5.0e0, i)) *
				    (imsl_fi_power (3.0e0, j))))) / log (2.0e0)); k <= nint ((log (((Mfloat) (nz2) /
				    (imsl_fi_power (5.0e0, i))) / (imsl_fi_power (3.0e0, j)))) / (log (2.0e0))); k++) {
		    valnew = (imsl_fi_power (2.0e0, k)) * (imsl_fi_power (3.0e0, j)) * (imsl_fi_power (5.0e0, i));
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
    if (*ipad == 1)
	n = valold;
    *nz = n;
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /* FILL AND PAD WORK VECTORS. */
    scopy (*nx, x, 1, xwk, 1);
    if (*nz > *nx)
	sset (*nz - *nx, 0.0e0, &xwk[*nx], 1);

    scopy (*ny, y, 1, ywk, 1);
    if (*nz > *ny)
	sset (*nz - *ny, 0.0e0, &ywk[*ny], 1);

    /*
     * REAL FORWARD TRANSFORM VECTORS INTO XWK & YWK.
     */
    iwfftc = 1;
    if ((*ido == 0) || (*ido == 1))
	l_fftri (nz, &wk[iwfftc - 1]);
    imsl_f2trf (nz, xwk, xwk, &wk[iwfftc - 1]);
    imsl_f2trf (nz, ywk, ywk, &wk[iwfftc - 1]);
    /*
     * COMPUTE PRE-REVERSE-TRANSFORM REAL VECTOR.  THIS LOOP MIMICKS THE THE
     * PROCESS OF TAKING THE COMPLEX FORWARD TRANSFORM OF THE ORIGINAL
     * VECTORS, MULTIPLYING THE ELEMENTS OF THE TRANSFORMED VECTORS, AND
     * CONVERTING THAT BACK TO THE REAL VECTOR TO BE REVERSE-TRANSFORMED.
     */
    zhat[0] = xwk[0] * ywk[0];
    for (i = 1; i <= ((*nz - 1) / 2); i++) {
	zhat[i * 2 - 1] = xwk[i * 2 - 1] * ywk[i * 2 - 1] - xwk[i * 2] * ywk[i * 2];
	zhat[i * 2] = xwk[i * 2] * ywk[i * 2 - 1] + xwk[i * 2 - 1] * ywk[i * 2];
    }
    if (mod (*nz, 2) == 0)
	zhat[*nz - 1] = xwk[*nz - 1] * ywk[*nz - 1];
    /* REVERSE TRANSFORM REAL VECTOR. */
    imsl_f2trb (nz, zhat, z, &wk[iwfftc - 1]);
    /* SCALE ANSWER */
    for (i = 1; i <= *nz; i++) {
	z[i - 1] /= *nz;
    }
    /* FILL REST OF Z WITH ZEROS */
    if (nzold > *nz)
	sset (nzold - *nz, 0.0e0, &z[*nz], 1);
L_9000:
    imsl_e1pop ("R2ONV ");
    return;
}				/* end of function */

/*Translated by FOR_C++, v0.1, on 06/17/91 at 13:47:46 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/17/91 at 13:47:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  FFTRI/DFFTRI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 10, 1991

    Purpose:    Compute parameters needed by FFTRF and FFTRB.

    Usage:      CALL FFTRI (N, WFFTR)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WFFTR  - Array of length 2N+15 containing parameters needed by
                FFTRF and FFTRB.  (Output)

    Remark:
       Different WFFTR arrays are needed for different values of N.

    Keywords:   Initialization; Transforms; Trigonometric

    GAMS:       J1a1

    Chapters:   MATH/LIBRARY Transforms
                STAT/LIBRARY Mathematical Support

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_fftri (Mint *n, Mfloat wfftr[])
#else
static void l_fftri (n, wfftr)
    Mint        *n;
    Mfloat       wfftr[];
#endif
{

    /* CHECK ARGUMENT N */
    if (*n < 1) {
	imsl_e1psh ("FFTRI ");
	imsl_e1sti (1, *n);
/*
	imsl_e1mes (5, 1, "The length of the sequence N = %(i1).  It must be at least 1.");
*/
	imsl_ermes(IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	imsl_e1pop ("FFTRI ");
	goto L_9000;
    }

    if (*n > 1) {
	l_f3tri (n, &wfftr[*n], &wfftr[*n * 2]);
    }

L_9000:
    ;
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F3TRI/DF3TRI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute parameters needed by RFFTF and RFFTB.

    Usage:      CALL F3TRI (N, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WA     - Real vector.  (Output)
       FAC    - Real vector.  (Output)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f3tri (Mint *n, Mfloat wa[], Mfloat imsl_fac[])
#else
static void l_f3tri (n, wa, imsl_fac)
    Mint        *n;
    Mfloat       wa[], imsl_fac[];
#endif
{
    Mint         i, ido, ii, ip, ipm, is, j, k1, l1, l2, ld, nf, nfm1, nl, nq,
                nr, ntry;
    Mfloat       arg, argh, argld, fi, tpi;
    static Mint  ntryh[4] = {4, 2, 3, 5};



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
    imsl_fac[nf + 1] = (Mfloat) (ntry);
    nl = nq;
    if (ntry != 2)
	goto L_60;
    if (nf == 1)
	goto L_60;
    imsl_scopy (nf - 1, &imsl_fac[2], -1, &imsl_fac[3], -1);
    imsl_fac[2] = 2.0;
L_60:
    if (nl != 1)
	goto L_40;
    imsl_fac[0] = (Mfloat) (*n);
    imsl_fac[1] = (Mfloat) (nf);
    tpi = 2.0 * 3.1415926535897932384626433831;
    argh = tpi / (Mfloat) (*n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0)
	return;
    for (k1 = 1; k1 <= nfm1; k1++) {
	ip = nint (imsl_fac[k1 + 1]);
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	for (j = 1; j <= ipm; j++) {
	    ld += l1;
	    i = is;
	    argld = (Mfloat) (ld) * argh;
	    fi = 0.0;
	    for (ii = 3; ii <= ido; ii += 2) {
		i += 2;
		fi += 1.0;
		arg = fi * argld;
		wa[i - 2] = cos (arg);
		wa[i - 1] = sin (arg);
	    }
	    is += ido;
	}
	l1 = l2;
    }
    return;
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 06/19/91 at 11:26:25 */

/* Structured by FOR_STRUCT, v0.2, on 06/19/91 at 11:26:21
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  R2ORL/DR2ORL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 11, 1990

    Purpose:    Compute the correlation of two real vectors.

    Usage:      CALL R2ORL (IDO, N, X, Y, IPAD, NZ, Z, ZHAT,
                            XWK, YWK, WK)

    Arguments:  (See RCORL)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2orl (Mint *ido, Mint *n, Mfloat x[], Mfloat y[],
		Mint *ipad, Mint *nz, Mfloat z[], Mfloat zhat[],
		Mfloat xwk[], Mfloat ywk[], Mfloat wk[])
#else
static void l_r2orl (ido, n, x, y, ipad, nz, z, zhat, xwk, ywk,
                wk)
    Mint        *ido, *n;
    Mfloat       x[], y[];
    Mint        *ipad, *nz;
    Mfloat       z[], zhat[], xwk[], ywk[], wk[];
#endif
{
    Mint         i, iwfftc, j, k, nz2, nzold;
    Mfloat       valnew, valold;


    imsl_e1psh ("R2ORL ");

    if (imsl_n1rty (0) != 0)
        goto L_9000;
    /* Compute correct NZ */
    nzold = *nz;
    nz2 = *n;
    if ((*ipad == 1) || (*ipad == 3)) {
        nz2 = 2 ** n - 1;
        valold = 2.0 * (nz2);
        for (i = 0; i <= nint (log ((Mfloat) (nz2)) / log (5.0e0)); i++) {
            for (j = 0; j <= nint ((log ((Mfloat) (nz2) / (imsl_fi_power (5.0e0, i)))) /
                    (log (3.0e0))); j++) {
                for (k = nint ((log ((Mfloat) (nz2) / ((imsl_fi_power (5.0e0, i)) *
                                    (imsl_fi_power (3.0e0, j))))) / log (2.0e0)); k <= nint ((log (((Mfloat) (nz2) /
                                    (imsl_fi_power (5.0e0, i))) / (imsl_fi_power (3.0e0, j)))) / (log (2.0e0))); k++) {
                    valnew = (imsl_fi_power (2.0e0, k)) * (imsl_fi_power (3.0e0, j)) * (imsl_fi_power (5.0e0, i));
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
    *nz = nz2;
    if (imsl_n1rty (0) != 0)
        goto L_9000;
    /*
     * If IPAD is zero or one, we have to transform both vectors.  If IPAD is
     * two or three we can compute the correlation using only X.
     */
    if ((*ipad == 0) || (*ipad == 1)) {
        /* Fill and pad work vectors. */
        scopy (*n, x, 1, xwk, 1);
        scopy (*n, y, 1, ywk, 1);
        sset (*nz - *n, 0.0e0, &xwk[*n], 1);
        sset (*nz - *n, 0.0e0, &ywk[*n], 1);

        /*
         * real forward transform vectors into XWK & YWK.
         */
        iwfftc = 1;
        if ((*ido == 0) || (*ido == 1))
            l_fftri (nz, &wk[iwfftc - 1]);
        imsl_f2trf (nz, xwk, xwk, &wk[iwfftc - 1]);
        imsl_f2trf (nz, ywk, ywk, &wk[iwfftc - 1]);
        /*
         * Compute pre-reverse-transform real vector.  This loop mimicks the
         * the process of taking the d_complex forward transform of the
         * original vectors, multiplying the ith element of the X-vector
         * times the conjugate of the ith element of the Y-vector, and then
         * converting the d_complex vector holding these products back to a
         * real vector to be backwards transformed by the real backwards FFT.
         */
        zhat[0] = xwk[0] * ywk[0];
        for (i = 1; i <= ((*nz - 1) / 2); i++) {
            zhat[i * 2 - 1] = xwk[i * 2 - 1] * ywk[i * 2 - 1] + xwk[i * 2] * ywk[i * 2];
            zhat[i * 2] = xwk[i * 2] * ywk[i * 2 - 1] - xwk[i * 2 - 1] * ywk[i * 2];
        }
        if (mod (*nz, 2) == 0)
            zhat[*nz - 1] = xwk[*nz - 1] * ywk[*nz - 1];

    }
    else {
        /*
         * This else-clause is executed if IPAD is two or three. Fill and pad
         * XWK.
         */
        scopy (*n, x, 1, xwk, 1);
        sset (*nz - *n, 0.0e0, &xwk[*n], 1);

        /*
         * Real forward transform vectors into XWK & YWK.
         */
        iwfftc = 1;
        if ((*ido == 0) || (*ido == 1))
            l_fftri (nz, &wk[iwfftc - 1]);
        imsl_f2trf (nz, xwk, xwk, &wk[iwfftc - 1]);
        /*
         * Compute pre-reverse-transform real vector.  This loop mimicks the
         * the process of taking the d_complex forward transform of the
         * original vectors, multiplying the ith element of the X-vector
         * times the conjugate of the ith element of the Y-vector, and then
         * converting the d_complex vector holding these products back to a
         * real vector to be backwards transformed by the real backwards FFT.
         */
        zhat[0] = xwk[0] * xwk[0];
        for (i = 1; i <= ((*nz - 1) / 2); i++) {
            zhat[i * 2 - 1] = xwk[i * 2 - 1] * xwk[i * 2 - 1] + xwk[i * 2] * xwk[i * 2];
            zhat[i * 2] = xwk[i * 2] * xwk[i * 2 - 1] - xwk[i * 2 - 1] * xwk[i * 2];
        }
        if (mod (*nz, 2) == 0)
            zhat[*nz - 1] = xwk[*nz - 1] * xwk[*nz - 1];
    }
    /* Reverse transform real vector. */
    imsl_f2trb (nz, zhat, z, &wk[iwfftc - 1]);
    /* Scale answer */
    for (i = 1; i <= *nz; i++) {
        z[i - 1] /= (Mfloat) (*nz);
    }
    /* fill rest of Z with zeros. */
    sset (nzold - *nz, 0.0e0, &z[*nz], 1);
L_9000:
    imsl_e1pop ("R2ORL ");
    return;
}                               /* end of function */
