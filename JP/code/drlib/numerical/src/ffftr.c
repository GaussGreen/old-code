#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_fft_real (Mint n, Mfloat *seq, va_list argptr);
void imsl_f2trf (Mint *n, Mfloat seq[], Mfloat coef[], Mfloat wfftr[]);
void imsl_f2trb (Mint *n, Mfloat *coef, Mfloat *seq, Mfloat *wfftr);
static void l_f3trf (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[],
                Mfloat fac[]);
    static void l_f3trb (Mint *n, Mfloat *c, Mfloat *ch, Mfloat *wa,
                Mfloat *fac);
    static void l_f4trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1);
    static void l_f5trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2);
    static void l_f5trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2);
    static void l_f6trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3);
    static void l_f7trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3, Mfloat *wa4);
    static void l_f6trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3);
    static void l_f7trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3, Mfloat *wa4);
    static void l_f7trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3, Mfloat *wa4);
    static void l_f8trf (Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat *wa);
    static void l_f8trb (Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat *wa);
#else
static VA_LIST_HACK l_fft_real ();
void imsl_f2trf ();
void imsl_f2trb ();
static void l_f3trf ();
static void l_f3trb ();
static void l_f4trb ();
static void l_f5trf ();
static void l_f5trb ();
static void l_f6trf ();
static void l_f6trb ();
static void l_f7trf ();
static void l_f7trb ();
static void l_f8trf ();
static void l_f8trb ();
#endif

static Mfloat *lv_coef;
#ifdef ANSI
Mfloat     *imsl_f_fft_real (Mint n, Mfloat *seq,...)
#else
Mfloat     *imsl_f_fft_real (n, seq, va_alist)
    Mint        n;
    Mfloat     *seq;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, seq);
    E1PSH ("imsl_f_fft_real", "imsl_d_fft_real");
    lv_coef = NULL;
    IMSL_CALL (l_fft_real (n, seq, argptr));
    va_end (argptr);
    E1POP ("imsl_f_fft_real", "imsl_d_fft_real");
    return lv_coef;
}


#ifdef ANSI
static VA_LIST_HACK l_fft_real (Mint n, Mfloat *seq, va_list argptr)
#else
static VA_LIST_HACK l_fft_real (n, seq, argptr)
    Mint        n;
    Mfloat     *seq;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 2;
    Mint        user_solve = 0;
    Mint        forward = 1;
    Mint        backward = 0;
    Mfloat      supplied_params = 0;
    Mfloat     *params = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_coef = va_arg (argptr, Mfloat *);
	    user_solve = 1;
	    arg_number++;
	    break;
	case IMSL_BACKWARD:
	    forward = 0;
	    backward = 1;
	    break;
	case IMSL_PARAMS:
	    supplied_params = 1;
	    params = va_arg (argptr, Mfloat *);
	    arg_number++;
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

    if (user_solve && (lv_coef == NULL)) {
	imsl_e1stl (1, "q");
	imsl_e1stl (2, "IMSL_RETURN_USER");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }

    if (supplied_params && (params == NULL)) {
	imsl_e1stl (1, "params");
	imsl_e1stl (2, "IMSL_PARAMS");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }

    if (seq == NULL) {
	imsl_e1stl (1, "seq");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
	goto RETURN;
    }

    /* CHECK ARGUMENT N */
    if (n < 1) {
	/* The length of the sequence n = %(I1). */
	/* It must be at least one. */
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	goto RETURN;
    }

    if (!user_solve) {
	lv_coef = (Mfloat *) imsl_malloc (n * sizeof (*lv_coef));
	if (lv_coef == NULL) {
	    /* Not enough memory, with %(L1) = %(I1). */
	    imsl_e1sti (1, n);
	    imsl_e1stl (1, "n");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }

    if (!supplied_params)
	params = imsl_f_fft_real_init (n);
    if (imsl_n1rty (0))
	goto RETURN;

    if (forward) {
	imsl_f2trf (&n, seq, lv_coef, params);
    }

    if (backward) {
	imsl_f2trb (&n, seq, lv_coef, params);
    }
FREE_SPACE:
    if (!supplied_params && (params != NULL)) {
	imsl_free (params);
    }
RETURN:
    if (imsl_n1rty (0) > 3) {
	if (lv_coef != NULL)
	    imsl_free (lv_coef);
	lv_coef = NULL;
    }
    return (argptr);
}


/*Translated by FOR_C++, v0.1, on 06/11/90 at 09:55:48 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 09:55:45
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2TRF/DF2TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence.

    Usage:      CALL F2TRF (N, SEQ, COEF, WFFTR)

    Arguments:  (See FFTRF)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f2trf (Mint *n, Mfloat seq[], Mfloat coef[], Mfloat wfftr[])
#else
void imsl_f2trf (n, seq, coef, wfftr)
    Mint       *n;
    Mfloat      seq[], coef[], wfftr[];
#endif
{


    /* COPY SEQ TO COEF */
    scopy (*n, seq, 1, coef, 1);

    if (*n > 1) {
	l_f3trf (n, coef, wfftr, &wfftr[*n], &wfftr[*n * 2]);
    }
L_9000:
    ;
    return;
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 06/11/90 at 09:57:08 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 09:57:05
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F3TRF/DF3TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 22, 1990

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence.

    Usage:      CALL F3TRF (N, C, CH, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       C      - Real array.  (Input/Output)
       CH     - Real array.  (Workspace)
       WA     - Real array from FFTRI containing the initialized values.
                (Input)
       FAC    - Real array from FFTRI containing the length of the
                sequence, N, in the first location, the number of
                factors, NF, in the second location and the factors in
                the remaining NF locations.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f3trf (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[], Mfloat fac[])
#else
static void l_f3trf (n, c, ch, wa, fac)
    Mint       *n;
    Mfloat      c[], ch[], wa[], fac[];
#endif
{
    Mint        i, ic, idl1, ido, idp2, ip, iw, ix2, ix3, ix4, k1, kh,
                l1, l2, na, nf;
    Mfloat      ti2, tr2;


    if (*n == 1)
	goto L_40;
    nf = nint (fac[1]);
    na = 1;
    l2 = *n;
    iw = *n;
    for (k1 = 1; k1 <= nf; k1++) {
	kh = nf - k1;
	ip = nint (fac[kh + 2]);
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw += -(ip - 1) * ido;
	na = 1 - na;
	if (ip == 4) {
	    ix2 = iw + ido;
	    ix3 = ix2 + ido;
	    if (na != 0) {
		l_f6trf (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	    else {
		l_f6trf (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	}
	else if (ip == 2) {
	    if (na != 0) {
		c[0] = ch[0] + ch[ido];
		c[ido * 2 - 1] = ch[0] - ch[ido];

		if (ido != 1) {
		    idp2 = ido + 2;
		    for (i = 3; i <= ido; i += 2) {
			ic = idp2 - i;
			tr2 = wa[iw + i - 4] * ch[ido + i - 2] + wa[iw + i - 3] *
			    ch[ido + i - 1];
			ti2 = wa[iw + i - 4] * ch[ido + i - 1] - wa[iw + i - 3] *
			    ch[ido + i - 2];
			c[i - 1] = ch[i - 1] + ti2;
			c[ic + ido - 1] = ti2 - ch[i - 1];
			c[i - 2] = ch[i - 2] + tr2;
			c[ic + ido - 2] = ch[i - 2] - tr2;
		    }

		    if (mod (ido, 2) != 1) {
			c[ido] = -ch[ido + ido - 1];
			c[ido - 1] = ch[ido - 1];
		    }
		}
	    }
	    else {
		ch[0] = c[0] + c[ido];
		ch[ido * 2 - 1] = c[0] - c[ido];

		if (ido != 1) {
		    idp2 = ido + 2;
		    for (i = 3; i <= ido; i += 2) {
			ic = idp2 - i;
			tr2 = wa[iw + i - 4] * c[ido + i - 2] + wa[iw + i - 3] *
			    c[ido + i - 1];
			ti2 = wa[iw + i - 4] * c[ido + i - 1] - wa[iw + i - 3] *
			    c[ido + i - 2];
			ch[i - 1] = c[i - 1] + ti2;
			ch[ic + ido - 1] = ti2 - c[i - 1];
			ch[i - 2] = c[i - 2] + tr2;
			ch[ic + ido - 2] = c[i - 2] - tr2;
		    }

		    if (mod (ido, 2) != 1) {
			ch[ido] = -c[ido + ido - 1];
			ch[ido - 1] = c[ido - 1];
		    }
		}
	    }
	}
	else if (ip == 3) {
	    ix2 = iw + ido;
	    if (na != 0) {
		l_f5trf (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	    else {
		l_f5trf (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	}
	else if (ip == 5) {
	    ix2 = iw + ido;
	    ix3 = ix2 + ido;
	    ix4 = ix3 + ido;
	    if (na != 0) {
		l_f7trf (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	    else {
		l_f7trf (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	}
	else {
	    if (ido == 1)
		na = 1 - na;
	    if (na != 0) {
		l_f8trf (&ido, &ip, &l1, &idl1, ch, ch, ch, c, c, &wa[iw - 1]);
		na = 0;
	    }
	    else {
		l_f8trf (&ido, &ip, &l1, &idl1, c, c, c, ch, ch, &wa[iw - 1]);
		na = 1;
	    }
	}
	l2 = l1;
    }
    if (na != 1)
	scopy (*n, ch, 1, c, 1);
L_40:
    return;
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 06/11/90 at 10:00:27 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 10:00:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F5TRF/DF5TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence when the length is divisible by three.

    Usage:      CALL F5TRF (IDO, L1, CC, CH, WA1, WA2)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f5trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat *wa1,
                Mfloat *wa2)
#else
static void l_f5trf (ido, l1, cc, ch, wa1, wa2)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(3)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, cr2, di2, di3, dr2, dr3, ti2, ti3, tr2, tr3;
    static Mfloat taur = -0.5e0;
    static Mfloat taui = 0.866025403784438646763723170753e0;



    for (k = 1; k <= sl1; k++) {
	cr2 = *CC (1, k - 1, 0) + *CC (2, k - 1, 0);
	*CH (k - 1, 0, 0) = *CC (0, k - 1, 0) + cr2;
	*CH (k - 1, 2, 0) = taui * (*CC (2, k - 1, 0) - *CC (1, k - 1, 0));
	*CH (k - 1, 1, sido - 1) = *CC (0, k - 1, 0) + taur * cr2;
    }

    if (sido != 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    dr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    di2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    dr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    di3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    cr2 = dr2 + dr3;
		    ci2 = di2 + di3;
		    *CH (k - 1, 0, i - 2) = *CC (0, k - 1, i - 2) + cr2;
		    *CH (k - 1, 0, i - 1) = *CC (0, k - 1, i - 1) + ci2;
		    tr2 = *CC (0, k - 1, i - 2) + taur * cr2;
		    ti2 = *CC (0, k - 1, i - 1) + taur * ci2;
		    tr3 = taui * (di2 - di3);
		    ti3 = taui * (dr3 - dr2);
		    *CH (k - 1, 2, i - 2) = tr2 + tr3;
		    *CH (k - 1, 1, ic - 2) = tr2 - tr3;
		    *CH (k - 1, 2, i - 1) = ti2 + ti3;
		    *CH (k - 1, 1, ic - 1) = ti3 - ti2;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    dr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    di2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    dr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    di3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    cr2 = dr2 + dr3;
		    ci2 = di2 + di3;
		    *CH (k - 1, 0, i - 2) = *CC (0, k - 1, i - 2) + cr2;
		    *CH (k - 1, 0, i - 1) = *CC (0, k - 1, i - 1) + ci2;
		    tr2 = *CC (0, k - 1, i - 2) + taur * cr2;
		    ti2 = *CC (0, k - 1, i - 1) + taur * ci2;
		    tr3 = taui * (di2 - di3);
		    ti3 = taui * (dr3 - dr2);
		    *CH (k - 1, 2, i - 2) = tr2 + tr3;
		    *CH (k - 1, 1, ic - 2) = tr2 - tr3;
		    *CH (k - 1, 2, i - 1) = ti2 + ti3;
		    *CH (k - 1, 1, ic - 1) = ti3 - ti2;
		}
	    }
	}
    }
    return;
}				/* end of function */

#undef CC
#undef CH



/*Translated by FOR_C++, v0.1, on 06/11/90 at 10:02:05 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 10:01:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F6TRF/DF6TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence when the length is divisible by four.

    Usage:      CALL F6TRF (IDO, L1, CC, CH, WA1, WA2, WA3)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)
       WA3    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f6trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3)
#else
static void l_f6trf (ido, l1, cc, ch, wa1, wa2, wa3)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(4)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1,
                tr2, tr3, tr4;
    static Mfloat hsqt2 = 0.707106781186547524400844362105e0;
    for (k = 1; k <= sl1; k++) {
	tr1 = *CC (1, k - 1, 0) + *CC (3, k - 1, 0);
	tr2 = *CC (0, k - 1, 0) + *CC (2, k - 1, 0);
	*CH (k - 1, 0, 0) = tr1 + tr2;
	*CH (k - 1, 3, sido - 1) = tr2 - tr1;
	*CH (k - 1, 1, sido - 1) = *CC (0, k - 1, 0) - *CC (2, k - 1, 0);
	*CH (k - 1, 2, 0) = *CC (3, k - 1, 0) - *CC (1, k - 1, 0);
    }

    if (sido > 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    cr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    ci2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    cr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    ci3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    cr4 = wa3[i - 3] ** CC (3, k - 1, i - 2) + wa3[i - 2] *
			*CC (3, k - 1, i - 1);
		    ci4 = wa3[i - 3] ** CC (3, k - 1, i - 1) - wa3[i - 2] *
			*CC (3, k - 1, i - 2);
		    tr1 = cr2 + cr4;
		    tr4 = cr4 - cr2;
		    ti1 = ci2 + ci4;
		    ti4 = ci2 - ci4;
		    ti2 = *CC (0, k - 1, i - 1) + ci3;
		    ti3 = *CC (0, k - 1, i - 1) - ci3;
		    tr2 = *CC (0, k - 1, i - 2) + cr3;
		    tr3 = *CC (0, k - 1, i - 2) - cr3;
		    *CH (k - 1, 0, i - 2) = tr1 + tr2;
		    *CH (k - 1, 3, ic - 2) = tr2 - tr1;
		    *CH (k - 1, 0, i - 1) = ti1 + ti2;
		    *CH (k - 1, 3, ic - 1) = ti1 - ti2;
		    *CH (k - 1, 2, i - 2) = ti4 + tr3;
		    *CH (k - 1, 1, ic - 2) = tr3 - ti4;
		    *CH (k - 1, 2, i - 1) = tr4 + ti3;
		    *CH (k - 1, 1, ic - 1) = tr4 - ti3;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    cr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    ci2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    cr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    ci3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    cr4 = wa3[i - 3] ** CC (3, k - 1, i - 2) + wa3[i - 2] *
			*CC (3, k - 1, i - 1);
		    ci4 = wa3[i - 3] ** CC (3, k - 1, i - 1) - wa3[i - 2] *
			*CC (3, k - 1, i - 2);
		    tr1 = cr2 + cr4;
		    tr4 = cr4 - cr2;
		    ti1 = ci2 + ci4;
		    ti4 = ci2 - ci4;
		    ti2 = *CC (0, k - 1, i - 1) + ci3;
		    ti3 = *CC (0, k - 1, i - 1) - ci3;
		    tr2 = *CC (0, k - 1, i - 2) + cr3;
		    tr3 = *CC (0, k - 1, i - 2) - cr3;
		    *CH (k - 1, 0, i - 2) = tr1 + tr2;
		    *CH (k - 1, 3, ic - 2) = tr2 - tr1;
		    *CH (k - 1, 0, i - 1) = ti1 + ti2;
		    *CH (k - 1, 3, ic - 1) = ti1 - ti2;
		    *CH (k - 1, 2, i - 2) = ti4 + tr3;
		    *CH (k - 1, 1, ic - 2) = tr3 - ti4;
		    *CH (k - 1, 2, i - 1) = tr4 + ti3;
		    *CH (k - 1, 1, ic - 1) = tr4 - ti3;
		}
	    }
	}

	if (mod (sido, 2) == 0) {
	    for (k = 1; k <= sl1; k++) {
		ti1 = -hsqt2 * (*CC (1, k - 1, sido - 1) + *CC (3, k - 1, sido - 1));
		tr1 = hsqt2 * (*CC (1, k - 1, sido - 1) - *CC (3, k - 1, sido - 1));
		*CH (k - 1, 0, sido - 1) = tr1 + *CC (0, k - 1, sido - 1);
		*CH (k - 1, 2, sido - 1) = *CC (0, k - 1, sido - 1) - tr1;
		*CH (k - 1, 1, 0) = ti1 - *CC (2, k - 1, sido - 1);
		*CH (k - 1, 3, 0) = ti1 + *CC (2, k - 1, sido - 1);
	    }
	}
    }
    return;
}				/* end of function */

#undef CC
#undef CH



/*Translated by FOR_C++, v0.1, on 06/11/90 at 10:04:35 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 10:04:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F7TRF/DF7TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence when the length is divisible by five.

    Usage:      CALL F7TRF (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)
       WA3    - Real vector.  (Input)
       WA4    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f7trf (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat *wa1, Mfloat *wa2, Mfloat *wa3, Mfloat *wa4)
#else
static void l_f7trf (ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[], wa4[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(5)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3, di4,
                di5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3,
                tr4, tr5;
    static Mfloat tr11 = 0.309016994374947424102293417183e0;
    static Mfloat ti11 = 0.951056516295153572116439333379e0;
    static Mfloat tr12 = -0.809016994374947424102293417183e0;
    static Mfloat ti12 = 0.587785252292473129168705954639e0;
    for (k = 1; k <= sl1; k++) {
	cr2 = *CC (4, k - 1, 0) + *CC (1, k - 1, 0);
	ci5 = *CC (4, k - 1, 0) - *CC (1, k - 1, 0);
	cr3 = *CC (3, k - 1, 0) + *CC (2, k - 1, 0);
	ci4 = *CC (3, k - 1, 0) - *CC (2, k - 1, 0);
	*CH (k - 1, 0, 0) = *CC (0, k - 1, 0) + cr2 + cr3;
	*CH (k - 1, 1, sido - 1) = *CC (0, k - 1, 0) + tr11 * cr2 + tr12 * cr3;
	*CH (k - 1, 2, 0) = ti11 * ci5 + ti12 * ci4;
	*CH (k - 1, 3, sido - 1) = *CC (0, k - 1, 0) + tr12 * cr2 + tr11 * cr3;
	*CH (k - 1, 4, 0) = ti12 * ci5 - ti11 * ci4;
    }

    if (sido != 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    dr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    di2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    dr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    di3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    dr4 = wa3[i - 3] ** CC (3, k - 1, i - 2) + wa3[i - 2] *
			*CC (3, k - 1, i - 1);
		    di4 = wa3[i - 3] ** CC (3, k - 1, i - 1) - wa3[i - 2] *
			*CC (3, k - 1, i - 2);
		    dr5 = wa4[i - 3] ** CC (4, k - 1, i - 2) + wa4[i - 2] *
			*CC (4, k - 1, i - 1);
		    di5 = wa4[i - 3] ** CC (4, k - 1, i - 1) - wa4[i - 2] *
			*CC (4, k - 1, i - 2);
		    cr2 = dr2 + dr5;
		    ci5 = dr5 - dr2;
		    cr5 = di2 - di5;
		    ci2 = di2 + di5;
		    cr3 = dr3 + dr4;
		    ci4 = dr4 - dr3;
		    cr4 = di3 - di4;
		    ci3 = di3 + di4;
		    *CH (k - 1, 0, i - 2) = *CC (0, k - 1, i - 2) + cr2 +
			cr3;
		    *CH (k - 1, 0, i - 1) = *CC (0, k - 1, i - 1) + ci2 +
			ci3;
		    tr2 = *CC (0, k - 1, i - 2) + tr11 * cr2 + tr12 * cr3;
		    ti2 = *CC (0, k - 1, i - 1) + tr11 * ci2 + tr12 * ci3;
		    tr3 = *CC (0, k - 1, i - 2) + tr12 * cr2 + tr11 * cr3;
		    ti3 = *CC (0, k - 1, i - 1) + tr12 * ci2 + tr11 * ci3;
		    tr5 = ti11 * cr5 + ti12 * cr4;
		    ti5 = ti11 * ci5 + ti12 * ci4;
		    tr4 = ti12 * cr5 - ti11 * cr4;
		    ti4 = ti12 * ci5 - ti11 * ci4;
		    *CH (k - 1, 2, i - 2) = tr2 + tr5;
		    *CH (k - 1, 1, ic - 2) = tr2 - tr5;
		    *CH (k - 1, 2, i - 1) = ti2 + ti5;
		    *CH (k - 1, 1, ic - 1) = ti5 - ti2;
		    *CH (k - 1, 4, i - 2) = tr3 + tr4;
		    *CH (k - 1, 3, ic - 2) = tr3 - tr4;
		    *CH (k - 1, 4, i - 1) = ti3 + ti4;
		    *CH (k - 1, 3, ic - 1) = ti4 - ti3;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    dr2 = wa1[i - 3] ** CC (1, k - 1, i - 2) + wa1[i - 2] *
			*CC (1, k - 1, i - 1);
		    di2 = wa1[i - 3] ** CC (1, k - 1, i - 1) - wa1[i - 2] *
			*CC (1, k - 1, i - 2);
		    dr3 = wa2[i - 3] ** CC (2, k - 1, i - 2) + wa2[i - 2] *
			*CC (2, k - 1, i - 1);
		    di3 = wa2[i - 3] ** CC (2, k - 1, i - 1) - wa2[i - 2] *
			*CC (2, k - 1, i - 2);
		    dr4 = wa3[i - 3] ** CC (3, k - 1, i - 2) + wa3[i - 2] *
			*CC (3, k - 1, i - 1);
		    di4 = wa3[i - 3] ** CC (3, k - 1, i - 1) - wa3[i - 2] *
			*CC (3, k - 1, i - 2);
		    dr5 = wa4[i - 3] ** CC (4, k - 1, i - 2) + wa4[i - 2] *
			*CC (4, k - 1, i - 1);
		    di5 = wa4[i - 3] ** CC (4, k - 1, i - 1) - wa4[i - 2] *
			*CC (4, k - 1, i - 2);
		    cr2 = dr2 + dr5;
		    ci5 = dr5 - dr2;
		    cr5 = di2 - di5;
		    ci2 = di2 + di5;
		    cr3 = dr3 + dr4;
		    ci4 = dr4 - dr3;
		    cr4 = di3 - di4;
		    ci3 = di3 + di4;
		    *CH (k - 1, 0, i - 2) = *CC (0, k - 1, i - 2) + cr2 +
			cr3;
		    *CH (k - 1, 0, i - 1) = *CC (0, k - 1, i - 1) + ci2 +
			ci3;
		    tr2 = *CC (0, k - 1, i - 2) + tr11 * cr2 + tr12 * cr3;
		    ti2 = *CC (0, k - 1, i - 1) + tr11 * ci2 + tr12 * ci3;
		    tr3 = *CC (0, k - 1, i - 2) + tr12 * cr2 + tr11 * cr3;
		    ti3 = *CC (0, k - 1, i - 1) + tr12 * ci2 + tr11 * ci3;
		    tr5 = ti11 * cr5 + ti12 * cr4;
		    ti5 = ti11 * ci5 + ti12 * ci4;
		    tr4 = ti12 * cr5 - ti11 * cr4;
		    ti4 = ti12 * ci5 - ti11 * ci4;
		    *CH (k - 1, 2, i - 2) = tr2 + tr5;
		    *CH (k - 1, 1, ic - 2) = tr2 - tr5;
		    *CH (k - 1, 2, i - 1) = ti2 + ti5;
		    *CH (k - 1, 1, ic - 1) = ti5 - ti2;
		    *CH (k - 1, 4, i - 2) = tr3 + tr4;
		    *CH (k - 1, 3, ic - 2) = tr3 - tr4;
		    *CH (k - 1, 4, i - 1) = ti3 + ti4;
		    *CH (k - 1, 3, ic - 1) = ti4 - ti3;
		}
	    }
	}
    }
    return;
}				/* end of function */

#undef CC
#undef CH




/*Translated by FOR_C++, v0.1, on 06/11/90 at 10:06:40 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 10:06:35
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F8TRF/DF8TRF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 29, 1990

    Purpose:    Compute the Fourier coefficients of a real periodic
                sequence when the length is divisible by an integer
                other than two, three, four or five.

    Usage:      CALL F8TRF (IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)

    Arguments:
       IDO    - Integer scalar.  (Input)
       IP     - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       IDL1   - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       C1     - Real array.  (Input/Output)
       C2     - Real array.  (Input/Output)
       CH     - Real array.  (Output)
       CH2    - Real array.  (Output)
       WA     - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	LSIZE	256
#ifdef ANSI
static void l_f8trf (Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat *wa)
#else
static void l_f8trf (ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    Mint       *ido, *ip, *l1, *idl1;
    Mfloat     *cc, *c1, *c2, *ch, *ch2, wa[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(sip)*(sido)+(J_)*(sido)+(K_))
#define C1(I_,J_,K_)	(c1+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define C2(I_,J_)	(c2+(I_)*(sidl1)+(J_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH2(I_,J_)	(ch2+(I_)*(sidl1)+(J_))
    Mint sl1 = *l1;
    Mint sidl1 = *idl1;
    Mint        sido = *ido;
    Mint        sip = *ip;
    Mint        i, ic, idp2, ik, ipp2, ipph, is, j, j2, j2s, j3s, jc, k,
                l, lc, nbd;
    Mfloat      ai1, ai2, ai21u[LSIZE], ar1, ar1h, ar2, ar21u[LSIZE], ar2h,
                arg, dc2, dc21u[LSIZE], dcp, ds2, ds21u[LSIZE], dsp, r1s,
                r1v[LSIZE], r2v[LSIZE];


    arg = F_TWO * 3.1415926535897932384626433831 / (Mfloat) (sip);
    dcp = cos (arg);
    dsp = sin (arg);
    ipph = (sip + 1) / 2;
    ipp2 = sip + 2;
    idp2 = sido + 2;
    nbd = (sido - 1) / 2;

    if (sido == 1) {
	*C2 (0, 0) = *CH2 (0, 0);
	if (*idl1 != 1) {
	    for (i = 2; i <= *idl1; i++) {
		*C2 (0, i - 1) = *CH2 (0, i - 1);
	    }
	}
    }
    else {
	*CH2 (0, 0) = *C2 (0, 0);
	if (*idl1 != 1) {
	    for (i = 2; i <= *idl1; i++) {
		*CH2 (0, i - 1) = *C2 (0, i - 1);
	    }
	}
	for (j = 2; j <= sip; j++) {
	    for (k = 1; k <= *l1; k++) {
		*CH (j - 1, k - 1, 0) = *C1 (j - 1, k - 1, 0);
	    }
	}
	if (nbd < *l1) {
	    is = -sido;
	    for (j = 1; j <= (sip - 1); j++) {
		for (i = 1; i <= ((sido - 1) / 2); i++) {
		    for (k = 1; k <= *l1; k++) {
			*CH (j, k - 1, i * 2 - 1) = wa[is + sido * j + (i - 1) * 2] *
			    *C1 (j, k - 1, i * 2 - 1) + wa[is + sido * j + i * 2 - 1] *
			    *C1 (j, k - 1, (i - 1) * 2 + 2);
			*CH (j, k - 1, (i - 1) * 2 + 2) = wa[is + sido * j + (i - 1) * 2] *
			    *C1 (j, k - 1, (i - 1) * 2 + 2) - wa[is + sido * j + i * 2 - 1] *
			    *C1 (j, k - 1, i * 2 - 1);
		    }
		}
	    }

	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (i = 3; i <= sido; i += 2) {
		    for (k = 1; k <= *l1; k++) {
			*C1 (j - 1, k - 1, i - 2) = *CH (j - 1, k - 1, i - 2) +
			    *CH (jc - 1, k - 1, i - 2);
			*C1 (jc - 1, k - 1, i - 2) = *CH (j - 1, k - 1, i - 1) -
			    *CH (jc - 1, k - 1, i - 1);
			*C1 (j - 1, k - 1, i - 1) = *CH (j - 1, k - 1, i - 1) +
			    *CH (jc - 1, k - 1, i - 1);
			*C1 (jc - 1, k - 1, i - 1) = *CH (jc - 1, k - 1, i - 2) -
			    *CH (j - 1, k - 1, i - 2);
		    }
		}
	    }
	}
	else {
	    is = -sido;
	    for (j = 1; j <= (sip - 1); j++) {
		for (k = 1; k <= *l1; k++) {
		    for (i = 1; i <= ((sido - 1) / 2); i++) {
			*CH (j, k - 1, i * 2 - 1) = wa[is + sido * j + i * 2 - 2] *
			    *C1 (j, k - 1, i * 2 - 1) + wa[is + sido * j + i * 2 - 1] *
			    *C1 (j, k - 1, i * 2);
			*CH (j, k - 1, i * 2) = wa[is + sido * j + i * 2 - 2] *
			    *C1 (j, k - 1, i * 2) - wa[is + sido * j + i * 2 - 1] *
			    *C1 (j, k - 1, i * 2 - 1);
		    }
		}
	    }

	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (k = 1; k <= *l1; k++) {
		    for (i = 3; i <= sido; i += 2) {
			*C1 (j - 1, k - 1, i - 2) = *CH (j - 1, k - 1, i - 2) +
			    *CH (jc - 1, k - 1, i - 2);
			*C1 (jc - 1, k - 1, i - 2) = *CH (j - 1, k - 1, i - 1) -
			    *CH (jc - 1, k - 1, i - 1);
			*C1 (j - 1, k - 1, i - 1) = *CH (j - 1, k - 1, i - 1) +
			    *CH (jc - 1, k - 1, i - 1);
			*C1 (jc - 1, k - 1, i - 1) = *CH (jc - 1, k - 1, i - 2) -
			    *CH (j - 1, k - 1, i - 2);
		    }
		}
	    }
	}
    }

    for (j = 2; j <= ipph; j++) {
	jc = ipp2 - j;
	for (k = 1; k <= *l1; k++) {
	    *C1 (j - 1, k - 1, 0) = *CH (j - 1, k - 1, 0) + *CH (jc - 1, k - 1, 0);
	    *C1 (jc - 1, k - 1, 0) = *CH (jc - 1, k - 1, 0) - *CH (j - 1, k - 1, 0);
	}
    }

    ar1 = F_ONE;
    ai1 = F_ZERO;
    if (*idl1 != 1) {
	for (l = 2; l <= ipph; l++) {
	    lc = ipp2 - l;
	    ar1h = dcp * ar1 - dsp * ai1;
	    ai1 = dcp * ai1 + dsp * ar1;
	    ar1 = ar1h;

	    for (ik = 1; ik <= *idl1; ik++) {
		*CH2 (l - 1, ik - 1) = *C2 (0, ik - 1) + ar1 ** C2 (1, ik - 1);
		*CH2 (lc - 1, ik - 1) = ai1 ** C2 (sip - 1, ik - 1);
	    }

	    dc2 = ar1;
	    ds2 = ai1;
	    ar2 = ar1;
	    ai2 = ai1;

	    for (j = 3; j <= ipph; j++) {
		jc = ipp2 - j;
		ar2h = dc2 * ar2 - ds2 * ai2;
		ai2 = dc2 * ai2 + ds2 * ar2;
		ar2 = ar2h;
		for (ik = 1; ik <= *idl1; ik++) {
		    *CH2 (l - 1, ik - 1) += ar2 ** C2 (j - 1, ik - 1);
		    *CH2 (lc - 1, ik - 1) += ai2 ** C2 (jc - 1, ik - 1);
		}
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    for (ik = 1; ik <= *idl1; ik++) {
		*CH2 (0, ik - 1) += *C2 (j - 1, ik - 1);
	    }
	}
    }
    else {
	for (j2s = 0; j2s <= (ipph - 2); j2s += LSIZE) {
	    j3s = imsl_i_min (ipph - 1 - j2s, LSIZE);
	    for (l = 1; l <= j3s; l++) {
		ar1h = dcp * ar1 - dsp * ai1;
		ai1 = dcp * ai1 + dsp * ar1;
		ar1 = ar1h;
		r1v[l - 1] = ar1;
		dc21u[l - 1] = ar1;
		ar21u[l - 1] = ar1;
		ds21u[l - 1] = ai1;
		r2v[l - 1] = ai1;
		ai21u[l - 1] = ai1;
	    }

	    for (l = 1; l <= j3s; l++) {
		*CH2 (j2s + l, 0) = *C2 (0, 0) + r1v[l - 1] ** C2 (1, 0);
		*CH2 (ipp2 - j2s - l - 2, 0) = r2v[l - 1] ** C2 (sip - 1, 0);
	    }

	    for (j = 1; j <= (ipph - 2); j++) {
		for (l = 1; l <= j3s; l++) {
		    r1s = ds21u[l - 1] * ar21u[l - 1];
		    ar2h = dc21u[l - 1] * ar21u[l - 1] - ds21u[l - 1] *
			ai21u[l - 1];
		    ar21u[l - 1] = ar2h;
		    *CH2 (j2s + l, 0) += ar21u[l - 1] ** C2 (j + 1, 0);
		    ai21u[l - 1] = dc21u[l - 1] * ai21u[l - 1] + r1s;
		    *CH2 (ipp2 - j2s - l - 2, 0) += ai21u[l - 1] ** C2 (ipp2 - (j - 1) - 4, 0);
		}
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    *CH2 (0, 0) += *C2 (j - 1, 0);
	}
    }

    if (sido >= *l1) {
	for (k = 1; k <= *l1; k++) {
	    for (i = 1; i <= sido; i++) {
		*CC (k - 1, 0, i - 1) = *CH (0, k - 1, i - 1);
	    }
	}
    }
    else {
	for (i = 1; i <= sido; i++) {
	    for (k = 1; k <= *l1; k++) {
		*CC (k - 1, 0, i - 1) = *CH (0, k - 1, i - 1);
	    }
	}
    }

    for (j = 2; j <= ipph; j++) {
	jc = ipp2 - j;
	j2 = j + j;
	for (k = 1; k <= *l1; k++) {
	    *CC (k - 1, j2 - 3, sido - 1) = *CH (j - 1, k - 1, 0);
	    *CC (k - 1, j2 - 2, 0) = *CH (jc - 1, k - 1, 0);
	}
    }

    if (sido != 1) {
	if (nbd >= *l1) {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		j2 = j + j;
		for (k = 1; k <= *l1; k++) {
		    for (i = 3; i <= sido; i += 2) {
			ic = idp2 - i;
			*CC (k - 1, j2 - 2, i - 2) = *CH (j - 1, k - 1, i - 2) +
			    *CH (jc - 1, k - 1, i - 2);
			*CC (k - 1, j2 - 3, ic - 2) = *CH (j - 1, k - 1, i - 2) -
			    *CH (jc - 1, k - 1, i - 2);
			*CC (k - 1, j2 - 2, i - 1) = *CH (j - 1, k - 1, i - 1) +
			    *CH (jc - 1, k - 1, i - 1);
			*CC (k - 1, j2 - 3, ic - 1) = *CH (jc - 1, k - 1, i - 1) -
			    *CH (j - 1, k - 1, i - 1);
		    }
		}
	    }
	}
	else {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		j2 = j + j;
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    for (k = 1; k <= *l1; k++) {
			*CC (k - 1, j2 - 2, i - 2) = *CH (j - 1, k - 1, i - 2) +
			    *CH (jc - 1, k - 1, i - 2);
			*CC (k - 1, j2 - 3, ic - 2) = *CH (j - 1, k - 1, i - 2) -
			    *CH (jc - 1, k - 1, i - 2);
			*CC (k - 1, j2 - 2, i - 1) = *CH (j - 1, k - 1, i - 1) +
			    *CH (jc - 1, k - 1, i - 1);
			*CC (k - 1, j2 - 3, ic - 1) = *CH (jc - 1, k - 1, i - 1) -
			    *CH (j - 1, k - 1, i - 1);
		    }
		}
	    }
	}
    }
    return;
}				/* end of function */

#undef CC
#undef C1
#undef C2
#undef CH
#undef CH2


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:49:19 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:49:17
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2TRB/DF2TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F2TRB (N, COEF, SEQ, WFFTR)

    Arguments:  (See FFTRB).

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f2trb (Mint *n, Mfloat *coef, Mfloat *seq, Mfloat *wfftr)
#else
void imsl_f2trb (n, coef, seq, wfftr)
    Mint       *n;
    Mfloat      coef[], seq[], wfftr[];
#endif
{


    /* CHECK ARGUMENT N */
    if (*n < 1) {
	imsl_e1psh ("l_f2trb");
	/* The length of the sequence n = %(I1). */
	/* It must be at least one. */
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	imsl_e1pop ("l_f2trb");
	goto L_9000;
    }
    /* COPY COEF TO SEQ */
    scopy (*n, coef, 1, seq, 1);

    if (*n > 1) {
	l_f3trb (n, seq, wfftr, &wfftr[*n], &wfftr[*n * 2]);
    }
L_9000:
    ;
    return;
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:50:33 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:50:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F3TRB/DF3TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F3TRB (N, C, CH, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       C      - Real array.  (Input/Output)
       CH     - Real array.  (Workspace)
       WA     - Real array from FFTRI.  (Input)
       FAC    - Real array from FFTRI.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f3trb (Mint *n, Mfloat *c, Mfloat *ch, Mfloat *wa, Mfloat *fac)
#else
static void l_f3trb (n, c, ch, wa, fac)
    Mint       *n;
    Mfloat      c[], ch[], wa[], fac[];
#endif
{
    Mint        i, ic, idl1, ido, idp2, ip, iw, ix2, ix3, ix4, k1, l1,
                l2, na, nf;
    Mfloat      ti2, tr2;


    if (*n > 1) {
	nf = nint (fac[1]);
	na = 0;
	l1 = 1;
	iw = 1;
	for (k1 = 1; k1 <= nf; k1++) {
	    ip = nint (fac[k1 + 1]);
	    l2 = ip * l1;
	    ido = *n / l2;
	    idl1 = ido * l1;
	    if (ip == 4) {
		ix2 = iw + ido;
		ix3 = ix2 + ido;
		if (na != 0) {
		    l_f6trb (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
			&wa[ix3 - 1]);
		}
		else {
		    l_f6trb (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
			&wa[ix3 - 1]);
		}
		na = 1 - na;
	    }
	    else if (ip == 2) {
		if (na != 0) {
		    /*
		     * CALL F4TRB (IDO, L1, CH, C, WA(IW))
		     */
		    c[0] = ch[0] + ch[ido * 2 - 1];
		    c[ido] = ch[0] - ch[ido * 2 - 1];
		    if (ido != 1) {
			idp2 = ido + 2;
			for (i = 3; i <= ido; i += 2) {
			    ic = idp2 - i;
			    c[i - 2] = ch[i - 2] + ch[ic + ido - 2];
			    tr2 = ch[i - 2] - ch[ic + ido - 2];
			    c[i - 1] = ch[i - 1] - ch[ic + ido - 1];
			    ti2 = ch[i - 1] + ch[ic + ido - 1];
			    c[i + ido - 2] = wa[iw + i - 4] * tr2 -
				wa[iw + i - 3] * ti2;
			    c[i + ido - 1] = wa[iw + i - 4] * ti2 +
				wa[iw + i - 3] * tr2;
			}
			if (mod (ido, 2) != 1) {
			    c[ido - 1] = ch[ido - 1] + ch[ido - 1];
			    c[ido * 2 - 1] = -(ch[ido] + ch[ido]);
			}
		    }
		}
		else {
		    /*
		     * CALL F4TRB (IDO, L1, C, CH, WA(IW))
		     */
		    ch[0] = c[0] + c[ido * 2 - 1];
		    ch[ido] = c[0] - c[ido * 2 - 1];
		    if (ido != 1) {
			idp2 = ido + 2;
			for (i = 3; i <= ido; i += 2) {
			    ic = idp2 - i;
			    ch[i - 2] = c[i - 2] + c[ic + ido - 2];
			    tr2 = c[i - 2] - c[ic + ido - 2];
			    ch[i - 1] = c[i - 1] - c[ic + ido - 1];
			    ti2 = c[i - 1] + c[ic + ido - 1];
			    ch[i + ido - 2] = wa[iw + i - 4] * tr2 -
				wa[iw + i - 3] * ti2;
			    ch[i + ido - 1] = wa[iw + i - 4] * ti2 +
				wa[iw + i - 3] * tr2;
			}
			if (mod (ido, 2) != 1) {
			    ch[ido - 1] = c[ido - 1] + c[ido - 1];
			    ch[ido * 2 - 1] = -(c[ido] + c[ido]);
			}
		    }
		}
		na = 1 - na;
	    }
	    else if (ip == 3) {
		ix2 = iw + ido;
		if (na != 0) {
		    l_f5trb (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1]);
		}
		else {
		    l_f5trb (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1]);
		}
		na = 1 - na;
	    }
	    else if (ip == 5) {
		ix2 = iw + ido;
		ix3 = ix2 + ido;
		ix4 = ix3 + ido;
		if (na != 0) {
		    l_f7trb (&ido, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
			&wa[ix3 - 1], &wa[ix4 - 1]);
		}
		else {
		    l_f7trb (&ido, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
			&wa[ix3 - 1], &wa[ix4 - 1]);
		}
		na = 1 - na;
	    }
	    else {
		if (na != 0) {
		    l_f8trb (&ido, &ip, &l1, &idl1, ch, ch, ch, c, c,
			&wa[iw - 1]);
		}
		else {
		    l_f8trb (&ido, &ip, &l1, &idl1, c, c, c, ch, ch,
			&wa[iw - 1]);
		}
		if (ido == 1)
		    na = 1 - na;
	    }

	    l1 = l2;
	    iw += (ip - 1) * ido;
	}
	if (na != 0)
	    scopy (*n, ch, 1, c, 1);
    }
    return;
}				/* end of function */




/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:51:51 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:51:49
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F4TRB/DF4TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F4TRB (IDO, L1, CC, CH, WA1)

    Arguments:
       IDO    - Real scalar.  (Input)
       L1     - Real scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f4trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat *wa1)
#else
static void l_f4trb (ido, l1, cc, ch, wa1)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(2)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ti2, tr2;


    for (k = 1; k <= sl1; k++) {
	*CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + *CC (k - 1, 1, sido - 1);
	*CH (1, k - 1, 0) = *CC (k - 1, 0, 0) - *CC (k - 1, 1, sido - 1);
    }
    if (sido - 2 < 0)
	goto L_70;
    if (sido - 2 > 0)
	goto L_20;
    goto L_50;
L_20:
    idp2 = sido + 2;
    for (k = 1; k <= sl1; k++) {
	for (i = 3; i <= sido; i += 2) {
	    ic = idp2 - i;
	    *CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + *CC (k - 1, 1, ic - 2);
	    tr2 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 1, ic - 2);
	    *CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) - *CC (k - 1, 1, ic - 1);
	    ti2 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 1, ic - 1);
	    *CH (1, k - 1, i - 2) = wa1[i - 3] * tr2 - wa1[i - 2] * ti2;
	    *CH (1, k - 1, i - 1) = wa1[i - 3] * ti2 + wa1[i - 2] * tr2;
	}
    }
    if (mod (sido, 2) == 1)
	goto L_70;
L_50:
    for (k = 1; k <= sl1; k++) {
	*CH (0, k - 1, sido - 1) = *CC (k - 1, 0, sido - 1) + *CC (k - 1, 0, sido - 1);
	*CH (1, k - 1, sido - 1) = -(*CC (k - 1, 1, 0) + *CC (k - 1, 1, 0));
    }
L_70:
    return;
}				/* end of function */

#undef CC
#undef CH


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:53:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:53:05
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F5TRB/DF5TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients when the length is divisible by three.

    Usage:      CALL F5TRB (IDO, L1, CC, CH, WA1, WA2)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f5trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat *wa1,
                Mfloat *wa2)
#else
static void l_f5trb (ido, l1, cc, ch, wa1, wa2)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(3)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, ti2, tr2;
    static Mfloat taur = -0.5e0;
    static Mfloat taui = 0.866025403784438646763723170753e0;



    for (k = 1; k <= sl1; k++) {
	tr2 = *CC (k - 1, 1, sido - 1) + *CC (k - 1, 1, sido - 1);
	cr2 = *CC (k - 1, 0, 0) + taur * tr2;
	*CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2;
	ci3 = taui * (*CC (k - 1, 2, 0) + *CC (k - 1, 2, 0));
	*CH (1, k - 1, 0) = cr2 - ci3;
	*CH (2, k - 1, 0) = cr2 + ci3;
    }

    if (sido != 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    tr2 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		    *CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		    ti2 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		    *CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		    cr3 = taui * (*CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2));
		    ci3 = taui * (*CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1));
		    dr2 = cr2 - ci3;
		    dr3 = cr2 + ci3;
		    di2 = ci2 + cr3;
		    di3 = ci2 - cr3;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * dr2 - wa1[i - 2] *
			di2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * di2 + wa1[i - 2] *
			dr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * dr3 - wa2[i - 2] *
			di3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * di3 + wa2[i - 2] *
			dr3;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    tr2 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		    *CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		    ti2 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		    *CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		    cr3 = taui * (*CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2));
		    ci3 = taui * (*CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1));
		    dr2 = cr2 - ci3;
		    dr3 = cr2 + ci3;
		    di2 = ci2 + cr3;
		    di3 = ci2 - cr3;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * dr2 - wa1[i - 2] *
			di2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * di2 + wa1[i - 2] *
			dr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * dr3 - wa2[i - 2] *
			di3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * di3 + wa2[i - 2] *
			dr3;
		}
	    }
	}
    }
    return;
}				/* end of function */


#undef CC
#undef CH


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:54:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:54:31
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F6TRB/DF6TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients when the length is divisible by four.

    Usage:      CALL F6TRB (IDO, L1, CC, CH, WA1, WA2, WA3)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)
       WA3    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f6trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat *wa1,
                Mfloat *wa2, Mfloat *wa3)
#else
static void l_f6trb (ido, l1, cc, ch, wa1, wa2, wa3)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(4)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1,
                tr2, tr3, tr4;
    static Mfloat sqrt2 = 1.41421356237309504880168872421e0;
    for (k = 1; k <= sl1; k++) {
	tr1 = *CC (k - 1, 0, 0) - *CC (k - 1, 3, sido - 1);
	tr2 = *CC (k - 1, 0, 0) + *CC (k - 1, 3, sido - 1);
	tr3 = *CC (k - 1, 1, sido - 1) + *CC (k - 1, 1, sido - 1);
	tr4 = *CC (k - 1, 2, 0) + *CC (k - 1, 2, 0);
	*CH (0, k - 1, 0) = tr2 + tr3;
	*CH (1, k - 1, 0) = tr1 - tr4;
	*CH (2, k - 1, 0) = tr2 - tr3;
	*CH (3, k - 1, 0) = tr1 + tr4;
    }

    if (sido > 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    ti1 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 3, ic - 1);
		    ti2 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 3, ic - 1);
		    ti3 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    tr4 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1);
		    tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 3, ic - 2);
		    tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 3, ic - 2);
		    ti4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2);
		    tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    *CH (0, k - 1, i - 2) = tr2 + tr3;
		    cr3 = tr2 - tr3;
		    *CH (0, k - 1, i - 1) = ti2 + ti3;
		    ci3 = ti2 - ti3;
		    cr2 = tr1 - tr4;
		    cr4 = tr1 + tr4;
		    ci2 = ti1 + ti4;
		    ci4 = ti1 - ti4;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * cr2 - wa1[i - 2] *
			ci2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * ci2 + wa1[i - 2] *
			cr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * cr3 - wa2[i - 2] *
			ci3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * ci3 + wa2[i - 2] *
			cr3;
		    *CH (3, k - 1, i - 2) = wa3[i - 3] * cr4 - wa3[i - 2] *
			ci4;
		    *CH (3, k - 1, i - 1) = wa3[i - 3] * ci4 + wa3[i - 2] *
			cr4;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    ti1 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 3, ic - 1);
		    ti2 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 3, ic - 1);
		    ti3 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    tr4 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1);
		    tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 3, ic - 2);
		    tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 3, ic - 2);
		    ti4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2);
		    tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    *CH (0, k - 1, i - 2) = tr2 + tr3;
		    cr3 = tr2 - tr3;
		    *CH (0, k - 1, i - 1) = ti2 + ti3;
		    ci3 = ti2 - ti3;
		    cr2 = tr1 - tr4;
		    cr4 = tr1 + tr4;
		    ci2 = ti1 + ti4;
		    ci4 = ti1 - ti4;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * cr2 - wa1[i - 2] *
			ci2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * ci2 + wa1[i - 2] *
			cr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * cr3 - wa2[i - 2] *
			ci3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * ci3 + wa2[i - 2] *
			cr3;
		    *CH (3, k - 1, i - 2) = wa3[i - 3] * cr4 - wa3[i - 2] *
			ci4;
		    *CH (3, k - 1, i - 1) = wa3[i - 3] * ci4 + wa3[i - 2] *
			cr4;
		}
	    }
	}

	if (mod (sido, 2) != 1) {
	    for (k = 1; k <= sl1; k++) {
		ti1 = *CC (k - 1, 1, 0) + *CC (k - 1, 3, 0);
		ti2 = *CC (k - 1, 3, 0) - *CC (k - 1, 1, 0);
		tr1 = *CC (k - 1, 0, sido - 1) - *CC (k - 1, 2, sido - 1);
		tr2 = *CC (k - 1, 0, sido - 1) + *CC (k - 1, 2, sido - 1);
		*CH (0, k - 1, sido - 1) = tr2 + tr2;
		*CH (1, k - 1, sido - 1) = sqrt2 * (tr1 - ti1);
		*CH (2, k - 1, sido - 1) = ti2 + ti2;
		*CH (3, k - 1, sido - 1) = -sqrt2 * (tr1 + ti1);
	    }
	}
    }
    return;
}				/* end of function */


#undef CC
#undef CH


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:56:01 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:55:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F7TRB/DF7TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 28, 1990

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients when the length is divisible by five.

    Usage:      CALL F7TRB (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)

    Arguments:
       IDO    - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       CH     - Real array.  (Input/Output)
       WA1    - Real vector.  (Input)
       WA2    - Real vector.  (Input)
       WA3    - Real vector.  (Input)
       WA4    - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_f7trb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat *wa1,
                Mfloat *wa2, Mfloat *wa3, Mfloat *wa4)
#else
static void l_f7trb (ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[], wa4[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(5)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, ic, idp2, k;
    Mfloat      ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3, di4;
    Mfloat      di5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3;
    Mfloat      tr4, tr5;
    static Mfloat tr11 = 0.309016994374947424102293417183e0;
    static Mfloat ti11 = 0.951056516295153572116439333379e0;
    static Mfloat tr12 = -0.809016994374947424102293417183e0;
    static Mfloat ti12 = 0.587785252292473129168705954639e0;
    for (k = 1; k <= sl1; k++) {
	ti5 = *CC (k - 1, 2, 0) + *CC (k - 1, 2, 0);
	ti4 = *CC (k - 1, 4, 0) + *CC (k - 1, 4, 0);
	tr2 = *CC (k - 1, 1, sido - 1) + *CC (k - 1, 1, sido - 1);
	tr3 = *CC (k - 1, 3, sido - 1) + *CC (k - 1, 3, sido - 1);
	*CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2 + tr3;
	cr2 = *CC (k - 1, 0, 0) + tr11 * tr2 + tr12 * tr3;
	cr3 = *CC (k - 1, 0, 0) + tr12 * tr2 + tr11 * tr3;
	ci5 = ti11 * ti5 + ti12 * ti4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	*CH (1, k - 1, 0) = cr2 - ci5;
	*CH (2, k - 1, 0) = cr3 - ci4;
	*CH (3, k - 1, 0) = cr3 + ci4;
	*CH (4, k - 1, 0) = cr2 + ci5;
    }

    if (sido > 1) {
	idp2 = sido + 2;
	if ((sido - 1) / 2 >= sl1) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    ti5 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1);
		    ti2 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    ti4 = *CC (k - 1, 4, i - 1) + *CC (k - 1, 3, ic - 1);
		    ti3 = *CC (k - 1, 4, i - 1) - *CC (k - 1, 3, ic - 1);
		    tr5 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2);
		    tr2 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    tr4 = *CC (k - 1, 4, i - 2) - *CC (k - 1, 3, ic - 2);
		    tr3 = *CC (k - 1, 4, i - 2) + *CC (k - 1, 3, ic - 2);
		    *CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 +
			tr3;
		    *CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 +
			ti3;
		    cr2 = *CC (k - 1, 0, i - 2) + tr11 * tr2 + tr12 * tr3;
		    ci2 = *CC (k - 1, 0, i - 1) + tr11 * ti2 + tr12 * ti3;
		    cr3 = *CC (k - 1, 0, i - 2) + tr12 * tr2 + tr11 * tr3;
		    ci3 = *CC (k - 1, 0, i - 1) + tr12 * ti2 + tr11 * ti3;
		    cr5 = ti11 * tr5 + ti12 * tr4;
		    ci5 = ti11 * ti5 + ti12 * ti4;
		    cr4 = ti12 * tr5 - ti11 * tr4;
		    ci4 = ti12 * ti5 - ti11 * ti4;
		    dr3 = cr3 - ci4;
		    dr4 = cr3 + ci4;
		    di3 = ci3 + cr4;
		    di4 = ci3 - cr4;
		    dr5 = cr2 + ci5;
		    dr2 = cr2 - ci5;
		    di5 = ci2 - cr5;
		    di2 = ci2 + cr5;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * dr2 - wa1[i - 2] *
			di2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * di2 + wa1[i - 2] *
			dr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * dr3 - wa2[i - 2] *
			di3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * di3 + wa2[i - 2] *
			dr3;
		    *CH (3, k - 1, i - 2) = wa3[i - 3] * dr4 - wa3[i - 2] *
			di4;
		    *CH (3, k - 1, i - 1) = wa3[i - 3] * di4 + wa3[i - 2] *
			dr4;
		    *CH (4, k - 1, i - 2) = wa4[i - 3] * dr5 - wa4[i - 2] *
			di5;
		    *CH (4, k - 1, i - 1) = wa4[i - 3] * di5 + wa4[i - 2] *
			dr5;
		}
	    }
	}
	else {
	    for (i = 3; i <= sido; i += 2) {
		for (k = 1; k <= sl1; k++) {
		    ic = idp2 - i;
		    ti5 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 1, ic - 1);
		    ti2 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 1, ic - 1);
		    ti4 = *CC (k - 1, 4, i - 1) + *CC (k - 1, 3, ic - 1);
		    ti3 = *CC (k - 1, 4, i - 1) - *CC (k - 1, 3, ic - 1);
		    tr5 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 1, ic - 2);
		    tr2 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 1, ic - 2);
		    tr4 = *CC (k - 1, 4, i - 2) - *CC (k - 1, 3, ic - 2);
		    tr3 = *CC (k - 1, 4, i - 2) + *CC (k - 1, 3, ic - 2);
		    *CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 +
			tr3;
		    *CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 +
			ti3;
		    cr2 = *CC (k - 1, 0, i - 2) + tr11 * tr2 + tr12 * tr3;
		    ci2 = *CC (k - 1, 0, i - 1) + tr11 * ti2 + tr12 * ti3;
		    cr3 = *CC (k - 1, 0, i - 2) + tr12 * tr2 + tr11 * tr3;
		    ci3 = *CC (k - 1, 0, i - 1) + tr12 * ti2 + tr11 * ti3;
		    cr5 = ti11 * tr5 + ti12 * tr4;
		    ci5 = ti11 * ti5 + ti12 * ti4;
		    cr4 = ti12 * tr5 - ti11 * tr4;
		    ci4 = ti12 * ti5 - ti11 * ti4;
		    dr3 = cr3 - ci4;
		    dr4 = cr3 + ci4;
		    di3 = ci3 + cr4;
		    di4 = ci3 - cr4;
		    dr5 = cr2 + ci5;
		    dr2 = cr2 - ci5;
		    di5 = ci2 - cr5;
		    di2 = ci2 + cr5;
		    *CH (1, k - 1, i - 2) = wa1[i - 3] * dr2 - wa1[i - 2] *
			di2;
		    *CH (1, k - 1, i - 1) = wa1[i - 3] * di2 + wa1[i - 2] *
			dr2;
		    *CH (2, k - 1, i - 2) = wa2[i - 3] * dr3 - wa2[i - 2] *
			di3;
		    *CH (2, k - 1, i - 1) = wa2[i - 3] * di3 + wa2[i - 2] *
			dr3;
		    *CH (3, k - 1, i - 2) = wa3[i - 3] * dr4 - wa3[i - 2] *
			di4;
		    *CH (3, k - 1, i - 1) = wa3[i - 3] * di4 + wa3[i - 2] *
			dr4;
		    *CH (4, k - 1, i - 2) = wa4[i - 3] * dr5 - wa4[i - 2] *
			di5;
		    *CH (4, k - 1, i - 1) = wa4[i - 3] * di5 + wa4[i - 2] *
			dr5;
		}
	    }
	}
    }
    return;
}				/* end of function */


#undef CC
#undef CH


/*Translated by FOR_C++, v0.1, on 06/11/90 at 13:57:45 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 13:57:41
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F8TRB/DF8TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 29, 1990

    Purpose:    Compute the real periodic sequence from its Fourier
                coefficients when the length is divisible by integers
                other than two, three, four or five.

    Usage:      CALL F8TRB (IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)

    Arguments:
       IDO    - Integer scalar.  (Input)
       IP     - Integer scalar.  (Input)
       L1     - Integer scalar.  (Input)
       IDL1   - Integer scalar.  (Input)
       CC     - Real array.  (Input/Output)
       C1     - Real array.  (Input/Output)
       C2     - Real array.  (Input/Output)
       CH     - Real array.  (Output)
       CH2    - Real array.  (Output)
       WA     - Real vector.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	LSIZE	256
#ifdef ANSI
static void l_f8trb (Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat *wa)
#else
static void l_f8trb (ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    Mint       *ido, *ip, *l1, *idl1;
    Mfloat     *cc, *c1, *c2, *ch, *ch2, wa[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(sip)*(sido)+(J_)*(sido)+(K_))
#define C1(I_,J_,K_)	(c1+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define C2(I_,J_)	(c2+(I_)*(sidl1)+(J_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH2(I_,J_)	(ch2+(I_)*(sidl1)+(J_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        sip = *ip;
    Mint        sidl1 = *idl1;
    Mint        i, ic, idp2, ik, ipp2, ipph, is, j, j2, j2s, j3s;
    Mint        jc, k, l, lc, nbd;
    Mfloat      ai1, ai2, ai21u[LSIZE], ar1, ar1h, ar2, ar21u[LSIZE];
    Mfloat      ar2h, arg, dc2, dc21u[LSIZE], dcp, ds2, ds21u[LSIZE];
    Mfloat      dsp, r1s, r1v[LSIZE], r2v[LSIZE];
    arg = F_TWO * 3.1415926535897932384626433831 / (Mfloat) (sip);
    dcp = cos (arg);
    dsp = sin (arg);
    idp2 = sido + 2;
    nbd = (sido - 1) / 2;
    ipp2 = sip + 2;
    ipph = (sip + 1) / 2;

    if (sido >= sl1) {
	for (k = 1; k <= sl1; k++) {
	    for (i = 1; i <= sido; i++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}
    }
    else {
	for (i = 1; i <= sido; i++) {
	    for (k = 1; k <= sl1; k++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}
    }

    for (j = 2; j <= ipph; j++) {
	jc = ipp2 - j;
	j2 = j + j;
	for (k = 1; k <= sl1; k++) {
	    *CH (j - 1, k - 1, 0) = *CC (k - 1, j2 - 3, sido - 1) + *CC (k - 1, j2 - 3, sido - 1);
	    *CH (jc - 1, k - 1, 0) = *CC (k - 1, j2 - 2, 0) + *CC (k - 1, j2 - 2, 0);
	}
    }

    if (sido != 1) {
	if (nbd >= sl1) {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (k = 1; k <= sl1; k++) {
		    for (i = 3; i <= sido; i += 2) {
			ic = idp2 - i;
			*CH (j - 1, k - 1, i - 2) = *CC (k - 1, j * 2 - 2, i - 2) +
			    *CC (k - 1, j * 2 - 3, ic - 2);
			*CH (jc - 1, k - 1, i - 2) = *CC (k - 1, j * 2 - 2, i - 2) -
			    *CC (k - 1, j * 2 - 3, ic - 2);
			*CH (j - 1, k - 1, i - 1) = *CC (k - 1, j * 2 - 2, i - 1) -
			    *CC (k - 1, j * 2 - 3, ic - 1);
			*CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j * 2 - 2, i - 1) +
			    *CC (k - 1, j * 2 - 3, ic - 1);
		    }
		}
	    }
	}
	else {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (i = 3; i <= sido; i += 2) {
		    ic = idp2 - i;
		    for (k = 1; k <= sl1; k++) {
			*CH (j - 1, k - 1, i - 2) = *CC (k - 1, j * 2 - 2, i - 2) +
			    *CC (k - 1, j * 2 - 3, ic - 2);
			*CH (jc - 1, k - 1, i - 2) = *CC (k - 1, j * 2 - 2, i - 2) -
			    *CC (k - 1, j * 2 - 3, ic - 2);
			*CH (j - 1, k - 1, i - 1) = *CC (k - 1, j * 2 - 2, i - 1) -
			    *CC (k - 1, j * 2 - 3, ic - 1);
			*CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j * 2 - 2, i - 1) +
			    *CC (k - 1, j * 2 - 3, ic - 1);
		    }
		}
	    }
	}
    }
    ar1 = F_ONE;
    ai1 = F_ZERO;
    if (sidl1 != 1) {
	for (l = 2; l <= ipph; l++) {
	    lc = ipp2 - l;
	    ar1h = dcp * ar1 - dsp * ai1;
	    ai1 = dcp * ai1 + dsp * ar1;
	    ar1 = ar1h;
	    for (ik = 1; ik <= sidl1; ik++) {
		*C2 (l - 1, ik - 1) = *CH2 (0, ik - 1) + ar1 ** CH2 (1, ik - 1);
		*C2 (lc - 1, ik - 1) = ai1 ** CH2 (sip - 1, ik - 1);
	    }

	    dc2 = ar1;
	    ds2 = ai1;
	    ar2 = ar1;
	    ai2 = ai1;

	    for (j = 3; j <= ipph; j++) {
		jc = ipp2 - j;
		ar2h = dc2 * ar2 - ds2 * ai2;
		ai2 = dc2 * ai2 + ds2 * ar2;
		ar2 = ar2h;
		for (ik = 1; ik <= sidl1; ik++) {
		    *C2 (l - 1, ik - 1) += ar2 ** CH2 (j - 1, ik - 1);
		    *C2 (lc - 1, ik - 1) += ai2 ** CH2 (jc - 1, ik - 1);
		}
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    for (ik = 1; ik <= sidl1; ik++) {
		*CH2 (0, ik - 1) += *CH2 (j - 1, ik - 1);
	    }
	}
    }
    else {
	for (j2s = 0; j2s <= (ipph - 2); j2s += LSIZE) {
	    j3s = imsl_i_min (ipph - 1 - j2s, LSIZE);
	    for (l = 1; l <= j3s; l++) {
		ar1h = dcp * ar1 - dsp * ai1;
		ai1 = dcp * ai1 + dsp * ar1;
		ar1 = ar1h;
		r1v[l - 1] = ar1;
		dc21u[l - 1] = ar1;
		ar21u[l - 1] = ar1;
		ds21u[l - 1] = ai1;
		r2v[l - 1] = ai1;
		ai21u[l - 1] = ai1;
	    }

	    for (l = 1; l <= j3s; l++) {
		*C2 (j2s + l, 0) = *CH2 (0, 0) + r1v[l - 1] ** CH2 (1, 0);
		*C2 (ipp2 - j2s - l - 2, 0) = r2v[l - 1] ** CH2 (sip - 1, 0);
	    }

	    for (j = 1; j <= (ipph - 2); j++) {
		for (l = 1; l <= j3s; l++) {
		    r1s = ds21u[l - 1] * ar21u[l - 1];
		    ar2h = dc21u[l - 1] * ar21u[l - 1] - ds21u[l - 1] *
			ai21u[l - 1];
		    ar21u[l - 1] = ar2h;
		    *C2 (j2s + l, 0) += ar21u[l - 1] ** CH2 (j + 1, 0);
		    ai21u[l - 1] = dc21u[l - 1] * ai21u[l - 1] + r1s;
		    *C2 (ipp2 - j2s - l - 2, 0) += ai21u[l - 1] ** CH2 (ipp2 - (j - 1) - 4, 0);
		}
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    *CH2 (0, 0) += *CH2 (j - 1, 0);
	}
    }

    for (j = 2; j <= ipph; j++) {
	jc = ipp2 - j;
	for (k = 1; k <= sl1; k++) {
	    *CH (j - 1, k - 1, 0) = *C1 (j - 1, k - 1, 0) - *C1 (jc - 1, k - 1, 0);
	    *CH (jc - 1, k - 1, 0) = *C1 (j - 1, k - 1, 0) + *C1 (jc - 1, k - 1, 0);
	}
    }

    if (sido != 1) {
	if (nbd >= sl1) {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (k = 1; k <= sl1; k++) {
		    for (i = 3; i <= sido; i += 2) {
			*CH (j - 1, k - 1, i - 2) = *C1 (j - 1, k - 1, i - 2) -
			    *C1 (jc - 1, k - 1, i - 1);
			*CH (jc - 1, k - 1, i - 2) = *C1 (j - 1, k - 1, i - 2) +
			    *C1 (jc - 1, k - 1, i - 1);
			*CH (j - 1, k - 1, i - 1) = *C1 (j - 1, k - 1, i - 1) +
			    *C1 (jc - 1, k - 1, i - 2);
			*CH (jc - 1, k - 1, i - 1) = *C1 (j - 1, k - 1, i - 1) -
			    *C1 (jc - 1, k - 1, i - 2);
		    }
		}
	    }
	}
	else {
	    for (j = 2; j <= ipph; j++) {
		jc = ipp2 - j;
		for (i = 3; i <= sido; i += 2) {
		    for (k = 1; k <= sl1; k++) {
			*CH (j - 1, k - 1, i - 2) = *C1 (j - 1, k - 1, i - 2) -
			    *C1 (jc - 1, k - 1, i - 1);
			*CH (jc - 1, k - 1, i - 2) = *C1 (j - 1, k - 1, i - 2) +
			    *C1 (jc - 1, k - 1, i - 1);
			*CH (j - 1, k - 1, i - 1) = *C1 (j - 1, k - 1, i - 1) +
			    *C1 (jc - 1, k - 1, i - 2);
			*CH (jc - 1, k - 1, i - 1) = *C1 (j - 1, k - 1, i - 1) -
			    *C1 (jc - 1, k - 1, i - 2);
		    }
		}
	    }
	}

	*C2 (0, 0) = *CH2 (0, 0);
	if (sidl1 != 1) {
	    for (i = 2; i <= sidl1; i++) {
		*C2 (0, i - 1) = *CH2 (0, i - 1);
	    }
	}
	for (j = 2; j <= sip; j++) {
	    for (k = 1; k <= sl1; k++) {
		*C1 (j - 1, k - 1, 0) = *CH (j - 1, k - 1, 0);
	    }
	}
	if (nbd < sl1) {
	    is = -sido;
	    for (j = 1; j <= (sip - 1); j++) {
		for (i = 1; i <= ((sido - 1) / 2); i++) {
		    for (k = 1; k <= sl1; k++) {
			*C1 (j, k - 1, i * 2 - 1) = wa[is + sido * j + (i - 1) * 2] *
			    *CH (j, k - 1, i * 2 - 1) - wa[is + sido * j + i * 2 - 1] *
			    *CH (j, k - 1, (i - 1) * 2 + 2);
			*C1 (j, k - 1, (i - 1) * 2 + 2) = wa[is + sido * j + (i - 1) * 2] *
			    *CH (j, k - 1, (i - 1) * 2 + 2) + wa[is + sido * j + i * 2 - 1] *
			    *CH (j, k - 1, i * 2 - 1);
		    }
		}
	    }
	}
	else {
	    is = -sido;
	    for (j = 1; j <= (sip - 1); j++) {
		for (k = 1; k <= sl1; k++) {
		    for (i = 1; i <= ((sido - 1) / 2); i++) {
			*C1 (j, k - 1, i * 2 - 1) = wa[is + sido * j + i * 2 - 2] *
			    *CH (j, k - 1, i * 2 - 1) - wa[is + sido * j + i * 2 - 1] *
			    *CH (j, k - 1, i * 2);
			*C1 (j, k - 1, i * 2) = wa[is + sido * j + i * 2 - 2] *
			    *CH (j, k - 1, i * 2) + wa[is + sido * j + i * 2 - 1] *
			    *CH (j, k - 1, i * 2 - 1);
		    }
		}
	    }
	}
    }
    return;
}				/* end of function */

#undef CC
#undef C1
#undef C2
#undef CH
#undef CH2
