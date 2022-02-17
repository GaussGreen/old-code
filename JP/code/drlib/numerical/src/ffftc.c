#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_fft_complex (Mint n, Mf_complex *seq, va_list argptr);
void imsl_f2tcf (Mint *n, Mf_complex seq[], Mf_complex coef[], Mfloat wfftc[], Mfloat cpy[]);
void imsl_f3tcf (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[], Mfloat
                imsl_fac[]);
void imsl_f5tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[]);
void imsl_f6tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[]);
void imsl_f7tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[], Mfloat wa4[]);
void imsl_f8tcf (Mint *nac, Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat wa[]);
void imsl_f2tcb (Mint *n, Mf_complex coef[], Mf_complex seq[], Mfloat wfftc[], Mfloat cpy[]);
void imsl_f3tcb (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[], Mfloat
                imsl_fac[]);
void imsl_f5tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat wa1[],
                Mfloat wa2[]);
void imsl_f6tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat wa1[], Mfloat wa2[], Mfloat wa3[]);
void imsl_f7tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[], Mfloat wa4[]);
void imsl_f8tcb (Mint *nac, Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat cc[], Mfloat c1[], Mfloat c2[], Mfloat ch[],
                Mfloat ch2[], Mfloat wa[]);
#else
static VA_LIST_HACK l_fft_complex ();
void imsl_f2tcf ();
void        imsl_f3tcf ();
void        imsl_f5tcf ();
void        imsl_f6tcf ();
void        imsl_f7tcf ();
void        imsl_f8tcf ();
void imsl_f2tcb ();
void        imsl_f3tcb ();
void        imsl_f5tcb ();
void        imsl_f6tcb ();
void        imsl_f7tcb ();
void        imsl_f8tcb ();
#endif

static Mf_complex *lv_coef;
#ifdef ANSI
Mf_complex *imsl_c_fft_complex (Mint n, Mf_complex *seq,...)
#else
Mf_complex *imsl_c_fft_complex (n, seq, va_alist)
    Mint        n;
    Mf_complex *seq;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, seq);
    E1PSH ("imsl_c_fft_complex", "imsl_z_fft_complex");
    lv_coef = NULL;
    IMSL_CALL (l_fft_complex (n, seq, argptr));
    va_end (argptr);
    E1POP ("imsl_c_fft_complex", "imsl_z_fft_complex");
    return lv_coef;
}


#ifdef ANSI
static VA_LIST_HACK l_fft_complex (Mint n, Mf_complex *seq, va_list argptr)
#else
static VA_LIST_HACK l_fft_complex (n, seq, argptr)
    Mint        n;
    Mf_complex *seq;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 2;
    Mint        forward = 1;
    Mint        backward = 0;
    Mint        user_solve = 0;
    Mint        supplied_params = 0;
    Mfloat     *params = NULL;
    Mfloat     *work = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_coef = va_arg (argptr, Mf_complex *);
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

    if (imsl_n1rty (0))
	goto RETURN;

    if (!supplied_params)
	params = imsl_c_fft_complex_init (n);
    if (imsl_n1rty (0))
	goto RETURN;

    work = (Mfloat *) imsl_malloc (2 * n * sizeof (*work));

    if (!user_solve) {
	lv_coef = (Mf_complex *) imsl_malloc (n * sizeof (*lv_coef));
    }

    if (lv_coef == NULL || work == NULL) {
	/* Not enough memory, with %(L1) = %(I1). */
	imsl_e1sti (1, n);
	imsl_e1stl (1, "n");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }

    if (forward) {
	imsl_f2tcf (&n, seq, lv_coef, params, work);
    }

    if (backward) {
	imsl_f2tcb (&n, seq, lv_coef, params, work);
    }
FREE_SPACE:
    if (!supplied_params && (params != NULL)) {
	imsl_free (params);
    }
    if (work != NULL) {
	imsl_free (work);
    }
RETURN:
    if (imsl_n1rty (0) > 3) {
	if (!user_solve) {
	    if (lv_coef != NULL)
		imsl_free (lv_coef);
	}
	lv_coef = NULL;
    }
    return (argptr);
}



/* Structured by FOR_STRUCT, v0.2, on 08/30/90 at 09:32:04
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2TCF/DF2TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F2TCF (N, SEQ, COEF, WFFTC, CPY)

    Arguments:  (See FFTCF)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f2tcf (Mint *n, Mf_complex seq[], Mf_complex coef[], Mfloat wfftc[], Mfloat cpy[])
#else
void imsl_f2tcf (n, seq, coef, wfftc, cpy)
    Mint       *n;
    Mf_complex  seq[], coef[];
    Mfloat      wfftc[], cpy[];
#endif
{
    Mint        i;
    /* CHECK ARGUMENT N */
    if (*n < 1) {
	imsl_e1psh ("l_f2tcf");
	imsl_e1sti (1, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	imsl_e1pop ("l_f2tcf");
	goto L_9000;
    }
    /* COPY SEQ TO CPY */
    for (i = 1; i <= *n; i++) {
	cpy[i * 2 - 2] = imsl_fc_convert (seq[i - 1]);
	cpy[i * 2 - 1] = imsl_c_aimag (seq[i - 1]);
    }

    if (*n > 1) {
	imsl_f3tcf (n, cpy, wfftc, &wfftc[*n * 2], &wfftc[*n * 4]);
    }
    for (i = 1; i <= *n; i++) {
	coef[i - 1] = imsl_cf_convert (cpy[i * 2 - 2], cpy[i * 2 - 1]);
    }

L_9000:
    ;
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F3TCF/DF3TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 22, 1990

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F3TCF (N, C, CH, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       C      - Real array of length 2*N, on input containing SEQ,
                on output containing the Fourier coefficients.
                (Input/Output)
       CH     - Real array of length 2*N, needed as workspace. Will
                contain the Fourier coefficients, which get copied
                into C.
       WA     - Real array of length 2*N from FFTCI containing powers of
                e**(2*PI*i/N).  (Input)
       FAC    - Real array of length 15 from FFTCI containing N (the
                length of the sequence), NF (the number of factors of
                N), and the imsl_prime factors of N.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f3tcf (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[], Mfloat imsl_fac[])
#else
void        imsl_f3tcf (n, c, ch, wa, imsl_fac)
    Mint       *n;
    Mfloat      c[], ch[], wa[], imsl_fac[];
#endif
{
    Mint        i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2,
                na, nac, nf;
    Mfloat      ti2, tr2;


    nf = nint (imsl_fac[1]);
    na = 0;
    l1 = 1;
    iw = 1;
    for (k1 = 1; k1 <= nf; k1++) {
	ip = nint (imsl_fac[k1 + 1]);
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip == 4) {
	    /* Case for N divisible by 4 */
	    ix2 = iw + idot;
	    ix3 = ix2 + idot;
	    if (na != 0) {
		imsl_f6tcf (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	    else {
		imsl_f6tcf (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	    na = 1 - na;
	}
	else if (ip == 2) {
	    /* Case for N divisible by 2 */
	    if (na != 0) {
		if (idot <= 2) {
		    c[0] = ch[0] + ch[idot];
		    c[idot] = ch[0] - ch[idot];
		    c[1] = ch[1] + ch[idot * 2 - 1];
		    c[idot * 2 - 1] = ch[1] - ch[idot * 2 - 1];
		}
		else {

		    /*
		     * **** Vector directive location ****
		     */
		    for (i = 2; i <= idot; i += 2) {
			c[i - 2] = ch[i - 2] + ch[idot + i - 2];
			tr2 = ch[i - 2] - ch[idot + i - 2];
			c[i - 1] = ch[i - 1] + ch[idot + i - 1];
			ti2 = ch[i - 1] - ch[idot + i - 1];
			c[i + idot - 1] = wa[iw + i - 3] * ti2 - wa[iw + i - 2] *
			    tr2;
			c[i + idot - 2] = wa[iw + i - 3] * tr2 + wa[iw + i - 2] *
			    ti2;
		    }
		}
	    }
	    else {
		if (ido <= 2) {
		    ch[0] = c[0] + c[idot];
		    ch[idot] = c[0] - c[idot];
		    ch[1] = c[1] + c[idot * 2 - 1];
		    ch[idot * 2 - 1] = c[1] - c[idot * 2 - 1];
		}
		else {

		    /*
		     * **** Vector directive location ****
		     */
		    for (i = 2; i <= idot; i += 2) {
			ch[i - 2] = c[i - 2] + c[idot + i - 2];
			tr2 = c[i - 2] - c[idot + i - 2];
			ch[i - 1] = c[i - 1] + c[idot + i - 1];
			ti2 = c[i - 1] - c[idot + i - 1];
			ch[i + idot - 1] = wa[iw + i - 3] * ti2 - wa[iw + i - 2] *
			    tr2;
			ch[i + idot - 2] = wa[iw + i - 3] * tr2 + wa[iw + i - 2] *
			    ti2;
		    }
		}
	    }
	    na = 1 - na;
	}
	else if (ip == 3) {
	    /* Case for N divisible by 3 */
	    ix2 = iw + idot;
	    if (na != 0) {
		imsl_f5tcf (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	    else {
		imsl_f5tcf (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	    na = 1 - na;
	}
	else if (ip == 5) {
	    /* Case for N divisible by 5 */
	    ix2 = iw + idot;
	    ix3 = ix2 + idot;
	    ix4 = ix3 + idot;
	    if (na != 0) {
		imsl_f7tcf (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	    else {
		imsl_f7tcf (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	    na = 1 - na;
	}
	else {
	    /* Case for remaining factors of N */
	    if (na != 0) {
		imsl_f8tcf (&nac, &idot, &ip, &l1, &idl1, ch, ch, ch, c,
		    c, &wa[iw - 1]);
	    }
	    else {
		imsl_f8tcf (&nac, &idot, &ip, &l1, &idl1, c, c, c, ch,
		    ch, &wa[iw - 1]);
	    }
	    if (nac != 0)
		na = 1 - na;
	}

	l1 = l2;
	iw += (ip - 1) * idot;
    }
    if (na != 0)
	scopy (2 ** n, ch, 1, c, 1);
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F5TCF/DF5TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 24, 1990

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F5TCF (IDO, L1, CC, CH, WA1, WA2)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCF.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing needed elements of array WA
                in F3TCF.  (Input)
       WA2    - Real array containing other needed elements of array
                WA in F3TCF.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f5tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[])
#else
void        imsl_f5tcf (ido, l1, cc, ch, wa1, wa2)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(3)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, ti2, tr2;
    static Mfloat taur = -0.5e0;
    static Mfloat taui = -0.866025403784438646763723170753e0;



    if (sido == 2) {

	/*
	 * **** Vector directive location ****
	 */
	for (k = 1; k <= sl1; k++) {
	    tr2 = *CC (k - 1, 1, 0) + *CC (k - 1, 2, 0);
	    cr2 = *CC (k - 1, 0, 0) + taur * tr2;
	    *CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2;
	    ti2 = *CC (k - 1, 1, 1) + *CC (k - 1, 2, 1);
	    ci2 = *CC (k - 1, 0, 1) + taur * ti2;
	    *CH (0, k - 1, 1) = *CC (k - 1, 0, 1) + ti2;
	    cr3 = taui * (*CC (k - 1, 1, 0) - *CC (k - 1, 2, 0));
	    ci3 = taui * (*CC (k - 1, 1, 1) - *CC (k - 1, 2, 1));
	    *CH (1, k - 1, 0) = cr2 - ci3;
	    *CH (2, k - 1, 0) = cr2 + ci3;
	    *CH (1, k - 1, 1) = ci2 + cr3;
	    *CH (2, k - 1, 1) = ci2 - cr3;
	}
    }
    else if ((sido - 1) / 2 > sl1) {
	for (k = 1; k <= sl1; k++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (i = 2; i <= sido; i += 2) {
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 2, i - 2);
		cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 2, i - 1);
		ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		cr3 = taui * (*CC (k - 1, 1, i - 2) - *CC (k - 1, 2, i - 2));
		ci3 = taui * (*CC (k - 1, 1, i - 1) - *CC (k - 1, 2, i - 1));
		dr2 = cr2 - ci3;
		dr3 = cr2 + ci3;
		di2 = ci2 + cr3;
		di3 = ci2 - cr3;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 - wa1[i - 1] * dr2;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 + wa1[i - 1] * di2;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 - wa2[i - 1] * dr3;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 + wa2[i - 1] * di3;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (k = 1; k <= sl1; k++) {
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 2, i - 2);
		cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 2, i - 1);
		ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		cr3 = taui * (*CC (k - 1, 1, i - 2) - *CC (k - 1, 2, i - 2));
		ci3 = taui * (*CC (k - 1, 1, i - 1) - *CC (k - 1, 2, i - 1));
		dr2 = cr2 - ci3;
		dr3 = cr2 + ci3;
		di2 = ci2 + cr3;
		di3 = ci2 - cr3;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 - wa1[i - 1] * dr2;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 + wa1[i - 1] * di2;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 - wa2[i - 1] * dr3;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 + wa2[i - 1] * di3;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F6TCF/DF6TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 24, 1990

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F6TCF (IDO, L1, CC, CH, WA1, WA2, WA3)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCF.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing the needed elements of array WA
                in F3TCF.  (Input)
       WA2    - Real array containing other needed elements of array
                WA in F3TCF.  (Input)
       WA3    - Real array containing other needed elements of array
                WA in F3TCF.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f6tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[])
#else
void        imsl_f6tcf (ido, l1, cc, ch, wa1, wa2, wa3)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(4)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1,
                tr2, tr3, tr4;


    if (sido == 2) {

	/*
	 * **** Vector directive location ****
	 */
	for (k = 1; k <= sl1; k++) {
	    ti1 = *CC (k - 1, 0, 1) - *CC (k - 1, 2, 1);
	    ti2 = *CC (k - 1, 0, 1) + *CC (k - 1, 2, 1);
	    tr4 = *CC (k - 1, 1, 1) - *CC (k - 1, 3, 1);
	    ti3 = *CC (k - 1, 1, 1) + *CC (k - 1, 3, 1);
	    tr1 = *CC (k - 1, 0, 0) - *CC (k - 1, 2, 0);
	    tr2 = *CC (k - 1, 0, 0) + *CC (k - 1, 2, 0);
	    ti4 = *CC (k - 1, 3, 0) - *CC (k - 1, 1, 0);
	    tr3 = *CC (k - 1, 1, 0) + *CC (k - 1, 3, 0);
	    *CH (0, k - 1, 0) = tr2 + tr3;
	    *CH (2, k - 1, 0) = tr2 - tr3;
	    *CH (0, k - 1, 1) = ti2 + ti3;
	    *CH (2, k - 1, 1) = ti2 - ti3;
	    *CH (1, k - 1, 0) = tr1 + tr4;
	    *CH (3, k - 1, 0) = tr1 - tr4;
	    *CH (1, k - 1, 1) = ti1 + ti4;
	    *CH (3, k - 1, 1) = ti1 - ti4;
	}
    }
    else if ((sido - 1) / 2 > sl1) {
	for (k = 1; k <= sl1; k++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (i = 2; i <= sido; i += 2) {
		ti1 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 2, i - 1);
		ti2 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 2, i - 1);
		ti3 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 3, i - 1);
		tr4 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 3, i - 1);
		tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 2, i - 2);
		tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 2, i - 2);
		ti4 = *CC (k - 1, 3, i - 2) - *CC (k - 1, 1, i - 2);
		tr3 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = tr2 + tr3;
		cr3 = tr2 - tr3;
		*CH (0, k - 1, i - 1) = ti2 + ti3;
		ci3 = ti2 - ti3;
		cr2 = tr1 + tr4;
		cr4 = tr1 - tr4;
		ci2 = ti1 + ti4;
		ci4 = ti1 - ti4;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * cr2 + wa1[i - 1] * ci2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * ci2 - wa1[i - 1] * cr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * cr3 + wa2[i - 1] * ci3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * ci3 - wa2[i - 1] * cr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * cr4 + wa3[i - 1] * ci4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * ci4 - wa3[i - 1] * cr4;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (k = 1; k <= sl1; k++) {
		ti1 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 2, i - 1);
		ti2 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 2, i - 1);
		ti3 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 3, i - 1);
		tr4 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 3, i - 1);
		tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 2, i - 2);
		tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 2, i - 2);
		ti4 = *CC (k - 1, 3, i - 2) - *CC (k - 1, 1, i - 2);
		tr3 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = tr2 + tr3;
		cr3 = tr2 - tr3;
		*CH (0, k - 1, i - 1) = ti2 + ti3;
		ci3 = ti2 - ti3;
		cr2 = tr1 + tr4;
		cr4 = tr1 - tr4;
		ci2 = ti1 + ti4;
		ci4 = ti1 - ti4;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * cr2 + wa1[i - 1] * ci2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * ci2 - wa1[i - 1] * cr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * cr3 + wa2[i - 1] * ci3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * ci3 - wa2[i - 1] * cr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * cr4 + wa3[i - 1] * ci4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * ci4 - wa3[i - 1] * cr4;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F7TCF/DF7TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 24, 1990

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F7TCF (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCF.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing the needed elements of array WA
                om F3TCF.  (Input)
       WA2    - Real array containing other needed elements of array WA
                om F3TCF.  (Input)
       WA3    - Real array containing other needed elements of array WA
                om F3TCF.  (Input)
       WA4    - Real array containing other needed elements of array WA
                om F3TCF.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f7tcf (Mint *ido, Mint *l1, Mfloat cc[], Mfloat ch[], Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[], Mfloat wa4[])
#else
void        imsl_f7tcf (ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[], wa4[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(5)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3, di4,
                di5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3,
                tr4, tr5;
    static Mfloat tr11 = 0.309016994374947424102293417183e0;
    static Mfloat ti11 = -0.951056516295153572116439333379e0;
    static Mfloat tr12 = -0.809016994374947424102293417183e0;
    static Mfloat ti12 = -0.587785252292473129168705954639e0;



    if (sido == 2) {

	/*
	 * **** Vector directive location ****
	 */
	for (k = 1; k <= sl1; k++) {
	    ti5 = *CC (k - 1, 1, 1) - *CC (k - 1, 4, 1);
	    ti2 = *CC (k - 1, 1, 1) + *CC (k - 1, 4, 1);
	    ti4 = *CC (k - 1, 2, 1) - *CC (k - 1, 3, 1);
	    ti3 = *CC (k - 1, 2, 1) + *CC (k - 1, 3, 1);
	    tr5 = *CC (k - 1, 1, 0) - *CC (k - 1, 4, 0);
	    tr2 = *CC (k - 1, 1, 0) + *CC (k - 1, 4, 0);
	    tr4 = *CC (k - 1, 2, 0) - *CC (k - 1, 3, 0);
	    tr3 = *CC (k - 1, 2, 0) + *CC (k - 1, 3, 0);
	    *CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2 + tr3;
	    *CH (0, k - 1, 1) = *CC (k - 1, 0, 1) + ti2 + ti3;
	    cr2 = *CC (k - 1, 0, 0) + tr11 * tr2 + tr12 * tr3;
	    ci2 = *CC (k - 1, 0, 1) + tr11 * ti2 + tr12 * ti3;
	    cr3 = *CC (k - 1, 0, 0) + tr12 * tr2 + tr11 * tr3;
	    ci3 = *CC (k - 1, 0, 1) + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    *CH (1, k - 1, 0) = cr2 - ci5;
	    *CH (4, k - 1, 0) = cr2 + ci5;
	    *CH (1, k - 1, 1) = ci2 + cr5;
	    *CH (2, k - 1, 1) = ci3 + cr4;
	    *CH (2, k - 1, 0) = cr3 - ci4;
	    *CH (3, k - 1, 0) = cr3 + ci4;
	    *CH (3, k - 1, 1) = ci3 - cr4;
	    *CH (4, k - 1, 1) = ci2 - cr5;
	}
    }
    else if ((sido - 1) / 2 > sl1) {
	for (k = 1; k <= sl1; k++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (i = 2; i <= sido; i += 2) {
		ti5 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 4, i - 1);
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 4, i - 1);
		ti4 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 3, i - 1);
		ti3 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 3, i - 1);
		tr5 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 4, i - 2);
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 4, i - 2);
		tr4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 + tr3;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 + ti3;
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
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 + wa1[i - 1] * di2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 - wa1[i - 1] * dr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 + wa2[i - 1] * di3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 - wa2[i - 1] * dr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * dr4 + wa3[i - 1] * di4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * di4 - wa3[i - 1] * dr4;
		*CH (4, k - 1, i - 2) = wa4[i - 2] * dr5 + wa4[i - 1] * di5;
		*CH (4, k - 1, i - 1) = wa4[i - 2] * di5 - wa4[i - 1] * dr5;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (k = 1; k <= sl1; k++) {
		ti5 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 4, i - 1);
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 4, i - 1);
		ti4 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 3, i - 1);
		ti3 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 3, i - 1);
		tr5 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 4, i - 2);
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 4, i - 2);
		tr4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 + tr3;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 + ti3;
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
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 + wa1[i - 1] * di2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 - wa1[i - 1] * dr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 + wa2[i - 1] * di3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 - wa2[i - 1] * dr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * dr4 + wa3[i - 1] * di4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * di4 - wa3[i - 1] * dr4;
		*CH (4, k - 1, i - 2) = wa4[i - 2] * dr5 + wa4[i - 1] * di5;
		*CH (4, k - 1, i - 1) = wa4[i - 2] * di5 - wa4[i - 1] * dr5;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F8TCF/DF8TCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 15, 1990

    Purpose:    Compute the Fourier coefficients of a d_complex periodic
                sequence.

    Usage:      CALL F8TCF (NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH,
                            CH2, WA)

    Arguments:
       NAC    - Integer set to 1 if IDO=2, set to 0 otherwise.  (Output)
       IDO    - 2*N divided by present and previous factors.  (Input)
       IP     - The factor now being used.  (Input)
       L1     - Product of previous factors used.  (Input)
       IDL1   - ID0*L1
       CC     - Matrix containing C or CH from F3TCF.  (Input/Output)
       C1     - Matrix containing C or CH from F3TCF.  (Input/Output)
       C2     - Matrix containing C or CH from F3TCF.  (Input/Output)
       CH     - Matrix that contains a better estimate for
                COEF.  (Output)
       CH2    - Matrix that contains a better estimate for
                COEF.  (Output)  Note:  The roles of CC and CH are
                reversed when NA=1.
       WA     - Real array containing the needed elements of array WA
                in F3TCF.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f8tcf (Mint *nac, Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat *cc, Mfloat *c1, Mfloat *c2, Mfloat *ch,
                Mfloat *ch2, Mfloat wa[])
#else
void        imsl_f8tcf (nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2,
                wa)
    Mint       *nac, *ido, *ip, *l1, *idl1;
    Mfloat     *cc, *c1, *c2, *ch, *ch2, wa[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(*ip)*(sido)+(J_)*(sido)+(K_))
#define C1(I_,J_,K_)	(c1+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define C2(I_,J_)	(c2+(I_)*(sidl1)+(J_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH2(I_,J_)	(ch2+(I_)*(sidl1)+(J_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        sidl1 = *idl1;
    Mint        i, idij, idj, idl, idlj, idlj1u[512], idot, idp, ik, ipp2,
                ipph, j, j2s, j3s, j4s, j5s, jc, k, l;
    Mfloat      sumi, sumr, wai, war;


    idot = sido / 2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    idp = *ip ** ido;

    if (sido > sl1) {
	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;
	    for (k = 1; k <= sl1; k++) {

		/*
		 * **** Vector directive location ****
		 */
		for (i = 1; i <= sido; i++) {
		    *CH (j - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) +
			*CC (k - 1, jc - 1, i - 1);
		    *CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) -
			*CC (k - 1, jc - 1, i - 1);
		}
	    }
	}

	for (k = 1; k <= sl1; k++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (i = 1; i <= sido; i++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}
    }
    else {
	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;
	    for (i = 1; i <= sido; i++) {

		/*
		 * **** Vector directive location ****
		 */
		for (k = 1; k <= sl1; k++) {
		    *CH (j - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) +
			*CC (k - 1, jc - 1, i - 1);
		    *CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) -
			*CC (k - 1, jc - 1, i - 1);
		}
	    }
	}

	for (i = 1; i <= sido; i++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (k = 1; k <= sl1; k++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}

	;
    }

    idl = 2 - sido;
    if (sidl1 > 2) {
	for (l = 1; l <= (ipph - 1); l++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (ik = 1; ik <= sidl1; ik++) {
		*C2 (l, ik - 1) = *CH2 (0, ik - 1) + wa[idl + sido * l - 2] *
		    *CH2 (1, ik - 1);
		*C2 (ipp2 - l - 2, ik - 1) = -wa[idl + sido * l - 1] ** CH2 (*ip - 1, ik - 1);
	    }
	    idlj = idl + sido * l;

	    for (j = 1; j <= (ipph - 2); j++) {
		idlj += sido * l;
		if (idlj > idp)
		    idlj -= idp;
		wai = -wa[idlj - 1];
		war = wa[idlj - 2];

		/*
		 * **** Vector directive location ****
		 */
		for (ik = 1; ik <= sidl1; ik++) {
		    *C2 (l, ik - 1) += war ** CH2 (j + 1, ik - 1);
		    *C2 (ipp2 - l - 2, ik - 1) += wai ** CH2 (ipp2 - j - 3, ik - 1);
		}
	    }
	}

	for (j = 2; j <= ipph; j++) {

	    /*
	     * **** Vector directive location ****
	     */
	    for (ik = 1; ik <= sidl1; ik++) {
		*CH2 (0, ik - 1) += *CH2 (j - 1, ik - 1);
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;

	    /*
	     * **** Vector directive location ****
	     */
	    for (ik = 2; ik <= sidl1; ik += 2) {
		*CH2 (j - 1, ik - 2) = *C2 (j - 1, ik - 2) - *C2 (jc - 1, ik - 1);
		*CH2 (jc - 1, ik - 2) = *C2 (j - 1, ik - 2) + *C2 (jc - 1, ik - 1);
		*CH2 (j - 1, ik - 1) = *C2 (j - 1, ik - 1) + *C2 (jc - 1, ik - 2);
		*CH2 (jc - 1, ik - 1) = *C2 (j - 1, ik - 1) - *C2 (jc - 1, ik - 2);
	    }
	}
    }
    else {
	for (j4s = 0; j4s <= (ipph - 2); j4s += 511) {
	    j5s = imsl_i_min (ipph - 1 - j4s, 511);
	    /*
	     * The local integer vector idlj1u was inserted by fpp on the
	     * Cray to imsl_aid in the vectorization of the code. DIR$
	     * IVDEP
	     */
	    for (l = 1; l <= j5s; l++) {
		*C2 (j4s + l, 0) = *CH2 (0, 0) + wa[idl + (j4s + l) ** ido - 2] *
		    *CH2 (1, 0);
		*C2 (ipp2 - j4s - l - 2, 0) = -wa[idl + (j4s + l) ** ido - 1] *
		    *CH2 (*ip - 1, 0);
		*C2 (j4s + l, 1) = *CH2 (0, 1) + wa[idl + (j4s + l) ** ido - 2] *
		    *CH2 (1, 1);
		*C2 (ipp2 - j4s - l - 2, 1) = -wa[idl + (j4s + l) ** ido - 1] *
		    *CH2 (*ip - 1, 1);
		idlj1u[l - 1] = idl + sido * (l + j4s);
	    }

	    for (j = 1; j <= (ipph - 2); j++) {
		/* DIR$          IVDEP */
		for (l = 1; l <= j5s; l++) {
		    idlj1u[l - 1] += sido * (l + j4s);
		    if (idlj1u[l - 1] > idp)
			idlj1u[l - 1] -= idp;
		    j2s = idlj1u[l - 1];
		    j3s = idlj1u[l - 1];
		    war = wa[j3s - 2];
		    wai = -wa[j2s - 1];

		    /*
		     * **** Vector directive location ****
		     */
		    *C2 (j4s + l, 0) += war ** CH2 (j + 1, 0);
		    *C2 (ipp2 - j4s - l - 2, 0) += wai ** CH2 (ipp2 - (j - 1) - 4, 0);
		    *C2 (j4s + l, 1) += war ** CH2 (j + 1, 1);
		    *C2 (ipp2 - j4s - l - 2, 1) += wai ** CH2 (ipp2 - (j - 1) - 4, 1);
		}
	    }
	}

	sumr = F_ZERO;
	sumi = F_ZERO;
	for (j = 1; j <= ipph; j++) {
	    sumr += *CH2 (j - 1, 0);
	    sumi += *CH2 (j - 1, 1);
	}
	*CH2 (0, 0) = sumr;
	*CH2 (0, 1) = sumi;

	/* DIR$    IVDEP */
	for (j = 1; j <= (ipph - 1); j++) {
	    *CH2 (j, 0) = *C2 (j, 0) - *C2 (ipp2 - j - 2, 1);
	    *CH2 (ipp2 - j - 2, 0) = *C2 (j, 0) + *C2 (ipp2 - j - 2, 1);
	    *CH2 (j, 1) = *C2 (j, 1) + *C2 (ipp2 - j - 2, 0);
	    *CH2 (ipp2 - j - 2, 1) = *C2 (j, 1) - *C2 (ipp2 - j - 2, 0);
	}
    }
    *nac = 1;
    if (sido == 2)
	goto L_350;
    *nac = 0;

    /*
     * **** Vector directive location ****
     */
    for (ik = 1; ik <= sidl1; ik++) {
	*C2 (0, ik - 1) = *CH2 (0, ik - 1);
    }

    for (j = 2; j <= *ip; j++) {

	/*
	 * **** Vector directive location ****
	 */
	for (k = 1; k <= sl1; k++) {
	    *C1 (j - 1, k - 1, 0) = *CH (j - 1, k - 1, 0);
	    *C1 (j - 1, k - 1, 1) = *CH (j - 1, k - 1, 1);
	}
    }

    if (idot < sl1) {
	idij = 0;
	for (j = 2; j <= *ip; j++) {
	    idij += 2;
	    for (i = 0; i <= ((sido - 2) / 2 - 1); i++) {

		/*
		 * **** Vector directive location ****
		 */
		for (k = 1; k <= sl1; k++) {
		    *C1 (j - 1, k - 1, i * 2 + 2) = wa[idij + i * 2] ** CH (j - 1, k - 1, i * 2 + 2) +
			wa[idij + (i + 1) * 2 - 1] ** CH (j - 1, k - 1, i * 2 + 3);
		    *C1 (j - 1, k - 1, i * 2 + 3) = wa[idij + i * 2] ** CH (j - 1, k - 1, i * 2 + 3) -
			wa[idij + (i + 1) * 2 - 1] ** CH (j - 1, k - 1, i * 2 + 2);
		}
	    }
	    if ((sido - 2) / 2 > 0)
		idij += (sido - 2) / 2 * 2;
	}
    }
    else {
	idj = 2 - sido;
	for (j = 1; j <= (*ip - 1); j++) {
	    idl = idj + sido * j;
	    for (k = 1; k <= sl1; k++) {

		/*
		 * **** Vector directive location ****
		 */
		for (i = 1; i <= ((sido - 2) / 2); i++) {
		    *C1 (j, k - 1, i * 2) = wa[i * 2 + idl - 2] ** CH (j, k - 1, i * 2) +
			wa[i * 2 + idl - 1] ** CH (j, k - 1, i * 2 + 1);
		    *C1 (j, k - 1, i * 2 + 1) = wa[i * 2 + idl - 2] ** CH (j, k - 1, i * 2 + 1) -
			wa[i * 2 + idl - 1] ** CH (j, k - 1, i * 2);
		}
	    }
	}
    }
L_350:
    return;
}				/* end of function */
#undef CC
#undef C1
#undef C2
#undef CH
#undef CH2
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F2TCB/DF2TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F2TCB (N, COEF, SEQ, WFFTC, CPY)

    Arguments:  (See FFTCB)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f2tcb (Mint *n, Mf_complex coef[], Mf_complex seq[], Mfloat wfftc[], Mfloat cpy[])
#else
void imsl_f2tcb (n, coef, seq, wfftc, cpy)
    Mint       *n;
    Mf_complex  coef[], seq[];
    Mfloat      wfftc[], cpy[];
#endif
{
    Mint        i;
    /* CHECK ARGUMENT N */
    if (*n < 1) {
	imsl_e1psh ("l_f2tcb");
	imsl_e1sti (1, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	imsl_e1pop ("l_f2tcb");
	goto L_9000;
    }
    /* COPY COEF TO CPY */
    for (i = 1; i <= *n; i++) {
	cpy[i * 2 - 2] = imsl_fc_convert (coef[i - 1]);
	cpy[i * 2 - 1] = imsl_c_aimag (coef[i - 1]);
    }

    if (*n > 1) {
	imsl_f3tcb (n, cpy, wfftc, &wfftc[*n * 2], &wfftc[*n * 4]);
    }
    /* COPY CPY TO SEQ */
    for (i = 1; i <= *n; i++) {
	seq[i - 1] = imsl_cf_convert (cpy[i * 2 - 2], cpy[i * 2 - 1]);
    }

L_9000:
    ;
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F3TCB/DF3TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 21, 1990

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F3TCB (N, C, CH, WA, FAC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       C      - Real array of length 2*N, on input containing SEQ,
                on output containing the Fourier coefficients.
                (Input/Output)
       CH     - Real array of length 2*N, needed as workspace.  Will
                contain the Fourier coefficients, which get copied
                into C.  (Workspace)
       WA     - Real array of length 2*N from FFTCI containing powers of
                e**(2*PI*i/N).  (Input)
       FAC    - Real array of length 15 from FFTCI containing N (the
                length of the sequence), NF (the number of factors of
                N), and the imsl_prime factors of N.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f3tcb (Mint *n, Mfloat c[], Mfloat ch[], Mfloat wa[], Mfloat imsl_fac[])
#else
void        imsl_f3tcb (n, c, ch, wa, imsl_fac)
    Mint       *n;
    Mfloat      c[], ch[], wa[], imsl_fac[];
#endif
{
    Mint        i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2,
                na, nac, nf;
    Mfloat      ti2, tr2;


    nf = nint (imsl_fac[1]);
    na = 0;
    l1 = 1;
    iw = 1;
    for (k1 = 1; k1 <= nf; k1++) {
	ip = nint (imsl_fac[k1 + 1]);
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;

	if (ip == 4) {
	    /* Case for N divisible by 4 */
	    ix2 = iw + idot;
	    ix3 = ix2 + idot;
	    if (na != 0) {
		imsl_f6tcb (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	    else {
		imsl_f6tcb (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1]);
	    }
	    na = 1 - na;
	}
	else if (ip == 2) {
	    /* Case for N divisible by 2 */
	    if (na != 0) {
		if (idot <= 2) {
		    c[0] = ch[0] + ch[idot];
		    c[idot] = ch[0] - ch[idot];
		    c[1] = ch[1] + ch[idot * 2 - 1];
		    c[idot * 2 - 1] = ch[1] - ch[idot * 2 - 1];
		}
		else {
		    for (i = 2; i <= idot; i += 2) {
			c[i - 2] = ch[i - 2] + ch[idot + i - 2];
			tr2 = ch[i - 2] - ch[idot + i - 2];
			c[i - 1] = ch[i - 1] + ch[idot + i - 1];
			ti2 = ch[i - 1] - ch[idot + i - 1];
			c[idot + i - 1] = wa[iw + i - 3] * ti2 + wa[iw + i - 2] *
			    tr2;
			c[idot + i - 2] = wa[iw + i - 3] * tr2 - wa[iw + i - 2] *
			    ti2;
		    }
		}
	    }
	    else {
		if (idot <= 2) {
		    ch[0] = c[0] + c[idot];
		    ch[idot] = c[0] - c[idot];
		    ch[1] = c[1] + c[idot * 2 - 1];
		    ch[idot * 2 - 1] = c[1] - c[idot * 2 - 1];
		}
		else {
		    for (i = 2; i <= idot; i += 2) {
			ch[i - 2] = c[i - 2] + c[idot + i - 2];
			tr2 = c[i - 2] - c[idot + i - 2];
			ch[i - 1] = c[i - 1] + c[idot + i - 1];
			ti2 = c[i - 1] - c[idot + i - 1];
			ch[idot + i - 1] = wa[iw + i - 3] * ti2 + wa[iw + i - 2] *
			    tr2;
			ch[idot + i - 2] = wa[iw + i - 3] * tr2 - wa[iw + i - 2] *
			    ti2;
		    }
		}
	    }
	    na = 1 - na;
	}
	else if (ip == 3) {
	    /* Case for N divisible by 3 */
	    ix2 = iw + idot;
	    if (na != 0) {
		imsl_f5tcb (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	    else {
		imsl_f5tcb (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1]);
	    }
	    na = 1 - na;
	}
	else if (ip == 5) {
	    /* Case for N divisible by 5 */
	    ix2 = iw + idot;
	    ix3 = ix2 + idot;
	    ix4 = ix3 + idot;
	    if (na != 0) {
		imsl_f7tcb (&idot, &l1, ch, c, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	    else {
		imsl_f7tcb (&idot, &l1, c, ch, &wa[iw - 1], &wa[ix2 - 1],
		    &wa[ix3 - 1], &wa[ix4 - 1]);
	    }
	    na = 1 - na;
	}
	else {
	    /* Case for remaining factors of N */
	    if (na != 0) {
		imsl_f8tcb (&nac, &idot, &ip, &l1, &idl1, ch, ch, ch, c,
		    c, &wa[iw - 1]);
	    }
	    else {
		imsl_f8tcb (&nac, &idot, &ip, &l1, &idl1, c, c, c, ch,
		    ch, &wa[iw - 1]);
	    }
	    if (nac != 0)
		na = 1 - na;
	}
	l1 = l2;
	iw += (ip - 1) * idot;
    }
    if (na != 0)
	scopy (2 ** n, ch, 1, c, 1);
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F5TCB/DF5TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 14, 1990

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F5TCB (IDO, L1, CC, CH, WA1, WA2)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCB.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing needed elements of array WA
                in F3TCB.  (Input)
       WA2    - Real array containing other needed elements of array
                WA in F3TCB.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f5tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat wa1[],
                Mfloat wa2[])
#else
void        imsl_f5tcb (ido, l1, cc, ch, wa1, wa2)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(3)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, ti2, tr2;
    static Mfloat taur = -0.5e0;
    static Mfloat taui = 0.866025403784438646763723170753e0;



    if (sido == 2) {
	for (k = 1; k <= sl1; k++) {
	    tr2 = *CC (k - 1, 1, 0) + *CC (k - 1, 2, 0);
	    cr2 = *CC (k - 1, 0, 0) + taur * tr2;
	    *CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2;
	    ti2 = *CC (k - 1, 1, 1) + *CC (k - 1, 2, 1);
	    ci2 = *CC (k - 1, 0, 1) + taur * ti2;
	    *CH (0, k - 1, 1) = *CC (k - 1, 0, 1) + ti2;
	    cr3 = taui * (*CC (k - 1, 1, 0) - *CC (k - 1, 2, 0));
	    ci3 = taui * (*CC (k - 1, 1, 1) - *CC (k - 1, 2, 1));
	    *CH (1, k - 1, 0) = cr2 - ci3;
	    *CH (2, k - 1, 0) = cr2 + ci3;
	    *CH (1, k - 1, 1) = ci2 + cr3;
	    *CH (2, k - 1, 1) = ci2 - cr3;
	}
    }
    else if (sl1 <= sido / 2) {
	for (k = 1; k <= sl1; k++) {
	    for (i = 2; i <= sido; i += 2) {
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 2, i - 2);
		cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 2, i - 1);
		ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		cr3 = taui * (*CC (k - 1, 1, i - 2) - *CC (k - 1, 2, i - 2));
		ci3 = taui * (*CC (k - 1, 1, i - 1) - *CC (k - 1, 2, i - 1));
		dr2 = cr2 - ci3;
		dr3 = cr2 + ci3;
		di2 = ci2 + cr3;
		di3 = ci2 - cr3;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {
	    for (k = 1; k <= sl1; k++) {
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 2, i - 2);
		cr2 = *CC (k - 1, 0, i - 2) + taur * tr2;
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2;
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 2, i - 1);
		ci2 = *CC (k - 1, 0, i - 1) + taur * ti2;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2;
		cr3 = taui * (*CC (k - 1, 1, i - 2) - *CC (k - 1, 2, i - 2));
		ci3 = taui * (*CC (k - 1, 1, i - 1) - *CC (k - 1, 2, i - 1));
		dr2 = cr2 - ci3;
		dr3 = cr2 + ci3;
		di2 = ci2 + cr3;
		di3 = ci2 - cr3;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/* Structured by FOR_STRUCT, v0.2, on 08/30/90 at 10:21:26
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F6TCB/DF6TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 14, 1990

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F6TCB (IDO, L1, CC, CH, WA1, WA2, WA3)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCB.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing the needed elements of array WA
                in F3TCB.  (Input)
       WA2    - Real array containing other needed elements of array
                WA in F3TCB.  (Input)
       WA3    - Real array containing other needed elements of array
                WA in F3TCB.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f6tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch,
                Mfloat wa1[], Mfloat wa2[], Mfloat wa3[])
#else
void        imsl_f6tcb (ido, l1, cc, ch, wa1, wa2, wa3)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(4)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1,
                tr2, tr3, tr4;


    if (sido == 2) {
	for (k = 1; k <= sl1; k++) {
	    ti1 = *CC (k - 1, 0, 1) - *CC (k - 1, 2, 1);
	    ti2 = *CC (k - 1, 0, 1) + *CC (k - 1, 2, 1);
	    tr4 = *CC (k - 1, 3, 1) - *CC (k - 1, 1, 1);
	    ti3 = *CC (k - 1, 1, 1) + *CC (k - 1, 3, 1);
	    tr1 = *CC (k - 1, 0, 0) - *CC (k - 1, 2, 0);
	    tr2 = *CC (k - 1, 0, 0) + *CC (k - 1, 2, 0);
	    ti4 = *CC (k - 1, 1, 0) - *CC (k - 1, 3, 0);
	    tr3 = *CC (k - 1, 1, 0) + *CC (k - 1, 3, 0);
	    *CH (0, k - 1, 0) = tr2 + tr3;
	    *CH (2, k - 1, 0) = tr2 - tr3;
	    *CH (0, k - 1, 1) = ti2 + ti3;
	    *CH (2, k - 1, 1) = ti2 - ti3;
	    *CH (1, k - 1, 0) = tr1 + tr4;
	    *CH (3, k - 1, 0) = tr1 - tr4;
	    *CH (1, k - 1, 1) = ti1 + ti4;
	    *CH (3, k - 1, 1) = ti1 - ti4;
	}
    }
    else if (sl1 <= sido / 2) {
	for (k = 1; k <= sl1; k++) {
	    for (i = 2; i <= sido; i += 2) {
		ti1 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 2, i - 1);
		ti2 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 2, i - 1);
		ti3 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 3, i - 1);
		tr4 = *CC (k - 1, 3, i - 1) - *CC (k - 1, 1, i - 1);
		tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 2, i - 2);
		tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 2, i - 2);
		ti4 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = tr2 + tr3;
		cr3 = tr2 - tr3;
		*CH (0, k - 1, i - 1) = ti2 + ti3;
		ci3 = ti2 - ti3;
		cr2 = tr1 + tr4;
		cr4 = tr1 - tr4;
		ci2 = ti1 + ti4;
		ci4 = ti1 - ti4;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * cr2 - wa1[i - 1] * ci2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * ci2 + wa1[i - 1] * cr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * cr3 - wa2[i - 1] * ci3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * ci3 + wa2[i - 1] * cr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * cr4 - wa3[i - 1] * ci4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * ci4 + wa3[i - 1] * cr4;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {
	    for (k = 1; k <= sl1; k++) {
		ti1 = *CC (k - 1, 0, i - 1) - *CC (k - 1, 2, i - 1);
		ti2 = *CC (k - 1, 0, i - 1) + *CC (k - 1, 2, i - 1);
		ti3 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 3, i - 1);
		tr4 = *CC (k - 1, 3, i - 1) - *CC (k - 1, 1, i - 1);
		tr1 = *CC (k - 1, 0, i - 2) - *CC (k - 1, 2, i - 2);
		tr2 = *CC (k - 1, 0, i - 2) + *CC (k - 1, 2, i - 2);
		ti4 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = tr2 + tr3;
		cr3 = tr2 - tr3;
		*CH (0, k - 1, i - 1) = ti2 + ti3;
		ci3 = ti2 - ti3;
		cr2 = tr1 + tr4;
		cr4 = tr1 - tr4;
		ci2 = ti1 + ti4;
		ci4 = ti1 - ti4;
		*CH (1, k - 1, i - 2) = wa1[i - 2] * cr2 - wa1[i - 1] * ci2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * ci2 + wa1[i - 1] * cr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * cr3 - wa2[i - 1] * ci3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * ci3 + wa2[i - 1] * cr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * cr4 - wa3[i - 1] * ci4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * ci4 + wa3[i - 1] * cr4;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F7TCB/DF7TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 14, 1990

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F7TCB (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)

    Arguments:
       IDO    - 2*N divided by present and previous factors.  (Input)
       L1     - Product of previous factors used.  (Input)
       CC     - SEQ or CH from F3TCB.  (Input)
       CH     - A better estimate for COEF.  (Output)  Note:  The roles
                of CC and CH are reversed when NA=1.
       WA1    - Real array containing the needed elements of array WA
                om F3TCB.  (Input)
       WA2    - Real array containing other needed elements of array WA
                om F3TCB.  (Input)
       WA3    - Real array containing other needed elements of array WA
                om F3TCB.  (Input)
       WA4    - Real array containing other needed elements of array WA
                om F3TCB.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f7tcb (Mint *ido, Mint *l1, Mfloat *cc, Mfloat *ch, Mfloat wa1[],
                Mfloat wa2[], Mfloat wa3[], Mfloat wa4[])
#else
void        imsl_f7tcb (ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    Mint       *ido, *l1;
    Mfloat     *cc, *ch, wa1[], wa2[], wa3[], wa4[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(5)*(sido)+(J_)*(sido)+(K_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        i, k;
    Mfloat      ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3, di4,
                di5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3,
                tr4, tr5;
    static Mfloat tr11 = 0.309016994374947424102293417183e0;
    static Mfloat ti11 = 0.951056516295153572116439333379e0;
    static Mfloat tr12 = -0.809016994374947424102293417183e0;
    static Mfloat ti12 = 0.587785252292473129168705954639e0;



    if (sido == 2) {
	for (k = 1; k <= sl1; k++) {
	    ti5 = *CC (k - 1, 1, 1) - *CC (k - 1, 4, 1);
	    ti2 = *CC (k - 1, 1, 1) + *CC (k - 1, 4, 1);
	    ti4 = *CC (k - 1, 2, 1) - *CC (k - 1, 3, 1);
	    ti3 = *CC (k - 1, 2, 1) + *CC (k - 1, 3, 1);
	    tr5 = *CC (k - 1, 1, 0) - *CC (k - 1, 4, 0);
	    tr2 = *CC (k - 1, 1, 0) + *CC (k - 1, 4, 0);
	    tr4 = *CC (k - 1, 2, 0) - *CC (k - 1, 3, 0);
	    tr3 = *CC (k - 1, 2, 0) + *CC (k - 1, 3, 0);
	    *CH (0, k - 1, 0) = *CC (k - 1, 0, 0) + tr2 + tr3;
	    *CH (0, k - 1, 1) = *CC (k - 1, 0, 1) + ti2 + ti3;
	    cr2 = *CC (k - 1, 0, 0) + tr11 * tr2 + tr12 * tr3;
	    ci2 = *CC (k - 1, 0, 1) + tr11 * ti2 + tr12 * ti3;
	    cr3 = *CC (k - 1, 0, 0) + tr12 * tr2 + tr11 * tr3;
	    ci3 = *CC (k - 1, 0, 1) + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    *CH (1, k - 1, 0) = cr2 - ci5;
	    *CH (4, k - 1, 0) = cr2 + ci5;
	    *CH (1, k - 1, 1) = ci2 + cr5;
	    *CH (2, k - 1, 1) = ci3 + cr4;
	    *CH (2, k - 1, 0) = cr3 - ci4;
	    *CH (3, k - 1, 0) = cr3 + ci4;
	    *CH (3, k - 1, 1) = ci3 - cr4;
	    *CH (4, k - 1, 1) = ci2 - cr5;
	}
    }
    else if (sl1 <= sido / 2) {
	for (k = 1; k <= sl1; k++) {
	    for (i = 2; i <= sido; i += 2) {
		ti5 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 4, i - 1);
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 4, i - 1);
		ti4 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 3, i - 1);
		ti3 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 3, i - 1);
		tr5 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 4, i - 2);
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 4, i - 2);
		tr4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 + tr3;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 + ti3;
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
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * dr4 - wa3[i - 1] * di4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * di4 + wa3[i - 1] * dr4;
		*CH (4, k - 1, i - 2) = wa4[i - 2] * dr5 - wa4[i - 1] * di5;
		*CH (4, k - 1, i - 1) = wa4[i - 2] * di5 + wa4[i - 1] * dr5;
	    }
	}
    }
    else {
	for (i = 2; i <= sido; i += 2) {
	    for (k = 1; k <= sl1; k++) {
		ti5 = *CC (k - 1, 1, i - 1) - *CC (k - 1, 4, i - 1);
		ti2 = *CC (k - 1, 1, i - 1) + *CC (k - 1, 4, i - 1);
		ti4 = *CC (k - 1, 2, i - 1) - *CC (k - 1, 3, i - 1);
		ti3 = *CC (k - 1, 2, i - 1) + *CC (k - 1, 3, i - 1);
		tr5 = *CC (k - 1, 1, i - 2) - *CC (k - 1, 4, i - 2);
		tr2 = *CC (k - 1, 1, i - 2) + *CC (k - 1, 4, i - 2);
		tr4 = *CC (k - 1, 2, i - 2) - *CC (k - 1, 3, i - 2);
		tr3 = *CC (k - 1, 2, i - 2) + *CC (k - 1, 3, i - 2);
		*CH (0, k - 1, i - 2) = *CC (k - 1, 0, i - 2) + tr2 + tr3;
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1) + ti2 + ti3;
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
		*CH (1, k - 1, i - 2) = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
		*CH (1, k - 1, i - 1) = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
		*CH (2, k - 1, i - 2) = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
		*CH (2, k - 1, i - 1) = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
		*CH (3, k - 1, i - 2) = wa3[i - 2] * dr4 - wa3[i - 1] * di4;
		*CH (3, k - 1, i - 1) = wa3[i - 2] * di4 + wa3[i - 1] * dr4;
		*CH (4, k - 1, i - 2) = wa4[i - 2] * dr5 - wa4[i - 1] * di5;
		*CH (4, k - 1, i - 1) = wa4[i - 2] * di5 + wa4[i - 1] * dr5;
	    }
	}
    }
    return;
}				/* end of function */
#undef  CC
#undef  CH
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F8TCB/DF8TCB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 15, 1990

    Purpose:    Compute the d_complex periodic sequence from its Fourier
                coefficients.

    Usage:      CALL F8TCB (NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2,
                            WA)

    Arguments:
       NAC    - Integer set to 1 if IDO=2, set to 0 otherwise.  (Output)
       IDO    - 2*N divided by present and previous factors.  (Input)
       IP     - The factor now being used.  (Input)
       L1     - Product of previous factors used.  (Input)
       IDL1   - ID0*L1
       CC     - Matrix containing C or CH from F3TCB.  (Input/Output)
       C1     - Matrix containing C or CH from F3TCB.  (Input/Output)
       C2     - Matrix containing C or CH from F3TCB.  (Input/Output)
       CH     - Matrix that contains a better estimate for
                COEF.  (Output)
       CH2    - Matrix that contains a better estimate for
                COEF.  (Output)  Note:  The roles of CC and CH are
                reversed when NA=1.
       WA     - Real array containing the needed elements of array WA
                in F3TCB.  (Input)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_f8tcb (Mint *nac, Mint *ido, Mint *ip, Mint *l1, Mint *idl1,
                Mfloat cc[], Mfloat c1[], Mfloat c2[], Mfloat ch[],
                Mfloat ch2[], Mfloat wa[])
#else
void        imsl_f8tcb (nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2,
                wa)
    Mint       *nac, *ido, *ip, *l1, *idl1;
    Mfloat     *cc, *c1, *c2, *ch, *ch2, wa[];
#endif
{
#define CC(I_,J_,K_)	(cc+(I_)*(*ip)*(sido)+(J_)*(sido)+(K_))
#define C1(I_,J_,K_)	(c1+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define C2(I_,J_)	(c2+(I_)*(sidl1)+(J_))
#define CH(I_,J_,K_)	(ch+(I_)*(sl1)*(sido)+(J_)*(sido)+(K_))
#define CH2(I_,J_)	(ch2+(I_)*(sidl1)+(J_))
    Mint        sido = *ido;
    Mint        sl1 = *l1;
    Mint        sidl1 = *idl1;
    Mint        i, idij, idj, idl, idlj, idlj1u[512], idot, idp, ik, inc,
                ipp2, ipph, j, j5s, j6s, j7s, j8s, jc, k, l;
    Mfloat      r1s, r2s, sumi, sumr, wai, war;


    idot = sido / 2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    idp = *ip ** ido;

    if (sido >= sl1) {
	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;
	    for (k = 1; k <= sl1; k++) {
		for (i = 1; i <= sido; i++) {
		    *CH (j - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) +
			*CC (k - 1, jc - 1, i - 1);
		    *CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) -
			*CC (k - 1, jc - 1, i - 1);
		}
	    }
	}

	for (k = 1; k <= sl1; k++) {
	    for (i = 1; i <= sido; i++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}
    }
    else {
	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;
	    for (i = 1; i <= sido; i++) {
		for (k = 1; k <= sl1; k++) {
		    *CH (j - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) +
			*CC (k - 1, jc - 1, i - 1);
		    *CH (jc - 1, k - 1, i - 1) = *CC (k - 1, j - 1, i - 1) -
			*CC (k - 1, jc - 1, i - 1);
		}
	    }
	}

	for (i = 1; i <= sido; i++) {
	    for (k = 1; k <= sl1; k++) {
		*CH (0, k - 1, i - 1) = *CC (k - 1, 0, i - 1);
	    }
	}
    }

    idl = 2 - sido;
    inc = 0;
    if (sidl1 > 2) {
	for (l = 1; l <= (ipph - 1); l++) {
	    for (ik = 1; ik <= sidl1; ik++) {
		*C2 (l, ik - 1) = *CH2 (0, ik - 1) + wa[idl + sido * l - 2] *
		    *CH2 (1, ik - 1);
		*C2 (ipp2 - (l - 1) - 3, ik - 1) = wa[idl + sido * l - 1] *
		    *CH2 (*ip - 1, ik - 1);
	    }

	    idlj = idl + sido + sido * (l - 1);
	    for (j = 1; j <= (ipph - 2); j++) {
		idlj += inc + sido + sido * (l - 1);
		if (idlj > idp)
		    idlj -= idp;
		wai = wa[idlj - 1];
		war = wa[idlj - 2];
		for (ik = 1; ik <= sidl1; ik++) {
		    r1s = war ** CH2 (j + 1, ik - 1);
		    *C2 (l, ik - 1) += r1s;
		    r2s = wai ** CH2 (ipp2 - (j - 1) - 4, ik - 1);
		    *C2 (ipp2 - (l - 1) - 3, ik - 1) += r2s;
		}
	    }
	}
	if (ipph - 1 > 0) {
	    idl += (ipph - 1) ** ido;
	    inc += (ipph - 1) ** ido;
	}
	for (j = 2; j <= ipph; j++) {
	    for (ik = 1; ik <= sidl1; ik++) {
		*CH2 (0, ik - 1) += *CH2 (j - 1, ik - 1);
	    }
	}

	for (j = 2; j <= ipph; j++) {
	    jc = ipp2 - j;
	    for (ik = 2; ik <= sidl1; ik += 2) {
		*CH2 (j - 1, ik - 2) = *C2 (j - 1, ik - 2) - *C2 (jc - 1, ik - 1);
		*CH2 (jc - 1, ik - 2) = *C2 (j - 1, ik - 2) + *C2 (jc - 1, ik - 1);
		*CH2 (j - 1, ik - 1) = *C2 (j - 1, ik - 1) + *C2 (jc - 1, ik - 2);
		*CH2 (jc - 1, ik - 1) = *C2 (j - 1, ik - 1) - *C2 (jc - 1, ik - 2);
	    }
	}
    }
    else {
	for (j7s = 0; j7s <= (ipph - 2); j7s += 511) {
	    j8s = imsl_i_min (ipph - 1 - j7s, 511);
	    /*
	     * The local integer vector idlj1u was inserted by fpp on the
	     * Cray to imsl_aid in the vectorization of the code. DIR$
	     * IVDEP
	     */
	    for (l = 1; l <= j8s; l++) {
		*C2 (j7s + l, 0) = *CH2 (0, 0) + wa[idl + (j7s + l) ** ido - 2] *
		    *CH2 (1, 0);
		*C2 (ipp2 - j7s - l - 2, 0) = wa[idl + (j7s + l) ** ido - 1] *
		    *CH2 (*ip - 1, 0);
		*C2 (j7s + l, 1) = *CH2 (0, 1) + wa[idl + (j7s + l) ** ido - 2] *
		    *CH2 (1, 1);
		*C2 (ipp2 - j7s - l - 2, 1) = wa[idl + (j7s + l) ** ido - 1] *
		    *CH2 (*ip - 1, 1);

		idlj1u[l - 1] = idl + sido + sido * (l - 1 + j7s);
	    }
	    for (j = 1; j <= (ipph - 2); j++) {
		/* DIR$          IVDEP */
		for (l = 1; l <= j8s; l++) {
		    idlj1u[l - 1] += inc + sido + sido * (l - 1 + j7s);
		    if (idlj1u[l - 1] > idp)
			idlj1u[l - 1] -= idp;
		    j5s = idlj1u[l - 1];
		    j6s = idlj1u[l - 1];
		    wai = wa[j6s - 1];
		    war = wa[j5s - 2];
		    *C2 (j7s + l, 0) += war ** CH2 (j + 1, 0);
		    *C2 (ipp2 - j7s - l - 2, 0) += wai ** CH2 (ipp2 - (j - 1) - 4, 0);
		    *C2 (j7s + l, 1) += war ** CH2 (j + 1, 1);
		    *C2 (ipp2 - j7s - l - 2, 1) += wai ** CH2 (ipp2 - (j - 1) - 4, 1);
		}
	    }
	}

	sumr = F_ZERO;
	sumi = F_ZERO;
	for (j = 1; j <= ipph; j++) {
	    sumr += *CH2 (j - 1, 0);
	    sumi += *CH2 (j - 1, 1);
	}
	*CH2 (0, 0) = sumr;
	*CH2 (0, 1) = sumi;

	/* DIR$    IVDEP */
	for (j = 1; j <= (ipph - 1); j++) {
	    *CH2 (j, 0) = *C2 (j, 0) - *C2 (ipp2 - j - 2, 1);
	    *CH2 (ipp2 - j - 2, 0) = *C2 (j, 0) + *C2 (ipp2 - j - 2, 1);
	    *CH2 (j, 1) = *C2 (j, 1) + *C2 (ipp2 - j - 2, 0);
	    *CH2 (ipp2 - j - 2, 1) = *C2 (j, 1) - *C2 (ipp2 - j - 2, 0);
	}
    }
    *nac = 1;

    if (sido == 2)
	goto L_360;
    *nac = 0;
    for (ik = 1; ik <= sidl1; ik++) {
	*C2 (0, ik - 1) = *CH2 (0, ik - 1);
    }

    if (*ip - 1 < sl1) {
	for (j = 2; j <= *ip; j++) {
	    for (k = 1; k <= sl1; k++) {
		*C1 (j - 1, k - 1, 0) = *CH (j - 1, k - 1, 0);
		*C1 (j - 1, k - 1, 1) = *CH (j - 1, k - 1, 1);
	    }
	}
    }
    else {
	for (k = 1; k <= sl1; k++) {
	    for (j = 2; j <= *ip; j++) {
		*C1 (j - 1, k - 1, 0) = *CH (j - 1, k - 1, 0);
		*C1 (j - 1, k - 1, 1) = *CH (j - 1, k - 1, 1);
	    }
	}
    }

    if (idot <= sl1) {
	idij = 0;
	for (j = 2; j <= *ip; j++) {
	    idij += 2;
	    for (i = 1; i <= ((sido - 2) / 2); i++) {
		for (k = 1; k <= sl1; k++) {
		    *C1 (j - 1, k - 1, (i - 1) * 2 + 2) = wa[idij + (i - 1) * 2] *
			*CH (j - 1, k - 1, (i - 1) * 2 + 2) - wa[idij + i * 2 - 1] *
			*CH (j - 1, k - 1, (i - 1) * 2 + 3);
		    *C1 (j - 1, k - 1, (i - 1) * 2 + 3) = wa[idij + (i - 1) * 2] *
			*CH (j - 1, k - 1, (i - 1) * 2 + 3) + wa[idij + i * 2 - 1] *
			*CH (j - 1, k - 1, (i - 1) * 2 + 2);
		}
	    }
	    if ((sido - 2) / 2 > 0)
		idij += (sido - 2) / 2 * 2;
	}
    }
    else {
	idj = 2 - sido;
	for (j = 1; j <= (*ip - 1); j++) {
	    for (k = 1; k <= sl1; k++) {
		for (i = 1; i <= ((sido - 2) / 2); i++) {
		    *C1 (j, k - 1, i * 2) = wa[idj + sido * j + i * 2 - 2] *
			*CH (j, k - 1, i * 2) - wa[idj + sido * j + i * 2 - 1] *
			*CH (j, k - 1, (i + 1) * 2 - 1);
		    *C1 (j, k - 1, (i + 1) * 2 - 1) = wa[idj + sido * j + i * 2 - 2] *
			*CH (j, k - 1, (i + 1) * 2 - 1) + wa[idj + sido * j + i * 2 - 1] *
			*CH (j, k - 1, i * 2);
		}
	    }
	}
    }
L_360:
    return;
}				/* end of function */
#undef CC
#undef C1
#undef C2
#undef CH
#undef CH2
