#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static void l_f3tci (Mint *n, Mfloat *wa, Mfloat *imsl_fac);
#else
static void l_f3tci ();
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/30/90 at 11:20:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  FFTCI/DFFTCI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute parameters needed by FFTCF and FFTCB.

    Usage:      CALL FFTCI (N, WFFTC)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WFFTC  - Array of length 4*N+15 containing parameters needed by
                FFTCF and FFTCB.  (Output)

    Remark:
       Different WFFTC arrays are needed for different values of N.

    Keywords:   Complex exponential FFT; Initialization

    GAMS:       J1a2

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat     *imsl_c_fft_complex_init (Mint n)
#else
Mfloat     *imsl_c_fft_complex_init (n)
    Mint        n;
#endif
{
    Mfloat     *wfftc = NULL;
    E1PSH ("imsl_c_fft_complex_init", "imsl_z_fft_complex_init");
    /* CHECK ARGUMENT N */
    if (n < 1) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	/* The length of the sequence N = %(i1). */
	/* It must be at least 1. */
	goto L_9000;
    }
    if (n > 1) {
	wfftc = (Mfloat *) imsl_malloc ((4 * n + 15) * sizeof (*wfftc));
	if (wfftc == NULL) {
	    /* Not enough memory, with %(L1) = %(I1). */
	    imsl_e1sti (1, n);
	    imsl_e1stl (1, "n");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	}
	else {
	    sset (4 * n + 15, F_ZERO, wfftc, 1);
	    l_f3tci (&n, &wfftc[n * 2], &wfftc[n * 4]);
	}
    }
L_9000:
    ;
    if (imsl_n1rty (0) > 3) {
	if (wfftc != NULL)
	    imsl_free (wfftc);
	wfftc = NULL;
    }
    E1POP ("imsl_c_fft_complex_init", "imsl_z_fft_complex_init");
    return wfftc;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F3TCI/DF3TCI (Single/Double precision version)

    Computer:   FORC/SINGLE

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
static void l_f3tci (Mint *n, Mfloat *wa, Mfloat *imsl_fac)
#else
static void l_f3tci (n, wa, imsl_fac)
    Mint       *n;
    Mfloat      wa[], imsl_fac[];
#endif
{
    Mint        i, i1, ido, idot, ii, ip, ipm, j, k1, l1, l2, ld, nf, nl,
                nq, nr, ntry;
    Mfloat      arg, argh, argld, fi, tpi;
    static Mint ntryh[4] = {3, 4, 2, 5};


    /*
     * FAC WILL BE A VECTOR OF LENGTH 14. FAC(1) WILL CONTAIN N, FAC(2) THRU
     * FAC(14) WILL HAVE THE PRIME FACTORS OF N, IN ASCENDING ORDER. WA WILL
     * BE A VECTOR OF LENGTH 2*N, CONTAINING THE REAL AND IMAGINARY PARTS OF
     * CEXP(2*PI*?/N)
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
    imsl_fac[nf + 1] = (Mfloat) (ntry);
    nl = nq;
    if (ntry != 2)
	goto L_60;
    if (nf == 1)
	goto L_60;
    imsl_scopy (nf - 1, &imsl_fac[2], -1, &imsl_fac[3], -1);
    imsl_fac[2] = F_TWO;
L_60:
    if (nl != 1)
	goto L_40;
    imsl_fac[0] = (Mfloat) (*n);
    imsl_fac[1] = (Mfloat) (nf);
    tpi = F_TWO * 3.1415926535897932384626433831;
    argh = tpi / (Mfloat) (*n);
    i = 2;
    l1 = 1;
    for (k1 = 1; k1 <= nf; k1++) {
	ip = nint (imsl_fac[k1 + 1]);
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	idot = ido + ido + 2;
	ipm = ip - 1;
	for (j = 1; j <= ipm; j++) {
	    i1 = i;
	    wa[i - 2] = F_ONE;
	    wa[i - 1] = F_ZERO;
	    ld += l1;
	    fi = F_ZERO;
	    argld = (Mfloat) (ld) * argh;
	    for (ii = 4; ii <= idot; ii += 2) {
		i += 2;
		fi += F_ONE;
		arg = fi * argld;
		wa[i - 2] = cos (arg);
		wa[i - 1] = sin (arg);
	    }
	    if (ip > 5) {
		wa[i1 - 2] = wa[i - 2];
		wa[i1 - 1] = wa[i - 1];
	    }
	}
	l1 = l2;
    }
    return;
}				/* end of function */
