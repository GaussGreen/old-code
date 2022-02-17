#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static void l_f3tri (Mint *n, Mfloat *wa, Mfloat *fac);
#else
static void l_f3tri ();
#endif
/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 10:09:44
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  FFTRI/DFFTRI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute parameters needed by FFTRF and FFTRB.

    Usage:      CALL FFTRI (N, WFFTR)

    Arguments:
       N      - Length of the sequence to be transformed.  (Input)
       WFFTR  - Array of length 2*N+15 containing parameters needed by
                FFTRF and FFTRB.  (Output)

    Remark:
       Different WFFTR arrays are needed for different values of N.

    Keywords:   Initialization; Transforms; Trigonometric

    GAMS:       J1a1

    Chapters:   MATH/LIBRARY Transforms
                STAT/LIBRARY Mathematical Support

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat     *imsl_f_fft_real_init (Mint n)
#else
Mfloat     *imsl_f_fft_real_init (n)
    Mint        n;
#endif
{
    Mfloat     *wfftr = NULL;
    E1PSH ("imsl_f_fft_real_init", "imsl_d_fft_real_init");
    /* CHECK ARGUMENT N */
    if (n < 1) {
	/* The length of the sequence n = %(I1). */
	/* It must be at least one. */
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_SEQUENCE_LENGTH);
	goto L_9000;
    }
    if (n > 1) {
	wfftr = (Mfloat *) imsl_malloc ((2 * n + 15) * sizeof (*wfftr));
	if (wfftr == NULL) {
	    /* Not enough memory, with %(L1) = %(I1). */
	    imsl_e1sti (1, n);
	    imsl_e1stl (1, "n");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	}
	else {
	    sset (2 * n + 15, F_ZERO, wfftr, 1);
	    l_f3tri (&n, &wfftr[n], &wfftr[n * 2]);
	}
    }
L_9000:
    ;
    if (imsl_n1rty (0) > 3) {
	if (wfftr != NULL)
	    imsl_free (wfftr);
	wfftr = NULL;
    }
    E1POP ("imsl_f_fft_real_init", "imsl_d_fft_real_init");
    return (wfftr);
}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 06/11/90 at 09:58:39
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F3TRI/DF3TRI (Single/Double precision version)

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
static void l_f3tri (n, wa, fac)
    Mint       *n;
    Mfloat      wa[], fac[];
{
    Mint        i, ido, ii, ip, ipm, is, j, k1, l1, l2, ld, nf;
    Mint        nfm1, nl, nq, nr, ntry;
    Mfloat      arg, argh, argld, fi, tpi;
    static Mint ntryh[4] = {4, 2, 3, 5};
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
    fac[nf + 1] = (Mfloat) (ntry);
    nl = nq;
    if (ntry != 2)
	goto L_60;
    if (nf == 1)
	goto L_60;
    imsl_scopy (nf - 1, &fac[2], -1, &fac[3], -1);
    fac[2] = F_TWO;
L_60:
    if (nl != 1)
	goto L_40;
    fac[0] = (Mfloat) (*n);
    fac[1] = (Mfloat) (nf);
    tpi = F_TWO * 3.1415926535897932384626433831;
    argh = tpi / (Mfloat) (*n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0)
	return;
    for (k1 = 1; k1 <= nfm1; k1++) {
	ip = nint (fac[k1 + 1]);
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	for (j = 1; j <= ipm; j++) {
	    ld += l1;
	    i = is;
	    argld = (Mfloat) (ld) * argh;
	    fi = F_ZERO;
	    for (ii = 3; ii <= ido; ii += 2) {
		i += 2;
		fi += F_ONE;
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
